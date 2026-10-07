#!/usr/bin/env Rscript
# Revision step 3: per-sample doublet detection (scDblFinder) and ambient-RNA correction
# (SoupX with the empty-droplet soup profile; DecontX with the empty-droplet background).
# Input: <work_dir>/<GSM>/{cells_counts.mtx, cells.csv, genes.csv, background.mtx} and
#        <work_dir>/soup_profiles_counts.csv written by scripts/11_sc_molecule_info_soup.py
# Output per sample: doublets.csv, soupx_counts.mtx, decontx_counts.mtx, decontx_contamination.csv
# and <work_dir>/ambient_summary.csv
suppressPackageStartupMessages({
  library(Matrix); library(SingleCellExperiment); library(scDblFinder)
  library(SoupX); library(celda)
})
args <- commandArgs(trailingOnly = TRUE)
work <- args[1]
soup_all <- read.csv(file.path(work, "soup_profiles_counts.csv"), row.names = 1, check.names = FALSE)
gsms <- list.dirs(work, full.names = FALSE, recursive = FALSE)
gsms <- gsms[grepl("^GSM", gsms)]
summ <- list()
for (g in gsms) {
  d <- file.path(work, g)
  M <- as(readMM(file.path(d, "cells_counts.mtx")), "CsparseMatrix")
  genes <- read.csv(file.path(d, "genes.csv"))$gene
  cells <- read.csv(file.path(d, "cells.csv"), colClasses = "character")
  rownames(M) <- make.unique(genes); colnames(M) <- cells$barcode
  cat(g, ":", ncol(M), "cells\n")

  # ---- doublets: scDblFinder (random artificial doublets; default rate ~1% per 1000 cells)
  set.seed(0)
  sce <- scDblFinder(SingleCellExperiment(list(counts = M)), verbose = FALSE)
  dbl <- data.frame(barcode = colnames(M), scDblFinder_class = sce$scDblFinder.class,
                    scDblFinder_score = sce$scDblFinder.score)
  write.csv(dbl, file.path(d, "doublets.csv"), row.names = FALSE)

  # ---- SoupX with the empty-droplet soup profile (droplets <=100 UMIs, non-cell)
  soup <- soup_all[[g]]; names(soup) <- rownames(soup_all)
  soup <- soup[rownames(M) %in% names(soup) | TRUE]
  sp <- setNames(rep(0, nrow(M)), rownames(M))
  common <- intersect(names(soup), rownames(M))
  sp[common] <- soup[common]
  sc <- SoupChannel(tod = M, toc = M, calcSoupProfile = FALSE)
  sc <- setSoupProfile(sc, data.frame(est = sp / sum(sp), counts = sp, row.names = rownames(M)))
  sc <- setClusters(sc, setNames(cells$leiden, cells$barcode))
  rho <- tryCatch({
    sc <- autoEstCont(sc, doPlot = FALSE, verbose = FALSE); sc$fit$rhoEst
  }, error = function(e) { cat("autoEstCont failed:", conditionMessage(e), "-> rho=0.10\n"); NA })
  if (is.na(rho)) sc <- setContaminationFraction(sc, 0.10)
  rho_used <- mean(sc$metaData$rho)
  sx <- adjustCounts(sc, roundToInt = TRUE, verbose = 0)
  writeMM(sx, file.path(d, "soupx_counts.mtx"))

  # ---- DecontX with the empty-droplet background (droplets 10-100 UMIs)
  B <- as(readMM(file.path(d, "background.mtx")), "CsparseMatrix")
  rownames(B) <- rownames(M)
  dx <- decontX(x = M, z = cells$leiden, background = B, seed = 0, verbose = FALSE)
  writeMM(round(dx$decontXcounts), file.path(d, "decontx_counts.mtx"))
  write.csv(data.frame(barcode = colnames(M), decontX_contamination = dx$contamination),
            file.path(d, "decontx_contamination.csv"), row.names = FALSE)

  summ[[g]] <- data.frame(gsm = g, n_cells = ncol(M),
                          n_doublets_scDblFinder = sum(dbl$scDblFinder_class == "doublet"),
                          soupx_rho_est = rho, soupx_rho_used = rho_used,
                          decontx_median_contamination = median(dx$contamination))
  print(summ[[g]])
}
write.csv(do.call(rbind, summ), file.path(work, "ambient_summary.csv"), row.names = FALSE)
cat("done\n")
