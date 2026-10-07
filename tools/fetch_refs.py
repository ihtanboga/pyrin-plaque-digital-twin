"""Fetch PubMed metadata for every cited PMID and write a Vancouver-style reference database
(tools/refs.json: key -> formatted reference). Keys are the citation keys used in the manuscript."""
import json
import re
import sys
import time
import urllib.request

REFS = {
    # original references
    "grebe2018": "29880500", "duewell2010": "20428172", "ridker2017": "28845751", "tardif2019": "31733140",
    "nidorf2020": "32865380", "heilig2018": "29148036", "xu2014": "24919149", "malik2020": "32721043",
    "yan2007": "17850490", "park2016": "27270401", "yu2007": "17964261", "alsaigh2022": "36224302",
    "mahmoud2019": "31339449", "edgar2002": "11752295", "wolf2018": "29409532", "luecken2019": "31217225",
    "traag2019": "30914743", "virtanen2020": "32015543", "shi2015": "26375003", "liu2016": "27383986",
    "depuydt2020": "32981416", "fernandez2019": "31591603",
    # added in revision: pyrin in disease (Reviewer 3)
    "frenchfmf1997": "9288094", "intlfmf1997": "9288758", "schnappauf2019": "31456795", "chae2011": "21600797",
    "masters2016": "27030597", "gao2016": "27482109", "vangorp2016": "27911804", "akula2016": "27270400",
    "centola2000": "10807793",
    # pyrin / FMF and atherosclerosis or coronary disease (Reviewer 3)
    "merashli2022": "35933450", "basar2017": "24702757", "bagheri2018": "29707173", "langevitz2001": "11344818",
    "erken2018": "29051974", "jolly2025": "39555823",
    # cytokine release vs lysis (Reviewer 2)
    "evavold2018": "29195811", "saeki2020": "32592165", "shi2014": "25119034", "kayagaki2015": "26375259",
    "wang2017": "28459430", "kayagaki2021": "33472215",
    # methods added in revision
    "germain2021": "35814628", "wolock2019": "30954476", "young2020": "33367645", "yang2020": "32138770",
    "korsunsky2019": "31740819", "marino2008": "18572196", "heumos2023": "37002403",
    "zernecke2023": "36190844", "willemsen2020": "32003464",
}


def fetch(pmids):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
           + ",".join(pmids))
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)["result"]


def fmt(r, pmid):
    au = [a["name"] for a in r.get("authors", []) if a.get("authtype", "Author") == "Author"]
    coll = [a["name"] for a in r.get("authors", []) if a.get("authtype") == "CollectiveName"]
    if not au and coll:
        authors = coll[0]
    elif len(au) > 6:
        authors = ", ".join(au[:6]) + ", et al"
    else:
        authors = ", ".join(au)
    title = re.sub(r"<[^>]+>", "", r["title"]).strip()
    if not title.endswith((".", "?", "!")):
        title += "."
    src = r.get("source", "").replace(".", "")
    year = r.get("pubdate", "")[:4]
    vol, iss, pages = r.get("volume", ""), r.get("issue", ""), r.get("pages", "")
    doi = next((a["value"] for a in r.get("articleids", []) if a["idtype"] == "doi"), "")
    s = f"{authors}. {title} {src}. {year}"
    if vol:
        s += f";{vol}"
    if pages:
        s += f":{pages}"
    s += "."
    if doi:
        s += f" doi:{doi}."
    s += f" PMID:{pmid}."
    return s


def main():
    out = {}
    keys = list(REFS)
    for i in range(0, len(keys), 40):
        chunk = keys[i:i + 40]
        res = fetch([REFS[k] for k in chunk])
        for k in chunk:
            out[k] = fmt(res[REFS[k]], REFS[k])
        time.sleep(0.5)
    out["geo_gse120521"] = ("National Center for Biotechnology Information. GSE120521: RNA-seq of stable and unstable "
                            "sections of human atherosclerotic plaques. Gene Expression Omnibus [Internet]. "
                            "Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120521.")
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "tools/refs.json", "w"), indent=1, ensure_ascii=False)
    for k, v in out.items():
        print(k, "|", v[:150])


if __name__ == "__main__":
    main()
