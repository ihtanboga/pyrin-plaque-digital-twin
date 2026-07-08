# Day-4 ODE Model Design — Minimal Pyrin-Adapted Pyroptosis Scaffold

## Scope and disclaimers
This is a **minimal baseline scaffold only**. It is **dimensionless, hypothesis-generating, NOT calibrated to patients, and NOT a clinical prediction model.** The NLRP3 pyroptosis ODE literature is used purely as a **structural modeling template**; adopting a similar ODE shape does **not** assert that pyrin and NLRP3 are mechanistically equivalent. Full scenarios and sensitivity analysis are Day-5 work.

## State variables (8)
| # | Variable | Meaning |
|---|----------|---------|
| 0 | P_inactive | inactive pyrin sensor pool |
| 1 | P_active | active pyrin sensor |
| 2 | ASC_complex | ASC speck / adaptor complex |
| 3 | CASP1_active | active caspase-1 |
| 4 | GSDMD_N | cleaved gasdermin-D N-terminal (pore) |
| 5 | IL1B_external | externalized IL-1β (proxy) |
| 6 | IL18_external | externalized IL-18 (proxy) |
| 7 | rupture_proxy | membrane rupture / pyroptosis proxy |

## Parameters
priming_input, pyrin_activation_threshold, pyrin_activation_rate, ASC_recruitment_rate, CASP1_activation_rate, GSDMD_cleavage_rate, cytokine_release_rate_IL1B, cytokine_release_rate_IL18, rupture_rate, inhibition_factor, decay, saturation. All in `config/ode_parameters.yaml`.

## Structure
- Pyrin activation is **gated**: drive = max(0, priming_input − pyrin_activation_threshold), modulated by (1 − inhibition_factor). This encodes pyrin's threshold/priming behavior distinct from NLRP3's NEK7/K+-efflux triggering.
- Production terms use a Hill-like saturation `x/(1+x/ceiling)` to stay bounded and nonnegative.
- Linear first-order `decay` on active species.
- Signal propagates strictly forward: P_active → ASC → CASP1 → GSDMD_N → {IL1B, IL18, rupture}.

## Numerics
LSODA (scipy solve_ivp), rtol 1e-8, atol 1e-10, t∈[0,50], 501 points. Tiny negative numerical excursions are clipped to 0 (documented in `simulate`).

## Baseline result
Sequential cascade peaks (P_active early, cytokines/rupture later), all nonnegative, no NaN/inf. 13/13 unit tests pass including an inhibition smoke check (higher inhibition_factor lowers peak IL1B).

## Day-5 (NOT run today)
pyrin_sensitized (lower activation threshold), intervention (inhibition_factor>0), parameter sensitivity/uncertainty sweep, and linking cytokine outputs to the evidence graph as testable IL-1β/IL-18 hypotheses.
