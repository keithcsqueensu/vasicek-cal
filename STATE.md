# State

Last updated: 2026-10-12

The current position only. Decisions and their rationale are in [DECISIONS.md](DECISIONS.md); the
detailed record of M0 to M3 and the earlier plans are in [docs/history.md](docs/history.md).

## Now

- **M6, the AR(1)-aware binomial MLE: in progress.**
  - Merged: the likelihood (forward–backward–forward filter, φ = 0 bit for bit; D-195), the search and
    3-D profile intervals with the S-35 M6 arm registered first (D-196, PR #50), the validation
    (likelihood goldens, a two-formulation scipy replication, the search's recovery tests,
    `docs/methodology/ar1_binomial_mle.md`; D-197), the tool `study_lrpd_serial_m6` (D-198), P-18 and
    the order of the performance pass (D-199), and the validation run's registration (W1–W9 and
    `score.py` committed before the run; D-200, PR #54).
  - **The M6 validation run: done (D-202).** 7 of 8 strict predictions held; **W3 NOT HELD**, and the
    verdict stands. D-202 explains the miss without changing it: at φ = 0 the cross-information is
    exactly zero, so W3's denominator, sqrt(v₃ − v₂), was Monte Carlo noise in the information file;
    the paired difference itself is second order, as theory says. No estimator fix is indicated. The
    arm's registration has no prediction with an information-based denominator, so it is unchanged.
    `rows.csv` is on the release `m6-validation-data-2026-10-10`.
  - **The S-35 M6 arm: done (D-208).** 9 of 11 strict predictions held, 3 of 4 qualitative; NOT HELD:
    L21, L28 and L25's qualitative part.
    - **The headline:** the AR(1)-aware MLE restores PD coverage: 0.946 and 0.942 at T = 100, against
      0.69–0.73 for the iid MLE at φ = 0.6, a gain of about +0.22 at φ = 0.6 and almost none at φ = 0.
    - **L21:** by the registration's own alternative, the AR(1) model moves the median PD toward 0,
      consistently but by small amounts.
    - **L28:** the percentile bootstrap undercovers, consistent with S-10 and D-176.
    - **L25:** the bound is reachable at T = 20, so there is no implementation finding.
    - **Runtime:** the projection landed in range again (34.94 h against 32.9–37.0 h).
    - **Data:** `rows.csv` is on the release `s35-m6-data-2026-10-12`.
- **S-48, the HAC SE: done (D-211).** 5 of 6 held; **H2 NOT HELD** (ρ's coverage below the reference
  at φ = 0.6 in groups B and C).
  - For PD the HAC interval is a partial correction that grows with T: at φ = 0.6, 0.721, 0.805 and 0.852
    at T = 20, 40 and 100, against the iid Wald's 0.65–0.72 and the AR(1) MLE's 0.91–0.94. It costs
    0.017–0.043 under independence.
  - The HAC SE stays a `native` diagnostic for PD only, reported beside the iid SE, not added to the C
    ABI. HAC is not recommended for ρ, where it covers worse than the iid Wald.
  - S-48's negative-dependence part is registered with S-36 (D-401, #60).
  - The next decision number is D-214 (P-7 is D-212, #67; the SIMD sweep took D-213, #68).
  - **Runtime:** 304 s, in the projected range. The projection's sample ran 84 s, under D-193's 120 s
    minimum: a deviation (D-211). `rows.csv` is on the release `s48-data-2026-10-12`.

## Next

1. **Copilot-audit fix PRs**, one each:
   - **A, correctness;**
   - **B, data integrity and docs accuracy** (the MANIFEST `sha256` check was split out and is done:
     D-201, PR #55);
   - **C, cosmetic,** plus converting the study scripts' engine flags to `tests/engine_flags.json`.
2. **The performance pass, in the owner's order** (D-199): P-7 (`erfcx`), P-13, P-6; P-19 with the next
   golden refresh.
3. **The rate MLE as a first-class estimator** (D-177): analytic score and Hessian with Newton solves
   (it still uses D-119's central differences, D-185); interval coverage; an n·PD adequacy flag from
   S-27; rate-vs-count disagreement as a data diagnostic; its own data-sufficiency table.
4. **`docs/methodology/estimator_selection.md`** (D-177): comparison matrix, tiers with reasons, a
   selection rule fixed in advance. Not started.
5. **Z_t extraction,** E[Z_t | d_t] with its uncertainty, in adverse/benign terms (S-19, pulled forward
   from M7; D-177). Not started.
6. **Consolidation:** tag v0.2 (D-168, D-177); the README already carries the study results.
7. **The rest of M6:** multi-grade calibration with a shared ρ (D-151; S-29, S-30, S-40), multi-segment
   factors, alternative mixing distributions, Pluto–Tasche, moving-block coverage under AR(1) data,
   R-1. Then M4 (D-168).

**Backlog** (`native` options, each judged against the pinned recovery verdicts; D-131, D-137):
- a bias-corrected ρ̂ from exact delete-one refits, reported beside ρ̂ (S-3, D-159);
- BCa intervals: better than the percentile interval for ρ, no gain for PD, no fix near a bound (S-5,
  D-159); not recommended over Jeffreys or profile (D-176);
- the stationary bootstrap and a data-driven block length (Politis–White), with M6's AR(1) data;
- P-11 part 2 (pinning runs write goldens directly; D-179); P-20 and the in-house `log`/`exp` for
  `log_phi`, each a decision of its own (D-199).

## Waiting on the owner

None open. Recorded for later, not blocking:
the data repository for S-47's public-data part needs its own D-entry (D-188).

## Milestones

| M | Scope | Status | Key decisions |
|---|---|---|---|
| M1 | Core math, `ref/`, DGP, cross-reference, recovery harness, scipy replication; binomial MLE on the CPU | **done**, tag `v0.1.0-parity-binomial` | D-025–D-127 |
| M2 | Profile intervals, resampling engine (W × L), C ABI | **done** | D-128–D-144 |
| M3 | Vasicek-rate MLE, method of moments, grid Bayesian; C ABI 0.3, then 0.4 (S-4b) | **done** | D-164–D-166, D-169, D-184 |
| M4 | CUDA backend, CPU↔GPU parity suite | after M6; nvcc compile-only CI exists | D-168, D-100 |
| M5 | Monitoring: backtests, run ledger, report data, threshold table | planned | D-152 |
| M6 | Benchmark models: AR(1) factor, multi-grade, multi-segment, mixing distributions, Pluto–Tasche | **in progress** (AR(1)-aware MLE) | D-151, D-195–D-200, D-208 |
| M7 | Macro pipeline: transforms, pre-screens, batched fits, scorers, ranker | planned | D-148, R-2 |
| M8 | Public-data demonstrations, scipy replication set, ABI 1.0 | planned | D-046, D-188 |

## Performance

A track separate from the studies (D-177): P-items change how results are computed, never what. Each has
a decision, a benchmark in `perf/` and an identity or tolerance test; result changes need a golden
refresh with a flip check. Every item with its benchmark and write-up: [perf/README.md](perf/README.md).

| ID | Item | Status | Decision | Changes results? |
|---|---|---|---|---|
| P-1 | Surface rows cached by (n, d) | done | D-167 | no |
| P-2 | Profile and q solves by safeguarded Newton | done | D-170 | within solver tolerance; 0 flips |
| P-3 | Study fits compute only scored arms | done | D-171 | no |
| P-4 | Bootstrap reduction over distinct rows | done | D-172 | within rounding; 0 flips |
| P-5 | Bounded argmax, full-grid fallback | done | D-172 | no |
| P-6 | SIMD reduction across grid points | planned, fourth (5–7% of the run) | D-199 | no |
| P-7 | SIMD `erfcx` (narrowed; `log_phi_pair` extension done) | next | D-199 | no |
| P-8 | GPU backend | M4 | D-168 | within tolerance |
| P-9 | Memoise identical panels | measured, not built (~1%) | D-189 | no |
| P-10 | Dynamic scheduling across threads | planned | D-186 | no |
| P-11 | Checkpoint and resume | part 1 done; part 2 open | D-179 | no |
| P-12 | Warm-started bootstrap search | **rejected** (competes with P-6; extra fallback path) | D-199 | — |
| P-13 | LTO/PGO for local study builds | planned, third | D-199 | no |
| P-14 | Calibrate's checks per distinct observation | done | D-178 | no |
| P-15 | Calibrate's SEs from the analytic Hessian | done | D-185 | within tolerance (a correction) |
| P-16 | Sampled, reported quadrature check on bootstrap fits | done | D-180 | no |
| P-17 | SIMD Philox for the DGP's Bernoulli draws | done | D-194 | no |
| P-18 | AR(1) transition sums on kernels with shared `det_exp` | done | D-199 | within rounding (≤ 7.2e-15) |
| P-19 | Analytic SE shared with calibrate (~13% of S-10) | future, with next golden refresh | D-199 | within rounding; flip check |
| P-20 | Warm-start S-4b's W* refits (40% of S-10) | future | D-199 | yes; own flip check |
| P-21 | One `erfcx` per node on the derivative path | proposed | D-213 | no |
| P-22 | SIMD derivative node functions | proposed (after P-7, P-21) | D-213 | no |
| P-23 | Hint and extent search across grid points | proposed | D-213 | no |

Also future (D-199): the in-house `log`/`exp` for `log_phi`, a decision of its own (within rounding;
would make results platform-independent). Rejected: sequential stopping of Monte Carlo runs (a design
change, not a performance item). Runtimes are projected from a timing sample before every full run
(D-187, D-193); the method is **validated** on the M6 validation run (3.499 h against 3.3–3.7 h; D-202),
and landed in range again on the S-35 M6 arm (34.94 h against 32.9–37.0 h; D-208).

## Studies

The index, registrations and results are in [studies/README.md](studies/README.md) (D-148, D-153).

- **Finished:** S-1, S-2, S-3, S-4, S-5, S-8, S-9, S-10, S-13 (targeted part), S-15, S-21, S-23, S-27,
  S-28, S-34, S-35 (main arms, D-192; M6 arm, D-208: 9 of 11 strict held, 3 of 4 qualitative; L21, L28
  and L25's qualitative part NOT HELD), S-48's positive-dependence and independence parts (D-211: 5 of 6
  held, H2 NOT HELD), the M6 validation run (`m6-validation`: 7 of 8 strict held, W3
  NOT HELD; D-200, D-202).
- **Registered, not yet run:** none.
- **Planned (not registered):** S-6, S-7, S-11, S-12, S-13 (full), S-14, S-16 to S-20, S-22, S-24 to
  S-26, S-29 to S-33, S-36, S-38 to S-48; S-37 is documentation, not a study (D-188).

**Roadmap addendum (D-188, merged):** references in [docs/references.md](docs/references.md).

| Study | Placement | Note |
|---|---|---|
| S-35 | now; M6 arm | finished (D-192; M6 arm D-208) |
| S-36 | M6 | new generator code; registered (D-401, #60): two momentum arms, Q1–Q13; not run |
| S-37 | M7, documentation | an identity (D-164's closed form); **documented** in the rate MLE note, §3 |
| S-38 | now (with S-19); M7 | cheap |
| S-39 | M6 with M7 | plug-in factor |
| S-40 | M6, in multi-grade | |
| S-41 | M6 with M7 | one experiment with S-47's synthetic part |
| S-42 | M7 | the ranker's test bed |
| S-43 | M6 | a note with a simulation check |
| S-44 | M7 | likely a numerical check |
| S-45 | after M7 | with the lifetime PD feature |
| S-46 | after M5 | lowest priority |
| S-47 | synthetic M6/M7; public data M8 | needs a data-repository decision |
| S-48 | done; the S-36 part with S-36 | finished (D-211): 5 of 6 held, H2 NOT HELD |

Features: the published-method track (M6 step 2, M7 steps 1 and 3) and the lifetime PD term structure
(after M7).

## Open items / risks

- **R-1:** GPU-side DGP, revisit at M6 (large AR(1) parametric bootstraps; D-043).
- **R-2:** Z sign convention and macro sign filters, revisit at M7 (S-18).
- **PR B (queued):** `tools/check_manifests.py` lists files from the index but reads them at HEAD, so a
  staged, uncommitted manifest crashes it; it should read what it lists. CI is unaffected.
- **Toolchain floors:** a portability goal: CUDA ≥ 11.8, GCC 11, MSVC 14.39, device-visible code in
  C++17 (D-094); the development machine needs a CUDA toolkit before M4.

## Where history lives

- [DECISIONS.md](DECISIONS.md): every decision and its rationale (D-nnn).
- `git log`: what changed and when.
- [docs/history.md](docs/history.md): the detailed M0–M3 status record, the earlier priority lists and
  the full milestone plan as it stood on 2026-10-10.
