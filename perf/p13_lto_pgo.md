# P-13: link-time and profile-guided optimisation for local study builds (D-403)

## Expected, written before the build (2026-10-13)

- **Scope (D-192 addendum, D-199):** local study builds only. Options, off by default; CI's presets, the shared
  library and every golden are unchanged.
- **Where LTO can act:** almost nowhere. `core/`, `engine/`, `resample/` and `backend_cpu/` are header-only, so
  each study tool's hot code (the quadrature, the surface fill, the searches) is already one translation unit.
  The only calls across units go to the compiled kernels (`vcal_core_kernels`: erfcx; `vcal_dgp`: the DGP's
  draws; `vcal_benchmark_models`: the AR(1) transition) and the test harness. They are called once per batch or
  per panel, not per node. **Expected: 0-3%.**
- **Where PGO can act:** branch layout and inlining choices in code with data-dependent branches, such as the
  integrand's one-sided and interior paths, the hint's bracket and Newton, the composite rule's extent bisection,
  the argmax's tile visits and Newton's polish. **Expected: 3-10%** on the study tools, from a training run on
  the recovery subset.
- **Stop rule:** if LTO + PGO together give under 3% on both measured tools, P-13 is recorded as measured and
  the options are not added.
- **Bit for bit:** neither changes floating-point semantics. Contraction stays off (no `-mfma`; GCC and Clang
  in ISO mode; MSVC `/fp:precise`, which does not contract), and no reassociation is allowed. So every output must
  be byte-identical to the plain Release build's. The check compares each measured tool's rows byte for byte, and
  runs the non-slow test suite under the LTO + PGO build.
- **Measured on:** this cloud machine with GCC and Clang. The owner's laptop builds with MSVC, so its LTO
  (`/GL`, `/LTCG`) and PGO (`/GENPROFILE`, `/USEPROFILE`) are built here but measured and checked there.
