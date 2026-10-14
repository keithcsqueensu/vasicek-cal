# References

Published papers that motivate studies or features of this library. They are cited as public motivation
only. A study here tests a question the paper raises, on synthetic or public data; it does not reproduce
or endorse any particular implementation, and nothing here describes any organisation as using a method.

Each entry gives:
- **the citation and its sources:** the full citation, the DOI and the page it was checked against;
- **a summary**, in this project's own words;
- **what it motivates:** the studies (`studies/README.md`) or features (`STATE.md`) it motivates.

**How the citations were checked (2026-10-09):**
- Every DOI was resolved at doi.org, and its publisher-deposited Crossref record read.
- The abstract and volume and issue were read on the publisher's article page.
- The SSRN pages refused automated access, so the two SSRN ids are confirmed through SSRN's own Crossref
  DOIs (10.2139/ssrn.*). Miu and Ozdemir (2009) was also checked against its full text (the SSRN version).
- **Verified:** confirmed on a primary page. Any detail confirmed only on a secondary page, or not at all,
  is marked.

---

## Miu and Ozdemir (2005)

Peter Miu and Bogie Ozdemir (2005). "Practical and theoretical challenges in validating Basel parameters:
key learnings from the experience of a Canadian bank". *Journal of Credit Risk* 1(4): 89–136.
DOI [10.21314/JCR.2005.026](https://doi.org/10.21314/JCR.2005.026). **Verified** (DOI, Crossref,
publisher's page).

**Summary.** The paper discusses how to validate the parameters of an internal-ratings framework and
proposes a back-testing scheme that includes a simulation-based calibration test. Its central argument is
that the rating philosophy matters for both quantities: whether grades are point-in-time (PIT, grade PDs
move with the cycle) or through-the-cycle (TTC, grade PDs are stable and default rates move instead)
changes the PD a grade should carry and the default correlation its history implies. The correlation
appropriate for a TTC system is therefore larger than for a PIT one. It also discusses adjustments to LGD.

**Motivates:** S-36 (fitted ρ as the rating philosophy moves from TTC to PIT).

## Miu and Ozdemir (2006)

Peter Miu and Bogie Ozdemir (2006). "Basel requirements of downturn loss given default: modeling and
estimating probability of default and loss given default correlations". *Journal of Credit Risk* 2(2):
43–68. DOI [10.21314/JCR.2006.037](https://doi.org/10.21314/JCR.2006.037). **Verified.**

**Summary.** The paper models default and loss severity as driven by a common systematic factor and
separates the components of the PD–LGD correlation. It shows that omitting the correlation can be offset by
a conservative, cyclical LGD, and, calibrating the model to historical default and loss data and running an
economic-capital simulation, estimates how much the mean LGD must be raised to compensate: about 35–41% for
a corporate portfolio and about 16% for a mid-market one.

**Motivates:** S-46 (tail-loss understatement with fixed rather than factor-linked LGD and EAD).

## Miu and Ozdemir (2008)

Peter Miu and Bogie Ozdemir (2008). "Estimating and validating long-run probability of default with
respect to Basel II requirements". *Journal of Risk Model Validation* 2(2): 3–41.
DOI [10.21314/JRMV.2008.021](https://doi.org/10.21314/JRMV.2008.021). SSRN 1026181 (DOI
[10.2139/ssrn.1026181](https://doi.org/10.2139/ssrn.1026181)). **Verified.**

**Summary.** The paper derives maximum-likelihood estimators of the long-run PD (LRPD) for a Vasicek-type
model with both cross-sectional and serial correlation in the systematic factor. It treats an infinitely
granular portfolio first and then finite portfolios. In simulations the MLEs outperform the alternatives,
including in small samples, and the simple average of annual default rates tends to understate the LRPD.
Confidence intervals built on the MLE make fewer testing errors, and the framework lets external default
data supplement a short internal history.

**Motivates:** S-35 (LRPD by simple averaging against MLE under serial correlation; the highest priority of
this addendum). Its finite-portfolio limit (the infinite-granularity approximation failing below about 100
borrowers, as cited in the 2009 paper) also bears on S-38.

## Miu and Ozdemir (2009)

Peter Miu and Bogie Ozdemir (2009). "Stress-testing probability of default and migration rate with respect
to Basel II requirements". *Journal of Risk Model Validation* 3(4): 3–38.
DOI [10.21314/JRMV.2009.048](https://doi.org/10.21314/JRMV.2009.048). SSRN 1365842 (DOI
[10.2139/ssrn.1365842](https://doi.org/10.2139/ssrn.1365842)). **Verified,** including against the full text
(SSRN version).

**Summary.** The paper gives a method for stressed PDs and rating-migration rates conditional on chosen
macroeconomic variables, consistent with the single-factor model and usable with an existing internal rating
scale. It works in three steps:

1. **Extract a systematic factor** from a long external series of default rates (aggregated speculative-grade
   rates). The factor is modelled as a function of macroeconomic explanatory variables plus a residual that
   follows an AR(1) process. This step assumes infinite granularity for computational efficiency, and notes
   that the approximation becomes inappropriate below about 100 borrowers.
2. **Estimate the internal grades' parameters** (the pairwise correlation and each grade's default point and
   migration thresholds) from the short internal history, taking the internal portfolio's systematic factor
   to be identical to the external one.
3. **Combine the two** to evaluate stressed PDs and migration rates under chosen scenarios.

Some details the studies here rely on:
- **Specification search:** 154 candidate specifications are ranked by McFadden's adjusted pseudo-R² and by
  the median RMSE over leave-one-out hold-outs, and the two ranks are summed.
- **Overlapping series:** overlapping series (quarterly-observed 12-month default rates) give a higher
  estimate of the AR(1) coefficient than non-overlapping ones.
- **Long-run PDs and transitions** are evaluated by simulation, bootstrapping the explanatory variables.
- **The worked example** is the oil and gas sector. In it, R² is fixed at 0.20 rather than estimated jointly
  with the other parameters (footnotes 25 and 27); the general method estimates it.
- **The Alberta remark:** the paper argues for calibrating to portfolios that resemble the one being stressed.
  Its illustration is that a model calibrated to all corporations in Canada is unlikely to suit
  corporations in the Province of Alberta, whose economy is more sensitive to commodity prices. This is a
  motivating remark, not a worked example.

**Motivates:**
- **Studies:** S-38 (granularity), S-39 (the plug-in factor), S-40 (fixed, estimated or grade-specific
  correlation), S-41 (external factor against internal factor), S-42 (specification-search rules), S-43
  (overlapping series), S-44 (LRPD by bootstrapping the macro variables) and S-47 (regional and sector
  factors).
- **Features:** the published-method track (the three steps as separately callable `native` steps).

## Miu and Ozdemir (2017)

Peter Miu and Bogie Ozdemir (2017). "Adapting the Basel II advanced internal-ratings-based models for
International Financial Reporting Standard 9". *Journal of Credit Risk* 13(2): 53–83.
DOI [10.21314/JCR.2017.224](https://doi.org/10.21314/JCR.2017.224). **Verified.**

**Summary.** The paper sets out the changes needed to turn regulatory-capital PD, LGD and EAD models into
forward-looking, lifetime expected-credit-loss (ECL) estimates, with worked examples for secured lending.
Among the practical issues it discusses, it proposes a convexity adjustment to the expected loss under a
base scenario, so that probabilities need not be assigned to multiple macroeconomic scenarios.

**Motivates:**
- **Study:** S-45 (lifetime PD: the convexity adjustment against few-scenario weighting against full factor
  integration).
- **Feature:** the lifetime PD term structure conditional on a scenario path.

## Miu and Ozdemir (2023)

Peter Miu and Bogie Ozdemir (2023). "A coherent economic framework to model correlations between PD, LGD and
EaD, and its applications in EaD modelling and IFRS-9". *Journal of Risk Management in Financial
Institutions* 16(1): 52–78. DOI [10.69554/QCVI3102](https://doi.org/10.69554/QCVI3102).
- **Verified:** the authors, title, journal, volume, issue and first page, on the publisher's DOI page and
  in Crossref.
- **Secondary only:** the page range 52–78 is from RePEc
  (<https://ideas.repec.org/a/aza/rmfi00/y2023v16i1p52-78.html>).

**Summary.** The paper treats exposure at default as stochastic and models the PD–LGD, PD–EaD and LGD–EaD
correlations together through common systematic factors. It shows how to estimate them from historical data
and to simulate the three jointly. It then measures how much tail risk in a credit loss distribution, and in
lifetime ECL, is understated when the joint correlation structure is ignored, including the extra drawdowns
in a downturn.

**Motivates:** S-46.

## Yang and Du (2015)

Bill Huajian Yang and Zunwei Du (2015). "Stress testing and modeling of rating migration under the Vasicek
model framework: empirical approaches and technical implementation". *Journal of Risk Model Validation*
9(2): 33–47. DOI [10.21314/JRMV.2015.137](https://doi.org/10.21314/JRMV.2015.137). **Verified.**
- **Details filled in:** the volume, issue and DOI were not given with the original reference, and the
  paper's subject is the stress testing of rating migration.
- **Not the same paper as** the same authors' 2016 paper in the same journal, DOI 10.21314/JRMV.2016.155.

**Summary.** The paper extends Miu and Ozdemir's (2009) transition-probability stress-testing framework by
giving each non-default rating its own asset correlation and its own factor model. It proposes two
Vasicek-type models per rating:
- a single latent factor, for the rating-level correlation;
- a multifactor model with a latent effect, for systematic downgrade risk.

Both are fitted as nonlinear mixed models, and the conditional transition probabilities are derived in
closed form. Modelling downgrades rather than defaults addresses the scarcity of defaults in high-quality
grades.

**Motivates:** S-40 (fixed, estimated or grade-specific correlation).

## Yang (2013)

Bill Huajian Yang (2013). "Estimating long-run probability of default, asset correlation and
portfolio-level probability of default using Vasicek models". *Journal of Risk Model Validation* 7(4):
3–19. DOI [10.21314/JRMV.2013.112](https://doi.org/10.21314/JRMV.2013.112). **Verified.**

**Summary.** The paper proposes Vasicek-type models of portfolio-level default rates in which the asset
correlation and the long-run PD have closed-form estimators, with room for macroeconomic or market
covariates. It shows that the usual asymptotic (infinitely granular) maximum-likelihood approach is
equivalent to least-squares regression of the probit-transformed default rates. The closed-form LRPD
shows explicitly how parameter uncertainty raises it, and a bootstrap is recommended for serial
correlation. The applications use published sector default rates and public delinquency-rate series.

**Motivates:** S-37, which the addendum recommends folding into documentation rather than running (D-188).
- **Without covariates, the equivalence is already this library's closed form:** the Vasicek-rate MLE is the
  sample mean and 1/T variance of the probit rates, which is least squares on an intercept
  (`docs/methodology/vasicek_rate_mle.md`, D-164).
- **With covariates,** it is the first step of the published-method track, documented there as a replication
  aid rather than run as a study.

## Chawla, Forest and Aguais (2017)

Gaurav Chawla, Lawrence R. Forest Jr. and Scott D. Aguais (2017). "Convexity and correlation effects in
expected credit loss calculations for IFRS9/CECL and stress testing". *Journal of Risk Management in
Financial Institutions* 10(1): 99–110. DOI [10.69554/VWXW7019](https://doi.org/10.69554/VWXW7019).
- **Verified:** the authors, title, journal, volume, issue and first page, on the publisher's DOI page and
  in Crossref.
- **Secondary only:** the page range 99–110 is from RePEc
  (<https://ideas.repec.org/a/aza/rmfi00/y2017v10i1p99-110.html>).

**Summary.** The paper shows that two effects skew the credit-loss distribution, and that ignoring them
understates expected credit losses, by 20% or more in the authors' examples:
- **Convexity:** the PD as a function of the systematic factors is convex. This effect grows with the
  curvature and with the dispersion of the factors.
- **Correlation:** PD, LGD and EAD co-move. This effect grows with how strongly they do.

It proposes credit-cycle indices that turn existing models into co-moving point-in-time estimates, for
unbiased scenario-conditional losses.

**Motivates:**
- **Study:** S-45.
- **Feature:** the lifetime PD term structure conditional on a scenario path.

## Pluto and Tasche (2005)

Katja Pluto and Dirk Tasche (2005). "Estimating probabilities of default for low default portfolios". Working
paper, arXiv:cond-mat/0411699v3 (4 April 2005; v1 28 November 2004),
<https://arxiv.org/abs/cond-mat/0411699>. Published as "Thinking positively", *Risk*, August 2005, pp. 72–78
(<https://www.risk.net/derivatives/credit-derivatives/1500251/thinking-positively>).
- **Verified (2026-10-13):** the working paper, read in full from arXiv: authors, title, versions and dates.
  The *Risk* article's title and authors, on Risk.net's page title (the article is paywalled).
- **Secondary only:** the *Risk* issue (August 2005) and pages 72–78, from reference listings (e.g. Scientific
  Research Publishing's reference record). The volume number is not confirmed.

**Summary.** For portfolios with no or very few defaults, the paper proposes the "most prudent estimation"
principle: estimate each grade's PD by an upper confidence bound, assuming only that the ordinal ranking of the
grades is correct. Grade k's bound pools grade k with every worse grade, which makes the bounds monotone across
grades.
- **Independent defaults:** the bound is the one-sided Clopper–Pearson bound.
- **Correlated defaults:** the bound uses the Basel one-factor model (ρ = 12% in its examples).
- **Extensions:** scaling the bounds to a portfolio's central tendency, and a multi-period cohort model whose
  yearly factors are correlated at ϑ^|s−t|.
- **Open questions it leaves:** the confidence level; the number of defaults beyond which ordinary estimation
  should take over; and what a non-monotone sequence of bounds says about the ranking. In a footnote it also
  leaves open what enlarged confidence level for the independent bound would allow for correlation.

`docs/methodology/pluto_tasche.md` §6 compares the engine with it and lists its errata.

**Motivates:** S-6 (coverage and conservatism of the bounds), and the open questions as S-49 to S-52.
