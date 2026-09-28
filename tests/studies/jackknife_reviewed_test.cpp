// SPDX-License-Identifier: Apache-2.0
//
// The shared jackknife run (D-155, D-159): its committed summary holds raw band classes, and its
// reviewed labels are kept apart (studies/jackknife-bias-rho/reviewed.csv), so a review never needs
// a refit. Every out-of-band class must carry a label on its own side, and every label must still
// match an out-of-band class: a verdict that moves needs a new review, as for the recovery lists.
#include <map>
#include <set>
#include <string>
#include <utility>

#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"

VCAL_TEST(jackknife_reviewed_labels_match_the_summary) {
    const auto summary = vcal::test::read_golden_csv("../../studies/jackknife-bias-rho/summary.csv");
    const auto reviewed = vcal::test::read_golden_csv("../../studies/jackknife-bias-rho/reviewed.csv");
    VCAL_REQUIRE(summary.rows.size() == 81u);
    const std::pair<const char*, const char*> families[] = {{"shifted_rho", "shifted_class"},
                                                             {"bca_pd", "pd_bca_class"},
                                                             {"bca_rho", "rho_bca_class"},
                                                             {"q_wald_a", "qa_class"},
                                                             {"q_wald_b", "qb_class"}};
    const std::set<std::string> diagnoses = {"near_bound",          "conservative_box",      "small_rho_resolution",
                                             "refinement_resolution", "boundary_breakdown",  "resampling_variance",
                                             "se_overstated",       "se_tracks_error",       "jackknife_se_noisy"};
    std::map<std::pair<std::string, std::string>, std::string> label;  // (scenario, family) -> side
    for (const auto& r : reviewed.rows) {
        VCAL_CHECK(diagnoses.count(r[reviewed.column("diagnosis")]) == 1);
        const auto key = std::make_pair(r[reviewed.column("scenario")], r[reviewed.column("family")]);
        VCAL_CHECK(label.count(key) == 0);  // one label per verdict
        label[key] = r[reviewed.column("side")];
    }
    std::string unreviewed, stale;
    std::size_t out_of_band = 0;
    for (const auto& r : summary.rows) {
        const std::string id = r[summary.column("scenario")];
        for (const auto& [family, column] : families) {
            const std::string c = r[summary.column(column)];
            const auto it = label.find({id, family});
            if (c == "BELOW" || c == "ABOVE") {
                ++out_of_band;
                if (it == label.end() || it->second != c) unreviewed += " " + id + "/" + family + " " + c;
            } else if (it != label.end()) {
                stale += " " + id + "/" + family;
            }
        }
    }
    vcal::test::note(std::to_string(out_of_band) + " out-of-band verdicts, " + std::to_string(label.size()) +
                     " reviewed labels");
    if (!unreviewed.empty()) vcal::test::note("unreviewed:" + unreviewed);
    if (!stale.empty()) vcal::test::note("labelled but inside the band:" + stale);
    VCAL_CHECK(unreviewed.empty());
    VCAL_CHECK(stale.empty());
    VCAL_CHECK_EQ(out_of_band, label.size());
}
