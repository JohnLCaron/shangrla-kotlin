package org.cryptobiotic.rla

import kotlin.math.ceil
import kotlin.math.max
import org.cryptobiotic.shangrla.core.*

fun makeTestAudit(max_cards: Int) = AuditSimple(
    seed = 1234567890,
    sim_seed = 314159265,
    cvr_file = "test/data/rla/cvrFile",
    manifest_file = "test/data/rla/OC_full_manifest.xlsx",
    quantile = 0.8,
    error_rate_1 = 0.0,
    error_rate_2 = 0.0,
    reps = 100,
    max_cards = max_cards,
    use_styles = true,
)

// audit = Audit.from_dict({
//         'seed':           12345678901234567890,
//         'sim_seed':       314159265,
//         'cvr_file':       '/Users/amanda/Downloads/oc_cvrs.zip',
//         #'cvr_file':       '/Users/Jake/Desktop/oc_cvrs.zip',
//         #'manifest_file':  'data/OC_mock_manifest_detailed.xlsx',
//         'manifest_file': 'data/OC_full_manifest.xlsx',
//         #'manifest_file': 'tests/data/Hart_manifest.xlsx',
//         'sample_file':    '',
//         'mvr_file':       '',
//         'log_file':       'data/OC_example_log.json',
//         'quantile':       0.8,
//         'error_rate_1': 0,
//         'error_rate_2': 0,
//         #'error_rate_1':   0.001,
//         #'error_rate_2':   0.0001,
//         'reps':           100,
//         'strata':         {'stratum_1': {'max_cards':   3094308,
//                                          'use_style':   False,
//                                          'replacement': False,
//                                          'audit_type':  Audit.AUDIT_TYPE.BALLOT_COMPARISON, // TODO ??
//                                          'test':        NonnegMean.alpha_mart,
//                                          'estimator':   NonnegMean.optimal_comparison,
//                                          'test_kwargs': {}
//                                         }
//                           }
//        })
// audit = Audit.from_dict({
//         'seed':           12345678901234567890,
//         'sim_seed':       314159265,
//         'cvr_file':       './data/SFDA2019/SFDA2019_PrelimReport12VBMJustDASheets.raire',
//         'manifest_file':  './data/SFDA2019/N19 ballot manifest with WH location for RLA Upload VBM 11-14.xlsx',
//         'sample_file':    './data/sample.csv',
//         'mvr_file':       './data/mvr.json',
//         'log_file':       './data/log.json',
//         'quantile':       0.8,
//         'error_rate_1':   0.001,
//         'error_rate_2':   0.0,
//         'reps':           100,
//         'strata':         {'stratum_1': {'max_cards':   293555,
//                                          'use_style':   True,
//                                          'replacement': False,
//                                          'audit_type':  Audit.AUDIT_TYPE.CARD_COMPARISON,
//                                          'test':        NonnegMean.alpha_mart,
//                                          'estimator':   NonnegMean.optimal_comparison,
//                                          'test_kwargs': {}
//                                         }
//                           }
//        })

data class AuditSimple(
    val cvr_file: String,
    val manifest_file: String,
    val max_cards: Int,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val seed: Int = 1234567890,
    val sim_seed: Int = 314159265,
    val quantile: Double = 0.8,
    val error_rate_1: Double = 0.001,
    val error_rate_2: Double = 0.00,
    val reps: Int = 100,
    val use_styles: Boolean = true,
) {

    fun find_sample_size(
        contests: List<ContestSimple>,
        cvrs: List<CvrSimple>? = null,
        mvr_sample: List<CvrSimple> = emptyList(),
        cvr_sample: List<CvrSimple> = emptyList(),
    ): Int {
        // unless style information is being used, the sample size is the same for every contest.
        val old = if (this.use_styles) 0 else mvr_sample.size
        val old_sizes: MutableMap<String, Int> =
            contests.associate { it.id to old }.toMutableMap()  // old_sizes = {c:old for c in contests.keys()}

        for (contest in contests) {
            val contestId = contest.id
            if (this.use_styles) {
                requireNotNull(cvrs)
                // old_sizes[c] = np.sum(np.array([cvr.sampled for cvr in cvrs if cvr.has_contest(c)]))
                // TODO summing booleans?? Maybe 0 or 1 ??
                old_sizes[contestId] = cvrs.filter { it.has_contest(contestId) }.map { if (it.sampled) 1 else 0 }.sum()
            }
            var new_size = 0
            for ((_, asn) in contest.assertions) {
                if (!asn.proved) {
                    if (mvr_sample.isNotEmpty()) { // use MVRs to estimate the next sample size. Set `prefix=True` to use data
                        val (data, _) = asn.mvrs_to_data(mvr_sample, cvr_sample)
                        new_size = max(
                            new_size,
                            asn.find_sample_size(
                                data = data, prefix = true,
                                reps = this.reps, quantile = this.quantile,
                                seed = this.sim_seed
                            )
                        )
                    } else {
                        new_size = max(
                            new_size, asn.find_sample_size(
                                data = null, rate1 = this.error_rate_1,
                                rate2 = this.error_rate_2,
                                reps = this.reps, quantile = this.quantile,
                                seed = this.sim_seed
                            )
                        )
                    }
                }
                contest.sample_size = new_size
            }
        }

        // TODO why are we setting p here ??
        var total_size: Int
        if (this.use_styles) {
            requireNotNull(cvrs)
            for (cvr in cvrs) {
                if (cvr.sampled) {
                    cvr.p = 1.0
                } else {
                    cvr.p = 0.0
                    for (con in contests) {
                        if (cvr.has_contest(con.id) && !cvr.sampled) {
                            val p1 = con.sample_size!! / (con.ncards - old_sizes[con.id]!!)
                            cvr.p = max(p1.toDouble(), cvr.p) // TODO nullability
                        }
                    }
                }
            }
            // total_size = ceil(np.sum([x.p for x in cvrs if !x.phantom))
            val summ: Double = cvrs.filter { !it.phantom }.map { it.p!! }.sum()
            total_size = ceil(summ).toInt()
        } else {
            // total_size = np.max(np.array([con.sample_size for con in contests.values()]))
            total_size = contests.map { it.sample_size!! }.max()
        }
        return total_size
    }

    fun check_audit_parameters(contests: List<ContestSimple>) {
        require(this.error_rate_1 >= 0) { "expected rate of 1-vote errors must be nonnegative" }
        require(this.error_rate_2 >= 0) { "expected rate of 2-vote errors must be nonnegative" }
        for (con in contests) {
            require(con.risk_limit > 0) { "risk limit {con.risk_limit} negative in contest {c}" }
            require(con.risk_limit <= .5) { "risk limit {con.risk_limit} exceeds 1/2 in contest {c}" }
            // require(con.choice_function in Contest.SOCIAL_CHOICE_FUNCTION.SOCIAL_CHOICE_FUNCTIONS) { "unsupported choice function {con.choice_function} in contest {c}" }
            require(con.n_winners <= con.candidates.size) { "more winners than candidates in contest {c}" }
            require(con.winners.size == con.n_winners) { "number of reported winners does not equal n_winners in contest {c}" }
            for (w in con.winners) {
                require(con.candidates.contains(w)) { "reported winner {w} is not a candidate in contest {c}" }
            }
            if (con.choice_function in listOf(SocialChoiceFunction.IRV, SocialChoiceFunction.SUPERMAJORITY)) {
                require(con.n_winners == 1) { "{con.choice_function} can have only 1 winner in contest {c}" }
                //if (con.choice_function == SocialChoiceFunction.IRV) {
                //    require(con.assertion_file,) { "IRV contest {c} requires an assertion file" }
                //}
            }
        }
    }
}
