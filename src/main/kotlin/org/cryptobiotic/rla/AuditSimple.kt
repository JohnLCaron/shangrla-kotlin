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
    val error_rate_1: Double = 0.0,
    val error_rate_2: Double = 0.0,
    val reps: Int = 100,
    val use_styles: Boolean = true,
) {

    fun find_sample_size(
        contests: Map<String, ContestSimple>,
        cvrs: List<CvrSimple>?,
        mvr_sample: List<CvrSimple>,
        cvr_sample: List<CvrSimple>
    ): Int {
        /*
        Estimate sample size for each contest && overall to allow the audit to complete.
        Uses simulations. For speed, uses the numpy.r&&om Mersenne Twister instead of cryptor&&om, a higher-quality
        PRNG used to select the actual audit sample.

        Parameters
        ----------
        contests: dict of dicts; the contest data structure. outer keys are contest identifiers; inner keys are assertions
            TODO wrong, this is now just dict of Contests
        cvrs: list of CVR objects; the full set of CVRs
        mvr_sample: list of CVR objects; manually ascertained votes
        cvr_sample: list of CVR objects; CVRs corresponding to the cards that were manually inspected

        Returns
        -------
        new_size: int
        new sample size

        Side effects
        ------------
        sets c.sample_size for each Contest in contests
        if use_style, sets cvr.p for each CVR
        */
        //         if len(self.strata) > 1:
        //            raise NotImplementedError('Stratified audits are not currently implemented.')
        //        stratum = next(iter(self.strata.values())) # the only stratum
        //        if stratum.use_style and cvrs is None:
        //            raise ValueError("stratum.use_style==True but cvrs were not provided.")
        //        # unless style information is being used, the sample size is the same for every contest.
        //        old = (0 if stratum.use_style
        //               else len(mvr_sample))
        //        old_sizes = {c:old for c in contests.keys()}


        // unless style information is being used, the sample size is the same for every contest.
        val old = if (this.use_styles) 0 else mvr_sample.size
        val old_sizes: MutableMap<String, Int> =
            contests.keys.associate { it to old }.toMutableMap()  // old_sizes = {c:old for c in contests.keys()}

//        for c, con in contests.items():
//            if stratum.use_style:
//                old_sizes[c] = np.sum(np.array([cvr.sampled for cvr in cvrs if cvr.has_contest(c)]))
//            new_size = 0
//            for a, asn in con.assertions.items():
//                if not asn.proved:
//                    if mvr_sample is not None: # use MVRs to estimate the next sample size. Set `prefix=True` to use data
//                        data, u =  asn.mvrs_to_data(mvr_sample, cvr_sample)
//                        new_size = max(new_size, asn.find_sample_size(data=data, prefix=True,
//                                                                  reps=self.reps, quantile=self.quantile,
//                                                                  seed=self.sim_seed))
//                    else:
//                        data=None
//                        new_size = max(new_size, asn.find_sample_size(data=data, rate_1=self.error_rate_1,
//                                                                  rate_2=self.error_rate_2,
//                                                                  reps=self.reps, quantile=self.quantile,
//                                                                  seed=self.sim_seed))
//            con.sample_size = new_size

        for ((contestId, contest) in contests) {
            if (this.use_styles) {
                requireNotNull(cvrs)
                // old_sizes[c] = np.sum(np.array([cvr.sampled for cvr in cvrs if cvr.has_contest(c)]))
                // TODO summing booleans?? Maybe 0 or 1 ??
                old_sizes[contestId] = cvrs.filter { it.has_contest(contestId) }.map { if (it.sampled) 1 else 0 }.sum()
            }
            var new_size = 0
            for ((_, asn) in contest.assertions) {
                if (!asn.proved) {
                    if (mvr_sample != null) { // use MVRs to estimate the next sample size. Set `prefix=True` to use data
                        val (data, _) = asn.mvrs_to_data(mvr_sample, cvr_sample)
                        new_size = max(
                            new_size,
                            asn.find_sample_size(
                                data = data, prefix = true,
                                reps = this.reps, quantile = this.quantile,
                                seed = this.sim_seed
                            )
                        )
                    }
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
                contest.sample_size = new_size
            }
        }

//        if stratum.use_style:
//            for cvr in cvrs:
//                if cvr.sampled:
//                    cvr.p=1
//                else:
//                    cvr.p=0
//                    for c, con in contests.items():
//                        if cvr.has_contest(c) and not cvr.sampled:
//                            cvr.p = max(con.sample_size/(con.cards - old_sizes[c]), cvr.p)
//            total_size = math.ceil(np.sum([x.p for x in cvrs if not x.phantom]))
//        else:
//            total_size = np.max(np.array([con.sample_size for con in contests.values()]))
//        return total_size

        // TODO why are we setting p here ??
        var total_size: Int
        if (this.use_styles) {
            requireNotNull(cvrs)
            for (cvr in cvrs) {
                if (cvr.sampled) {
                    cvr.p = 1.0
                } else {
                    cvr.p = 0.0
                    for ((c, con) in contests) {
                        if (cvr.has_contest(c) && !cvr.sampled) {
                            val p1 = con.sample_size!! / (con.ncards - old_sizes[c]!!)
                            cvr.p = max(p1.toDouble(), cvr.p!!) // TODO nullability
                        }
                    }
                }
            }
            // total_size = ceil(np.sum([x.p for x in cvrs if !x.phantom))
            val summ: Double = cvrs.filter { !it.phantom }.map { it.p!! }.sum()
            total_size = ceil(summ).toInt()
        } else {
            // total_size = np.max(np.array([con.sample_size for con in contests.values()]))
            total_size = contests.values.map { it.sample_size!! }.max()
        }
        return total_size
    }
}

//        'seed':           12345678901234567890,
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
