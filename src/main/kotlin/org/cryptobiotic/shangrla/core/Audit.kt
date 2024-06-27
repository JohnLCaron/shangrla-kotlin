package org.cryptobiotic.shangrla.core

import kotlin.math.ceil
import kotlin.math.max

class Stratum(
    id: String? = null,
    val max_cards: Int? = null,
    val use_style: Boolean,
    replacement: Boolean,
    audit_type: AuditType,
    //test:  callable=None,
    //estimator:  callable=None,
    //bet: callable=None,
    //test_kwargs: dict=None):
)

enum class AuditType { POLLING, CARD_COMPARISON, ONEAUDIT }

class Audit(
    seed: Object? = null,
    val sim_seed: Int,
    cvr_file: String? = null,
    manifest_file: String? = null,
    sample_file: String? = null,
    mvr_file: String? = null,
    log_file: String? = null,
    val quantile: Double,
    val error_rate_1: Double = 0.001,
    val error_rate_2: Double = 0.00,
    val reps: Int = 100,
    max_cards: Int? = null,
    val strata: Map<String, Stratum>,
) {
    init {
        require(error_rate_1 >= 0) { "expected rate of 1-vote errors must be nonnegative" }
        require(error_rate_2 >= 0) { "expected rate of 2-vote errors must be nonnegative" }
    }

    fun find_sample_size(
        contests: Map<String, Contest>,
        cvrs: List<Cvr>? = null,
        mvr_sample: List<Cvr> = emptyList(),
        cvr_sample: List<Cvr> = emptyList(),
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


        // Currently, only unstratified audits are supported
        if (this.strata.size > 1) {
            throw NotImplementedError("Stratified audits are not currently implemented.")
        }
        val stratum = strata.values.first() // the only stratum
        if (stratum.use_style && (cvrs == null)) {
            throw Exception("stratum.use_style==True but cvrs were not provided.") // TODO
        }
        // unless style information is being used, the sample size is the same for every contest.
        val old = if (stratum.use_style) 0 else mvr_sample.size
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
            if (stratum.use_style) {
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
        var total_size: Int
        if (stratum.use_style) {
            requireNotNull(cvrs)
            for (cvr in cvrs) {
                if (cvr.sampled) {
                    cvr.p = 1.0
                } else {
                    cvr.p = 0.0
                    for ((c, con) in contests) {
                        if (cvr.has_contest(c) && !cvr.sampled) {
                            val p1 = con.sample_size!!.toDouble() / (con.ncards - old_sizes[c]!!)
                            cvr.p = max(p1, cvr.p!!) // TODO nullability
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



