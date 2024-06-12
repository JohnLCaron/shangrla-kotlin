package org.cryptobiotic.shangrla.core

import kotlin.math.max

enum class SocialChoiceFunction{ APPROVAL, PLURALITY, SUPERMAJORITY, IRV }

enum class Candidates { ALL, ALL_OTHERS, WRITE_IN, NO_CANDIDATE }

class Contest(
    val id: String,
    val name: String,
    val risk_limit: Double = 0.05,
    val cards: Int = 0,
    val choice_function: SocialChoiceFunction,
    val n_winners: Int = 1,
    val share_to_win: Double,
    val candidates: List<String>,
    val reported_winners: List<String>, // TODO was called winners?
    val assertion_file: String,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val testFn: TestFn?,
    val g: Double = 0.1,
    val estimFn: EstimatorFn?,
    val bet: BetFn?,
    val use_style: Boolean = true,
    var assertions: Map<String, Assertion>, // TODO why is this a map?
    var tally: MutableMap<String, Int>,
    var sample_size: Int,
    val sample_threshold: Double,
) {

    fun find_margins_from_tally() {
        /*
        Use the `Contest.tally` attribute to set the margins of the contest's assorters.

        Appropriate only for the social choice functions
        Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY,
        Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY,
        Contest.SOCIAL_CHOICE_FUNCTION.APPROVAL

        Side effects
        ------------
        sets Assertion.margin for all Assertions in the Contest
        */
        for ((_, assn) in this.assertions) {
            assn.find_margin_from_tally()
        }
    }

    fun find_sample_size(audit: Audit, mvr_sample: List<CVR>, cvr_sample: List<CVR>): Int {
        /*
        Estimate the sample size required to confirm the contest at its risk limit.

        This function can be used with or without data, for Audit.AUDIT_TYPE.POLLING,
        Audit.AUDIT_TYPE.CARD_COMPARISON, and Audit.AUDIT_TYPE.ONEAUDIT audits.

        The simulations in this implementation are inefficient because the randomization happens separately
        for every assorter, rather than in parallel.

        Parameters
        ----------
        cvr_sample: list of CVRs; data (or simulated data) to base the sample size estimates on
        mvr_sample: list of MVRs (CVR objects); manually read votes to base the sample size estimates on, if data are available.

        Returns
        -------
        estimated sample size

        Side effects
        ------------
        sets this.sample_size to the estimated sample size
        */

        this.sample_size = 0
        for (a in this.assertions.values) {
            val data = null
            if (mvr_sample != null) {  // process the MVRs/CVRs to get data appropriate to each assertion
                val (data, u) = a.mvrs_to_data(mvr_sample, cvr_sample)
            }
            this.sample_size = max(
                this.sample_size,
                a.find_sample_size(
                    data = data, rate1 = audit.error_rate_1, rate2 = audit.error_rate_2,
                    reps = audit.reps, quantile = audit.quantile, seed = audit.sim_seed
                )
            )
        }
        return this.sample_size
    }

    companion object {

        fun tally(contests: List<Contest>, cvr_list: List<CVR>) {
            /*
            Tally the votes in the contests from a collection of CVRs.
            Only tallies plurality, multi-winner plurality, supermajority, and approval contests

            Parameters
            ----------
            con_dict: dict; dict of Contest objects to find tallies for
            cvr_list: list[CVR]; list of CVRs containing the votes to tally

            Returns
            -------

            Side Effects
            ------------
            Sets the `tally` dict for the contests in con_list, if their social choice function is appropriate
            */
            //         tallies = {}
            //        cons = []
            //        for id, c in con_dict.items():
            //            if c.choice_function in [Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY,
            //                                     Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY,
            //                                     Contest.SOCIAL_CHOICE_FUNCTION.APPROVAL]:
            //                cons.append(c)
            //                c.tally = defaultdict(int)
            //            else:
            //                warnings.warn(f'contest {c.id} ({c.name}) has social choice function ' +
            //                              f'{c.choice_function}: not tabulated')
            //        for cvr in cvr_list:
            //            for c in cons:
            //                if cvr.has_contest(c.id):
            //                    for candidate, vote in cvr.votes[c.id].items():
            //                        if candidate:
            //                            c.tally[candidate] += int(bool(vote))
            val wantContests = mutableListOf<Contest>()
            for (contest in contests) {
                if (contest.choice_function in listOf(
                        SocialChoiceFunction.PLURALITY,
                        SocialChoiceFunction.SUPERMAJORITY,
                        SocialChoiceFunction.APPROVAL
                    )
                ) {
                    wantContests.add(contest)
                    contest.tally = mutableMapOf<String, Int>()
                } else {
                    println("contest ${contest.id} (${contest.name}) has social choice function ${contest.choice_function}: not tabulated")
                }
            }
            for (cvr in cvr_list) {
                for (contest in wantContests) {
                    if (cvr.has_contest(contest.id)) {
                        for ((candidate, vote) in cvr.votes[contest.id]!!) {
                            if (candidate != null && (vote > 0)) {
                                // contest.tally[candidate] += int(bool(vote)) // TODO int(bool(vote)) ??
                                val accum = contest.tally[candidate]!!
                                contest.tally[candidate] = accum + 1
                            }
                        }
                    }
                }
            }
        }

    }
}