package org.cryptobiotic.shangrla.core

import kotlin.math.max

enum class SocialChoiceFunction{ APPROVAL, PLURALITY, SUPERMAJORITY, IRV }

enum class Candidates { ALL, ALL_OTHERS, WRITE_IN, NO_CANDIDATE }

class Contest(
    val id: String,
    val name: String = id,
    val candidates: List<String>,
    val choice_function: SocialChoiceFunction,
    var assertions: MutableMap<String, Assertion> = mutableMapOf(),// key = winr + " v " + losr
    val assertion_file: String? = null,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val betFn: BetFn? = null,
    val estimFn: EstimatorFn? = null,
    val testFn: TestFn? = null,
    val g: Double = 0.1,
    val n_winners: Int = 1,
    val reported_winners: List<String> = emptyList(), // TODO was called winners?
    val risk_limit: Double = 0.05,
    var ncvrs: Int = 0,
    var ncards: Int = 0,
    var sample_size: Int? = null,
    var sample_threshold: Double? = null,
    val share_to_win: Double = 0.5,
    val tally: MutableMap<String, Int> = mutableMapOf(), // candidate name -> vote count
    val use_style: Boolean = true,
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

    fun find_sample_size(audit: Audit, mvr_sample: List<Cvr>, cvr_sample: List<Cvr>): Int {
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
                this.sample_size ?: 0,
                a.find_sample_size(
                    data = data, rate1 = audit.error_rate_1, rate2 = audit.error_rate_2,
                    reps = audit.reps, quantile = audit.quantile, seed = audit.sim_seed
                )
            )
        }
        return this.sample_size!!
    }

    companion object {

        /** Tally the votes in the contests from a collection of CVRs.
         * Only tallies plurality, multi-winner plurality, supermajority, and approval contests
         * Sets the `tally` dict for the contests in con_list, if their social choice function is appropriate
         */
        fun tally(contests: List<Contest>, cvr_list: List<Cvr>) {
            //        tallies = {}
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

            val wantContests = mutableListOf<Contest>()
            for (contest in contests) {
                if (contest.choice_function in listOf(
                        SocialChoiceFunction.PLURALITY,
                        SocialChoiceFunction.SUPERMAJORITY,
                        SocialChoiceFunction.APPROVAL
                    )
                ) {
                    wantContests.add(contest)
                } else {
                    println("contest ${contest.id} (${contest.name}) has social choice function ${contest.choice_function}: not tabulated")
                }
            }
            //        for cvr in cvr_list:
            //            for c in cons:
            //                if cvr.has_contest(c.id):
            //                    for candidate, vote in cvr.votes[c.id].items():
            //                        if candidate:
            //                            c.tally[candidate] += int(bool(vote))
            for (cvr in cvr_list) {
                for (contest in wantContests) {
                    if (cvr.has_contest(contest.id)) {
                        val candMap = cvr.votes[contest.id]!!
                        for ((candidate, vote) in candMap) {
                            val accum = contest.tally.getOrPut(candidate) { 0 }
                            val voteN = if (python_bool(vote)) 1 else 0
                            contest.tally[candidate] = accum + voteN
                        }
                    }
                }
            }
        }

        // val cvrs : List<Cvr> = emptyList()
        //val votes: Map<String, Map<String, Int>> = Cvr.tabulate_votes(cvrs)
        //val styles: Map<Set<String>, Int> = Cvr.tabulate_styles(cvrs)
        //val cards: Map<String, Int> = Cvr.tabulate_cards_contests(cvrs)

        //         if len(audit.strata) > 1:
        //            raise NotImplementedError("stratified audits not implemented")
        //        stratum = next(iter(audit.strata.values()))
        //        use_style = stratum.use_style
        //        max_cards = stratum.max_cards
        //        contest_dict = {}
        //        for key in votes:
        //            contest_name = str(key)
        //            cards_with_contest = cards[key]
        //            options = np.array(list(votes[key].keys()), dtype="str")
        //            tallies = np.array(list(votes[key].values()))
        //
        //            reported_winner = options[np.argmax(tallies)]
        //
        //            contest_dict[contest_name] = {
        //                "name": contest_name,
        //                "cards": cards_with_contest if use_style else max_cards,
        //                "choice_function": Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY,
        //                "n_winners": 1,
        //                "risk_limit": 0.05,
        //                "candidates": list(options),
        //                "winner": [reported_winner],
        //                "assertion_file": None,
        //                "audit_type": Audit.AUDIT_TYPE.CARD_COMPARISON,
        //                "test": NonnegMean.alpha_mart,
        //                "estim": NonnegMean.optimal_comparison,
        //                "bet": NonnegMean.fixed_bet,
        //            }
        //        contests = Contest.from_dict_of_dicts(contest_dict)
        //        return contests
        fun fromVotes(audit: Audit, votes: Map<String, Map<String, Int>>, cards: Map<String, Int>): List<Contest> {
            /*
            Create a contest dict containing all contests in cvr_list.
            Every contest is single-winner plurality by default, audited by ballot comparison
            */
            if (audit.strata.size > 1) {
                throw NotImplementedError ("stratified audits not implemented")
            }
            val stratum = audit.strata.values.first()
            val use_style = stratum.use_style
            val max_cards = stratum.max_cards!!
            val contests = mutableListOf<Contest>()

            for ((key, candidateMap) in votes) {
                /* ugh
                val options = np.array(list(votes[key].keys()), dtype = "str")
                val tallies = np.array(list(votes[key].values()))
                val reported_winner = options[np.argmax(tallies)]
                 */

                contests.add( Contest(
                    id = key,
                        name= key,
                        ncards = if (use_style)  cards[key]!! else max_cards, // LOOK ??
                        choice_function = SocialChoiceFunction.PLURALITY,
                        n_winners = 1,
                        risk_limit = 0.05,
                        candidates = candidateMap.keys.toList(),
                        reported_winners = emptyList<String>(), // listOf(reported_winner),
                        //assertion_file = None,
                        audit_type = AuditType.CARD_COMPARISON,
                        //test = NonnegMean.alpha_mart,
                        //estim = NonnegMean.optimal_comparison,
                        //bet = NonnegMean.fixed_bet,
                ))
            }
            return contests
        }

    }
}