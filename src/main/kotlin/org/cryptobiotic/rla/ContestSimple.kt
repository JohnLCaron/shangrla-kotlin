package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.*

data class ContestSimple (
    val id: String,
    val name: String = id,
    val choice_function: SocialChoiceFunction = SocialChoiceFunction.PLURALITY,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val n_winners: Int = 1,
    val candidates: List<String>,

    val betFn: BetFn? = null,
    val estimFn: EstimatorFn? = null,
    val testFn: TestFn? = null,
    val g: Double = 0.1,
    val winners: List<String>,
    val risk_limit: Double = 0.05,
    val share_to_win: Double = 0.5,
    val use_style: Boolean = true,

    // why are these mutable ??
    var assertions: MutableMap<String, AssertionSimple> = mutableMapOf(),// key = winr + " v " + losr
    var ncvrs: Int = 0,
    var ncards: Int = 0,
    var sample_size: Int? = null,
    var sample_threshold: Double? = null,
    val tally: MutableMap<String, Int> = mutableMapOf(), // candidate name -> vote count
) {
    companion object {
        fun fromVotes(audit: AuditSimple, votes: Map<String, Map<String, Int>>, cards: Map<String, Int>): List<ContestSimple> {
            /*
            Create a contest dict containing all contests in cvr_list.
            Every contest is single-winner plurality by default, audited by ballot comparison
            */
            val use_style = audit.use_styles
            val max_cards = audit.max_cards!!
            val contests = mutableListOf<ContestSimple>()

            for ((key, candidateMap) in votes) {
                /* ugh
                val options = np.array(list(votes[key].keys()), dtype = "str")
                val tallies = np.array(list(votes[key].values()))
                val reported_winner = options[np.argmax(tallies)]
                 */

                contests.add( ContestSimple(
                    id = key,
                    name= key,
                    ncards = if (use_style)  cards[key]!! else max_cards, // LOOK ??
                    choice_function = SocialChoiceFunction.PLURALITY,
                    n_winners = 1,
                    risk_limit = 0.05,
                    candidates = candidateMap.keys.toList(),
                    winners = emptyList<String>(), // listOf(reported_winner),
                    //assertion_file = None,
                    audit_type = AuditType.CARD_COMPARISON,
                    //test = NonnegMean.alpha_mart,
                    //estim = NonnegMean.optimal_comparison,
                    //bet = NonnegMean.fixed_bet,
                )
                )
            }
            return contests
        }
    }
}