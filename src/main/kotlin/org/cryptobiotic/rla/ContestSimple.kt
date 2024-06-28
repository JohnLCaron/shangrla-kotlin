package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.*

data class ContestSimple (
    val id: String,
    val name: String = id,
    val choice_function: SocialChoiceFunction = SocialChoiceFunction.PLURALITY,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val n_winners: Int = 1,
    var candidates: List<String>,
    val ncards: Int,
    val winners: List<String>,
    val risk_limit: Double = 0.05,

    var assertions : Map<String, AssertionSimple> = emptyMap(),
    var sample_threshold: Int = 0,
    var sample_size: Int = 0,
    var testOverride: NonnegMean? = null,
    var estimOverride: EstimFnType? = null,
    ) {

    fun getTestFn(): NonnegMean {
        val result = testOverride ?: NonnegMean(u = 1.0, N = this.ncards, t = .5, estimFnType =  this.estimOverride)
        return result
    }

    companion object {
        fun fromVotes(audit: AuditSimple,
                      votes: Map<String, Map<String, Int>>,  // contestId -> candidate -> votes
                      cards: Map<String, Int>, // contestId -> ncards
                      winnerOverride: List<String>? = null,
        ): List<ContestSimple> {

            val use_style = audit.use_styles
            val max_cards = audit.max_cards
            val contests = mutableListOf<ContestSimple>()

            var count = 0
            for ((key, candidateMap) in votes) {
                val winner = if (winnerOverride == null) candidateMap.maxBy { it.value }.key else winnerOverride[count++]

                contests.add( ContestSimple(
                    id = key,
                    name= key,
                    ncards = if (use_style)  cards[key]!! else max_cards, // LOOK ??
                    choice_function = SocialChoiceFunction.PLURALITY,
                    n_winners = 1,
                    candidates = candidateMap.keys.toList(),
                    winners = listOf(winner),
                    audit_type = AuditType.CARD_COMPARISON,
                )
                )
            }
            return contests
        }
    }
}