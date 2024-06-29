package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.AuditType
import org.cryptobiotic.shangrla.core.SocialChoiceFunction

data class Contest(
    val id: String,
    val n_winners: Int = 1,
    val candidates: List<String>,
    val choice_function: SocialChoiceFunction = SocialChoiceFunction.PLURALITY,
    val audit_type: AuditType = AuditType.CARD_COMPARISON,
    val risk_limit: Double = 0.05,
    val share_to_win: Double = .5, // for supermajority

    val ncards: Int,
    val winners: List<String>,
    val tally: Map<String, Int>,
) {
    val assertions = mutableMapOf<String, Assertion>() // does it really want to be a map ??
    var sample_size: Int = 0
    var sample_threshold: Int = 0

    init {
        require (n_winners > 0 && n_winners < candidates.size)
        require (share_to_win >= .5) // could allow this
    }

    fun addAssertion(a: Assertion) = assertions.put(a.name(), a)

    companion object {

        fun fromVotes(
            audit: Audit,
            votes: Map<String, Map<String, Int>>,  // contestId -> candidate -> votes
            cards: Map<String, Int>, // contestId -> ncards
        ): List<Contest> {

            val contests = mutableListOf<Contest>()

            for ((contestId, candidateMap) in votes) {
                val winner = candidateMap.maxBy { it.value }.key

                contests.add(
                    Contest(
                        id = contestId,
     //  ncards = if (use_style)  cards[key]!! else max_cards, // LOOK ??
                        ncards = cards[contestId]!!,
                        n_winners = 1,
                        candidates = candidateMap.keys.toList(),
                        winners = listOf(winner),
                        tally = candidateMap,
                    )
                )
            }

            return contests
        }
    }
}