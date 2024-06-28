package org.cryptobiotic.start


data class Contest(
    val id: String,
    val n_winners: Int = 1,
    var candidates: List<String>,
    val ncards: Int,
    val winners: List<String>,
    val risk_limit: Double = 0.05,
) {
    val assertions = mutableMapOf<String, Assertion>() // does it really want to be a map ??
    var sample_size: Int = 0
    var sample_threshold: Int = 0

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
                        ncards = cards[contestId]!!,
                        n_winners = 1,
                        candidates = candidateMap.keys.toList(),
                        winners = listOf(winner),
                    )
                )
            }

            return contests
        }
    }
}