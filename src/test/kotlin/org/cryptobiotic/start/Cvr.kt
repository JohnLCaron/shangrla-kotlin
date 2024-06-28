package org.cryptobiotic.start

data class Cvr(
    val id: String,
    val votes: Map<String, Map<String, Int>>, // contest : candidate : vote
    val phantom: Boolean = false
) {
    var sample_num: Int = 0
    var sampled: Boolean = false
    var p: Double = Double.POSITIVE_INFINITY

    fun has_contest(contest_id: String): Boolean = votes[contest_id] != null

    fun get_vote_for(contestId: String, candidate: String): Int {
        return if (votes[contestId] == null) 0 else votes[contestId]!![candidate] ?: 0
    }

    companion object {
        fun tabulateVotes(cvrs: List<Cvr>): Map<String, Map<String, Int>> {
            val r = mutableMapOf<String, MutableMap<String, Int>>()
            for (cvr in cvrs) {
                for ((con, conVotes) in cvr.votes) {
                    val accumVotes = r.getOrPut(con) { mutableMapOf() }
                    for ((cand, vote) in conVotes) {
                        val accum = accumVotes.getOrPut(cand) { 0 }
                        accumVotes[cand] = accum + vote
                    }
                }
            }
            return r
        }

        // Number of cards in each contest, return contestId -> ncards
        fun cardsPerContest(cvrs: List<Cvr>): Map<String, Int> {
            val d = mutableMapOf<String, Int>()
            for (cvr in cvrs) {
                for (con in cvr.votes.keys) {
                    val accum = d.getOrPut(con) { 0 }
                    d[con] = accum + 1
                }
            }
            return d
        }
    }
}