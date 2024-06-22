package org.cryptobiotic.rla

import org.cryptobiotic.rla.CvrSimple.Companion.makeSimple
import java.security.SecureRandom

// always one plurality contest, always two candidates "A" and "B"
data class CvrSimple(
    val id: String,
    val votes : Map<String, Map<String, Int>>, // contest : candidate : vote
    val phantom : Boolean = false,
    val sampled: Boolean = false,   // what does this mean ??
    //  tally_pool  # what tallying pool of cards does this CVR belong to (used by ONEAudit)?
    //  pool  # pool votes on this CVR within its tally_pool?
) {
    var p: Double = if (sampled) 1.0 else 0.0
    val sample_num = random.nextDouble()

    fun votedForA(): Boolean {
        return votes.values.first()["A"]!! == 1
    }

    fun get_vote_for(contestId: String, candidate: String) : Int {
        return votes.values.first()[candidate]!!
    }

    fun has_contest(contest_id: String): Boolean = votes[contest_id] != null

    companion object {
        private val random = SecureRandom.getInstanceStrong()
        private var nextId = 0

        // always one plurality contest, always two candidates "A" and "B", randomly choose which one gets a vote
        fun makeSimple(): CvrSimple {
            val which = random.nextBoolean()
            val avote = if (which) 1 else 0
            val bvote = if (which) 0 else 1
            return CvrSimple(
                id = "ballot${nextId++}",
                votes = mapOf("contest" to mapOf("A" to avote, "B" to bvote))
            )
        }

        /**
         * Find unique styles and count them.
         * A style is a unique set of contest ids.
         */
        fun tabulate_styles(cvr_list: List<CvrSimple>): Map<Set<String>, Int> {
            val style_counts = mutableMapOf<Set<String>, Int>()
            for (cvr in cvr_list) {
                val style: Set<String> = cvr.votes.keys.toSet()
                val accum = style_counts.getOrPut(style) { 0 }
                style_counts[style] = accum + 1
            }
            return style_counts

            //val style_counts = defaultdict(int)
            //for (cvr in cvr_list) {
            //     style_counts[frozenset(cvr.votes.keys)] += 1
            //}
            //return style_counts
        }

        /**
         * Tabulate total votes for each candidate in each contest in cvr_list.
         * For plurality, supermajority, and approval. Not useful for ranked-choice voting.
         * return contestId -> candidate -> votes
         */
        fun tabulate_votes(cvr_list: List<CvrSimple>): Map<String, Map<String, Int>> {
            val r = mutableMapOf<String, MutableMap<String, Int>>()
            for (cvr in cvr_list) {
                for ((con, conVotes) in cvr.votes) {
                    val accumVotes = r.getOrPut(con) { mutableMapOf() }
                    for ((cand, vote) in conVotes) {
                        val accum = accumVotes.getOrPut(cand) { 0 }
                        accumVotes[cand] = accum + vote
                        // d[con][cand] += CVR.as_vote(c.get_vote_for(con, cand)) // TODO normalized to 0, 1 ?
                    }
                }
            }
            return r
        }

        // Number of cards in each contest, return contestId -> ncards
        fun tabulate_cards_contests(cvr_list: List<CvrSimple>): Map<String, Int> {
            val d = mutableMapOf<String, Int>()
            for (cvr in cvr_list) {
                for (con in cvr.votes.keys) {
                    val accum = d.getOrPut(con) { 0 }
                    d[con] = accum + 1
                }
            }
            return d
        }
    }
}

fun makeSimpleCvrs(n: Int): List<CvrSimple> {
    val result = mutableListOf<CvrSimple>()
    repeat(n) { result.add(makeSimple()) }
    return result
}

