package org.cryptobiotic.simple

import org.cryptobiotic.simple.CvrSimple.Companion.makeSimple
import java.security.SecureRandom

// always one plurality contest, always two candidates "A" and "B"
data class CvrSimple(
    val id: String,
    val votes : Map<String, Map<String, Int>>, // contest : candidate : vote
    val phantom : Boolean = false,
    var sampled: Boolean = false,
    var p: Double = if (sampled) 1.0 else 0.0
    //  tally_pool  # what tallying pool of cards does this CVR belong to (used by ONEAudit)?
    //  pool  # pool votes on this CVR within its tally_pool?
) {
    var sample_num = random.nextInt()

    fun votedForA(): Boolean {
        return votes.values.first()["A"]!! == 1
    }

    //     def get_vote_for(self, contest_id: str, candidate: str):
    //        return (
    //            False
    //            if (contest_id not in self.votes or candidate not in self.votes[contest_id])
    //            else self.votes[contest_id][candidate]
    //        )
    fun get_vote_for(contestId: String, candidate: String): Int {
        return votes.values.first()[candidate] ?: 0
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
         * Find unique styles and count how many cvrs use them.
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


        fun assign_sample_nums(cvr_list: List<CvrSimple>) {
            //         prng = SHA256(audit.seed)  // presumable returns 32 random bytes
            //     def assign_sample_nums(cls, cvr_list: list["CVR"], prng: "np.RandomState") -> bool:
            //         for cvr in cvr_list:
            //            cvr.sample_num = int_from_hash(prng.nextRandom()) // TODO wtf does int_from_hash do ?

            cvr_list.forEach {
                it.sample_num = random.nextInt() // TODO modify the cvr!
            } // so sample_num is a random integer, uses all 32 bits.
        }

        fun consistent_sampling(
            cvr_list: List<CvrSimple>,
            contests: List<ContestSimple>,
            sampled_cvr_indices_starting: List<Int> = emptyList()
        ): List<Int> {
            /*
        Sample CVR ids for contests to attain sample sizes in contests, a dict of Contest objects

        Assumes that phantoms have already been generated and sample_num has been assigned
        to every CVR, including phantoms

        Parameters
        ----------
        cvr_list: list of CVR objects
        contests: Contest objects. Contest sample sizes must be set before calling this function.
        sampled_cvr_indices: indices of cvrs already in the sample

        Returns sampled_cvr_indices: indices of CVRs to sample (0-indexed)
        */

            //        current_sizes = defaultdict(int)
            //        contest_in_progress = lambda c: (current_sizes[c.id] < c.sample_size)
            //        if sampled_cvr_indices is None:
            //            sampled_cvr_indices = []
            //        else:
            //            for sam in sampled_cvr_indices:
            //                for c, con in contests.items():
            //                    current_sizes[c] += 1 if cvr_list[sam].has_contest(con.id) else 0

            val current_sizes = mutableMapOf<String, Int>() // contest: Accum vote over sampled cvrs
            for (sam in sampled_cvr_indices_starting) {
                for (con in contests) {
                    if (cvr_list[sam].has_contest(con.id)) {
                        val accum = current_sizes.getOrPut(con.id) { 0 }
                        current_sizes[con.id] = accum + 1
                    }
                }
            }
            val sampled_cvr_indices = mutableListOf<Int>()
            sampled_cvr_indices.addAll(sampled_cvr_indices_starting)

            //        sorted_cvr_indices = [
            //            i for i, cv in sorted(enumerate(cvr_list), key=lambda x: x[1].sample_num)
            //        ]
            // enumerate returns Pair(index, element)

            // map of cvr index to random "sample_num", sorted by sample num? why sorted? why store sample_num in cvr?
            val sorted_cvr_indices: Map<Int, Int> =
                cvr_list.mapIndexed { idx, it -> Pair(idx, it.sample_num) }.sortedBy { it.second }.toMap()

            // val contests_in_progress = lambda c: (current_sizes[c.id] < c.sample_size)
            //        inx = len(sampled_cvr_indices)
            //        while any([contest_in_progress(con) for c, con in contests.items()]):
            //            if any(
            //                [
            //                    (
            //                        contest_in_progress(con)
            //                        and cvr_list[sorted_cvr_indices[inx]].has_contest(con.id)
            //                    )
            //                    for c, con in contests.items()
            //                ]
            //            ):
            //                sampled_cvr_indices.append(sorted_cvr_indices[inx])
            //                for c, con in contests.items():
            //                    if cvr_list[sorted_cvr_indices[inx]].has_contest(
            //                        con.id
            //                    ) and contest_in_progress(con):
            //                        con.sample_threshold = cvr_list[
            //                            sorted_cvr_indices[inx]
            //                        ].sample_num
            //                        current_sizes[c] += 1
            //            inx += 1


            val contest_in_progress = { con: ContestSimple -> current_sizes[con.id]!! < con.sample_size }

            var inx = sampled_cvr_indices.size
            while (contests.firstOrNull { contest_in_progress(it) } != null) {
                val index = sorted_cvr_indices[inx]!! // why non null?
                if (contests.firstOrNull { contest -> contest_in_progress(contest) && cvr_list[index].has_contest(contest.id) } != null) {
                    sampled_cvr_indices.add(index)
                    contests.forEach { con ->
                        if (cvr_list[index].has_contest(con.id) && contest_in_progress(con)) {
                            con.sample_threshold = cvr_list[index].sample_num // TODO modify the contest!
                            val accum = current_sizes.getOrPut(con.id) { 0 }
                            current_sizes[con.id]  = accum + 1
                        }
                    }
                    inx += 1
                }
            }

            //        for i in range(len(cvr_list)):
            //            if i in sampled_cvr_indices:
            //                cvr_list[i].sampled = True
            //        return sampled_cvr_indices
            for (i in 0..cvr_list.size) {
                if (i in sampled_cvr_indices) {
                    cvr_list[i].sampled = true // TODO modify the cvr!
                }
            }
            return sampled_cvr_indices
        }

    } // companion object
}


fun makeSimpleCvrs(n: Int): List<CvrSimple> {
    val result = mutableListOf<CvrSimple>()
    repeat(n) { result.add(makeSimple()) }
    return result
}

