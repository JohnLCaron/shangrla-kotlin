package org.cryptobiotic.shangrla

/*
    Generic class for cast-vote records.

    The CVR class DOES NOT IMPOSE VOTING RULES. For instance, the social choice
    function might consider a CVR that contains two votes in a contest to be an overvote.

    Rather, a CVR is supposed to reflect what the ballot shows, even if the ballot does not
    contain a valid vote in one or more contests.

    Class method get_vote_for returns the vote for a given candidate if the candidate is a
    key in the CVR, or False if the candidate is not in the CVR.

    This allows very flexible representation of votes, including ranked voting.

    For instance, in a plurality contest with four candidates, a vote for Alice (and only Alice)
    in a mayoral contest could be represented by any of the following:
            {"id": "A-001-01", "pool": False, "pool_group": "ABC", "phantom:: False"votes": {"mayor": {"Alice": True}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": "marked"}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": 5}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": 1, "Bob": 0, "Candy": 0, "Dan": ""}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": True, "Bob": False}}}
    A CVR that contains a vote for Alice for "mayor" and a vote for Bob for "DA" could be represented as
            {"id": "A-001-01", "votes": {"mayor": {"Alice": True}, "DA": {"Bob": True}}}

    NOTE: some methods distinguish between a CVR that contains a particular contest, but no valid
    vote in that contest, and a CVR that does not contain that contest at all. Thus, the following
    are not equivalent:
            {"id": "A-001-01", "votes": {"mayor": {}} }
            and
            {"id": "A-001-01", "votes": {} }

    Ranked votes also have simple representation, e.g., if the CVR is
            {"id": "A-001-01", "votes": {"mayor": {"Alice": 1, "Bob": 2, "Candy": 3, "Dan": ''}}}
    Then int(vote_for("Candy","mayor"))=3, Candy's rank in the "mayor" contest.

    CVRs can be flagged as `phantoms` to account for ballot cards not listed in the manifest using the boolean
    `phantom` attribute.

    CVRs can be assigned to a `tally_pool`, useful for the ONEAudit method or batch-level comparison audits
    using the `batch` attribute (batch-level comparison audits are not currently implemented)

    CVRs can be flagged for use in ONEAudit "pool" assorter means. When a CVR is flagged this way, the
    value of the assorter applied to the MVR is compared to the mean value of the assorter applied to the
    CVRs in the tally batch the CVR belongs to.

    CVRs can include sampling probabilities `p` and sample numbers `sample_num` (pseudo-random numbers
    to facilitate consistent sampling)

    CVRs can include a sequence number to facilitate ordering, sorting, and permuting

    Methods:
    --------

    get_vote_for:
         get_vote_for(candidate, contest_id) returns the value in the votes dict for the key `candidate`, or
         False if the candidate did not get a vote or the contest_id is not in the CVR
    has_contest: returns bool
         does the CVR have the contest?
    cvrs_to_json:
         represent CVR list as json
    from_dict: create a CVR from a dict
    from_dict_of_dicts:
         create dict of CVRs from a list of dicts
    from_raire:
         create CVRs from the RAIRE representation
 */
class CVR (
    val id: String,
    staringVotes: Map<String, Map<String, Int>>, // contest/vote dict
    var phantom: Boolean,
    var tally_pool: Any, // what tallying pool of cards does this CVR belong to (used by ONEAudit)? TODO String?
    var pool: Boolean,      // pool votes on this CVR within its tally_pool?
    val sample_num: Float,  // pseudorandom number used for consistent sampling
    var p: Float,           // sampling probability
    val sampled: Boolean,   // is this CVR in the sample?
) {
    val votes = mutableMapOf<String, MutableMap<String, Int>>() // Map(contestId, Map(selectionId, vote))


    fun get_vote_for(contest_id: String, candidate: String): Int? {
        val contest: Map<String, Int>? = votes[contest_id]
        return if (contest == null) null else contest[candidate]
    }

    fun has_contest(contest_id: String): Boolean = votes[contest_id] != null

    // add votes
    fun update_votes(votes: Map<String, Map<String, Int>>): Boolean {

        /*
        Update the votes for any contests the CVR already contains; add any contests and votes not already contained

        Parameters
        ----------
        votes: dict of dict of dicts
        key is a contest id; value is a dict of votes--keys and values

        Returns
        -------
        added: bool
        True if the contest was already present; else false
        TODO true if any contests were added.

        Side effects
        ------------
        updates the CVR to add the contest if it was not already present and to update the votes
        */
        //         added = False
        //        for c, v in votes.items():
        //            if self.has_contest(c):
        //                self.votes[c].update(v)
        //            else:
        //                self.votes[c] = v
        //                added = True
        //        return added

        var added = false
        for ((contest, contestVotes) in votes) {
            if (this.has_contest(contest)) {
                val myContestVotes = this.votes[contest]!!
                myContestVotes += contestVotes
            } else {
                this.votes[contest] = contestVotes.toMutableMap()
                added = true
            }
        }
        return added
    }

    fun has_one_vote(contest_id: String, candidates: List<String>): Boolean {
        /*
        Is there exactly one vote among the candidates in the contest `contest_id`?

        Parameters:
        -----------
        contest_id: identifier of contest
        candidates: list of identifiers of candidates

        Returns:
        ----------
        True if there is exactly one vote among those candidates in that contest, where a
        vote means that the value for that key casts as boolean True.
    */
        //         v = np.sum([0 if c not in self.votes[contest_id] else bool(self.votes[contest_id][c]) \
        //                    for c in candidates])
        //        return True if v==1 else False

        val contestVotes = this.votes[contest_id]
        if (contestVotes == null) return false

        val totalVotes = candidates.map{ contestVotes[it] ?: 0 }.sum()
        return (totalVotes == 1)
    }

    fun rcv_lfunc_wo(contest_id: String, winner: String, loser: String): Int {
        /*
            Check whether vote is a vote for the loser with respect to a 'winner only'
            assertion between the given 'winner' and 'loser'.Parameters:
            -----------
            contest_id: identifier for the contest
                    winner: string
            identifier for winning candidate
                    loser: string
            identifier for losing candidate
                    cvr: CVR object

            Returns:
            --------
            1 if the given vote is a vote for 'loser' and 0 otherwise
        */

        val rank_winner = this.get_vote_for(contest_id, winner)
        val rank_loser = this.get_vote_for(contest_id, loser)

        //         if not bool(rank_winner) and bool(rank_loser):
        if ((rank_winner == null) && (rank_loser != null)) return 1
        //         elif bool(rank_winner) and bool(rank_loser) and rank_loser < rank_winner:
        return if ((rank_winner != null) && (rank_loser != null) && rank_loser < rank_winner) 1 else 0
    }

    fun rcv_votefor_cand(contest_id: String, cand: String, remaining: List<String>): Int {
        /*
        Check whether 'vote' is a vote for the given candidate in the context
        where only candidates in 'remaining' remain standing.

        Parameters:
        -----------
        contest_id: string
        identifier of the contest used in the CVRs
        cand: string
        identifier for candidate
        remaining: list
        list of identifiers of candidates still standing

        vote: dict of dicts

        Returns:
        --------
        1 if the given vote for the contest counts as a vote for 'cand' and 0 otherwise.
        Essentially, if you reduce the ballot down to only those candidates in 'remaining',
        and 'cand' is the first preference, return 1; otherwise return 0.
        */
        if (!remaining.contains(cand)) return 0

        val rank_cand = this.get_vote_for(contest_id, cand)
        if (rank_cand == null) return 0

        for (altc in remaining) {
            if (altc == cand) continue
            val rank_altc = this.get_vote_for(contest_id, altc)
            if ((rank_altc != null) && rank_altc <= rank_cand) return 0
        }
        return 1
    }

    companion object {
        /*
        fun from_dict(cvr_dict: List<Any>): List<CVR> {
            /*
            Construct a list of CVR objects from a list of dicts containing cvr data

            Parameters:
            -----------
            cvr_dict: a list of dicts, one per cvr

            Returns:
            ---------
            list of CVR objects
            */
            val cvr_list = mutableListOf<CVR>()
            for (c in cvr_dict) {
                val phantom = False if ("phantom") not in c . keys () else c["phantom"]
                val pool = None if "pool" not in c . keys () else c["pool"]
                val tally_pool = None if "tally_pool" not in c . keys () else c["tally_pool"]
                val sample_num = None if "sample_num" not in c . keys () else c["sample_num"]
                val p = None if "p" not in c . keys () else c["p"]
                val sampled = None if "sampled" not in c . keys () else c["sampled"]
                cvr_list.append(
                    CVR(
                        id = c["id"],
                        startingVotes = c["votes"],
                        phantom = phantom,
                        pool = pool,
                        tally_pool = tally_pool,
                        sample_num = sample_num,
                        p = p,
                        sampled = sampled
                    )
                )
            }
            return cvr_list
        }

        fun from_raire(raire: List<Any>, phantom: Boolean): Pair<List<Any>, Int> {
            /*
            Create a list of CVR objects from a list of cvrs in RAIRE format

            Parameters:
            -----------
            raire: list of comma-separated values
            source in RAIRE format. From the RAIRE documentation:
            The RAIRE format (for later processing) is a CSV file.
            First line: number of contests.
            Next, a line for each contest
            Contest,id,N,C1,C2,C3 ...
            id is the contest_id
            N is the number of candidates in that contest
            and C1, ... are the candidate id's relevant to that contest.
            Then a line for every ranking that appears on a ballot:
            Contest id,Ballot id,R1,R2,R3,...
            where the Ri's are the unique candidate ids.

            The CVR file is assumed to have been read using csv.reader(), so each row has been split.

            Returns:
            --------
            list of CVR objects corresponding to the RAIRE cvrs, merged number of CVRs read (before merging)
            */
            val skip = int(raire[0][0])
            cvr_list = []
            for c in raire[(skip+1):]:
            contest_id = c[0]
            id = c[1]
            votes =
            {}
            for j in range(2, len(c)):
            votes[str(c[j])] = j-1
            cvr_list.append(CVR.from_vote(votes, id=id, contest_id=contest_id, phantom=phantom))
            return CVR.merge_cvrs(cvr_list), len(raire)-skip
        } */

        fun merge_cvrs(cvr_list: List<CVR>): List<CVR> {
            /*
            Takes a list of CVRs that might contain duplicated ballot ids and merges the votes
            so that each identifier is listed only once, and votes from different records for that
            identifier are merged.
            The merge is in the order of the list: if a later mention of a ballot id has votes
            for the same contest as a previous mention, the votes in that contest are updated
            per the later mention.

            If any of the CVRs has phantom==False, sets phantom=False in the result.
            If only one of a multiple has `tally_pool`, set the tally_pool to that value; if they disagree, throw an error
            Set `pool=True` if any CVR with the ID has `pool=True`


            Parameters:
            -----------
            cvr_list: list of CVRs

            Returns:
            -----------
            list of merged CVRs
            */

            //         od = OrderedDict()
            //        for c in cvr_list:
            //            if c.id not in od:
            //                od[c.id] = c
            //            else:
            //                od[c.id].votes = {**od[c.id].votes, **c.votes}
            //                od[c.id].phantom = (c.phantom and od[c.id].phantom)
            //                od[c.id].pool = (c.pool or od[c.id].pool)
            //                if (
            //                    (od[c.id].tally_pool is None and c.tally_pool is None)
            //                    or (od[c.id].tally_pool is not None and c.tally_pool is None)
            //                    or (od[c.id].tally_pool == c.tally_pool)
            //                   ):
            //                    pass
            //                elif od[c.id].tally_pool is None and c.tally_pool is not None:
            //                    od[c.id].tally_pool = c.tally_pool
            //                else:
            //                    throw ValueError(f'two CVRs with the same ID have different tally_pools: \n{str(od)=}\n{str(c)=}')
            //                od[c.id].pool = od[c.id] or c.pool
            //        return [v for v in od.values()]

        val od = sortedMapOf<String, CVR>()
        for (cvr: CVR in cvr_list) {
            if (!od.contains(cvr.id))
                od[cvr.id] = cvr
            else {
                val useCvr = od[cvr.id]!!
                // TODO useCvr.votes = { **od[c.id].votes, **c.votes }
                useCvr.phantom = (cvr.phantom and useCvr.phantom)
                useCvr.pool = (cvr.pool or useCvr.pool)
                if (
                    (useCvr.tally_pool == null && cvr.tally_pool == null) ||
                    (useCvr.tally_pool != null && cvr.tally_pool == null) ||
                    (useCvr.tally_pool == cvr.tally_pool)
                ) {
                    ;
                } else if ((useCvr.tally_pool == null) && (cvr.tally_pool != null)) {
                    useCvr.tally_pool = cvr.tally_pool
                } else {
                    throw Exception("two CVRs with the same ID have different tally_pools: \n{str(od)=}\n{str(c)=}")
                }
                // od[c.id].pool = od[c.id] or c.pool // TODO error in srla?
                useCvr.pool = useCvr.pool || cvr.pool
            }
        }
            return od.values.toList()
        }
    }

}