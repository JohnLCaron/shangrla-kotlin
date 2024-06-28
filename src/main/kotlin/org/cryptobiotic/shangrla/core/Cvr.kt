package org.cryptobiotic.shangrla.core

import kotlin.random.Random

/*
    Generic class for cast-vote records.

    The CVR class DOES NOT IMPOSE VOTING RULES. For instance, the social choice
    function might consider a CVR that contains two votes in a contest to be an overvote.

    Rather, a CVR is supposed to reflect what the ballot shows, even if the ballot does not
    contain a valid vote in one or more contests.

    Class method get_vote_for returns the vote for a given candidate if the candidate is a
    key in the CVR, or false if the candidate is not in the CVR.

    This allows very flexible representation of votes, including ranked voting.

    For instance, in a plurality contest with four candidates, a vote for Alice (and only Alice)
    in a mayoral contest could be represented by any of the following:
            {"id": "A-001-01", "pool": false, "pool_group": "ABC", "phantom:: false"votes": {"mayor": {"Alice": true}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": "marked"}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": 5}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": 1, "Bob": 0, "Candy": 0, "Dan": ""}}}
            {"id": "A-001-01", "votes": {"mayor": {"Alice": true, "Bob": false}}}
    A CVR that contains a vote for Alice for "mayor" and a vote for Bob for "DA" could be represented as
            {"id": "A-001-01", "votes": {"mayor": {"Alice": true}, "DA": {"Bob": true}}}

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
         false if the candidate did not get a vote or the contest_id is not in the CVR
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
class Cvr (
    val id: String,  // aka ballot id?
    var phantom: Boolean,
    val votes : MutableMap<String, MutableMap<String, Int>>, // contest : candidate : vote
    var tally_pool: String? = null,    // what tallying pool of cards does this CVR belong to (used by ONEAudit)? TODO String?
    var pool: Boolean = false,      // pool votes on this CVR within its tally_pool?
    var sample_num: Int? = null,  // pseudorandom number used for consistent sampling
    var p: Double? = null,           // sampling probability
    var sampled: Boolean = false,   // is this CVR in the sample?
) {

    // TODO not distinguishing "not voted" vs "voted no on this candidate"
    fun get_vote_for(contest_id: String, candidate: String): Int {
        val contest: Map<String, Int>? = votes[contest_id]
        return if (contest == null) 0 else contest[candidate] ?: 0

        //         return (
        //            False
        //            if (contest_id not in self.votes or candidate not in self.votes[contest_id])
        //            else self.votes[contest_id][candidate]
        //        )
    }

    fun has_contest(contest_id: String): Boolean = votes[contest_id] != null

    fun update_votes(cvr: Cvr): Boolean {
        /*
        Update the votes for any contests the CVR already contains; add any contests and votes not already contained

        Parameters
        ----------
        votes: key is a contest id; value is a dict of votes--keys and values

        Returns
        -------
        added: True if the contest was already present; else false
        added: True if any contests were added

        Side effects
        ------------
        updates the CVR to add the contest if it was not already present and to update the votes
        */
        var added = false
        for ((c, v) in cvr.votes) {
           if (this.has_contest(c)) {
               val wtf: MutableMap<String, Int> = this.votes[c]!!
               this.votes[c] = (wtf + v).toMutableMap()
            } else {
                this.votes[c] = v
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
        true if there is exactly one vote among those candidates in that contest, where a
        vote means that the value for that key casts as boolean true.
    */
        //         v = np.sum([0 if c not in self.votes[contest_id] else bool(self.votes[contest_id][c]) \
        //                    for c in candidates])
        //        return true if v==1 else false

        // [0 if c not in self.votes[contest_id] else bool(self.votes[contest_id][c]) for c in candidates]
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
        //            return 1
        //        elif bool(rank_winner) and bool(rank_loser) and rank_loser < rank_winner:
        //            return 1
        //        else:
        //            return 0

        if (!python_bool(rank_winner) && python_bool(rank_loser)) return 1
        return if (python_bool(rank_winner) && python_bool(rank_loser) && (rank_loser!! < rank_winner!!)) 1 else 0
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
        if (!python_bool(rank_cand)) return 0

        for (altc in remaining) {
            if (altc == cand) continue
            val rank_altc = this.get_vote_for(contest_id, altc)
            if (python_bool(rank_altc) && rank_altc!! <= rank_cand!!) return 0
        }
        return 1

        //         if not cand in remaining:
        //            return 0
        //
        //        if not bool(rank_cand := self.get_vote_for(contest_id, cand)):
        //            return 0
        //        else:
        //            for altc in remaining:
        //                if altc == cand:
        //                    continue
        //                rank_altc = self.get_vote_for(contest_id, altc)
        //                if bool(rank_altc) and rank_altc <= rank_cand:
        //                    return 0
        //            return 1
    }

    companion object {
        fun as_vote(vote: Int?): Int {
            return if (vote == null || vote == 0) 0 else vote
        }

        /*

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

        fun merge_cvrs(cvr_list: List<Cvr>): List<Cvr> {
            /*
            Takes a list of CVRs that might contain duplicated ballot ids and merges the votes
            so that each identifier is listed only once, and votes from different records for that
            identifier are merged.
            The merge is in the order of the list: if a later mention of a ballot id has votes
            for the same contest as a previous mention, the votes in that contest are updated
            per the later mention.

            If any of the CVRs has phantom==false, sets phantom=false in the result.
            If only one of a multiple has `tally_pool`, set the tally_pool to that value; if they disagree, throw an error
            Set `pool=true` if any CVR with the ID has `pool=true`


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

            val od = sortedMapOf<String, Cvr>()
            for (cvr: Cvr in cvr_list) {
                if (!od.contains(cvr.id))
                    od[cvr.id] = cvr
                else {
                    val useCvr = od[cvr.id]!!
                    // TODO useCvr.votes = { **od[c.id].votes, **c.votes }
                    useCvr.phantom = (cvr.phantom and useCvr.phantom)
                    useCvr.pool = (cvr.pool || useCvr.pool)
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

        fun make_phantoms(max_cards: Int, cvr_list: List<Cvr>, contests: List<Contest>, use_style: Boolean=true, prefix: String = "")
        : Pair<List<Cvr>, Int> {
            /*
            Make phantom CVRs as needed for phantom cards; set contest parameters `cards` (if not set) and `cvrs`

            If use_style, phantoms are "per contest": each contest needs enough to account for the difference between
            the number of cards that might contain the contest and the number of CVRs that contain the contest. This can
            result in having more cards in all (manifest and phantoms) than max_cards, the maximum cast.

            If not use_style, phantoms are for the election as a whole: need enough to account for the difference
            between the number of cards in the manifest and the number of CVRs that contain the contest. Then, the total
            number of cards (manifest plus phantoms) equals max_cards.

            Parameters
            ----------
            max_cards : int; upper bound on the number of ballot cards
            cvr_list : list of CVR objects; the reported CVRs
            contests : dict of contests; information about each contest under audit
            prefix : String; prefix for ids for phantom CVRs to be added
            use_style : Boolean
            does the sampling use style information?

            Returns
            -------
            cvr_list : list of CVR objects; the reported CVRs and the phantom CVRs
            n_phantoms : int; number of phantom cards added

            Side effects
            ------------
            for each contest in `contests`, sets `cards` to max_cards if not specified by the user
            for each contest in `contests`, set `cvrs` to be the number of (real) CVRs that contain the contest
            */
            //         phantom_vrs = []
            //        n_cvrs = len(cvr_list)
            //        for c, v in contests.items():  # set contest parameters
            //            v['cvrs'] = np.sum([cvr.has_contest(c) for cvr in cvr_list if not cvr.is_phantom()])
            //            v['cards'] = max_cards if v['cards'] is None else v['cards'] // upper bound on cards cast in the contest
            //        if not use_style:              #  make (max_cards - len(cvr_list)) phantoms
            //            phantoms = max_cards - n_cvrs
            //            for i in range(phantoms):
            //                phantom_vrs.append(CVR(id=prefix+str(i+1), votes={}, phantom=true))
            //        else:                          # create phantom CVRs as needed for each contest
            //            for c, v in contests.items():
            //                phantoms_needed = v['cards']-v['cvrs']
            //                while len(phantom_vrs) < phantoms_needed:
            //                    phantom_vrs.append(CVR(id=prefix+str(len(phantom_vrs)+1), votes={}, phantom=true))
            //                for i in range(phantoms_needed):
            //                    phantom_vrs[i].votes[c]={}  # list contest c on the phantom CVR
            //            phantoms = len(phantom_vrs)
            //        cvr_list = cvr_list + phantom_vrs
            //        return cvr_list, phantoms

            val phantom_vrs = mutableListOf<Cvr>()
            var n_phantoms: Int
            val n_cvrs = cvr_list.size
            for (contest in contests) { // } set contest parameters
                contest.ncvrs = cvr_list.filter { it.has_contest(contest.id) && !it.phantom }.count()
                contest.ncards =  if (contest.ncards == null) max_cards else contest.ncards // upper bound on cards cast in the contest TODO
            }
            if (!use_style) {              //  make (max_cards-len(cvr_list)) phantoms
                n_phantoms = max_cards - n_cvrs
                repeat (n_phantoms) {
                    phantom_vrs.add(Cvr(id = "${prefix}${it + 1}",true, mutableMapOf())) // TODO wrong
                }
            } else {                         // create phantom CVRs as needed for each contest
                for (contest in contests) {
                    val phantoms_needed = contest.ncards - contest.ncvrs
                    while (phantom_vrs.size < phantoms_needed) {
                        phantom_vrs.add(Cvr(id = "${prefix}${phantom_vrs.size + 1}", true, mutableMapOf())) // TODO wrong
                    }
                    phantom_vrs.forEach {
                        it.votes[contest.id] = mutableMapOf()  // list contest c on the phantom CVR
                    }
                }
                n_phantoms = phantom_vrs.size
            }
            val result = cvr_list + phantom_vrs
            return Pair(result, n_phantoms)
        }

        /*
        fun check_tally_pools(cvr_list: List<CVR>, force: Boolean = true): Pair<List<CVR>, Boolean> {
            /*
            Checks whether every CVR in each tally_pool has the same value of `pool`.
            If `force==true`, set them all to true if any of them is true

            Parameters:
            -----------
            cvr_list: collection of CVRs to be merged
            force: set pool equal to the logical union of the pool values for each tally group

            Returns:
            -----------
            tuple: list of CVRs and a bool. The bool is true if `force==true` and any value of `pool` was changed.
            */
    //         for c in cvr_list:
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
    //                    oc[c.id].tally_pool = c.tally_pool
    //                else:
    //                    raise ValueError(f'two CVRs with the same ID have different tally_pools: \n{str(od)=}\n{str(c)=}')
    //                oc[c.id].pool = oc[c.id] or c.pool
    //        return [v for v in od.values()]
            val od = mutableMapOf<String, CVR>()
            for (c in cvr_list) {
                if (c.id not in od) { // TODO od not defined
                    od[c.id] = c
                } else {
                    od[c.id].votes = { **od[c.id].votes, **c.votes }
                    od[c.id].phantom = (c.phantom && od[c.id].phantom)
                    od[c.id].pool = (c.pool || od[c.id].pool)
                    if ((od[c.id].tally_pool == null && c.tally_pool == null)
                        || (od[c.id].tally_pool != null && c.tally_pool == null)
                        || (od[c.id].tally_pool == c.tally_pool)) {
                        ;
                    } else if (od [c.id].tally_pool == null && c.tally_pool != null) {
                        oc[c.id].tally_pool = c.tally_pool
                    } else {
                        throw Exception ("two CVRs with the same ID have different tally_pools: \n{str(od)=}\n{str(c)=}")
                    }
                    oc[c.id].pool = oc[c.id] or c.pool
                }
                return [v for v in od.values()]
        }

        fun pool_contests(cvrs: List<CVR>): Set<String> {
            /*
            create a set containing all contest ids in the list of CVRs

            Parameters
            ----------
            cvrs : list of CVR objects
            the set to collect contests from

            Returns
            -------
            a set containing the ID of every contest mentioned in the CVR list
            */
            val contests = mutableSetOf<String>()
            for (c in cvrs) {
                contests.addAll(c.votes.keys)
            }
            return contests
        }

        fun add_pool_contests(cvrs: List<Cvr>, tally_pools: dict): Boolean {
            /*
            for each tally_pool, ensure every CVR in that pool has every contest in that pool

            Parameters
            ----------
            cvrs : list of CVR objects
            the set to update with additional contests as needed

            tally_pools : dict
            keys are tally_pool ids, values are sets of contests every CVR in that pool should have

            Returns
            -------
            bool : true if any contest is added to any CVR
            */
            var added = false
            for (cvr in cvrs) {
                // added = c.update_votes({ con: {} for con in tally_pools[c.tally_pool] }) || added
            } // note: order of terms matters!
            return added
        }
         */

        // Assigns a pseudo-random sample number to each cvr in cvr_list
        fun assign_sample_nums(cvr_list: List<Cvr>) {
             for (cvr in cvr_list) {
                cvr.sample_num = Random.nextInt() // TODO SecureRandom ?
            }
        }

        /*

        fun prep_comparison_sample(mvr_sample: List<CVR>, cvr_sample: List<CVR>, sample_order: list)
        {
            /*
            prepare the MVRs and CVRs for comparison by putting them into the same (random) order
            in which the CVRs were selected

            conduct data integrity checks.

            Side-effects: sorts the mvr sample into the same order as the cvr sample

            Parameters
            ----------
            mvr_sample: the manually determined votes for the audited cards
            cvr_sample: the electronic vote record for the audited cards
            sample_order: dict to look up selection order of the cards. Keys are card identifiers. Values are dicts
                containing "selection_order" (which draw yielded the card) and "serial" (the card's original position)

            Side effects
            ------------
            sorts the mvr sample into the same order as the cvr sample
            */
            mvr_sample.sort(key = lambda x : sample_order [x.id]["selection_order"])
            cvr_sample.sort(key = lambda x : sample_order [x.id]["selection_order"])
            require(cvr_sample.size == mvr_sample.size) {
                "Number of cvrs (${cvr_sample.size}) and number of mvrs (${mvr_sample.size}) differ" }
            for (i in range(len(cvr_sample))) {
                require (mvr_sample [i].id == cvr_sample[i].id)
                { "Mismatch between id of cvr (${cvr_sample[i].id}) and mvr (${mvr_sample[i].id})" }
        }

        fun prep_polling_sample(mvr_sample: List<CVR>, sample_order: dict) {
            /*
            Put the mvr sample back into the random selection order.
            Only about the side effects.

            Parameters
            ----------
            mvr_sample: list of CVR objects
            sample_order: dict to look up selection order of the cards. Keys are card identifiers. Values are dicts
              containing "selection_order" (which draw yielded the card) and "serial" (the card's original position)

            Side Effects
            -------------
            mvr_sample is reordered into the random selection order
            */
            mvr_sample.sort(key = lambda x : sample_order [x.id]["selection_order"])
        }

        fun sort_cvr_sample_num(cvr_list: List<CVR>) {
            /*
            Sort cvr_list by sample_num
            Only about the side effects.

            Parameters
            ----------
            cvr_list: list
            list of CVR objects
            Side effects
            ------------
            cvr_list is sorted by sample_num
            */
            cvr_list.sort(key = lambda x : x . sample_num)
        }

        fun consistent_sampling(cvr_list: List<CVR>, contests: List<Contest>, sampled_cvr_indices: List<Int>): List<Int> {
            /*
            Sample CVR ids for contests to attain sample sizes in contests, a dict of Contest objects

            Assumes that phantoms have already been generated and sample_num has been assigned
            to every CVR, including phantoms

            Parameters
            ----------
            cvr_list: list of CVR objects
            contests: dict of Contest objects. Contest sample sizes must be set before calling this function.
            sampled_cvr_indices: indices of cvrs already in the sample

            Returns
            -------
            sampled_cvr_indices: indices of CVRs to sample (0-indexed)
            */
    //         current_sizes = defaultdict(int)
    //        contest_in_progress = lambda c: (current_sizes[c.id] < c.sample_size)
    //        if sampled_cvr_indices is None:
    //            sampled_cvr_indices = []
    //        else:
    //            for sam in sampled_cvr_indices:
    //                for c, con in contests.items():
    //                    current_sizes[c] += (1 if cvr_list[sam].has_contest(con.id) else 0)
    //        sorted_cvr_indices = [i for i, cv in sorted(enumerate(cvr_list), key = lambda x: x[1].sample_num)]
    //        inx = len(sampled_cvr_indices)
    //        while any([contest_in_progress(con) for c, con in contests.items()]):
    //            if any([(contest_in_progress(con) and cvr_list[sorted_cvr_indices[inx]].has_contest(con.id))
    //                     for c, con in contests.items()]):
    //                sampled_cvr_indices.append(sorted_cvr_indices[inx])
    //                for c, con in contests.items():
    //                    if cvr_list[sorted_cvr_indices[inx]].has_contest(con.id) and contest_in_progress(con):
    //                        con.sample_threshold = cvr_list[sorted_cvr_indices[inx]].sample_num
    //                        current_sizes[c] += 1
    //            inx += 1
    //        for i in range(len(cvr_list)):
    //            if i in sampled_cvr_indices:
    //                cvr_list[i].sampled = true
    //        return sampled_cvr_indices

            val current_sizes = defaultdict(int)
            val contest_in_progress = lambda c :(current_sizes[c.id] < c.sample_size)
            if (sampled_cvr_indices == null) {
                sampled_cvr_indices = []
            } else {
            for (sam in sampled_cvr_indices) {
                for (con in contests) {
                    current_sizes[c] += (1 if cvr_list[sam].has_contest(con.id) else 0)
                }
            }
            val sorted_cvr_indices = [i for i, cv in sorted(enumerate(cvr_list), key = lambda x: x[1].sample_num)]
            var inx = sampled_cvr_indices.size
            while (any([contest_in_progress(con) for c, con in contests.items()]) {
                if (any(
                    (contest_in_progress(con) and cvr_list[sorted_cvr_indices[inx]].has_contest(con.id)
                    for c, con in contests.items()])) {
                        sampled_cvr_indices.add(sorted_cvr_indices[inx])
                        for (con in contests) {
                            if (cvr_list[sorted_cvr_indices[inx]].has_contest(con.id) && contest_in_progress(con)) {
                                con.sample_threshold = cvr_list[sorted_cvr_indices[inx]].sample_num
                                current_sizes[c] += 1
                            }
                        }
                    }
                inx += 1
                }
            for (i in range(len(cvr_list))) {
                if (i in sampled_cvr_indices) {
                    cvr_list[i].sampled = true
                }
            }
            return sampled_cvr_indices
        }
         */

        /**
         * Find unique styles and count them.
        * A style is a unique set of contest ids.
        */
        fun tabulate_styles(cvr_list: List<Cvr>): Map<Set<String>, Int> {
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
        fun tabulate_votes(cvr_list: List<Cvr>): Map<String, Map<String, Int>>  {
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
        fun tabulate_cards_contests(cvr_list: List<Cvr>): Map<String, Int> {
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