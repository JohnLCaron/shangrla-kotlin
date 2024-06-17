package org.cryptobiotic.shangrla.raire

import org.cryptobiotic.shangrla.core.Cvr
import kotlin.math.max


data class Contest(
    val name: String,
    val candidates: List<String>,
    val winner: String,
    val tot_ballots: Int,
    val outcome: IntArray
)

/*
fun load_contests_from_txt(path: String) {
    /*
        Format:
        First line is a comma separated list of candidate identifiers, either
        ending with the winner expressed as ",winner,winner identifier" or
        addditionally specifying the full outcome with ",order,sequence"

        Second line is party identifiers for each candidate
        Third line is a separator (eg. -----)

        Each subsequent line has the form:
        (Comma separated list of candidate identifiers) : Number of ballots

        Each line defines a ballot signature, a preference ordering over
        candidates, and the number of ballots that have been cast with that
        signature.

        Use default contest name of "1".
    */

    contests = []
    cvrs = {}

    var tot_auditable_ballots = 0

    with(open(path, "r") as data) {
        val lines = data.readlines()

        val toks = [line.strip() for line in lines[0].strip().split(',')]
        val windx = toks.index("winner")
        val winner = toks[windx + 1]
        val cands = toks[:windx]

        val order = []
        if ("order" in toks) {
            order = toks[windx + 2:]
        }

        val bcntr = 0

        for (l in range(3, len(lines))) {
            toks = [line.strip() for line in lines[l].strip().split(':')]
        }

        val num = int(toks[1])

        val prefs = [p.strip() for p in toks[0][1:-1].split(',')]

        tot_auditable_ballots += num

        if (prefs == []) {
            continue
        }

        for (i in range(num)) {
            val ballot = {}
            for (c in cands) {
                if (c in prefs) {
                    idx = prefs.index(c)
                    ballot[c] = idx
                }

                cvrs[bcntr] = { 1 : ballot }
                bcntr += 1
            }
        }

        return [Contest(
            1, cands, winner, this.tot_ballots,
            order = order
        )], cvrs
    }
}

// Data file in .raire format.
fun load_contests_from_raire(path: String) {
    val contests = []

//  A map between ballot id and the relevant CVR.
    val cvrs = {}
    with(open(path, "r") as data) {
        val lines = data.readlines()

//  Total number of contests described in data file
        val ncontests = int(lines[0])

//  Map between contest id and number of ballots involving that contest
        val num_ballots = {}

//  Map between contest id and the candidates & winner of that contest.
        val contest_info = {}

        for (i in range(ncontests)) {
            toks = [line.strip() for line in lines[1 + i].strip().split(',')]
        }

//  Get contest id and number of candidates in that contest
        val cid = toks[1]
        val ncands = int(toks[2])

//  Get list of candidate identifiers
        val cands = []

        for (j in range(ncands)) {
            cands.append(toks[3 + j])
        }

        val windx = toks.index("winner")
        val winner = toks[windx + 1]

        var informal = 0
        var inf_index = null
        if ("informal" in toks) {
            inf_index = toks.index("informal")
            informal = int(toks[inf_index + 1])
        }

        val order = []
        if ("order" in toks) {
            order = toks[windx + 2:inf_index] if inf_index != null  else
            toks[windx + 2:]
        }

        contest_info[cid] = (cands, winner, order)
        num_ballots[cid] = informal

        for (l in range(ncontests + 1, len(lines))) {
            toks = [line.strip() for line in lines[l].strip().split(',')]


            val cid = toks[0]
            val bid = toks[1]
            val prefs = toks[2:]

            val ballot = {}
            for (c in contest_info[cid][0]) {
                if (c in prefs) {
                    idx = prefs.index(c)
                    ballot[c] = idx
                }
            }
            num_ballots[cid] += 1

            if (!bid in cvrs) {
                cvrs[bid] = { cid: ballot }
            } else {
                cvrs[bid][cid] = ballot
            }
        }

        for (cid, (cands, winner, order) in contest_info.items()) {
        con = Contest(cid, cands, winner, num_ballots[cid], order = order)
        contests.append(con)
    }
    }

    return Pair(contests, cvrs)
}

fun index_of(cand, list_of_cand) {
    /*
    Returns position of given candidate 'cand' in the list of candidates
    'list_of_cand'. Returns -1 if 'cand' is not in the given list.

    Input:
    cand : string       - Identifier of candidate we are looking for
    list_of_cand : list - List of candidate identifiers

    Output:
    Index (starting at 0) of 'cand' in 'list_of_cand', and -1 if
    'cand' is not in 'list_of_cand'.
    */
    for (i in range(len(list_of_cand))) {
        if (list_of_cand[i] == cand) {
            return i
        }
    }

    return -1
}

 */

fun ranking(cand: String, ballot: Map<String, Int>): Int {
    /*
    Input:
    cand           -   identifier for candidate
    ballot         -   mapping between candidate name and their position in the ranking for a relevant contest on a given ballot.

    Output:
    Returns the position of candidate 'cand' in the ranking of the
    given ballot 'ballot'. Returns -1 if 'cand' is not preferenced on the ballot.
    */
    return ballot[cand] ?: -1
}


fun vote_for_cand(cand: String, eliminated: List<String>, ballot: Map<String, Int>): Int {
    /*
    Input:
    cand - identifier for candidate
    eliminated : list-identifiers of eliminated candidates
    ballot - mapping between candidate name and their
            position in the ranking for a relevant contest on a given ballot .
    Output:
    Returns 1 if the given 'ballot' is a vote for the given candidate 'cand'
    in the context where candidates in 'eliminated' have been eliminated.
    Otherwise, return 0 as the 'ballot' is not a vote for 'cand'.
    */

//  If 'cand' is not in the set of candidates assumed still standing,
//  'cand' does not get this vote.
    if (cand in eliminated) return 0

//  If 'cand' does not appear on the ballot, they do not get this vote.
    val c_idx = ranking(cand, ballot)
    if (c_idx == -1) return 0

    for ((alt_c, a_idx) in ballot) {
        if (alt_c == cand) continue
        if (alt_c in eliminated) continue

        if (a_idx < c_idx) return 0
    }

    return 1
}


open class RaireAssertion(val contest_name: String, val winner: String, val loser: String): Comparable<RaireAssertion> {
    /*
        Initializes a RAIRE assertion involving a comparison between
        the tallies of a candidate labelled 'winner' and a candidate
        labelled 'loser'. This assertion 'asserts' that the tally of
        the winner is larger than the tally of the loser in some context.

        Each assertion will have an estimated 'difficulty' related to
        the anticipated number of ballot checks required to audit it.

        Each assertion will have a margin defined as the difference in
        tallies ascribed to 'winner' and 'loser'
        */

    var votes_for_winner = 0
    var votes_for_loser = 0

    var margin = -1
    var difficulty = Double.POSITIVE_INFINITY

    val rules_out = mutableSetOf<RaireAssertion>()

    open fun is_vote_for_winner(cvr: Cvr): Int {
        /*
        Input:
            cvr - cast vote record

        Output:
            Returns 1 if the given cvr represents a vote for the assertions
            winner, and 0 otherwise.
        */
        return 0 // TODO
    }

    open fun is_vote_for_loser(cvr: Cvr): Int {
        /*
        Input:
            cvr - cast vote record

        Output:
            Returns 1 if the given cvr represents a vote for the assertions
            loser, and 0 otherwise.
        */
        return 0 // TODO
    }

    open fun subsumes(other: RaireAssertion): Boolean {
        /*
        Returns true if this assertion 'subsumes' the input assertion 'other'.
        An assertion 'A' subsumes assertion 'B' if the alternate outcomes
        ruled out by 'B' is a subset of those ruled out by 'A'. If we include
        'A' in an audit, we don't need to include 'B'.

        Input:
        other : RaireAssertion   - Assertion 'B'

        Output:
        Returns true if this assertion subsumes assertion 'other'.
        */
        return true // TODO
    }

    open fun same_as(other: RaireAssertion): Boolean {
        /*
        Returns true  if this assertion is equal to 'other' (i.e., they
        are the same assertion), and false  otherwise.
        */
        return true // TODO
    }

    /*  Assertions are ordered in terms of how many alternate outcomes that they are able to rule out.
    fun __lt__(other: RaireAssertion) {
        val self_rule_out = if (!this.rules_out) -1 else min([len(ro) for ro in this.rules_out])

        other_rule_out = -1 if not other . rules_out else
        min([len(ro) for ro in other.rules_out])

        return self_rule_out < other_rule_out
    }

    fun __gt__(other: RaireAssertion) {
        self_rule_out = -1 if not self . rules_out else
        min([len(ro) for ro in this.rules_out])

        other_rule_out = -1 if not other . rules_out else
        min([len(ro) for ro in other.rules_out])

        return self_rule_out > other_rule_out
    }

     */

    open fun to_str(): String {
        return "TODO"
    }

    override fun compareTo(other: RaireAssertion): Int {
        //         val self_rule_out = if (!this.rules_out) -1 else min([len(ro) for ro in this.rules_out])
        val self_rule_out = if (this.rules_out.isEmpty()) -1 else 0 // this.rules_out.map { it.size }.min() TODO
        val other_rule_out = if (other.rules_out.isEmpty()) -1 else 0 // other.rules_out.map { it.size }.min()

        return self_rule_out - other_rule_out
    }
}

class NEBAssertion(contest_name: String, winner: String, loser: String) : RaireAssertion(contest_name, winner, loser) {

    /*
    A Not-Eliminated-Before (NEB) assertion between a candidate 'winner' and
    a candidate 'loser' compares the minimum possible tally 'winner' could
    have (their first preference tally) with the maximum possible tally
    candidate 'loser' could have while 'winner' is still standing.

    We give 'winner' only those votes that rank 'winner' first.

    We give 'loser' ALL votes in which 'loser' appears in the ranking and
    'winner' does not, or 'loser' is ranked higher than 'winner'.

    This assertion "asserts" that the tally of 'winner' is larger than the
    tally of the 'loser'. This means that 'winner' could never be eliminated
    prior to 'loser'.
    */

    override fun is_vote_for_winner(cvr: Cvr): Int {
        if (!cvr.has_contest(this.contest_name)) return 0

        return if (ranking(this.winner, cvr.votes[this.contest_name]!!) == 0) 1 else 0
    }

    override fun is_vote_for_loser(cvr: Cvr): Int {
        if (!cvr.has_contest(this.contest_name)) return 0

        val w_idx = ranking(this.winner, cvr.votes[this.contest_name]!!)
        val l_idx = ranking(this.loser, cvr.votes[this.contest_name]!!)

        // return 1 if l_idx != -1 and (w_idx == -1 or (w_idx != -1 and l_idx < w_idx)) else 0
        return if (l_idx != -1 && (w_idx == -1 || (w_idx != -1 && l_idx < w_idx))) 1 else 0
    }

    override fun same_as(other: RaireAssertion): Boolean {
        return this.contest_name == other.contest_name && this.winner == other.winner
                && this.loser == other.loser
    }

    override fun subsumes(other: RaireAssertion): Boolean {
        /*
        An NEBAssertion 'A' subsumes an assertion 'other' if:
        - 'other' is not an NEBAssertion
        - Both assertions have the same winner & loser
        - 'other' rules out an outcome with the tail 'Tail' and either the
        winner of this NEBAssertion assertion appears before the loser in
        'Tail' or the loser appears and the winner does not.
        */

        //         if type(other) == NEBAssertion:
        //            return False
        //
        //        if self.winner == other.winner and self.loser == other.loser:
        //            return True
        //
        //        if self.winner == other.winner and not(self.loser in \
        //            other.eliminated):
        //            return True
        //
        //        elif self.winner in other.eliminated and not(self.loser in \
        //            other.eliminated):
        //            return True
        //
        //        else:
        //            # For all outcomes that 'other' is ruling out, this NEB
        //            # rules them all out.
        //            for ro in other.rules_out:
        //                idxw = -1 if not self.winner in ro else ro.index(self.winner)
        //                idxl = -1 if not self.loser in ro else ro.index(self.loser)
        //
        //                if idxw == idxl or (idxl < idxw):
        //                    return False
        //
        //            return True
        //
        //        return False

        if (other is NEBAssertion) return false

        if (this.winner == other.winner && this.loser == other.loser) return true

        /* TODO
        if (this.winner == other.winner && !(this.loser in other.eliminated)) {
            return true

        } else if (this.winner in other.eliminated && !(this.loser in other.eliminated)) {
            return true

        } else {
            //  For all outcomes that 'other' is ruling out, this NEB rules them all out.
            for (ro in other.rules_out) {
                val idxw: Int = if (ro.contains(this.winner)) -1 else ro.indexOf(this.winner)
                val idxl: Int = if (ro.contains(this.loser)) -1 else ro.indexOf(this.loser)

                if (idxw == idxl || (idxl < idxw)) {
                    return false
                }
            }
            return true
        } */
        return false
    }

    override fun to_str(): String {
        return "NEB Winner ${this.winner},Loser ${this.loser},diff est ${this.difficulty}"
    }
}

// Returns true if listb = some_list + lista
fun is_suffix(lista: List<Any>, listb: List<Any>): Boolean {
    val len_lista = lista.size
    val len_listb = listb.size

    if (len_listb < len_lista) return false

    val ss = listb.subList(len_listb - len_lista, listb.size)
    return ss == lista
}


class NENAssertion(contest_name: String, winner: String, loser: String, val eliminated: List<String>) :
    RaireAssertion(contest_name, winner, loser) {
    /*
    A Not-Eliminated-Next (NEN) assertion between a candidate 'winner' and
    a candidate 'loser' compares the tally of the two candidates in the
    context where a given set of candidates have been eliminated.

    We give 'winner' all votes in which they are preferenced first AFTER
    the candidates in 'eliminated' are removed from the ranking.

    We give 'loser' all votes in which they are preferenced first AFTER
    the candidates in 'eliminated' are removed from the ranking.

    This assertion "asserts" that the tally of the 'winner' in this context,
    where the specified candidates have been eliminated, is larger than that
    of 'loser'.
    */

    override fun is_vote_for_winner(cvr: Cvr): Int {
        if (!cvr.has_contest(this.contest_name)) return 0
        return vote_for_cand(this.winner, this.eliminated, cvr.votes[this.contest_name]!!)
    }

    override fun is_vote_for_loser(cvr: Cvr): Int {
        if (!cvr.has_contest(this.contest_name)) return 0
        return vote_for_cand(this.loser, this.eliminated, cvr.votes[this.contest_name]!!)
    }

    override fun same_as(other: RaireAssertion): Boolean {
        if (other !is NENAssertion) return false
        return this.contest_name == other.contest_name
                && this.winner == other.winner
                && this.loser == other.loser
                && this.eliminated == other.eliminated
    }

    override fun subsumes(other: RaireAssertion): Boolean {
        /*
        An NENAssertion 'A' subsumes an assertion 'other' if 'other' is
        not an NEBAssertion, the outcomes that 'A' rules out are suffixes of
        the outcomes that 'B' rules out.
        */

        //         if type(other) == NEBAssertion:
        //            return False
        //
        //        other_ro = set(other.rules_out)
        //
        //        for ro in self.rules_out:
        //            other_ro = [o for o in other_ro if not(is_suffix(ro, o))]
        //
        //        return other_ro == []

        if (other !is NENAssertion) return false

        // TODO not sure whats correct
        val other_ro = mutableListOf(other.rules_out)

        /* TODO
        // for ro in this.rules_out: other_ro = [o for o in other_ro if not(is_suffix(ro, o))]
        this.rules_out.forEach { ro ->
            other_ro.forEach { oro ->
                if (!is_suffix(ro, oro)) other_ro.add(oro)
            }
        }

         */

        return other_ro.isEmpty()
    }

/*
    override fun to_str() = buildString {
        append( "NEN Winner= ${winner} Loser= ${loser} Eliminated= ")

        for (cand in this.eliminated:
        result += ",{}".format(cand)

        result += ",diff est {}, rules out: {}".format(this.difficulty, \
        this.rules_out)
        return result
    }

 */
}


class RaireNode(val tail: List<String>) {
    //  Tail of an "imagined" elimination sequence representing the
    //  outcome of an IRV election. The last candidate in the tail is
    //  the "imagined" winner of the election.

    //  Lowest cost assertion that, if true, can rule out any election
    //  outcome that *ends* with the given tail.
    var best_assertion: RaireAssertion? = null

    //  An "ancestor" of this node is a node whose tail equals the latter
    //  part of val tail (i.e., if this.tail is ["A", "B", "C"], the node
    //  will have an ancestor with tail ["B", "C"].
    var best_ancestor: RaireNode? = null

    //  If there are candidates not mentioned in this.tail, this node is not a leaf and it can be expanded.
    var expandable = true

    //  Estimate of difficulty of ruling out the outcome this node represents.
    var estimate = Double.NaN

    //  Record of the children of this node that have already been
    //  considered (for example, through diving). These children are
    //  represented by the candidate that was added to the front of
    //  val tail when the child was created.
    var explored = mutableListOf<String>()

    //  Flag to indicate if node was created as part of a dive.
    var dive_node = false

    fun is_descendent_of(node: RaireNode): Boolean {
        /*
        Determines if the given 'node' is an ancestor of this node in a
        tree of possible election outcomes. A node with a tail equal to
        [a,b,c,d] has ancestors with tails [b,c,d], [c,d], and [d].

        Input:
        node: RaireNode     -  Potential ancestor

        Output:
        Returns true  if the input 'node' is an ancestor of this node, and
        false  otherwise.
        */
        val l1 = this.tail.size
        val l2 = node.tail.size

        if (l1 <= l2) return false
        return (this.tail.subList(l1 - l2, l1) == node.tail)
    }

    /*
    fun display(self, stream = sys.stdout) {
        print("{} | ".format(this.tail[0]), file = stream, end = '')

        for i in range(1, len(this.tail)):
        print("{} ".format(this.tail[i]), file = stream, end = '')

        print("[{}]".format(this.estimate), file = stream, end = '')

        if this.best_ancestor != null :
        print(
            " (Best Ancestor {} | ".format(this.best_ancestor.tail[0]),
            file = stream, end = ''
        )

        for i in range(1, len(this.best_ancestor.tail)):
        print(
            "{} ".format(this.best_ancestor.tail[i]), file = stream,
            end = ''
        )
        print(
            "[{}])".format(this.best_ancestor.estimate), file = stream,
            end = ''
        )

        print("")
    }

     */
}


class RaireFrontier {
    val nodes = mutableListOf<RaireNode>()

    fun replace_descendents(parent: RaireNode, log: Boolean) {
        /*
        Remove all descendents of the input 'node' from the frontier, and
        insert 'node' to the frontier in the appropriate position.

        If 'log' is true, print logging statements to given 'stream'.
        */
        val descendents = mutableListOf<Int>()

        /*
        if (log) {
            print("Replacing descendents of ", file = stream, end = '')
            node.display(stream = stream)
        }
         */

        this.nodes.forEachIndexed { index, node ->
            //  Is node a descendent of the given node?
            if (node.is_descendent_of(parent)) {
                descendents.add(index)
            }
        }

        descendents.reverse()
        for (i in descendents) {
            this.nodes.removeAt(i)
        }

        this.insert_node(parent)
    }


    fun insert_node(node: RaireNode) {
        /*
        Insert given node into the frontier in the right position. Nodes
        that are not associated with an "invalidating" assertion are placed
        at the front of the frontier. After these nodes, nodes in frontier
        are ordered from most difficult to invalidate to easiest to
        invalidate. Leaf nodes -- nodes whose "tail" contains all candidates
        -- are placed at the end of the frontier.

        Input:
        node: RaireNode   - node, representing an alternate election outcome, to add to the frontier.
        */
        if (!node.expandable) {
            this.nodes.add(node)
        } else if (node.estimate == Double.POSITIVE_INFINITY) {
            this.nodes.add(0, node)
        } else {
            var i = 0
            while (i < this.nodes.size) {
                val n_est = this.nodes[i].estimate
                if (n_est <= node.estimate) break
                i += 1
            }
            this.nodes.add(i, node)
        }
    }
}

fun find_best_audit(
    contest: Contest,
    ballots: List<Map<String, Int>>,
    neb_matrix: Map<String, MutableMap<String, NEBAssertion?>>,
    node: RaireNode,
    asn_func: EstimatorFn
) {
    /*
    Input:
    node: RaireNode    -  A node in the tree of alternate election outcomes.
                          The node represents an election outcome that ends in the sequence node.tail.

    contest: Contest   -  Contest being audited.

    ballots            - this is the list of cvrs.votes for this contest, candidateID -> rank
                        "Details of reported ballots for this contest"

    neb_matrix         -  |Candidates| x |Candidates| dictionary where
                        neb_matrix[c1][c2] returns a NEBAssertion stating
                        that c1 cannot be eliminated before c2 (if one
                        exists) and null  otherwise.

    asn_func: Callable -  Function that takes an assertion margin and
        returns an estimate of how "difficult" it will
        be to audit that assertion.

    Output:
    Finds the least cost assertion that can be used to rule out all election
    outcomes that end with the sequence node.tail, and assigns that assertion
    to node.best_assertion. If no such assertion can be found, node.assertion
    will equal null  after this function is called.
    */

    val ntail = node.tail.size
    val first_in_tail = node.tail[0]

    var best_asrtn: RaireAssertion? = null

    //  We first consider if we can invalidate this outcome by showing that
    //  'first_in_tail' can not-be-eliminated-before a candidate that
    //  appears later in tail.
    for (later_cand in node.tail.subList(1, ntail)) {
        //  Can we show that the candidate 'later_cand' must come before
        //  candidate 'first_in_tail' in the elimination sequence?
        val neb = neb_matrix[first_in_tail]!![later_cand]
        if (neb != null && (best_asrtn == null || neb.difficulty < best_asrtn.difficulty)) best_asrtn = neb
    }

    //  'eliminated' is the list of candidates that are not mentioned in 'tail'.
    //   eliminated = [c for c in contest.candidates if not c in node . tail]
    val eliminated = contest.candidates.filter { !node.tail.contains(it) }

    //  We now look at whether there is a candidate not mentioned in
    //  'tail' (this means they are assumed to be eliminated at some prior
    //  point in the elimination sequence), that can not-be-eliminated-before
    //  'first_in_tail'.
    for (cand in eliminated) {
        for (cand_in_tail in node.tail) {
            val neb = neb_matrix[cand]!![cand_in_tail]
            if (neb != null && (best_asrtn == null || neb.difficulty < best_asrtn.difficulty)) best_asrtn = neb
        }
    }

    //  We now consider whether we can find a better NEN assertion. We
    //  want to show that at the point where all the candidates in 'tail'
    //  remain, 'first_in_tail' is not the candidate with the least number
    //  of votes. This means that 'first_in_tail' should not be eliminated next.

    //  Tally of the candidate 'first_in_tail'
    // val tally_first_in_tail = sum([vote_for_cand(first_in_tail, eliminated, blt) for blt in ballots])
    val tally_first_in_tail = ballots.map { vote_for_cand(first_in_tail, eliminated, it) }.sum()

    for (later_cand in node.tail.subList(1, node.tail.size)) {
        // tally_later_cand = sum([vote_for_cand(later_cand, eliminated, blt) for blt in ballots])
        val tally_later_cand = ballots.map { vote_for_cand(later_cand, eliminated, it) }.sum()

        if (tally_first_in_tail > tally_later_cand) {
            //  We can create a NEN assertion that says "first_in_cand"
            //  should not be eliminated next, after "eliminated" are
            //  eliminated, because "later_cand" actually has less votes at this point.
            val estimate = asn_func(
                tally_first_in_tail,
                tally_later_cand,
                contest.tot_ballots - (tally_first_in_tail + tally_later_cand),
                contest.tot_ballots
            )

            if (best_asrtn == null || estimate < best_asrtn.difficulty) {
                val nen = NENAssertion(contest.name, first_in_tail, later_cand, eliminated)

                // AFAICT, this is the only place that rules_out is added to?
                // rules_out has assewrtions, node.tail has Strings
                // nen.rules_out.add(tuple(node.tail))
                // nen.rules_out.addAll(node.tail) TODO
                nen.difficulty = estimate

                nen.votes_for_winner = tally_first_in_tail
                nen.votes_for_loser = tally_later_cand

                best_asrtn = nen
            }
        }
    }

    node.best_assertion = best_asrtn

    if (best_asrtn != null) {
        node.estimate = best_asrtn.difficulty
    }
}

data class ManageNodeReturn(val audit_not_possible: Boolean, val lowerBound: Double, val terminus: Boolean)

fun manage_node(newn: RaireNode, frontier: RaireFrontier, lowerbound: Double, log: Boolean): ManageNodeReturn {

    /*
Input:

newn: RaireNode    -  A node in the tree of alternate election outcomes that
has just been created and evaluated, but not yet
added to our frontier. We need to determine what this
node's evaluation means for our frontier.

frontier           -  Current frontier of our set of alternate outcome
trees.

lowerbound         -  Current lower bound on audit difficulty.

log                -  Flag indicating if logging statements should
be printed during the algorithm.

stream             -  Stream to which logging statements should
be printed.


Output:

Returns a triple:
audit_not_possible (Boolean), new lower bound, terminus (Boolean)

The first element of this triple is a boolean indicating whether or not
we have established that the audit is not possible. If so, this boolean
will be true , otherwise it will be false .

The second element indicates the new lower bound on audit difficulty
as a result of the node's evaluation (note it may not have changed from
the prior lower bound).

The third element indicates whether or not we will need to continue
exploring children of this node. The boolean 'terminus' will be set to
true  if we do not need to continue to explore children of this node, and
false  otherwise.

*/

    //     if not newn.expandable:
    //        # 'newn' is a leaf.
    //        if newn.estimate == np.inf and newn.best_ancestor.estimate == np.inf:
    //
    //            if log:
    //                print("Found branch that cannot be pruned.", file=stream)
    //
    //            return True, np.inf, True
    //
    //        if newn.best_ancestor.estimate <= newn.estimate:
    //            next_lowerbound = max(lowerbound, newn.best_ancestor.estimate)
    //            frontier.replace_descendents(newn.best_ancestor,log,stream=stream)
    //
    //            return False, next_lowerbound, True
    //
    //        else:
    //            next_lowerbound = max(lowerbound, newn.estimate)
    //            frontier.insert_node(newn)
    //
    //            if log:
    //                print("    Best audit ", file=stream, end='')
    //                newn.best_assertion.display(stream=stream)
    //
    //            return False, next_lowerbound, True
    //    else:
    //        frontier.insert_node(newn)
    //
    //        if log:
    //            if newn.best_assertion != None:
    //                print("    Best audit ", file=stream, end='')
    //                newn.best_assertion.display(stream=stream)
    //            else:
    //                print("    Cannot be disproved", file=stream)
    //
    //        return False, lowerbound, False

    if (!newn.expandable) {
        //  'newn' is a leaf.
        if (newn.estimate == Double.POSITIVE_INFINITY && newn.best_ancestor!!.estimate == Double.POSITIVE_INFINITY) {
            if (log) println("Found branch that cannot be pruned.")
            return ManageNodeReturn(true, Double.POSITIVE_INFINITY, true)
        }

        if (newn.best_ancestor!!.estimate <= newn.estimate) {
            val next_lowerbound = max(lowerbound, newn.best_ancestor!!.estimate)
            frontier.replace_descendents(newn.best_ancestor!!, log)
            return ManageNodeReturn(false, next_lowerbound, true)
        } else {
            val next_lowerbound = max(lowerbound, newn.estimate)
            frontier.insert_node(newn)
            if (log) {
                print("    Best audit ")
                //newn.best_assertion.display(stream = stream)
            }
            return ManageNodeReturn(false, next_lowerbound, true)
        }

    } else {
        frontier.insert_node(newn)

        if (log) {
            if (newn.best_assertion != null) println("    Best audit ") else println("    Cannot be disproved")
        }

        return ManageNodeReturn(false, lowerbound, false)
    }
}

fun perform_dive(
    node: RaireNode,
    contest: Contest,
    ballots: List<Map<String, Int>>,
    neb_matrix: Map<String, MutableMap<String, NEBAssertion?>>,
    asn_func: EstimatorFn,
    lower_bound: Double,
    frontier: RaireFrontier,
    log: Boolean
): Double {

    /*
    Input:
    node: RaireNode    -  A node in the tree of alternate election outcomes.
    Starting point of dive to a leaf.

    contest: Contest   -  Contest being audited.

    ballots:           -  Details of reported ballots for this contest.

    neb_matrix         -  |Candidates| x |Candidates| dictionary where
    neb_matrix[c1][c2] returns a NEBAssertion stating
    that c1 cannot be eliminated before c2 (if one
    exists) and null  otherwise.

    asn_func: Callable -  Function that takes an assertion margin and
    returns an estimate of how "difficult" it will
    be to audit that assertion.

    lower_bound        -  Current lower bound on audit difficulty.

    frontier           -  Current frontier of our set of alternate outcome
    trees.

    log                -  Flag indicating if logging statements should
    be printed during the algorithm.


    Output:
    Returns the difficulty estimate of the least-difficult-to-audit
    assertion that can be used to rule out at least one of the branches
    starting at the input 'node'. As this function dives from the given 'node'
    it will add nodes to the current frontier of our set of alternate outcome
    trees.
    */

//    ncands = len(contest.candidates)
//
//    rem_cands = [c for c in contest.candidates if not c in node.tail]
//
//    # sort rem_cands by position in contest.order if it is defined
//    next_cand = rem_cands[0]
//    if contest.outcome != []:
//        npos = contest.outcome.index(next_cand)
//
//        for i in range(1, len(rem_cands)):
//            c = rem_cands[i]
//            ipos = contest.outcome.index(c)
//
//            if ipos > npos:
//                next_cand = c
//                npos = ipos
//
//    newn = RaireNode([next_cand] + node.tail)
//    newn.expandable = False if len(newn.tail) == ncands else True
//    newn.dive_node = True
//
//    node.explored.append(next_cand)
//
//    # Assign a 'best ancestor' to the new node.
//    newn.best_ancestor = node.best_ancestor if \
//        node.best_ancestor != None and node.best_ancestor.estimate <= \
//        node.estimate else node
//
//    find_best_audit(contest, ballots, neb_matrix, newn, asn_func)
//
//    if log:
//        print("DIVE TESTED ", file=stream, end='')
//        newn.display(stream=stream)
//
//    audit_not_possible, next_lowerbound, dive_complete = manage_node(newn, \
//        frontier, lower_bound, log, stream=stream)
//
//    if audit_not_possible:
//        return np.inf
//
//    if dive_complete:
//        return next_lowerbound
//
//    return perform_dive(newn, contest, ballots, neb_matrix, asn_func, \
//            next_lowerbound, frontier, log, stream=stream)
    val ncands = contest.candidates.size

    // rem_cands = [c for c in contest.candidates if not c in node.tail]
    val rem_cands = contest.candidates.filter { !node.tail.contains(it) }

    //  sort rem_cands by position in contest.order if it is defined
    var next_cand = rem_cands[0] // TODO empty ?
    /* if (contest.outcome.isNotEmpty()) {
        var npos = contest.outcome.index(next_cand)

        for (i in 1..rem_cands.size) { // TODO check
            val c = rem_cands[i]
            val ipos = contest.outcome.indexOf(c)

            if (ipos > npos) {
                next_cand = c
                npos = ipos
            }
        }
    } */ // TODO

    val newn = RaireNode(listOf(next_cand) + node.tail)
    newn.expandable = (newn.tail.size != ncands)
    newn.dive_node = true

    node.explored.add(next_cand)

    //  Assign a 'best ancestor' to the new node.
    newn.best_ancestor = if (node.best_ancestor != null && node.best_ancestor!!.estimate <= node.estimate) node.best_ancestor else node

    find_best_audit(contest, ballots, neb_matrix, newn, asn_func)

    if (log) {
        println("DIVE TESTED ")
        //newn.display(stream = stream)
    }

    val (audit_not_possible, next_lowerbound, dive_complete) = manage_node(newn, frontier, lower_bound, log)

    if (audit_not_possible) {
        return Double.POSITIVE_INFINITY
    }

    if (dive_complete) {
        return next_lowerbound
    }

    return perform_dive(
        newn, contest, ballots, neb_matrix, asn_func,
        next_lowerbound, frontier, log,
    )
}
