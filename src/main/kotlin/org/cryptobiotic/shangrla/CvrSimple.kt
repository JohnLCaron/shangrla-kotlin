package org.cryptobiotic.shangrla

class CvrSimple(
    val id: String,
    val votes: Map<String, Map<String, Int>>, // contest/vote dict
    val phantom: Boolean = false,
) {

    fun has_contest(contestId: String): Boolean {
        return votes[contestId] != null
    }
}

fun makeCvrSimple(cvrId: String, phantom: Boolean, vararg votes: Vote): CvrSimple {
    val voteMap = mutableMapOf<String, MutableMap<String, Int>>() // contest/vote dict
    votes.forEach {
        val contestVotes = voteMap.getOrPut(it.contestId) { mutableMapOf() }
        val candidateVote = contestVotes.getOrPut(it.candidateId) { 0 }
        contestVotes[it.candidateId] = candidateVote + it.vote
    }
    return CvrSimple(cvrId, voteMap, phantom)
}

data class Vote(val contestId: String, val candidateId: String, val vote: Int = 1)


fun make_phantoms(max_cards: Int, cvr_list: List<CvrSimple>, contests: List<ContestSimple>,
                  use_style: Boolean=true, prefix: String = ""): Pair<List<CvrSimple>, Int> {
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
    use_style : Boolean; does the sampling use style information?

    Returns
    -------
    cvr_list : list of CVR objects; the reported CVRs and the phantom CVRs
    n_phantoms : int; number of phantom cards added

    Side effects
    ------------
    for each contest in `contests`, sets `cards` to max_cards if not specified by the user
    for each contest in `contests`, set `cvrs` to be the number of (real) CVRs that contain the contest
    */
    //        phantom_vrs = []
    //        n_cvrs = len(cvr_list)
    //        for c, v in contests.items():  # set contest parameters
    //            v['cvrs'] = np.sum([cvr.has_contest(c) for cvr in cvr_list if not cvr.is_phantom()])
    //            v['cards'] = max_cards if v['cards'] is None else v['cards'] // upper bound on cards cast in the contest
    //        if not use_style:              #  make (max_cards - len(cvr_list)) phantoms
    //            phantoms = max_cards - n_cvrs
    //            for i in range(phantoms):
    //                phantom_vrs.append(CVR(id=prefix+str(i+1), votes={}, phantom=True))
    //        else:                          # create phantom CVRs as needed for each contest
    //            for c, v in contests.items():
    //                phantoms_needed = v['cards']-v['cvrs']
    //                while len(phantom_vrs) < phantoms_needed:
    //                    phantom_vrs.append(CVR(id=prefix+str(len(phantom_vrs)+1), votes={}, phantom=True))
    //                for i in range(phantoms_needed):
    //                    phantom_vrs[i].votes[c]={}  # list contest c on the phantom CVR
    //            phantoms = len(phantom_vrs)
    //        cvr_list = cvr_list + phantom_vrs
    //        return cvr_list, phantoms

    val phantom_vrs = mutableListOf<CvrSimple>()
    var n_phantoms: Int
    val n_cvrs = cvr_list.size
    for (contest in contests) { // } set contest parameters
        // TODO these are intended to be set on the contest
        val cvrs = cvr_list.filter{ !it.phantom && it.has_contest(contest.id) }
        val cards = contest.cards ?: max_cards // upper bound on cards cast in the contest
    }

    if (!use_style) {              //  make (max_cards-len(cvr_list)) phantoms
        n_phantoms = max_cards - n_cvrs
        repeat(n_phantoms) {
            phantom_vrs.add( CvrSimple("$prefix${it + 1}", votes = emptyMap(), phantom = true))
        }
    } else {                         // create phantom CVRs as needed for each contest
        for (contest in contests) {
            val phantoms_needed = contest.cards!! // TODO - contest.cvrs
            repeat(phantoms_needed) {
                val votes = mutableMapOf<String, Map<String, Int>>()
                votes[contest.id] = emptyMap() // list contest c on the phantom CVR
                phantom_vrs.add( CvrSimple("$prefix${it + 1}", votes = votes, phantom = true))
            }
        }
        n_phantoms = phantom_vrs.size
    }
    val result = cvr_list + phantom_vrs
    return Pair(result, n_phantoms)
}

