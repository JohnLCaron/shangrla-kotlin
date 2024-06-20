package org.cryptobiotic.shangrla.reader

import org.cryptobiotic.shangrla.raire.RaireContest
import java.io.File

// Data file in .raire format.
fun load_contests_from_raire(fileName: String): Pair<List<RaireContest>, Map<String, MutableMap<String, MutableMap<String, Int>>>> {

//  A map between ballot id and the relevant CVR.
    val lines = File(fileName).bufferedReader().readLines()

//  Total number of contests described in data file
    var lineIndex = 0
    val ncontests = lines[lineIndex++].toInt()

//  Map between contest id and number of ballots involving that contest
    val num_ballots = mutableMapOf<String, Int>()

//  Map between contest id and the candidates & winner & order of that contest.
    val contest_info = mutableMapOf<String, ContestInfo>()

    repeat(ncontests) { contestIdx ->
        // toks = [line.strip() for line in lines[1 + i].strip().split(',')]
        // Contest,1,5,1,2,3,4,5,winner,4
        val toks: List<String> = lines[lineIndex++].split(",")

    //  Get contest id and number of candidates in that contest
        val cid = toks[1] // contest id
        val ncands = toks[2].toInt() // ncandidates

    //  Get list of candidate identifiers
        val cands = mutableListOf<String>()
        repeat(ncands) { cands.add(toks[3 + it]) }

        val windx = toks.indexOf("winner")
        val winner = toks[windx + 1]

        val inf_index = toks.indexOf("informal")
        val informal = if (inf_index < 0) 0 else toks[inf_index + 1].toInt()

        val order = mutableListOf<Int>()
        if (toks.contains("order")) {
            if (inf_index < 0) for (tokidx in windx + 2 until toks.size) order.add(toks[tokidx].toInt())
            else for (tokidx in windx + 2 until inf_index) order.add(toks[tokidx].toInt())
        }

        contest_info[cid] = ContestInfo(cands, winner, order)
        num_ballots[cid] = informal
    }

    // TODO I dont see how this code works for more than one contest

    //         for l in range(ncontests+1,len(lines)):
    //            toks = [line.strip() for line in lines[l].strip().split(',')]
    //
    //            cid = toks[0]
    //            bid = toks[1]
    //            prefs = toks[2:]
    //
    //            ballot = {}
    //            for c in contest_info[cid][0]:
    //                if c in prefs:
    //                    idx = prefs.index(c)
    //                    ballot[c] = idx
    //
    //            num_ballots[cid] += 1
    //
    //            if not bid in cvrs:
    //                cvrs[bid] = {cid: ballot}
    //            else:
    //                cvrs[bid][cid] = ballot

    val cvrs = mutableMapOf<String, MutableMap<String, MutableMap<String, Int>>>() //  ballotId -> contestId -> candidate -> Preference rank
    while (lineIndex < lines.size) {
        // 1,1,4,1,2,3
        // 1,2,4,1,2,3
        // ...
        val toks: List<String> = lines[lineIndex++].split(",")
        // toks = [line.strip() for line in lines[l].strip().split(',')]

        val cid = toks[0] // contest id
        val bid = toks[1] // ballot id
        val prefs = toks.subList(2, toks.size) // remaining "order"

        val ballot = mutableMapOf<String, Int>() // candidate to preference rank
        for (candidate in contest_info[cid]!!.candidates) {
            if (candidate in prefs) {
                ballot[candidate] = prefs.indexOf(candidate)
            }
        }
        val ncid = num_ballots[cid]!!
        num_ballots[cid] = ncid + 1

        val conBallotMap = cvrs.getOrPut(bid) { mutableMapOf() }
        conBallotMap[cid] = ballot
    }

    val raireContests = mutableListOf<RaireContest>()

    for ((cid, trip) in contest_info) {
            val (cands, winner, order) = trip
            val con = RaireContest(cid, cands, winner, num_ballots[cid]!!, outcome = IntArray(order.size) { order[it] })
            raireContests.add(con)
        }

    // TODO translate cvrs into immutables
    return Pair(raireContests, cvrs)
}

data class ContestInfo(val candidates: List<String>, val winner: String, val order: List<Int>)