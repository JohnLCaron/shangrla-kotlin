package org.cryptobiotic.simple

fun makeCvrsFromRaireBallots(
    ballots: Map<String, MutableMap<String, MutableMap<String, Int>>>,
    limit: Int = Integer.MAX_VALUE
): List<CvrSimple> {
    val result = mutableListOf<CvrSimple>()
    // Map<String, MutableMap<String, MutableMap<String, Int>>> ballotId: contestId : candidateId : count
    var ballotIter = ballots.entries.iterator()
    for (idx in 0..limit) {
        if (!ballotIter.hasNext()) break
        val (ballotId, cvr) = ballotIter.next()
        result.add(CvrSimple(ballotId, cvr))
    }
    return result
}

// change from IRV to plurality contest, for testing
fun makeCvrsFromRaireBallotsPlurality(
    ballots: Map<String, MutableMap<String, MutableMap<String, Int>>>, // ballotId : contestId: candidateId : order
    limit: Int = Integer.MAX_VALUE
): List<CvrSimple> {
    val result = mutableListOf<CvrSimple>()
    var ballotIter = ballots.entries.iterator()
    for (idx in 0..limit) {
        if (!ballotIter.hasNext()) break
        val (ballotId: String, cvr) = ballotIter.next()
        val cvrPlurality = mutableMapOf<String, MutableMap<String, Int>>()

        cvr.entries.forEach { (contestId, v1) ->
            val votesPlurality = mutableMapOf<String, Int>()
            v1.forEach { (candId, vote) ->
                if (vote == 1) votesPlurality[candId] = 1
            }
            cvrPlurality[contestId] = votesPlurality
        }

        result.add(CvrSimple(ballotId, cvrPlurality))
    }
    return result
}