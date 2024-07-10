package org.cryptobiotic.shangrla

import kotlin.math.ln

/*
fun from_vote(vote: String, id: String, contest_id: String, phantom: Boolean = false): CVR {
    /*
    Wraps a vote and creates a CVR, for unit tests

    Parameters:
    ----------
    vote: dict of votes in one contest
    id: str
    CVR id
    contest_id: str
    identifier of the contest

    Returns:
    --------
    CVR containing that vote in the contest "AvB", with CVR id=1.
    */
    return CVR(id = id, votes = { contest_id: vote }, phantom = phantom)
}

 */

class Bernoulli(p: Double) {
    val log_q = ln(1.0 - p)
    val n = 1.0

    fun get(): Double {
        var x = 0.0
        var sum = 0.0
        while (true) {
            val wtf = ln( Math.random()) / (n - x)
            sum += wtf
            if (sum < log_q) {
                return x
            }
            x++
        }
    }
}