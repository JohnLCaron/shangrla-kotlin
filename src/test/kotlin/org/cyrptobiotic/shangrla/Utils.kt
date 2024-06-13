package org.cyrptobiotic.shangrla

import org.cryptobiotic.shangrla.core.CVR

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