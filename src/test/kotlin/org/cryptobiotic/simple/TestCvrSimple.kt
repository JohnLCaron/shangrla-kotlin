package org.cryptobiotic.simple

import kotlin.test.*

class TestCvrSimple {

    @Test
    fun testGetVote() {
        val votes = mapOf( "id" to mapOf("cand1" to 1, "cand2" to 0),)
        val cvr = CvrSimple("id", votes)
        assertEquals(1, cvr.get_vote_for("id", "cand1"))
        assertEquals(0, cvr.get_vote_for("id", "cand2"))
        assertEquals(0, cvr.get_vote_for("id", "cand3"))
    }
}