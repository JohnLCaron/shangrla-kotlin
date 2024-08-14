package org.cryptobiotic.simulation

import org.cryptobiotic.shangrla.core.AuditType
import org.cryptobiotic.start.*
import kotlin.test.Test
import kotlin.test.assertEquals

class TestSimulation {

    @Test
    fun testWorkflow() {
        val audit = Audit(auditType = AuditType.CARD_COMPARISON)

        val cvrs = makeCvrsByCount(1000, 599)
        println("ncvrs = ${cvrs.size}")

        val votes: Map<String, Map<String, Int>> = Cvr.tabulateVotes(cvrs) // contest -> candidate -> count
        votes.forEach { key, cands ->
            println("contest ${key} ")
            cands.forEach { println("  ${it}") }
        }

        // make contests from cvrs
        val contests: List<Contest> = Contest.fromVotes(audit, votes, Cvr.cardsPerContest(cvrs))
        println("Contests")
        contests.forEach { println("  ${it}") }

        //// Create Assertions for every Contest, including an Assorter and NonnegMean for every Assertion
        audit.makeAssertions(contests)

        // Create CVRs for phantom cards
        // skip for now, no phantoms

        // sets margins on the assertions
        audit.set_all_margins_from_cvrs(contests, cvrs)
        // println("minimum assorter margin = ${min_margin}")

        // Set up for sampling
        val sample_size = audit.find_sample_size(contests, cvrs=cvrs)
        //println("sample_size = ${sample_size}")
        // val sample_size = 100

        val samples = audit.assign_sample_nums(cvrs, sample_size).toList()

        // Tst 1. suppose there are no errors, so that mvr == cvr Compute the p values
        val p_max = Assertion.set_p_values(contests=contests, mvr_sample=samples, cvr_sample=samples)
        println("p_max = ${p_max}")

        assertEquals(29, sample_size)
    }
}