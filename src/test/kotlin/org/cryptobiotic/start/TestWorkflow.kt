package org.cryptobiotic.start

import kotlin.test.Test
import kotlin.test.assertEquals

class TestWorkflow {

    @Test
    fun testWorkflow() {
        val audit = Audit()

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
        val min_margin = audit.set_all_margins_from_cvrs(contests, cvrs)
        // println("minimum assorter margin = ${min_margin}")

        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }

        ///* Set up for sampling
        val sample_size = audit.find_sample_size(contests, cvrs=cvrs)
        println("sample_size = ${sample_size}")

        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }

        val samples = audit.assign_sample_nums(cvrs, sample_size).toList()

        // Tst 1. suppose there are no errors, so that mvr == cvr
        // Compute the p values
        val p_max = Assertion.set_p_values(contests=contests, mvr_sample=samples, cvr_sample=samples)
        println("p_max = ${p_max}")

        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }

        assertEquals(29, sample_size)
    }
}