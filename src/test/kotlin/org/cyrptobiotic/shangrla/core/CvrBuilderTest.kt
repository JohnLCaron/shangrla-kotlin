package org.cyrptobiotic.shangrla.core

import org.junit.jupiter.api.Assertions.assertFalse
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CvrBuilderTest {

    @Test
    fun test_set_tally_pool_means() {

        val cvrBuilders = CvrBuilders()
            .add(
                id = "1", tally_pool = "1",
                ContestVotes("AvB", "Alice", 1),
                ContestVotes("CvD", "Candy", true)
            )
            .add(
                "2", "1",
                ContestVotes.add("CvD", Vote("Elvis", true), Vote("Candy", false)),
                ContestVotes("EvF")
            )
            .add(id = "3", tally_pool = "1", ContestVotes("GvH"))
            .add(
                "4", "2",
                ContestVotes("AvB", "Bob", 1),
                ContestVotes("CvD", "Candy", true)
            )
            .add(
                "5", "2",
                ContestVotes.add("CvD", Vote("Elvis", true), Vote("Candy", false)),
                ContestVotes("EvF")
            )

        val b1 = cvrBuilders.builders[0]
        assertTrue(b1.has_contest("AvB"))
        assertFalse(b1.has_contest("EvF"))
        val b1c1 = b1.contests["AvB"]!! // Map<String, Int>
        assertEquals(1, b1c1.size)
        assertEquals(1, b1c1["Alice"])
        val b1c2 = b1.contests["CvD"]!!
        assertEquals(1, b1c2.size)
        assertEquals(1, b1c2["Candy"])

        val b2 = cvrBuilders.builders[1]
        assertTrue(b2.has_contest("CvD"))
        assertFalse(b2.has_contest("AvB"))
        val b2c1 : Map<String, Int> = b2.contests["CvD"]!!
        assertEquals(2, b2c1.size)
        assertEquals(1, b2c1["Elvis"])
        assertEquals(0, b2c1["Candy"])

        println("1 ${cvrBuilders.show()}")

        val tally_pool = mutableMapOf<String, Set<String>>()
        for (p in cvrBuilders.poolSet()) {
            // tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
            tally_pool[p] = cvrBuilders.poolContests(p)
        }

        cvrBuilders.add_pool_contests(tally_pool)
        println("2 ${cvrBuilders.show()}")
    }
}