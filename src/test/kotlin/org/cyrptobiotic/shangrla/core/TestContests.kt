package org.cyrptobiotic.shangrla.core

import org.cryptobiotic.shangrla.core.AuditType
import org.cryptobiotic.shangrla.core.Contest
import org.cryptobiotic.shangrla.core.SocialChoiceFunction
import org.junit.jupiter.api.Assertions.assertEquals
import kotlin.test.Test
import kotlin.test.assertTrue

class TestContests {
    fun makeContests():  List<Contest> {
        val contests = listOf(
            Contest( id = "AvB",
                name = "contest_1",
                candidates = listOf("alice","bob","carol","dave","erin"),
                choice_function = SocialChoiceFunction.PLURALITY,
                audit_type = AuditType.CARD_COMPARISON,
                risk_limit = 0.05,
                ncards = 10000,
                n_winners = 2,
                reported_winners = listOf("alice","bob"),
                use_style = true,
            ),
            Contest( id = "CvD",
                name = "contest_2",
                candidates = listOf("alice","bob","carol","dave"),
                choice_function = SocialChoiceFunction.SUPERMAJORITY,
                audit_type = AuditType.POLLING,
                risk_limit = 0.05,
                ncards = 10000,
                n_winners = 1,
                reported_winners = listOf("alice"),
                use_style = false,
            ),
            Contest( id = "EvF",
                name = "contest_3",
                candidates = listOf("alice","bob","carol","dave"),
                choice_function = SocialChoiceFunction.IRV,
                audit_type = AuditType.CARD_COMPARISON,
                risk_limit = 0.05,
                ncards = 10000,
                n_winners = 1,
                reported_winners = listOf("alice"),
                use_style = false,
            ),
        )
        return contests
    }

    @Test
    fun test_tally() {
        val cvrs = CvrBuilders()
            .add(id = "1", ContestVotes("AvB", "Alice"), ContestVotes("CvD", "Candy"))
            .add(id = "2", ContestVotes("AvB", "Bob"), ContestVotes("CvD", Vote("Elvis"), Vote("Candy", 0)))
            .add(id = "3", ContestVotes("EvF",  Vote("Bob"), Vote("Edie", 2)),
                ContestVotes("CvD",  Vote("Elvis", 0), Vote("Candy")))
            .add(id = "4", ContestVotes("AvB", "Alice"), ContestVotes("CvD", "Candy"))
            .add(id = "5", ContestVotes("AvB", "Bob"), ContestVotes("CvD", Vote("Elvis"), Vote("Candy", 0)))
            .add(id = "6", ContestVotes("EvF", Vote("Bob", 2), Vote("Edie")), ContestVotes("CvD", Vote("Elvis", 0), Vote("Candy")))
            .add(id = "7", ContestVotes("AvB", Vote("Alice", 2)), ContestVotes("CvD", Vote("Elvis", 0), Vote("Candy")))
            .build()

//        cvr_dict = [{ "id": 1, "votes": { "AvB": { "Alice":True }, "CvD": { "Candy":True } } },
//            { "id": 2, "votes": { "AvB": { "Bob":True }, "CvD": { "Elvis":True, "Candy":False } } },
//            { "id": 3, "votes": { "EvF": { "Bob":1, "Edie":2 }, "CvD": { "Elvis":False, "Candy":True } } },
//            { "id": 4, "votes": { "AvB": { "Alice":1 }, "CvD": { "Candy":"yes" } } },
//            { "id": 5, "votes": { "AvB": { "Bob":True }, "CvD": { "Elvis":True, "Candy":False } } },
//            { "id": 6, "votes": { "EvF": { "Bob":2, "Edie":1 }, "CvD": { "Elvis":False, "Candy":True } } },
//            { "id": 7, "votes": { "AvB": { "Alice":2 }, "CvD": { "Elvis":False, "Candy":True } } }]
//
// TODO what does alice has 2 votes mean?. Note its normalized to non-zero in Contest.tally
//    Do the same in Cvr.tabulate_* ??

        val contests = makeContests()
        Contest.tally(contests, cvrs)
        assertEquals(mapOf("Alice" to 3, "Bob" to 2), contests[0].tally)
        assertEquals(mapOf("Candy" to 5, "Elvis" to 2), contests[1].tally)
        assertTrue(contests[2].tally.isEmpty())
    }
}