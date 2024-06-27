package org.cyrptobiotic.shangrla.core

import org.cryptobiotic.shangrla.core.Cvr
import kotlin.test.*

class TestCvr {

    @Test
    fun test_rcv_lfunc_wo() {
        val votes = ContestVotes("AvB", Vote("Alice"), Vote("Bob", 2), Vote("Candy", 3), Vote("Dan", 0))
        val cvr = CvrBuilder("id").addContestVotes(votes).build()
        assertEquals(1, cvr.rcv_lfunc_wo("AvB", "Bob", "Alice"))
        assertEquals(0, cvr.rcv_lfunc_wo("AvB", "Alice", "Candy"))
        assertEquals(1, cvr.rcv_lfunc_wo("AvB", "Dan", "Candy"))

        //     def test_rcv_lfunc_wo(self):
        //        votes = CVR.from_vote({"Alice": 1, "Bob": 2, "Candy": 3, "Dan": ''})
        //        assert votes.rcv_lfunc_wo("AvB", "Bob", "Alice") == 1
        //        assert votes.rcv_lfunc_wo("AvB", "Alice", "Candy") == 0
        //        assert votes.rcv_lfunc_wo("AvB", "Dan", "Candy") == 1
    }

    @Test
    fun test_rcv_votefor_cand() {
        val votes = ContestVotes("AvB", Vote("Alice"), Vote("Bob", 2), Vote("Candy", 3),
            Vote("Dan", 0), Vote("Ross", 4), Vote("Aaron", 5))
        val cvr = CvrBuilder("id").addContestVotes(votes).build()

        var remaining = listOf("Bob", "Dan", "Aaron", "Candy")
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Candy", remaining))
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Alice", remaining))
        assertEquals(1, cvr.rcv_votefor_cand("AvB", "Bob", remaining))
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Aaron", remaining))

        remaining = listOf("Dan", "Aaron", "Candy")
        assertEquals(1, cvr.rcv_votefor_cand("AvB", "Candy", remaining))
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Alice", remaining))
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Bob", remaining))
        assertEquals(0, cvr.rcv_votefor_cand("AvB", "Aaron", remaining))

        //         votes = CVR.from_vote({"Alice": 1, "Bob": 2, "Candy": 3, "Dan": '', "Ross": 4, "Aaron": 5})
        //        remaining = ["Bob","Dan","Aaron","Candy"]
        //        assert votes.rcv_votefor_cand("AvB", "Candy", remaining) == 0
        //        assert votes.rcv_votefor_cand("AvB", "Alice", remaining) == 0
        //        assert votes.rcv_votefor_cand("AvB", "Bob", remaining) == 1
        //        assert votes.rcv_votefor_cand("AvB", "Aaron", remaining) == 0
        //
        //        remaining = ["Dan","Aaron","Candy"]
        //        assert votes.rcv_votefor_cand("AvB", "Candy", remaining) == 1
        //        assert votes.rcv_votefor_cand("AvB", "Alice", remaining) == 0
        //        assert votes.rcv_votefor_cand("AvB", "Bob", remaining) == 0
        //        assert votes.rcv_votefor_cand("AvB", "Aaron", remaining) == 0
    }

    @Test
    fun test_cvr_from_dict() {
        //         cvr_dict = [{'id': 1, 'pool': True, 'tally_pool': 1, 'votes': {'AvB': {'Alice':True}, 'CvD': {'Candy':True}}},
        //                    {'id': 2, 'sample_num': 0.2, 'p': 0.5, 'sampled': True,
        //                              'votes': {'AvB': {'Bob':True}, 'CvD': {'Elvis':True, 'Candy':False}}},
        //                    {'id': 3, 'tally_pool': 'abc', 'votes': {'EvF': {'Bob':1, 'Edie':2}, 'CvD': {'Elvis':False, 'Candy':True}}}]
        //        cvr_list = CVR.from_dict(cvr_dict)

        val cvrs = CvrBuilders()
            .add(id = "1", tally_pool = "1", pool = true)
                .setContestVotes("1", ContestVotes("AvB", "Alice"), ContestVotes("CvD", "Candy"))
            .add(id = "2", sample_num = 2, p = 0.5, sampled = true)
                .setContestVotes("2", ContestVotes("AvB", "Bob"), ContestVotes("CvD", Vote("Elvis", 2), Vote("Candy", 0)))
            .add(id = "3", tally_pool = "abc")
                .setContestVotes("3", ContestVotes("EvF", Vote("Bob"), Vote("Edie", 2)), ContestVotes("CvD", Vote("Elvis", 0), Vote("Candy", 1)))
            .build()

        assertEquals(3, cvrs.size)
        assertEquals("1", cvrs[0].id)
        assertEquals("2", cvrs[1].id)
        assertEquals("3", cvrs[2].id)

        assertEquals(1, cvrs[0].get_vote_for("AvB", "Alice"))
        assertEquals(1, cvrs[0].get_vote_for("CvD", "Candy"))
        assertEquals(0, cvrs[0].get_vote_for("AvB", "Bob"))
        assertEquals(0, cvrs[0].get_vote_for("EvF", "Bob"))

        assertEquals(0, cvrs[1].get_vote_for("AvB", "Alice"))
        assertEquals(0, cvrs[1].get_vote_for("CvD", "Candy"))
        assertEquals(2, cvrs[1].get_vote_for("CvD", "Elvis"))
        assertEquals(0, cvrs[1].get_vote_for("CvD", "Candy"))
        assertEquals(0, cvrs[1].get_vote_for("CvD", "Edie"))
        assertEquals(1, cvrs[1].get_vote_for("AvB", "Bob"))
        assertEquals(0, cvrs[1].get_vote_for("EvF", "Bob"))

        assertEquals(0, cvrs[2].get_vote_for("AvB", "Alice"))
        assertEquals(1, cvrs[2].get_vote_for("CvD", "Candy"))
        assertEquals(0, cvrs[2].get_vote_for("CvD", "Edie"))
        assertEquals(0, cvrs[2].get_vote_for("AvB", "Bob"))
        assertEquals(1, cvrs[2].get_vote_for("EvF", "Bob"))
        assertEquals(2, cvrs[2].get_vote_for("EvF", "Edie"))
        assertEquals(0, cvrs[2].get_vote_for("EvF", "Alice"))

        assertTrue(cvrs[0].pool)
        assertEquals("1", cvrs[0].tally_pool)
        assertEquals(2, cvrs[1].sample_num)
        assertEquals(0.5, cvrs[1].p)
        assertTrue(cvrs[1].sampled)
        assertEquals("abc", cvrs[2].tally_pool)
    }

    @Test
    fun test_cvr_has_contest() {
        val cvrs = CvrBuilders()
            .add(id = "1", tally_pool = "1", ContestVotes("AvB"), ContestVotes("CvD", "Candy"))
            .add(id = "2", tally_pool = "1", ContestVotes("CvD", Vote("Elvis", 1), Vote("Candy", 0)))
            .build()

        //val cvr_dict = [{ "id": 1, "votes": { "AvB": {}, "CvD": { "Candy":True } } },
        //    { "id": 2, "votes": { "CvD": { "Elvis":True, "Candy":false } } }]

        assertTrue(cvrs[0].has_contest("AvB"))
        assertTrue(cvrs[0].has_contest("CvD"))
        assertFalse(cvrs[0].has_contest("EvF"))

        assertFalse(cvrs[1].has_contest("AvB"))
        assertTrue(cvrs[1].has_contest("CvD"))
        assertFalse(cvrs[1].has_contest("EvF"))
    }

    @Test
    fun test_cvr_add_votes() {
        val cvrs = CvrBuilders()
            .add(id = "1", tally_pool = "1", ContestVotes("AvB"), ContestVotes("CvD", "Candy"))
            .add(id = "2", tally_pool = "1", ContestVotes("CvD", Vote("Elvis", 1), Vote("Candy", 0)))
            .build()

        //val cvr_dicts = [{ "id": 1, "votes": { "AvB": {}, "CvD": { "Candy":True } } },
        //    { "id": 2, "votes": { "CvD": { "Elvis":True, "Candy":false } } }]

        assertFalse(cvrs[0].has_contest("QvR"))
        assertFalse(cvrs[1].has_contest("AvB"))
        assertTrue(cvrs[0].update_votes(CvrBuilder("id").addContest("QvR").build())) // new contest added
        assertNotNull(cvrs[1].get_vote_for("CvD", "Elvis"))
        assertFalse(cvrs[0].update_votes(CvrBuilder("id").addVote("CvD", "Dan", 7).build()))
        assertTrue(cvrs[1].update_votes(CvrBuilder("id")
            .addContest("QvR")
            .addContestVotes(ContestVotes("CvD", Vote("Dan", 7), Vote("Elvis", 0), Vote("Candy", 1)))
            .build()))
        //        assert not cvr_list[0].has_contest('QvR')
        //        assert not cvr_list[1].has_contest('AvB')
        //        assert cvr_list[0].update_votes({'QvR': {}})
        //        assert cvr_list[1].get_vote_for('CvD', 'Elvis')
        //        assert not cvr_list[0].update_votes({'CvD': {'Dan':7}})
        //        assert cvr_list[1].update_votes({'QvR': {}, 'CvD': {'Dan':7, 'Elvis':False, 'Candy': True}})

        for (c in cvrs) {
            assertTrue(c.has_contest("QvR"))
        }
        assertEquals(0, cvrs[0].get_vote_for("QvR", "Dan"))
        assertEquals(7, cvrs[0].get_vote_for("CvD", "Dan"))
        assertNotNull(cvrs[0].get_vote_for("CvD", "Candy"))
        assertEquals(0, cvrs[0].get_vote_for("CvD", "Elvis"))
        assertEquals(7, cvrs[1].get_vote_for("CvD", "Dan"))
        assertNotNull(cvrs[1].get_vote_for("CvD", "Candy"))
        assertEquals(0, cvrs[1].get_vote_for("CvD", "Elvis"))

        //        for c in cvr_list:
        //            assert c.has_contest('QvR')
        //        assert not cvr_list[0].get_vote_for('QvR', 'Dan')
        //        assert cvr_list[0].get_vote_for('CvD', 'Dan') == 7
        //        assert cvr_list[0].get_vote_for('CvD', 'Candy')
        //        assert not cvr_list[0].get_vote_for('CvD', 'Elvis')
        //        assert cvr_list[1].get_vote_for('CvD', 'Dan') == 7
        //        assert cvr_list[1].get_vote_for('CvD', 'Candy')
        //        assert not cvr_list[1].get_vote_for('CvD', 'Elvis')
    }

    @Test
    fun test_cvr_pool_contests() {
        val cvrbs = CvrBuilders()
            .add(id = "1", sample_num = 1).setContestVotes("1", ContestVotes("AvB"), ContestVotes("CvD", "Candy"))
            .add(id = "2", p = 0.5).setContestVotes("2", ContestVotes("CvD", Vote("Elvis", 1), Vote("Candy", 0)), ContestVotes("EvF"))
            .add(id = "3", tally_pool = "abc", sampled = true).setContestVotes("3", ContestVotes("GvH"))

        //val cvr_dicts = [{ "id": 1, "sample_num": 1, "votes": { "AvB": {}, "CvD": { "Candy":True } } },
        //    { "id": 2, "p": 0.5, "votes": { "CvD": { "Elvis":True, "Candy":false }, "EvF": {} } },
        //    { "id": 3, "tally_pool": "abc", "sampled": True, "votes": { "GvH": {} } } ]

        assertEquals(setOf("AvB", "CvD", "EvF", "GvH"), cvrbs.poolContests())
        // assert CVR.pool_contests(cvr_list) == {'AvB', 'CvD', 'EvF', 'GvH'}
    }

    @Test
    fun test_add_pool_contests() {
        val cvrbs = CvrBuilders()
            .add(id = "1", tally_pool = "1", ContestVotes("AvB"), ContestVotes("CvD", "Candy"))
            .add(id = "2",  tally_pool = "1", ContestVotes("CvD", Vote("Elvis", 1), Vote("Candy", 0)), ContestVotes("EvF"))
            .add(id = "3",  tally_pool = "1", ContestVotes("GvH"))
            .add(id = "4",  tally_pool = "2", ContestVotes("AvB"), ContestVotes("CvD", "Candy"))
            .add(id = "5",  tally_pool = "2", ContestVotes("CvD", Vote("Elvis", 1), Vote("Candy", 0)), ContestVotes("EvF"))

        //val cvr_dicts = [{ "id": 1, "tally_pool": 1, "votes": { "AvB": {}, "CvD": { "Candy":True } } },
        //    { "id": 2, "tally_pool": 1, "votes": { "CvD": { "Elvis":True, "Candy":false }, "EvF": {} } },
        //    { "id": 3, "tally_pool": 1, "votes": { "GvH": {} } },
        //    { "id": 4, "tally_pool": 2, "votes": { "AvB": {}, "CvD": { "Candy":True } } },
        //    { "id": 5, "tally_pool": 2, "votes": { "CvD": { "Elvis":True, "Candy":false }, "EvF": {} } } ]

        val tally_pools = mutableMapOf<String, Set<String>>()
        for (p in cvrbs.poolSet()) {
            // tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
            tally_pools[p] = cvrbs.poolContests(p)
        }
        // assertEquals(CVR.add_pool_contests(cvr_list, tally_pools)
        assertTrue(cvrbs.add_pool_contests(tally_pools))

//        for i in range(3):
//            assert set(cvr_list[i].votes.keys()) == {'AvB', 'CvD', 'EvF', 'GvH'}
//        for i in range(3,5):
//            assert set(cvr_list[i].votes.keys()) == {'AvB', 'CvD', 'EvF'}
//        assert not CVR.add_pool_contests(cvr_list, tally_pools)
        val cvrs = cvrbs.build()
        for (i in 0 until 3) {
            assertEquals(setOf("AvB", "CvD", "EvF", "GvH"), cvrs[i].votes.keys)
        }
        for (i in 3 until 5) {
            assertEquals(setOf("AvB", "CvD", "EvF"), cvrs[i].votes.keys)
        }
        assertFalse(cvrbs.add_pool_contests(tally_pools))
    }

/*
    fun test_cvr_from_raire() {
        val raire_cvrs = [["1"],
            ["Contest", "339", "5", "15", "16", "17", "18", "45"],
            ["339", "99813_1_1", "17"],
            ["339", "99813_1_3", "16"],
            ["339", "99813_1_6", "18", "17", "15", "16"],
            ["3", "99813_1_6", "2"]
        ]
        val (c, n) = CVR.from_raire(raire_cvrs)
        assertEquals(3, c.size)
        assertEquals("99813_1_1", c[0].id)
        assertEquals(c[0].votes == { "339": { "17":1 } })
        assertEquals("99813_1_6", c[2].id)
        assertEquals(c[2].votes == { "339": { "18":1, "17":2, "15":3, "16":4 }, "3": { "2":1 } }) // merges votes?
    }

    fun test_make_phantoms() {
        val audit = Audit.from_dict({
            "strata": {
            "stratum_1": {
            "max_cards":   8,
            "use_style":   True,
            "replacement": false
        }
        }
        })
        val contests = Contest.from_dict_of_dicts({
            "city_council": {
            "risk_limit":0.05,
            "id": "city_council",
            "cards": None,
            "choice_function":"plurality",
            "n_winners":3,
            "candidates":["Doug", "Emily", "Frank", "Gail", "Harry"],
            "winner": ["Doug", "Emily", "Frank"]
        },
            "measure_1":   {
            "risk_limit":0.05,
            "id": "measure_1",
            "cards": 5,
            "choice_function":"supermajority",
            "share_to_win":2/3,
            "n_winners":1,
            "candidates":["yes", "no"],
            "winner": ["yes"]
        }
        })
        val cvrs =
            [CVR(id = "1", votes = { "city_council": { "Alice": 1 }, "measure_1": { "yes": 1 } }, phantom = false),
                CVR(id = "2", votes = { "city_council": { "Bob": 1 }, "measure_1": { "yes": 1 } }, phantom = false),
                CVR(id = "3", votes = { "city_council": { "Bob": 1 }, "measure_1": { "no": 1 } }, phantom = false),
                CVR(id = "4", votes = { "city_council": { "Charlie": 1 } }, phantom = false),
                CVR(id = "5", votes = { "city_council": { "Doug": 1 } }, phantom = false),
                CVR(id = "6", votes = { "measure_1": { "no": 1 } }, phantom = false)
            ]
        val prefix = "phantom-"

        val (cvr_list, phantoms) = CVR.make_phantoms(
            audit = audit,
            contests = contests,
            cvr_list = cvrs,
            prefix = "phantom-"
        )
        assertEquals(9, cvr_list.size)
        assertEquals(3, phantoms)
        assertEquals(5, contests["city_council"].cvrs)
        assertEquals(4, contests["measure_1"].cvrs)
        assertEquals(8, contests["city_council"].cards)
        assertEquals(5, contests["measure_1"].cards)
        assertEquals(
            np.sum(
                [c.has_contest("city_council")
                for c in cvr_list]) == 8, np.sum([c.has_contest("city_council") for c in cvr_list]))
        assertEquals(
            np.sum(
                [c.has_contest("measure_1")
                for c in cvr_list]) == 5, np.sum([c.has_contest("measure_1") for c in cvr_list]))
        assertEquals(5, np.sum([c.has_contest("city_council") and !c.phantom for c in cvr_list]))
        assertEquals(4, np.sum([c.has_contest("measure_1") and !c.phantom for c in cvr_list]))

        audit.strata["stratum_1"].use_style = false
        val (cvr_list, phantoms) = CVR.make_phantoms(audit, contests, cvrs, prefix = "phantom-")
        assertEquals(8, cvr_list.size)
        assertEquals(2, phantoms)
        assertEquals(5, contests["city_council"].cvrs)
        assertEquals(4, contests["measure_1"].cvrs)
        assertEquals(8, contests["city_council"].cards)
        assertEquals(8, contests["measure_1"].cards)
        assertEquals(5, np.sum([c.has_contest("city_council") for c in cvr_list]))
        // np.sum([c.has_contest("city_council") for c in cvr_list])
        assertEquals(4, np.sum([c.has_contest("measure_1") for c in cvr_list]))
        // np.sum([c.has_contest("measure_1") for c in cvr_list])
        assertEquals(5, np.sum([c.has_contest("city_council") and !c.phantom for c in cvr_list]))
        assertEquals(4, np.sum([c.has_contest("measure_1") and !c.phantom for c in cvr_list]))
    }

    fun test_assign_sample_nums() {
    val cvrbs = CvrBuilders()
        .add(id = "1",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)))
        .add(id = "2",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("yes", 1)))
        .add(id = "3",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("no", 1)))
        .add(id = "4",  ContestVotes("city_council", Vote("Charlie", 1)))
        .add(id = "5",  ContestVotes("city_council", Vote("Doug", 1)))
        .add(id = "6",  ContestVotes("measure_1", Vote("no", 1)))

        //val cvrs =
        //    [CVR(id = "1", votes = { "city_council": { "Alice": 1 }, "measure_1": { "yes": 1 } }, phantom = false),
        //        CVR(id = "2", votes = { "city_council": { "Bob": 1 }, "measure_1": { "yes": 1 } }, phantom = false),
        //        CVR(id = "3", votes = { "city_council": { "Bob": 1 }, "measure_1": { "no": 1 } }, phantom = false),
        //        CVR(id = "4", votes = { "city_council": { "Charlie": 1 } }, phantom = false),
        //        CVR(id = "5", votes = { "city_council": { "Doug": 1 } }, phantom = false),
        //        CVR(id = "6", votes = { "measure_1": { "no": 1 } }, phantom = false) ]

        val prng = SHA256(1234567890)
        CVR.assign_sample_nums(cvrs, prng)
        assertEquals(cvrs[0].sample_num, 100208482908198438057700745423243738999845662853049614266130533283921761365671)
        assertEquals(cvrs[5].sample_num, 93838330019164869717966768063938259297046489853954854934402443181124696542865)
    }

    fun test_consistent_sampling() {
        val cvrs = CvrBuilders()
            .add(id = "1",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "2",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "3",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("no", 1)))
            .add(id = "4",  ContestVotes("city_council", Vote("Charlie", 1)))
            .add(id = "5",  ContestVotes("city_council", Vote("Doug", 1)))
            .add(id = "6",  ContestVotes("measure_1", Vote("no", 1)))
            .build()

        cvrs.forEachIndexed { idx, it ->
            it.sample_num = idx.toDouble()
        }

//            val contests = {
//                "city_council": {
//                "risk_limit":0.05,
//                "id": "city_council",
//                "cards": None,
//                "choice_function":"plurality",
//                "n_winners":3,
//                "candidates":["Doug", "Emily", "Frank", "Gail", "Harry"],
//                "winner": ["Doug", "Emily", "Frank"],
//                "sample_size": 3
//            },
//                "measure_1":   {
//                "risk_limit":0.05,
//                "id": "measure_1",
//                "cards": 5,
//                "choice_function":"supermajority",
//                "share_to_win":2/3,
//                "n_winners":1,
//                "candidates":["yes", "no"],
//                "winner": ["yes"],
//                "sample_size": 4
//            }
        val con_tests = Contest.from_dict_of_dicts(contests)
        val sample_cvr_indices = CVR.consistent_sampling(cvrs, con_tests)
        assertEquals(sample_cvr_indices == [0, 1, 2, 5])
        assertEquals(2, con_tests["city_council"].sample_threshold)
        assertEquals(5, con_tests["measure_1"].sample_threshold)

    }

 */

@Test
fun test_tabulate_styles() {
    val cvrs = CvrBuilders()
        .add(id = "1",  ContestVotes("city_council", "Alice"), ContestVotes("measure_1", "yes"))
        .add(id = "2",  ContestVotes("city_council", "Bob"), ContestVotes("measure_1", "yes"))
        .add(id = "3",  ContestVotes("city_council", "Bob"), ContestVotes("measure_1", "no"))
        .add(id = "4",  ContestVotes("city_council", "Charlie"))
        .add(id = "5",  ContestVotes("city_council", "Doug"))
        .add(id = "6",  ContestVotes("measure_1", "no"))
        .add(id = "7",  ContestVotes("city_council", "Alice"), ContestVotes("measure_1", "yes"),
            ContestVotes("measure_2", "yes"))
        .add(id = "8",  ContestVotes("measure_1", "no"), ContestVotes("measure_2", "yes"))
        .add(id = "9",  ContestVotes("measure_1", "no"), ContestVotes("measure_3", "yes"))
        .build()

        //         cvrs = [CVR(id="1", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}}, phantom=False),
    //                CVR(id="2", votes={"city_council": {"Bob": 1}, "measure_1": {"yes": 1}}, phantom=False),
    //                CVR(id="3", votes={"city_council": {"Bob": 1}, "measure_1": {"no": 1}}, phantom=False),
    //                CVR(id="4", votes={"city_council": {"Charlie": 1}}, phantom=False),
    //                CVR(id="5", votes={"city_council": {"Doug": 1}}, phantom=False),
    //                CVR(id="6", votes={"measure_1": {"no": 1}}, phantom=False),
    //                CVR(id="7", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}, "measure_2": {"no":1}},
    //                          phantom=False),
    //                CVR(id="8", votes={"measure_1": {"no": 1}, "measure_2": {"yes": 1}}, phantom=False),
    //                CVR(id="9", votes={"measure_1": {"no": 1}, "measure_3": {"yes": 1}}, phantom=False),
    //            ]

    val t = Cvr.tabulate_styles(cvrs)
    println(t)
    assertEquals(6, t.size)
        assertEquals(1, t[setOf("measure_1", "measure_2")])
        assertEquals(1, t[setOf("measure_1", "measure_3")])
        assertEquals(2, t[setOf("city_council")])
        assertEquals(3, t[setOf("city_council", "measure_1")])
        assertEquals(1, t[setOf("measure_1")])
        assertEquals(1, t[setOf("city_council", "measure_1", "measure_2")])
    }

    @Test
    fun test_tabulate_votes() {
        val cvrs = CvrBuilders()
            .add(id = "1",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "2",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "3",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("no", 1)))
            .add(id = "4",  ContestVotes("city_council", Vote("Charlie", 1)))
            .add(id = "5",  ContestVotes("city_council", Vote("Doug", 1)))
            .add(id = "6",  ContestVotes("measure_1", Vote("no", 1)))
            .add(id = "7",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)),
                ContestVotes("measure_2", Vote("yes", 1)))
            .add(id = "8",  ContestVotes("measure_1", "no"), ContestVotes("measure_2", "yes"))
            .add(id = "9",  ContestVotes("measure_1", "no"), ContestVotes("measure_3", "yes"))
            .build()

        //         cvrs = [CVR(id="1", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}}, phantom=False),
        //                CVR(id="2", votes={"city_council": {"Bob": 1}, "measure_1": {"yes": 1}}, phantom=False),
        //                CVR(id="3", votes={"city_council": {"Bob": 1}, "measure_1": {"no": 1}}, phantom=False),
        //                CVR(id="4", votes={"city_council": {"Charlie": 1}}, phantom=False),
        //                CVR(id="5", votes={"city_council": {"Doug": 1}}, phantom=False),
        //                CVR(id="6", votes={"measure_1": {"no": 1}}, phantom=False),
        //                CVR(id="7", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}, "measure_2": {"no":1}},
        //                          phantom=False),
        //                CVR(id="8", votes={"measure_1": {"no": 1}, "measure_2": {"yes": 1}}, phantom=False),
        //                CVR(id="9", votes={"measure_1": {"no": 1}, "measure_3": {"yes": 1}}, phantom=False),
        //            ]
        //        d = CVR.tabulate_votes(cvrs)
        //        assert d['city_council']['Alice'] == 2
        //        assert d['city_council']['Bob'] == 2
        //        assert d['city_council']['Doug'] == 1
        //        assert d['measure_1']['no'] == 4

        val d = Cvr.tabulate_votes(cvrs)
        println(d)

        assertEquals(2, d["city_council"]!!["Alice"]!!)
        assertEquals(2, d["city_council"]!!["Bob"]!!)
        assertEquals(1, d["city_council"]!!["Charlie"]!!)
        assertEquals(1, d["city_council"]!!["Doug"]!!)
        assertEquals(3, d["measure_1"]!!["yes"]!!)
        assertEquals(4, d["measure_1"]!!["no"]!!)
        assertEquals(2, d["measure_2"]!!["yes"]!!)
        assertEquals(1, d["measure_3"]!!["yes"]!!)
    }

    @Test
    fun test_tabulate_cards() {
        val cvrs = CvrBuilders()
            .add(id = "1",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "2",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("yes", 1)))
            .add(id = "3",  ContestVotes("city_council", Vote("Bob", 1)), ContestVotes("measure_1", Vote("no", 1)))
            .add(id = "4",  ContestVotes("city_council", Vote("Charlie", 1)))
            .add(id = "5",  ContestVotes("city_council", Vote("Doug", 1)))
            .add(id = "6",  ContestVotes("measure_1", Vote("no", 1)))
            .add(id = "7",  ContestVotes("city_council", Vote("Alice", 1)), ContestVotes("measure_1", Vote("yes", 1)),
                ContestVotes("measure_2", Vote("yes", 1)))
            .add(id = "8",  ContestVotes("measure_1", "no"), ContestVotes("measure_2", "yes"))
            .add(id = "9",  ContestVotes("measure_1", "no"), ContestVotes("measure_3", "yes"))
            .build()

        //         cvrs = [CVR(id="1", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}}, phantom=False),
        //                CVR(id="2", votes={"city_council": {"Bob": 1}, "measure_1": {"yes": 1}}, phantom=False),
        //                CVR(id="3", votes={"city_council": {"Bob": 1}, "measure_1": {"no": 1}}, phantom=False),
        //                CVR(id="4", votes={"city_council": {"Charlie": 1}}, phantom=False),
        //                CVR(id="5", votes={"city_council": {"Doug": 1}}, phantom=False),
        //                CVR(id="6", votes={"measure_1": {"no": 1}}, phantom=False),
        //                CVR(id="7", votes={"city_council": {"Alice": 1}, "measure_1": {"yes": 1}, "measure_2": {"no":1}},
        //                          phantom=False),
        //                CVR(id="8", votes={"measure_1": {"no": 1}, "measure_2": {"yes": 1}}, phantom=False),
        //                CVR(id="9", votes={"measure_1": {"no": 1}, "measure_3": {"yes": 1}}, phantom=False),
        //            ]
        //        d = CVR.tabulate_votes(cvrs)
        //        assert d['city_council']['Alice'] == 2
        //        assert d['city_council']['Bob'] == 2
        //        assert d['city_council']['Doug'] == 1
        //        assert d['measure_1']['no'] == 4

        val d = Cvr.tabulate_cards_contests(cvrs)
        println(d)

        assertEquals(6, d["city_council"]!!)
        assertEquals(7, d["measure_1"]!!)
        assertEquals(2, d["measure_2"]!!)
        assertEquals(1, d["measure_3"]!!)
    }
}
