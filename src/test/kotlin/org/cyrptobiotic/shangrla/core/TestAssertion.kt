package org.cyrptobiotic.shangrla.core

import org.cryptobiotic.shangrla.core.*
import kotlin.test.Test
import kotlin.test.assertEquals

class TestAssertion {
    val plur_cvr_list: List<Cvr>
    val plur_con_test: Contest
    val raw_AvB_asrtn: Assertion

    init {
        // objects used to test many Assertion functions in plurality contests
        plur_cvr_list = CvrBuilders()
            .add(id = "1_1", tally_pool = "1", ContestVotes("AvB", "Alice", 1))
            .add(id = "1_2", tally_pool = "1", ContestVotes("AvB", "Bob", true))
            .add(id = "1_3", tally_pool = "2", ContestVotes("AvB", "Alice", true))
            .add(id = "1_4", tally_pool = "2", ContestVotes("AvB", "Alice", true))
            .build()

        plur_con_test = Contest(
            id = "AvB",
            name = "AvB",
            ncards = 4,
            choice_function = SocialChoiceFunction.PLURALITY,
            // candidates = 3,
            candidates = listOf("Alice", "Bob", "Candy"),
            reported_winners = listOf("Alice"),
            audit_type = AuditType.CARD_COMPARISON,
            estimFn = { (x) -> doubleArrayOf(NonnegMean().optimal_comparison(0.2)) },
        )

        // assertion without a margin
        raw_AvB_asrtn = Assertion(
            contest = plur_con_test,
            winner = "Alice",
            loser = "Bob",
            test = NonnegMean(
                estimFnOverride = plur_con_test.estimFn,
                u = 1.0, N = plur_con_test.ncards, t = .5
            ),
            assorter = Assorter(
                contest = plur_con_test,
                assort = { cvr ->
                    (Cvr.as_vote(cvr.get_vote_for("AvB", "Alice")) -
                            Cvr.as_vote(cvr.get_vote_for("AvB", "Bob")) + 1) * 0.5
                },
                upper_bound = 1.0,
            )
        )
    }

    //    def test_min_p(self):
    //        asrtn1 = Assertion(p_value = 0.1, p_history = [1, 0.5, 0.1, 0.01, 0.1])
    //        asrtn2 = Assertion(p_value = 0.05, p_history = [0.05])
    //        assert Assertion.min_p(asrtn1) == 0.01
    //        assert Assertion.min_p(asrtn2) == 0.05
    @Test
    fun test_min_p() {
        val asrtn1 = raw_AvB_asrtn.copy(p_value = 0.1, p_history = listOf(1.0, 0.5, 0.1, 0.01, 0.1))
        val asrtn2 = raw_AvB_asrtn.copy(p_value = 0.05, p_history = listOf(0.05))
        assertEquals(0.01, asrtn1.min_p() )
        assertEquals(0.05, asrtn2.min_p() )
    }

    @Test
    fun test_margin() {
        assertEquals(0.5, raw_AvB_asrtn.margin( plur_cvr_list, true) )
    }

    //     def test_set_tally_pool_means(this) in test_Assertions.py
    //        cvr_dicts = [{"id": 1, "tally_pool": "1", "votes": {"AvB": {"Alice": 1}, "CvD": {"Candy":true}}},
    //                     {"id": 2, "tally_pool": "1", "votes": {"CvD": {"Elvis":true, "Candy":false}, "EvF": {}}},
    //                     {"id": 3, "tally_pool": "1", "votes": {"GvH": {}}},
    //                     {"id": 4, "tally_pool": "2", "votes": {"AvB": {"Bob": 1}, "CvD": {"Candy":true}}},
    //                     {"id": 5, "tally_pool": "2", "votes": {"CvD": {"Elvis":true, "Candy":false}, "EvF": {}}}
    //                   ]
    //        cvr_list = CVR.from_dict(cvr_dicts)
    //        pool_set = set(c.tally_pool for c in cvr_list)
    //        tally_pool = {}
    //        for p in pool_set:
    //            tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
    //        assert CVR.add_pool_contests(cvr_list, tally_pool)
    //        //
    //        // without use_style
    //        this.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=false)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1+1/2+1/2)/3)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0+1/2)/2)
    //        //
    //        // with use_style, but contests have already been added to every CVR in each pool
    //        this.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=true)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1+1/2+1/2)/3)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0+1/2)/2)
    //        //
    //        // with use_style, without adding contests to every CVR in each pool
    //        cvr_dicts = [{"id": 1, "tally_pool": "1", "votes": {"AvB": {"Alice": 1}, "CvD": {"Candy":true}}},
    //                     {"id": 2, "tally_pool": "1", "votes": {"CvD": {"Elvis":true, "Candy":false}, "EvF": {}}},
    //                     {"id": 3, "tally_pool": "1", "votes": {"GvH": {}}},
    //                     {"id": 4, "tally_pool": "2", "votes": {"AvB": {"Bob": 1}, "CvD": {"Candy":true}}},
    //                     {"id": 5, "tally_pool": "2", "votes": {"CvD": {"Elvis":true, "Candy":false}, "EvF": {}}}
    //                   ]
    //        cvr_list = CVR.from_dict(cvr_dicts)
    //        print(f"{list([str(c) for c in cvr_list])}")
    //        this.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=true)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["1"], 1)
    //        np.testing.assert_almost_equal(this.raw_AvB_asrtn.assorter.tally_pool_means["2"], 0)

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

        val tally_pool = mutableMapOf<String, Set<String>>()
        for (p in cvrBuilders.poolSet()) {
            // tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
            tally_pool[p] = cvrBuilders.poolContests(p)
        }

        cvrBuilders.add_pool_contests(tally_pool)
        val cvr_list = cvrBuilders.build()

        // without use_style
        // this.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=false)
        this.raw_AvB_asrtn.assorter.set_tally_pool_means(
            cvr_list = cvr_list,
            tally_pool = cvrBuilders.poolSet(),
            use_style = false
        )
        assertEquals((1.0 + (1.0 / 2) + (1.0 / 2)) / 3.0, this.raw_AvB_asrtn.assorter.tally_pool_means!!["1"]!!)
        assertEquals((0.0 + (1.0 / 2)) / 2.0, this.raw_AvB_asrtn.assorter.tally_pool_means!!["2"]!!)

        // with use_style, but contests have already been added to every CVR in each pool
        this.raw_AvB_asrtn.assorter.set_tally_pool_means(
            cvr_list = cvr_list,
            tally_pool = cvrBuilders.poolSet(),
            use_style = true
        )
        assertEquals((1.0 + 1.0 / 2 + 1.0 / 2) / 3, this.raw_AvB_asrtn.assorter.tally_pool_means!!["1"]!!)
        assertEquals((0.0 + 1.0 / 2) / 2, this.raw_AvB_asrtn.assorter.tally_pool_means!!["2"]!!)
        println("ok")
    }

    @Test
    fun test_set_tally_pool_means2() {
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

        val cvr_list = cvrBuilders.build()

        // with use_style, without adding contests to every CVR in each pool
        this.raw_AvB_asrtn.assorter.set_tally_pool_means(
            cvr_list = cvr_list,
            tally_pool = cvrBuilders.poolSet(),
            use_style = true
        )

        assertEquals(1.0, this.raw_AvB_asrtn.assorter.tally_pool_means!!["1"]!!)
        assertEquals(0.0, this.raw_AvB_asrtn.assorter.tally_pool_means!!["2"]!!)
        println("ok")
    }
}