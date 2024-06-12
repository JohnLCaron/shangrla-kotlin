package org.cyrptobiotic.shangrla

import org.cryptobiotic.shangrla.*
import org.cryptobiotic.shangrla.core.*
import kotlin.test.Test

class TestAssertion {
    val plur_con_test : Contest
    val raw_AvB_asrtn : Assertion

    init {
        plur_con_test = Contest(
            id = "AvB",
            name = "AvB",
            cards = 4,
            choice_function = SocialChoiceFunction.PLURALITY,
            // candidates = 3,
            candidates = listOf("Alice","Bob","Candy"),
            reported_winners = listOf("Alice"),
            audit_type = AuditType.CARD_COMPARISON,
            test = NonnegMean.alpha_mart,
            estim = NonnegMean.optimal_comparison,
        )
        
        // assertion without a margin
        raw_AvB_asrtn = Assertion(
            contest = plur_con_test,
            winner = "Alice",
            loser = "Bob",
            test = NonnegMean(test=plur_con_test.test, estim=plur_con_test.estim, u=1.0f,
                N=plur_con_test.cards, t=.5f, random_order=true),
            assorter = Assorter(
                contest = plur_con_test,
                assort = { cvr: CVR ->
                    (CVR.as_vote(cvr.get_vote_for("AvB", "Alice")) -
                     CVR.as_vote(cvr.get_vote_for("AvB", "Bob")) + 1)/2 },
                upper_bound = 1
            )
        )
    }

    //     def test_set_tally_pool_means(self):
    //        cvr_dicts = [{"id": 1, "tally_pool": "1", "votes": {"AvB": {"Alice": 1}, "CvD": {"Candy":True}}},
    //                     {"id": 2, "tally_pool": "1", "votes": {"CvD": {"Elvis":True, "Candy":False}, "EvF": {}}},
    //                     {"id": 3, "tally_pool": "1", "votes": {"GvH": {}}},
    //                     {"id": 4, "tally_pool": "2", "votes": {"AvB": {"Bob": 1}, "CvD": {"Candy":True}}},
    //                     {"id": 5, "tally_pool": "2", "votes": {"CvD": {"Elvis":True, "Candy":False}, "EvF": {}}}
    //                   ]
    //        cvr_list = CVR.from_dict(cvr_dicts)
    //        pool_set = set(c.tally_pool for c in cvr_list)
    //        tally_pool = {}
    //        for p in pool_set:
    //            tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
    //        assert CVR.add_pool_contests(cvr_list, tally_pool)
    //        //
    //        // without use_style
    //        self.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=False)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1+1/2+1/2)/3)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0+1/2)/2)
    //        //
    //        // with use_style, but contests have already been added to every CVR in each pool
    //        self.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=True)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1+1/2+1/2)/3)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0+1/2)/2)
    //        //
    //        // with use_style, without adding contests to every CVR in each pool
    //        cvr_dicts = [{"id": 1, "tally_pool": "1", "votes": {"AvB": {"Alice": 1}, "CvD": {"Candy":True}}},
    //                     {"id": 2, "tally_pool": "1", "votes": {"CvD": {"Elvis":True, "Candy":False}, "EvF": {}}},
    //                     {"id": 3, "tally_pool": "1", "votes": {"GvH": {}}},
    //                     {"id": 4, "tally_pool": "2", "votes": {"AvB": {"Bob": 1}, "CvD": {"Candy":True}}},
    //                     {"id": 5, "tally_pool": "2", "votes": {"CvD": {"Elvis":True, "Candy":False}, "EvF": {}}}
    //                   ]
    //        cvr_list = CVR.from_dict(cvr_dicts)
    //        print(f"{list([str(c) for c in cvr_list])}")
    //        self.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list=cvr_list, tally_pool=tally_pool, use_style=True)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], 1)
    //        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], 0)
    @Test
    fun test_set_tally_pool_means() {
        var cvr_dicts = [{ "id": 1, "tally_pool": "1", "votes": { "AvB": { "Alice": 1 }, "CvD": { "Candy":True } } },
            { "id": 2, "tally_pool": "1", "votes": { "CvD": { "Elvis":True, "Candy":False }, "EvF": {} } },
            { "id": 3, "tally_pool": "1", "votes": { "GvH": {} } },
            { "id": 4, "tally_pool": "2", "votes": { "AvB": { "Bob": 1 }, "CvD": { "Candy":True } } },
            { "id": 5, "tally_pool": "2", "votes": { "CvD": { "Elvis":True, "Candy":False }, "EvF": {} } }
        ]
        val cvr_list = CVR.from_dict(cvr_dicts)
        val pool_set = set(c.tally_pool for c in cvr_list)
        val tally_pool = {}
        for (p in pool_set) {
            tally_pool[p] = CVR.pool_contests(list([c for c in cvr_list if c.tally_pool == p]))
        }
        assert CVR . add_pool_contests (cvr_list, tally_pool)
        //
        // without use_style
        self.raw_AvB_asrtn.assorter.set_tally_pool_means(
            cvr_list = cvr_list,
            tally_pool = tally_pool,
            use_style = False
        )
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1 + 1 / 2 + 1 / 2) / 3)
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0 + 1 / 2) / 2)
        //
        // with use_style, but contests have already been added to every CVR in each pool
        self.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list = cvr_list, tally_pool = tally_pool, use_style = True)
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], (1 + 1 / 2 + 1 / 2) / 3)
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], (0 + 1 / 2) / 2)
        //
        // with use_style, without adding contests to every CVR in each pool
        cvr_dicts = [{ "id": 1, "tally_pool": "1", "votes": { "AvB": { "Alice": 1 }, "CvD": { "Candy":True } } },
            { "id": 2, "tally_pool": "1", "votes": { "CvD": { "Elvis":True, "Candy":False }, "EvF": {} } },
            { "id": 3, "tally_pool": "1", "votes": { "GvH": {} } },
            { "id": 4, "tally_pool": "2", "votes": { "AvB": { "Bob": 1 }, "CvD": { "Candy":True } } },
            { "id": 5, "tally_pool": "2", "votes": { "CvD": { "Elvis":True, "Candy":False }, "EvF": {} } }
        ]
        cvr_list = CVR.from_dict(cvr_dicts)
        print(f"{list([str(c) for c in cvr_list])}")
        self.raw_AvB_asrtn.assorter.set_tally_pool_means(cvr_list = cvr_list, tally_pool = tally_pool, use_style = True)
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["1"], 1)
        np.testing.assert_almost_equal(self.raw_AvB_asrtn.assorter.tally_pool_means["2"], 0)
    }
}