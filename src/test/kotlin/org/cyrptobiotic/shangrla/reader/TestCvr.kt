package org.cyrptobiotic.shangrla.reader

import org.cryptobiotic.shangrla.*
import org.cryptobiotic.shangrla.core.*
import org.cryptobiotic.shangrla.reader.sigfig
import org.cyrptobiotic.shangrla.AssertionUtils
import org.cyrptobiotic.shangrla.AssertionUtils.Companion.find_margins
import org.junit.jupiter.api.Test
import kotlin.random.Random

class TestCvr {

// def generate_fake_cvrs(contest_dict, style_dict):
//    fake_cvr_list = []
//    // loop through each style
//    for style in style_dict.keys():
//        // loop through the number of cards of that style
//        for i in range(style_dict[style]["cards"]):
//            // loop through the contests in that style and generate CVR
//            cvr = CVR(id = None, votes = {}, phantom=false, sample_num=None, p=None)
//            for contest in style_dict[style]["contests"]:
//                // randomly choose vote for that contest based on contest probabilities
//                cvr.set_votes({contest : {choice(contest_dict[contest]["candidates"],
//                                                        1, p = contest_dict[contest]["p"])[0] : true}})
//            // add cvr to list
//            fake_cvr_list.append(cvr)
//    // return the list of CVRs generated
//    return fake_cvr_list

    fun generate_fake_cvrs(contest_dict: Map<String, List<Candidate>>, style_dict: Map<String, Style>): List<TestCVR> {
        val fake_cvr_list = mutableListOf<TestCVR>()
        for ((name, style) in style_dict) {
            // loop through the number of cards of that style
            for (i in 0..style.cards-1) {
                // loop through the contests in that style and generate CVR
                val cvr = TestCVR()
                for (contestName in style.contests) {
                    val candidates: List<Candidate> = contest_dict[contestName]!!
                    val choice = choose(candidates)
                    // randomly choose vote for that contest based on contest probabilities
                    cvr.addVote(contestName, choice)
                }
                fake_cvr_list.add(cvr)
            }
        }
        // return the list of CVRs generated
        return fake_cvr_list
    }

    class TestCVR {
        val votes = mutableMapOf<String, MutableMap<String, Int>>() // Map(contestId, Map(selectionId, vote))

        fun addVote(contestName: String, candidateName: String) {
            val contest = votes.getOrPut(contestName) { mutableMapOf() }
            val vote: Int = contest.getOrPut(candidateName) { 0 }
            contest[candidateName] = vote + 1
        }

        fun hasContest(contestId: String): Boolean {
            return votes[contestId] != null
        }
    }

    val contest_dict = mapOf(
        "Contest 1" to makeCandidates(listOf("Candidate A", "Candidate B"), listOf(0.55, 0.45)),
        "Contest 2" to makeCandidates(listOf("Candidate A", "Candidate B"), listOf(0.7, 0.3)),
        "Contest 3" to makeCandidates(listOf("Candidate A", "Candidate B"), listOf(0.6, 0.4)),
        "Contest 4" to makeCandidates(listOf("Candidate A", "Candidate B"), listOf(0.2, 0.8)),
        "Contest 5" to makeCandidates(listOf("Candidate A", "Candidate B"), listOf(0.34, 0.66)),
    )

    val style_dict = mapOf(
        "style_1" to Style( listOf("Contest 1", "Contest 2"), 100),
        "style_2" to Style( listOf("Contest 3", "Contest 4", "Contest 5"), 200),
        "style_3" to Style( listOf("Contest 1", "Contest 2", "Contest 3", "Contest 4", "Contest 5"), 500),
    )

    fun makeCandidates(names: List<String>, ps: List<Double>): List<Candidate> {
        return names.zip(ps).map{ (name, p) -> Candidate(name, p) }
    }

    data class Candidate(val name: String, val p: Double)

    data class Style(val contests: List<String>, val cards: Int)

    fun choose(candidates: List<Candidate>): String {
        val random = Random.nextDouble(1.0)
        var accum = 0.0
        for (i in 0..candidates.size-1) {
            accum += candidates[i].p
            if (accum >= random) return candidates[i].name
        }
        throw IllegalStateException()
    }

    @Test
    fun testCounts() {
        val fake_cvrs = generate_fake_cvrs(contest_dict, style_dict)

        for (contest_name in contest_dict.keys) {
            println(contest_name)
            // Check vote counts
            val counts: MutableMap<String, Int> = mutableMapOf("Candidate A" to 0, "Candidate B" to 0)

            for (cvr in fake_cvrs) {
                if (cvr.hasContest(contest_name)) {
                    for ((candidate, vote) in cvr.votes[contest_name]!!) {
                        val accum: Int = counts[candidate]!!
                        counts[candidate] = accum + vote
                    }
                }
            }

            val cA = counts["Candidate A"]!!
            val cB = counts["Candidate B"]!!
            val total = cA + cB
            val percA = cA.toDouble() / total
            val percB = cB.toDouble() / total
            val candidates = contest_dict[contest_name]!!

            println("  ${counts} total=$total")
            println("  A ${percA.sigfig(2)} expect ${candidates[0].p.sigfig(2)}; B ${percB.sigfig(2)} expect ${candidates[1].p.sigfig(2)}")
        }
    }

    /*
    @Test
    fun testAudit() {
        val fake_cvrs = generate_fake_cvrs(contest_dict, style_dict)
        
        // set values
        val seed = 1234567890  // use, e.g., 20 rolls of a 10-sided die. Seed doesn"t have to be numeric
        val replacement = false
        val risk_function = "alpha_mart"

        // because comparison audit, may want to add f parameter to bias alpha towards u
        val risk_fn = { (x: DoubleArray, m, N) -> NonnegMean.alpha_mart(x, eta=(m+1)/2 , N=N, f=.1) }
        val g = 0.1
        val max_cards = 800
        val error_rate = 0.002

        // class Contest(
        //    val id: String,
        //    name: String,
        //    val risk_limit: Double,
        //    val cards: Int,
        //    val choice_function: SocialChoiceFunction,
        //    n_winners: Int,
        //    val share_to_win: Double,
        //    val candidates: List<Candidates>,
        //    winner: List<Any>,
        //    assertion_file: String,
        //    val audit_type: AuditType,
        //    //test: callable=None,
        //    g: Double,
        //    //estim: callable=None,
        //    //bet: callable=None,
        //    val use_style: Boolean,
        //    val assertions: Map<String, Assertion>, // TODO why is this a map?
        //    val tally: Map<Candidates, Int>,
        //    var sample_size: Int,
        //    val sample_threshold: Double,
        //)

        // Audit contest 2
        val contests: List<Contest> = listOf(
            Contest("Contest 1",
                cards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = "Candidate A",   // TODO no such field
            ),
            Contest("Contest 2",
                cards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = "Candidate A",
            ),
            Contest("Contest 1",
                cards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = "Candidate A",
            ),
        )

        // make assertions
        val all_assertions = Assertion.make_all_assertions(contests)

        val (cvr_list, phantom_vrs) = CVR.make_phantoms(max_cards, fake_cvrs, contests, use_style=true, prefix="phantom-1-")
        print("Created ${phantom_vrs} phantom records")

        // assign random sample nums including phantoms
        CVR.assign_sample_nums(cvr_list, prng=SHA256(32))

        // Find smallest margin
        val min_margin = find_margins(contests, cvr_list, use_style=true)
        println(min_margin)

        // Check audit parameters
        check_audit_parameters(risk_function, g, error_rate, contests)

        // find initial sample size
        val rf = lambda x,m,N: risk_fn(x,m,N)[1]   // p_history is the second returned value
        val ss_fn = lambda m, r, N: TestNonnegMean.initial_sample_size(risk_function=rf, N=N, margin=m, polling=false, error_rate=error_rate, alpha=r, reps=10) // change for comparison audits

        /*
        val (total_sample_size, sample_size_contests) =
            AssertionUtils.find_sample_size(contests, sample_size_function=ss_fn, use_style = true, cvr_list = cvr_list)

        println(sample_size_contests)
        println(total_sample_size)
        */

        val total_sample_size = AssertionUtils.find_sample_size(contests, sample_size_function=ss_fn)
        println(total_sample_size)

        // Created 0 phantom records
        //0.1200000000000001
        //{"Contest 1": 48, "Contest 2": 17, "Contest 4": 9}
        //50.99999999999999
    }

     */

    /*
    @Test
    fun testMakePhantoms() {
        val contests =  listOf(
            ContestSimple("city_council",
                risk_limit = 0.05,
                cards = 0,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 3,
                candidates = listOf("Doug","Emily","Frank","Gail","Harry"),
                reported_winners = listOf("Doug", "Emily", "Frank")
                ),
            ContestSimple("measure_1",
                risk_limit = 0.05,
                cards = 5,
                choice_function = SocialChoiceFunction.SUPERMAJORITY,
                share_to_win=2.0/3.0,
                n_winners = 1,
                candidates = listOf("yes","no"),
                reported_winners = listOf("yes")
            ),
        )


        val cvrs = listOf(
            makeCvrSimple("1", false, Vote("city_council", "Alice"), Vote("measure_1", "yes")),
            makeCvrSimple("2", false, Vote("city_council", "Bob"),   Vote("measure_1", "yes")),
            makeCvrSimple("3", false, Vote("city_council", "Bob"),   Vote("measure_1", "no")),
            makeCvrSimple("4", false, Vote("city_council", "Charlie")),
            makeCvrSimple("5", false, Vote("city_council", "Doug")),
            makeCvrSimple("6", false, Vote("measure_1", "no"))
        )
        val max_cards = 8
        val prefix = "phantom-"

        val (cvr_list, phantoms) = CvrSimple.make_phantoms(max_cards, cvrs, contests)

        /*
        assert len(cvr_list) == 9
        assert phantoms == 3
        assert contests["city_council"]["cvrs"] == 5
        assert contests["measure_1"]["cvrs"] == 4
        assert contests["city_council"]["cards"] == 8
        assert contests["measure_1"]["cards"] == 5
        assert np.sum([c.has_contest("city_council") for c in cvr_list]) == 8, \
        np.sum([c.has_contest("city_council") for c in cvr_list])
        assert np.sum([c.has_contest("measure_1") for c in cvr_list]) == 5, \
        np.sum([c.has_contest("measure_1") for c in cvr_list])
        assert np.sum([c.has_contest("city_council") and not c.is_phantom() for c in cvr_list]) ==  5
        assert np.sum([c.has_contest("measure_1") and not c.is_phantom() for c in cvr_list]) == 4

        cvr_list, phantoms = CVR.make_phantoms(max_cards, cvrs, contests, use_style=false, prefix='')
        assert len(cvr_list) == 8
        assert phantoms == 2
        assert contests["city_council"]["cvrs"] == 5
        assert contests["measure_1"]["cvrs"] == 4
        assert contests["city_council"]["cards"] == 8
        assert contests["measure_1"]["cards"] == 5
        assert np.sum([c.has_contest("city_council") for c in cvr_list]) == 5, \
        np.sum([c.has_contest("city_council") for c in cvr_list])
        assert np.sum([c.has_contest("measure_1") for c in cvr_list]) == 4, \
        np.sum([c.has_contest("measure_1") for c in cvr_list])
        assert np.sum([c.has_contest("city_council") and not c.is_phantom() for c in cvr_list]) ==  5
        assert np.sum([c.has_contest("measure_1") and not c.is_phantom() for c in cvr_list]) == 4

         */
    }

     */

}