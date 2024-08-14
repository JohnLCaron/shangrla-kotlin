package org.cryptobiotic.shangrla.core

import org.cryptobiotic.shangrla.reader.sigfig
import org.cryptobiotic.shangrla.core.AssertionUtils.Companion.find_margins
import org.junit.jupiter.api.Test
import kotlin.random.Random
import kotlin.test.assertEquals

class TestRaireContest {

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
        val fake_cvrs = generate_cvrs(contest_dict, style_dict)

        for (contest_name in contest_dict.keys) {
            println(contest_name)
            // Check vote counts
            val counts: MutableMap<String, Int> = mutableMapOf("Candidate A" to 0, "Candidate B" to 0)

            for (cvr in fake_cvrs) {
                if (cvr.has_contest(contest_name)) {
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

    @Test
    fun testAudit() {
        val cvrs = generate_cvrs(contest_dict, style_dict)
        
        // set values
        val seed = 1234567890  // use, e.g., 20 rolls of a 10-sided die. Seed doesn"t have to be numeric
        val replacement = false
        val risk_function = "alpha_mart"

        // because comparison audit, may want to add f parameter to bias alpha towards u
        // val risk_fn = { (x: DoubleArray, m, N) -> NonnegMean.alpha_mart(x, eta=(m+1)/2 , N=N, f=.1) }
        val g = 0.1
        val max_cards = 800
        val error_rate = 0.002

        // Audit contest 2
        val contests: List<Contest> = listOf(
            Contest("Contest 1",
                ncards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = listOf("Candidate A"),
            ),
            Contest("Contest 2",
                ncards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = listOf("Candidate A"),
            ),
            Contest("Contest 1",
                ncards = 600,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 1,
                candidates = listOf("Candidate A", "Candidate B"),
                reported_winners = listOf("Candidate A"),
            ),
        )

        // make assertions
        val all_assertions = Assertion.make_all_assertions(contests)

        val (cvr_list, phantom_vrs) = Cvr.make_phantoms(max_cards, cvrs, contests, use_style=true, prefix="phantom-1-")
        print("Created ${phantom_vrs} phantom records")

        // assign random sample nums including phantoms
        Cvr.assign_sample_nums(cvr_list)

        // Find smallest margin
        val min_margin = find_margins(contests, cvr_list, use_style=true)
        println(min_margin)

        // Check audit parameters
        // AssertionUtils.check_audit_parameters(risk_function, g, error_rate, contests)

        // find initial sample size
        //val rf = { (x,m,N) -> risk_fn(x,m,N)[1]   // p_history is the second returned value
        //val ss_fn = { (m, r, N) -> TestNonnegMean.initial_sample_size(risk_function=rf, N=N, margin=m, polling=false, error_rate=error_rate, alpha=r, reps=10) // change for comparison audits

        /*
        val (total_sample_size, sample_size_contests) =
            AssertionUtils.find_sample_size(contests, sample_size_function=ss_fn, use_style = true, cvr_list = cvr_list)

        println(sample_size_contests)
        println(total_sample_size)
        */

        //val total_sample_size = AssertionUtils.find_sample_size(contests, sample_size_function=ss_fn)
        //println(total_sample_size)

        // Created 0 phantom records
        //0.1200000000000001
        //{"Contest 1": 48, "Contest 2": 17, "Contest 4": 9}
        //50.99999999999999
    }

    @Test
    fun test_make_phantoms() {
        // def test_make_phantoms(self):
        //        audit = Audit.from_dict({'strata': {'stratum_1': {'max_cards':   8,
        //                                          'use_style':   True,
        //                                          'replacement': False
        //                                         }
        //                                      }})
        //        contests =  Contest.from_dict_of_dicts({'city_council': {'risk_limit':0.05,
        //                                     'id': 'city_council',
        //                                     'cards': None,
        //                                     'choice_function':'plurality',
        //                                     'n_winners':3,
        //                                     'candidates':['Doug','Emily','Frank','Gail','Harry'],
        //                                     'winner': ['Doug', 'Emily', 'Frank']
        //                                    },
        //                     'measure_1':   {'risk_limit':0.05,
        //                                     'id': 'measure_1',
        //                                     'cards': 5,
        //                                     'choice_function':'supermajority',
        //                                     'share_to_win':2/3,
        //                                     'n_winners':1,
        //                                     'candidates':['yes','no'],
        //                                     'winner': ['yes']
        //                                    }
        //                    })
        val contests =  listOf(
            ContestBuilder("city_council",
                risk_limit = 0.05,
                cards = null,
                choice_function = SocialChoiceFunction.PLURALITY,
                n_winners = 3,
                candidates = listOf("Doug","Emily","Frank","Gail","Harry"),
                reported_winners = listOf("Doug", "Emily", "Frank")
                ),
            ContestBuilder("measure_1",
                risk_limit = 0.05,
                cards = 5,
                choice_function = SocialChoiceFunction.SUPERMAJORITY,
                share_to_win=2.0/3.0,
                n_winners = 1,
                candidates = listOf("yes","no"),
                reported_winners = listOf("yes")
            ),
        )

        val cvrb = CvrBuilders()
            .add(id = "1", ContestVotes("city_council", "Alice"), ContestVotes("measure_1", "yes") )
            .add(id = "2", ContestVotes("city_council", "Bob"), ContestVotes("measure_1", "yes") )
            .add(id = "3", ContestVotes("city_council", "Bob"), ContestVotes("measure_1", "no") )
            .add(id = "4", ContestVotes("city_council", "Charlie"))
            .add(id = "5", ContestVotes("city_council", "Doug"))
            .add(id = "6", ContestVotes("measure_1", "no"))

        //        cvrs = [CVR(id="1", votes={"city_council": {"Alice": 1},     "measure_1": {"yes": 1}}, phantom=False),
        //                    CVR(id="2", votes={"city_council": {"Bob": 1},   "measure_1": {"yes": 1}}, phantom=False),
        //                    CVR(id="3", votes={"city_council": {"Bob": 1},   "measure_1": {"no": 1}}, phantom=False),
        //                    CVR(id="4", votes={"city_council": {"Charlie": 1}}, phantom=False),
        //                    CVR(id="5", votes={"city_council": {"Doug": 1}}, phantom=False),
        //                    CVR(id="6", votes={"measure_1": {"no": 1}}, phantom=False)
        //                ]

        val max_cards = 8
        val prefix = "phantom-"
        val (cvr_list, nphantoms) = make_phantoms(max_cards, cvrb.builders, contests)

        assertEquals(3, nphantoms)
        println(cvr_list)
        assertEquals(9, cvr_list.size)

        assertEquals(5, contests[0].ncvrs)
        assertEquals(4, contests[1].ncvrs)
        assertEquals(8, contests[0].cards)
        assertEquals(5, contests[1].cards)
        assertEquals(8, cvr_list.filter { it.has_contest("city_council") }.count())
        assertEquals(5, cvr_list.filter { it.has_contest("measure_1") }.count())
        assertEquals(5, cvr_list.filter { !it.phantom && it.has_contest("city_council") }.count())
        assertEquals(4, cvr_list.filter { !it.phantom && it.has_contest("measure_1") }.count())

        // TODO
        //         audit.strata['stratum_1'].use_style = False
        //        cvr_list, phantoms = CVR.make_phantoms(audit, contests, cvrs, prefix='phantom-')
        //        assert len(cvr_list) == 8
        //        assert phantoms == 2
        //        assert contests['city_council'].cvrs == 5
        //        assert contests['measure_1'].cvrs == 4
        //        assert contests['city_council'].cards == 8
        //        assert contests['measure_1'].cards == 8
        //        assert np.sum([c.has_contest('city_council') for c in cvr_list]) == 5, \
        //                       np.sum([c.has_contest('city_council') for c in cvr_list])
        //        assert np.sum([c.has_contest('measure_1') for c in cvr_list]) == 4, \
        //                       np.sum([c.has_contest('measure_1') for c in cvr_list])
        //        assert np.sum([c.has_contest('city_council') and not c.phantom for c in cvr_list]) ==  5
        //        assert np.sum([c.has_contest('measure_1') and not c.phantom for c in cvr_list]) == 4
    }

    fun generate_cvrs(contest_dict: Map<String, List<Candidate>>, style_dict: Map<String, Style>): List<Cvr> {
        val fake_cvr_list = mutableListOf<Cvr>()
        for ((name, style) in style_dict) {
            // loop through the number of cards of that style
            for (i in 0..style.cards-1) {
                // loop through the contests in that style and generate CVR
                val builder = CvrBuilder(name, false)
                for (contestName in style.contests) {
                    val candidates: List<Candidate> = contest_dict[contestName]!!
                    val choice = choose(candidates)
                    // randomly choose vote for that contest based on contest probabilities
                    builder.addVote(contestName, choice)
                }
                fake_cvr_list.add(builder.build())
            }
        }
        // return the list of CVRs generated
        return fake_cvr_list
    }

}