package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.AuditType
import org.cryptobiotic.shangrla.core.Candidates
import org.cryptobiotic.shangrla.core.SocialChoiceFunction
import org.cryptobiotic.start.NonnegMean.Companion.makeAlphaMart

import java.security.SecureRandom
import kotlin.math.ceil
import kotlin.math.max
import kotlin.math.min


class Audit(
    val quantile: Double = 0.8,
    val error_rate_1: Double = 0.001, // rate of 1-vote overstatement errors
    val error_rate_2: Double = 0.0, // rate of 2-vote overstatement errors
    val reps: Int = 10,
    val auditType: AuditType = AuditType.CARD_COMPARISON,
    val use_styles: Boolean = true,
    ) {

    init {
        require(quantile >= 0.0 && quantile <= 1.0) { "quantile must be between 0 and 1" }
        require(error_rate_1 >= 0) { "expected rate of 1-vote errors must be nonnegative" }
        require(error_rate_2 >= 0) { "expected rate of 2-vote errors must be nonnegative" }
    }

    fun makeAssertions(contests: List<Contest>) {
        for (contest in contests) {
            val losers = contest.candidates.filter { !contest.winners.contains(it) }
            when (contest.choice_function) {
                SocialChoiceFunction.PLURALITY -> {
                    make_plurality_assertions(contest, contest.winners, losers)
                }
                SocialChoiceFunction.SUPERMAJORITY -> {
                    require (contest.winners.size == 1)
                    make_supermajority_assertion(contest, contest.winners[0])
                }
                else -> {
                    throw RuntimeException("${contest.choice_function} not supported")
                }
            }
        }
    }

    // highest candidate wins
    fun make_plurality_assertions(contest: Contest, winners: List<String>, losers: List<String>) {

        // _test = NonnegMean(test=test, estim=estim, bet=bet, g=contest.g, u=1, N=contest.cards, t=1/2, random_order=true)
        for (winr in winners) {
            for (losr in losers) {
                val assertion = Assertion(
                    contest,
                    winner = winr,
                    loser = losr,
                    assorter = PluralityAssorter(contest, winr, losr),
                    initialTest = makeAlphaMart(u = 1.0, N = contest.ncards, t = .5)
                )
                contest.addAssertion(assertion)
            }
        }

    }

    // highest candidate with share_to_win percentage. share_to_win >= 1/2
    fun make_supermajority_assertion(contest: Contest, winner: String) {
        /*
        Construct assertion that winner got >= share_to_win \in (0,1) of the valid votes

        **TO DO: This method assumes there was a winner. To audit that there was no winner requires
        flipping things.**

        An equivalent condition is:

           (votes for winner)/(2*share_to_win) + (invalid votes)/2 > 1/2.

        Thus the correctness of a super-majority outcome (where share_to_win >= 1/2) can be checked with a single assertion.

        Note: share_to_win < 1/2 might be useful for some social choice functions, including
        primaries where candidates who receive less than some threshold share are eliminated.

        A CVR with a mark for more than one candidate in the contest is considered an invalid vote.
         */

        val asn = Assertion(contest, winner = winner, loser = Candidates.ALL_OTHERS.name,
            assorter = SupermajorityAssorter(
                contest = contest,
                upperBound = 1 / (2 * contest.share_to_win),
                winner = winner,
                loser = Candidates.ALL_OTHERS.name,
            ),
            initialTest = makeAlphaMart(u = 1 / (2 * contest.share_to_win), N = contest.ncards, t = 0.5)
        )
        contest.addAssertion(asn)

    }

    // TODO do this before you make assertions ?? while you make assertions ??
    fun set_all_margins_from_cvrs(contests: List<Contest>, cvrs: List<Cvr>): Double {
        var min_margin = Double.POSITIVE_INFINITY
        for (con in contests) {
            for ((_,asn) in con.assertions) {
                val margin = asn.set_margin_from_cvrs( cvrs)
                min_margin = min(min_margin, margin)
            }
        }
        return min_margin
    }

    fun set_all_margins_from_tallies(contests: List<Contest>): Double {
        var min_margin = Double.POSITIVE_INFINITY
        for (con in contests) {
            for ((_,asn) in con.assertions) {
                val margin = asn.set_margin_from_tally()
                min_margin = min(min_margin, margin)
            }
        }
        return min_margin
    }

    // use_styles or mvr_sample, from previous?
    fun find_sample_size(
        contests: List<Contest>,
        cvrs: List<Cvr>,
        mvr_sample: List<Cvr> = emptyList(),
        cvr_sample: List<Cvr> = emptyList(),
    ): Int {
        // unless style information is being used, the sample size is the same for every contest.
        val old_sizes: MutableMap<String, Int> =
            contests.associate { it.id to 0 }.toMutableMap()

        for (contest in contests) {
            val contestId = contest.id
            old_sizes[contestId] = cvrs.filter { it.has_contest(contestId) }.map { if (it.sampled) 1 else 0 }.sum()

            var new_size = 0 // max sample size over all assertions
            for ((_, asn) in contest.assertions) {
                if (!asn.proved) {
                    if (mvr_sample.isNotEmpty()) { // use MVRs to estimate the next sample size. Set `prefix=True` to use data
                        val (data, _) = asn.mvrs_to_data(mvr_sample, cvr_sample)
                        new_size = max(
                            new_size,
                            asn.find_sample_size_with_data(
                                data = data, prefix = true,
                                reps = this.reps, quantile = this.quantile,
                            )
                        )
                    } else {
                        new_size = max(
                            new_size, asn.estimateSampleSize(
                                rate1 = this.error_rate_1,
                                rate2 = this.error_rate_2,
                                reps = this.reps, quantile = this.quantile,
                            )
                        )
                    }
                }
            }
            contest.sample_size = new_size
        }

        //        for c, con in contests.items():
        //            if stratum.use_style:
        //                old_sizes[c] = np.sum(
        //                    np.array([cvr.sampled for cvr in cvrs if cvr.has_contest(c)])
        //                )
        //            new_size = 0
        //            for a, asn in con.assertions.items():
        //                if not asn.proved:
        //                    if (
        //                        mvr_sample is not None
        //                    ):  # use MVRs to estimate the next sample size. Set `prefix=True` to use data
        //                        data, u = asn.mvrs_to_data(mvr_sample, cvr_sample)
        //                        new_size = max(
        //                            new_size,
        //                            asn.find_sample_size(
        //                                data=data,
        //                                prefix=True,
        //                                reps=self.reps,
        //                                quantile=self.quantile,
        //                                seed=self.sim_seed,
        //                            ),
        //                        )
        //                    else:
        //                        data = None
        //                        new_size = max(
        //                            new_size,
        //                            asn.find_sample_size(
        //                                data=data,
        //                                rate_1=self.error_rate_1,
        //                                rate_2=self.error_rate_2,
        //                                reps=self.reps,
        //                                quantile=self.quantile,
        //                                seed=self.sim_seed,
        //                            ),
        //                        )
        //            con.sample_size = new_size
        //        if stratum.use_style:
        //            for cvr in cvrs:
        //                if cvr.sampled:
        //                    cvr.p = 1
        //                else:
        //                    cvr.p = 0
        //                    for c, con in contests.items():
        //                        if cvr.has_contest(c) and not cvr.sampled:
        //                            cvr.p = max(
        //                                con.sample_size / (con.cards - old_sizes[c]), cvr.p
        //                            )
        //            total_size = math.ceil(np.sum([x.p for x in cvrs if not x.phantom]))
        //        else:
        //            total_size = np.max(
        //                np.array([con.sample_size for con in contests.values()])
        //            )
        //        return total_size

        for (cvr in cvrs) {
            if (cvr.sampled) {
                cvr.p = 1.0
            } else {
                cvr.p = 0.0
                for (con in contests) {
                    if (cvr.has_contest(con.id) && !cvr.sampled) {
                        val p1 = con.sample_size.toDouble() / (con.ncards - old_sizes[con.id]!!)
                        cvr.p = max(p1, cvr.p) // TODO nullability
                    }
                }
            }
        }
        // when old_sizes == 0, total_size should be con.sample_size (61); python has roundoff to get 62
        // total_size = ceil(np.sum([x.p for x in cvrs if !x.phantom))
        val summ: Double = cvrs.filter { !it.phantom }.map { it.p }.sum()
        val total_size = ceil(summ).toInt()
        return total_size
    }

    // just pick sampleSize unique samples. do fancy stuff later
    fun assign_sample_nums(cvrs: List<Cvr>, sampleSize: Int): Set<Cvr> {
        val result = mutableSetOf<Cvr>()
        while (result.size < sampleSize) {
            val sample = random.nextInt(cvrs.size)
            result.add(cvrs[sample])
        }
        return result
    }

    companion object {
        val random = SecureRandom.getInstanceStrong()
    }

}