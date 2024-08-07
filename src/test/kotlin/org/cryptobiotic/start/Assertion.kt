package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.Assertion.Companion.interleave_values
import org.cryptobiotic.shangrla.core.AuditType
import org.cryptobiotic.shangrla.core.Candidates
import org.cryptobiotic.shangrla.core.SocialChoiceFunction
import org.cryptobiotic.shangrla.core.numpy_arange
import kotlin.math.max

class Assertion(
    val contest: Contest,
    val winner: String,
    val loser: String,
    val assorter: Assorter,
    initialTest: NonnegMean,
) {

    var p_value: Double? = null
    var p_history: List<Double>? = null
    var proved: Boolean = false
    var margin: Double = Double.POSITIVE_INFINITY   // "reported assorter margin"
    var sample_size: Int? = null
    private var test: NonnegMean = initialTest

    fun name() = "$winner v $loser"

    override fun toString(): String {
        return "Assertion(contest=${contest.id}, upperBound=${assorter.upperBound}, winner='$winner', loser='$loser', " +
                "margin=$margin, p_value=$p_value, p_history=${p_history?.size}, proved=$proved, sample_size=$sample_size)"
    }

    fun testFn(x: DoubleArray) = test.testFn.test(x)

    // set margin, test.mean
    fun set_margin_from_cvrs(cvrs: List<Cvr>): Double {
        val amean = this.assorter.mean(cvrs)
        if (amean < .5) {
            println("assertion $this not satisfied by CVRs: mean value is ${amean}")
        }

        // Define v ≡ 2Āc − 1, the reported assorter margin. In a two-candidate plurality contest, v
        // is the fraction of ballot cards with valid votes for the reported winner, minus the fraction
        // with valid votes for the reported loser. This is the diluted margin of [22,12]. (Margins are
        // traditionally calculated as the difference in votes divided by the number of valid votes.
        // Diluted refers to the fact that the denominator is the number of ballot cards, which is
        // greater than or equal to the number of valid votes.)
        this.margin = 2 * amean - 1 // eq 3

        val u = when (this.contest.audit_type) {
            AuditType.POLLING -> this.assorter.upperBound
            AuditType.CARD_COMPARISON, AuditType.ONEAUDIT -> 2 / (2 - this.margin / this.assorter.upperBound)
            else -> throw NotImplementedError("audit type ${this.contest.audit_type} not supported")
        }

        // Tests of the hypothesis that the mean of a population of values in [0, u] is less than or equal to t
        // so u must be the max value of the population.
        // the math seems to be on page 10, if you take tau = 2.
        this.test = this.test.changeMean(newu=u)
        return this.margin
    }

    fun set_margin_from_tally(tallyInput: Map<String, Int>? = null): Double { // cand -> count.
        /*
        find the assorter margin between implied by a tally.
        Generally useful only for approval, plurality, and supermajority contests.
        Assumes the number of cards containing the contest has been set.

        Parameters
        ----------
        tallyInput: dict of tallies for the candidates in the contest. Keys are candidates as listed
            in Contest.candidates. If `tally is None` tries to use the contest.tally.

        The margin for a supermajority contest with a winner is (see SHANRGLA section 2.3)
            2(pq/(2f) + (1 − q)/2 - 1/2) = q(p/f-1)
         where:
            q is the fraction of cards that have valid votes
            p is the fraction of cards that have votes for the winner
            f is the fraction of valid votes required to win.

        Side effects
        ------------
        sets this.margin

        */
        //        tally = tally if tally else self.contest.tally
        //        if self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY \
        //             or self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.APPROVAL:
        //            self.margin = (tally[self.winner]-tally[self.loser])/self.contest.cards
        //        elif self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY:
        //            if self.winner == Candidates.NO_CANDIDATE or self.loser != Candidates.ALL_OTHERS:
        //                raise NotImplementedError(f'TO DO: currently only support super-majority with a winner')
        //            else:
        //                q = np.sum([tally[c] for c in self.contest.candidates])/self.contest.cards
        //                p = tally[self.winner]/self.contest.cards
        //                self.margin = q*(p/self.contest.share_to_win - 1)
        //        else:
        //            raise NotImplementedError(f'social choice function {self.contest.choice_function} not supported')

        val tally = tallyInput ?: this.contest.tally
        if (this.contest.choice_function in listOf(SocialChoiceFunction.PLURALITY, SocialChoiceFunction.APPROVAL)) {
            this.margin = (tally[this.winner]!! - tally[this.loser]!!).toDouble() / this.contest.ncards // // TODO check nullable
        } else if (this.contest.choice_function == SocialChoiceFunction.SUPERMAJORITY) {
            if (this.winner == Candidates.NO_CANDIDATE.name || this.loser != Candidates.ALL_OTHERS.name) {
                throw NotImplementedError("TO DO: currently only support super-majority with a winner")
            } else {
                // val q = np.sum([tally[c] for c in this.contest.candidates])/this.contest.cards
                val q = this.contest.candidates.map { tally[it]!! }.sum() / this.contest.ncards // LOOK check nullable
                val p = tally[this.winner]!! / this.contest.ncards // LOOK check nullable
                this.margin = q * (p / this.contest.share_to_win - 1)
            }
        } else {
            throw NotImplementedError("social choice function ${this.contest.choice_function} not supported")
        }
        return this.margin
    }

    /*
    Process mvrs (and, for comparison audits, cvrs) to create data for the assertion"s test
    and for sample size simulations.

    Creates assorter values for the mvrs, or overstatement assorter values using the mvrs and cvrs,
    according to whether the audit uses ballot polling or card-level comparison

    The Assertion.margin must be set before calling this function.

    Parameters
    ----------
    mvr_sample: corresponding MVRs
    cvr_sample: sampled CVRs (for comparison audits)

    mvr_sample and cvr_sample should be ordered using CVR.prep_comparison_sample() or CVR.prep_polling_sample()
    before calling this routine

    Returns
    -------
    d: DoubleArray; either assorter values or overstatement assorter values, depending on the audit method
    u: upper bound for the test
*/
    fun mvrs_to_data(mvr_sample: List<Cvr>, cvr_sample: List<Cvr>): Pair<DoubleArray, Double> {
        require(this.margin != Double.POSITIVE_INFINITY) { "Margin is not set"  }
        require(this.margin > 0.0) { "Margin ${this.margin} is nonpositive" }
        val upper_bound = this.assorter.upperBound
        val con = this.contest

        var u: Double
        require(mvr_sample.size == cvr_sample.size)
        val cvr2: List<Pair<Cvr, Cvr>> = mvr_sample.zip(cvr_sample)
        val d = cvr2.filter { (_, cvr) -> (cvr.has_contest(con.id) && cvr.sample_num <= con.sample_threshold!!) }
            .map { (mvr, cvr) -> this.overstatement_assorter(mvr, cvr) }
        u = 2 / (2 - this.margin / upper_bound)

        // convert to double array
        val fa = DoubleArray(d.size) { d[it] }
        return Pair(fa, u)
    }

    /*
        Estimate sample size needed to reject the null hypothesis that the assorter mean is <=1/2,
        for the specified risk function, given:
        - for comparison audits, the assorter margin and assumptions about the rate of overstatement errors
        - for polling audits, either a set of assorter values, or the assumption that the reported tallies are correct

        If `data is not None`, uses data to make the estimate. There are three strategies:
            1. if `reps is None`, tile the data to make a list of length N
            2. if `reps is not None and not prefix`, sample from the data with replacement to make `reps` lists of length N
            3. if `reps is not None and prefix`, start with `data`, then draw N-len(data) times from data with
               replacement to make `reps` lists of length N

        If `data is None`, constructs values from scratch.
            - For polling audits, values are inferred from the reported tallies. Since contest.tally only reports
                actual candidate totals, not IRV/RAIRE pseudo-candidates, this is not implemented for IRV.
            - For comparison audits, there are two strategies to construct the values:
                1. Systematically interleave small and large values, starting with a small value (`reps is None`)
                2. Sample randomly from a set of such values
                The rate of small values is `rate_1` if `rate_1 is not None`. If `rate is None`, for POLLING audits, gets
                the rate of small values from the margin.
                For Audit.AUDIT_TYPE.POLLING audits, the small values are 0 and the large values are `u`; the rest are 1/2.
                For Audit.AUDIT_TYPE.CARD_COMPARISON audits, the small values are the overstatement assorter for an
                overstatement of `u/2` and the large values are the overstatement assorter for an overstatement of 0.

        **Assumes that this.test.u has been set appropriately for the audit type (polling or comparison).**
        **Thus, for comparison audits, the assorter margin should be be set before calling this function.**

        Parameters
        ----------
        data: DoubleArray; observations on which to base the calculation. If `data is not None`, uses them in a bootstrap
            approach, rather than simulating errors.
            If `this.contest.audit_type==Audit.POLLING`, the data should be (simulated or actual) values of
            the raw assorter.
            If `this.contest.audit_type==Audit.CARD_COMPARISON`, the data should be (simulated or actual)
            values of the overstatement assorter.
        prefix: bool; prefix the data, then sample or tile to produce the remaining values
        rate_1: float; assumed rate of "small" values for simulations (1-vote overstatements). Ignored if `data is not None`
            If `rate_1 is None and this.contest.audit_type==Audit.POLLING` the rate of small values is inferred
            from the margin
        rate_2: float; assumed rate of 0s for simulations (2-vote overstatements).
        reps: int; if `reps is None`, builds the data systematically
            if `reps is not None`, performs `reps` simulations to estimate the `quantile` quantile of sample size.
        quantile: float; if `reps is not None`, quantile of the distribution of sample sizes to return
            if `reps is None`, ignored
        seed: int; if `reps is not None`, use `seed` as the seed in numpy.random to estimate the quantile

        Returns
        -------
        sample_size: estimated to be sufficient to confirm the outcome if data are generated according to the assumptions

        Side effects
        ------------
        sets the sample_size attribute of the assertion
        */

    // break these two cases apart `data is not None` vs `data is None`
    fun find_sample_size_with_data(
        data: DoubleArray,
        prefix: Boolean = false,
        reps: Int,
        quantile: Double = 0.5, // quantile of the distribution of sample sizes to return
        seed: Int = 1234567890 // seed in numpy.random to estimate the quantile
    ): Int {
        return this.test.estimateSampleSize(
            x = data, alpha = this.contest.risk_limit, reps = reps,
            prefix = prefix, quantile = quantile, // seed = seed
        )
    }

    fun estimateSampleSize(
        prefix: Boolean = false,
        rate1: Double? = null,
        rate2: Double? = null,
        reps: Int,
        quantile: Double, // quantile of the distribution of sample sizes to return
    ): Int {
        require(this.margin != Double.POSITIVE_INFINITY) { "Margin is not set"  }
        require(this.margin > 0.0) { "Margin ${this.margin} is nonpositive" }

        return if (this.contest.audit_type == AuditType.POLLING) estimateSampleSizePolling(prefix, reps, quantile)
        else estimateSampleSizeComparision(prefix, rate1, rate2, reps, quantile)
    }

    fun estimateSampleSizePolling(
        prefix: Boolean = false,
        reps: Int,
        quantile: Double = 0.5, // quantile of the distribution of sample sizes to return
    ): Int {

        val big = this.assorter.upperBound

        val n_0 = this.contest.tally[this.loser]!! // LOOK nullability
        val n_big = this.contest.tally[this.winner]!! // LOOK nullability
        val n_half = this.test.N - n_0 - n_big
        //         fun interleave_values(n_small: Int, n_med: Int, n_big: Int, small: Double, med: Double, big: Double): DoubleArray {
        val x = interleave_values(n_0, n_half, n_big, big = big)

        val sample_size = this.test.estimateSampleSize(
            x = x,
            alpha = this.contest.risk_limit,
            reps = reps,
            prefix = prefix,
            quantile = quantile, // seed = seed
        )
        this.sample_size = sample_size
        return sample_size
    }

    fun estimateSampleSizeComparision(
        prefix: Boolean = false,
        rate1: Double? = null,
        rate2: Double? = null,
        reps: Int,
        quantile: Double = 0.5, // quantile of the distribution of sample sizes to return
    ): Int {
        val big = this.make_overstatement(overs = 0.0)
        val small = this.make_overstatement(overs = 0.5)
        println("  Assertion ${name()} big = ${big} small = ${small}")

        val rate_1 = rate1 ?: ((1 - this.margin) / 2)   // rate of small values
        val x = DoubleArray(this.test.N) { big } // array N floats, all equal to big

        val rate_1_i = numpy_arange(0, this.test.N, step = (1.0 / rate_1).toInt())
        val rate_2_i = if (rate2 != null  && rate2 != 0.0) numpy_arange(0, this.test.N, step = (1.0 / rate2).toInt()) else IntArray(0)
        rate_1_i.forEach { x[it] = small }
        rate_2_i.forEach { x[it] = 0.0 }

        val sample_size = this.test.estimateSampleSize(
            x = x,
            alpha = this.contest.risk_limit,
            reps = reps,
            prefix = prefix,
            quantile = quantile, // seed = seed
        )
        this.sample_size = sample_size
        return sample_size
    }


    fun make_overstatement(overs: Double): Double {
        // TODO reference in paper
        val result =  (1 - overs / this.assorter.upperBound) / (2 - this.margin / this.assorter.upperBound)
        return result
    }

    //  assorter that corresponds to normalized overstatement error for an assertion
    fun overstatement_assorter(mvr: Cvr, cvr: Cvr, use_style: Boolean = true): Double {
        return (1 - this.assorter.overstatement(mvr, cvr, use_style) / this.assorter.upperBound) /
                (2 - this.margin / this.assorter.upperBound)
    }


    companion object {

        fun set_p_values(contests: List<Contest>, mvr_sample: List<Cvr>, cvr_sample: List<Cvr>): Double {
            require(mvr_sample.size == cvr_sample.size) {"unequal numbers of cvrs and mvrs"}

            var p_max = 0.0
            for (con in contests) {
                for ((_, asn) in con.assertions) {
                    val (d, u) = asn.mvrs_to_data(mvr_sample, cvr_sample)
                    asn.test = asn.test.changeMean(newu=u)
                    // set upper bound for the test for each assorter
                    val (p_value, p_history) = asn.testFn(d)
                    asn.p_value = p_value
                    asn.p_history = p_history.toList()
                    asn.proved = (p_value <= con.risk_limit) || asn.proved
                    p_max = max(p_max, p_value)
                }
            }
            return p_max
        }
    }

}