package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.*
import org.cryptobiotic.shangrla.core.numpy_arange
import kotlin.math.min

data class AssertionSimple(
    val contest: ContestSimple,
    val assorter: AssorterFn,
    val winner: String,
    val loser: String,
    var test: NonnegMean,

    var margin: Double? = null,
    val p_value: Double? = null,
    val p_history: List<Double>? = null,
    var proved: Boolean = false,
    var sample_size: Int? = null,
    // var tally_pool_means: MutableMap<String, Double>? = null,
) {

    fun find_sample_size(
        data: DoubleArray? = null,
        prefix: Boolean = false,
        rate1: Double? = null,
        rate2: Double? = null,
        reps: Int? = null,
        quantile: Double = 0.5,
        seed: Int = 1234567890
    ): Int {

        if (data != null) {
            return this.test.sample_size(
                data, alpha = this.contest.risk_limit, reps = reps,
                prefix = prefix, quantile = quantile, // seed = seed TODO
            )
        }

        require(margin != null) { "Margin as not set"  }
        val amargin = margin!!
        require(amargin > 0.0) { "Margin ${amargin} is nonpositive" }

        val big = if (this.contest.audit_type == AuditType.POLLING) this.assorter.upper_bound
        else this.make_overstatement(overs = 0.0)
        val small = if (this.contest.audit_type == AuditType.POLLING) 0.0
                    else this.make_overstatement(overs = 0.5)
        println("  big = $big small = $small")

        val rate_1 = rate1 ?: ((1 - amargin) / 2)   // rate of small values
        var x = DoubleArray(this.test.N) { big } // array N floats, all equal to big
        // println("  x = ${x.contentToString()}")


     /*   if (this.contest.audit_type == AuditType.POLLING) {
            if (this.contest.choice_function == SocialChoiceFunction.IRV) {
                throw NotImplementedError("data must be provided to estimate sample sizes for IRV assertions")
            } else { // get tally
                if (this.contest.tally != null) {
                    val n_0 = this.contest.tally[this.loser]!! // LOOK nulllability
                    val n_big = this.contest.tally[this.winner]!! // LOOK nulllability
                    val n_half = this.test.N - n_0 - n_big
                    //         fun interleave_values(n_small: Int, n_med: Int, n_big: Int, small: Double, med: Double, big: Double): DoubleArray {
                    x = interleave_values(n_0, n_half, n_big, big = big)
                } else {
                    throw Exception("contest ${this.contest} tally required but not defined") // LOOK ValueError
                }
            }

        } else */
        if (this.contest.audit_type == AuditType.CARD_COMPARISON) { // # comparison audit
            //     arange([start,] stop[, step,], dtype=None, *, like=None) : Return evenly spaced values within a given interval.
            val rate_1_i = numpy_arange(0, this.test.N, step = (1.0 / rate_1).toInt())
            val rate_2_i = if (rate2 != null  && rate2 != 0.0) numpy_arange(0, this.test.N, step = (1.0 / rate2).toInt()) else IntArray(0)
            rate_1_i.forEach { x[it] = small }
            rate_2_i.forEach { x[it] = 0.0 }
        } else {
            throw NotImplementedError("audit type ${this.contest.audit_type} for contest ${this.contest} not implemented")
        }

        val sample_size = this.test.sample_size(
            x, alpha = this.contest.risk_limit, reps = reps,
            prefix = prefix, quantile = quantile, // seed = seed
        )
        this.sample_size = sample_size
        return sample_size
    }

    fun mvrs_to_data(mvr_sample: List<CvrSimple>, cvr_sample: List<CvrSimple>): Pair<DoubleArray, Double> {
        requireNotNull(this.margin)
        val margin = this.margin!!
        val upper_bound = this.assorter.upper_bound
        val con = this.contest
        val use_style = true // con.use_style

        var d: List<Double>
        var u: Double
        if (con.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            require(mvr_sample.size == cvr_sample.size)
            val cvr2: List<Pair<CvrSimple, CvrSimple>> = mvr_sample.zip(cvr_sample)
            d = cvr2.filter { (mvr, cvr) -> !use_style || (cvr.has_contest(con.id) && cvr.sample_num!! <= con.sample_threshold!!) }
                .map { (mvr, cvr) -> this.overstatement_assorter(mvr, cvr, use_style) }
            u = 2 / (2 - margin / upper_bound)
        } else if (con.audit_type == AuditType.POLLING) {  // Assume style information is irrelevant
            // d = np.array([this.assorter.assort(mvr_sample[i]) for i in range(len(mvr_sample))])
            d = mvr_sample.map { this.assorter.assort(it) }
            u = upper_bound
        } else {
            throw NotImplementedError("audit type ${con.audit_type} not implemented")
        }
        val fa = DoubleArray(d.size) { d[it]}
        return Pair(fa, u)
    }

    fun make_overstatement(overs: Double, use_style: Boolean = false): Double {
        val result =  (1 - overs / this.assorter.upper_bound) / (2 - this.margin!! / this.assorter.upper_bound)
        println("make_overstatement = $result this.margin = ${this.margin}")
        return result
    }

    //  assorter that corresponds to normalized overstatement error for an assertion
    fun overstatement_assorter(mvr: CvrSimple, cvr: CvrSimple, use_style: Boolean): Double {
        return (1 - this.assorter.overstatement(mvr, cvr, use_style) / this.assorter.upper_bound) /
                (2 - this.margin!! / this.assorter.upper_bound)
    }

    fun set_margin_from_cvrs(audit: AuditSimple, cvrs: List<CvrSimple>): AssertionSimple {
        val use_style = audit.use_styles
        val amean = this.assorter.mean(cvrs, use_style = use_style)
        if (amean < .5) {
            println("assertion $this not satisfied by CVRs: mean value is ${amean}")
        }
        this.margin = 2 * amean - 1
        println("assertion '${this.winner}v${this.loser}' amean=$amean, margin=${this.margin}")

        if (this.contest.audit_type == AuditType.POLLING) {
            this.test.u = this.assorter.upper_bound
        } else if (this.contest.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            this.test.u = 2 / (2 - this.margin!! / this.assorter.upper_bound)
        } else {
            throw NotImplementedError("audit type {this.contest.audit_type} not supported")
        }
        return this
    }

    override fun toString(): String {
        return "AssertionSimple(contest=${contest.id}, upper_bound=${assorter.upper_bound}, winner='$winner', loser='$loser', " +
                "margin=$margin, p_value=$p_value, p_history=$p_history, proved=$proved, sample_size=$sample_size)"
    }

}

fun make_all_assertions(contests: List<ContestSimple>) {
    for (contest in contests) {
        val scf = contest.choice_function
        val losers = contest.candidates.filter { !contest.winners.contains(it) }
        if (scf == SocialChoiceFunction.PLURALITY) {
            contest.assertions = // TODO mutable state
                make_plurality_assertions(contest = contest, winners = contest.winners, losers = losers)

        /* }
            if (scf == SocialChoiceFunction.SUPERMAJORITY) {
            contest.assertions.putAll(
                make_supermajority_assertion(
                contest = contest, winner = winrs[0],
                loser = losrs, /* share_to_win = contest.share_to_win, */
                test = test, estim = estim, bet = bet
            )
            )
                           } else if (scf == SocialChoiceFunction.IRV) {
                                // TODO like magic, contest.assertion_json appears from thin air!
                                // Assumption: contests[c].assertion_json yields list assertions in JSON format.
                                contest.assertions = make_assertions_from_json(
                                    contest = contest,
                                    candidates = contest.candidates,
                                    json_assertions = contest.assertion_json,
                                    test = test, estim = estim, bet = bet
                                ) */
        } else {
            throw Exception("Social choice function ${scf} is not implemented.")
        }
    }
}

fun make_plurality_assertions(contest: ContestSimple, winners: List<String>, losers: List<String>): Map<String, AssertionSimple> {

    val assertions = mutableMapOf<String, AssertionSimple>()

    for (winr in winners) {
        for (losr in losers) {
            val wl_pair = winr + " v " + losr

            assertions[wl_pair] = AssertionSimple(
                contest,
                winner = winr,
                loser = losr,
                assorter = AssorterFn(
                    contest = contest,

// (CVR.as_vote(c.get_vote_for(contest.id, winr)) - CVR.as_vote(c.get_vote_for(contest.id, losr)) + 1 ) / 2,
// CVR.as_vote: def as_vote(cls, v) -> int:
//        return int(bool(v))

                    assort =  { cvr: CvrSimple ->
                        val w = cvr.get_vote_for(contest.id, winr)
                        val l = cvr.get_vote_for(contest.id, losr)
                        val calc = (w - l + 1) * 0.5
                        calc },
                    upper_bound = 1.0,
                ),
                test = contest.getTestFn(),
            )
        }
    }

    return assertions
}

fun set_all_margins_from_cvrs(audit: AuditSimple, contests: List<ContestSimple>, cvr_list: List<CvrSimple>): Double {

    var min_margin = Double.POSITIVE_INFINITY
    for (con in contests) {
        val margins = mutableMapOf<String, Double>()
        for ((a,asn) in con.assertions) {
            asn.set_margin_from_cvrs(audit, cvr_list)
            val margin = asn.margin!!
            margins[a] = margin
            val u = if (con.audit_type == AuditType.POLLING) {
                asn.assorter.upper_bound
            } else if (con.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
                2 / (2 - margin / asn.assorter.upper_bound)
            } else {
                throw Exception ("audit type ${con.audit_type} not implemented")
            }
            asn.test.u = u
            min_margin = min(min_margin, margin)
        }
    }
    return min_margin
}