package org.cryptobiotic.shangrla.core

import kotlin.math.max
import kotlin.math.min


/*
    Parameters
        ----------
        contest: contest to which the assorter is relevant
        winner: identifier for the nominal "winner" for this assertion. Can be an element of this.contest.candidates,
            an element of Candidates, or an arbitrary label.
            Using an element of this.contest.candidates or an element of Candidates can be useful for
            setting the margin in approval, plurality, and supermajority contests.
        loser: identifier for the nominal "loser" for this assertion. Can be an element of this.contest.candidates,
            an element of Candidates, or an arbitrary label.
            Using an element of this.contest.candidates or an element of Candidates can be useful for
            setting the margin in approval, plurality, and supermajority contests.
        assorter: the assorter for the assertion
        margin: the assorter margin. Generally this will not be known when the assertion is created, but will be set later.
        test: NonnegMean; the function to find the p-value of the hypothesis that the assertion is true,
            i.e., that the assorter mean is <=1/2
        p_value: the current p-value for the complementary null hypothesis that the assertion is false
        p_history: the history of p-values, sample by sample. Generally, it is valid only for sequential risk-measuring
            functions.
        proved: has the complementary null hypothesis been rejected?
        sample_size: estimated total sample size to complete the audit of this assertion
        tally_pool_means: dict of reported assorter means for each `tally_pool`, for ONEAudit
 */
data class Assertion(
    val contest: Contest,
    val assorter: Assorter,
    val winner: String,
    val loser: String,
    var test: NonnegMean,
    var margin: Double? = null,
    var p_value: Double? = null,
    var p_history: List<Double>? = null,
    var proved: Boolean = false,
    var sample_size: Int? = null,
    var tally_pool_means: MutableMap<String, Double>? = null,
) {
    init {
        this.assorter.tally_pool_means = tally_pool_means
    }

    fun min_p(): Double {
        return java.util.Collections.min(p_history!!)
    }

    //    The margin for a list of CVRs.
    //    By definition, the margin is twice the mean of the assorter, minus 1.}
    fun margin(cvr_list: List<Cvr>, use_style: Boolean): Double {
        return 2 * this.assorter.mean(cvr_list, use_style = use_style) - 1
    }

    fun overstatement_assorter_margin(error_rate_1: Double = 0.0, error_rate_2: Double = 0.0): Double {
        /*
        find the overstatement assorter margin corresponding to an assumed rate of 1-vote and 2-vote overstatements

        Parameters
        ----------
        error_rate_1: the assumed rate of one-vote overstatement errors in the CVRs
        error_rate_2: the assumed rate of two-vote overstatement errors in the CVRs

        Returns
        -------
        the overstatement assorter margin implied by the reported margin and the assumed rates of overstatements
        */
        return (1 - (error_rate_2 + error_rate_1 / 2) * this.assorter.upper_bound / this.margin!!) /
                (2 * this.assorter.upper_bound / this.margin!! - 1)
    }


    fun overstatement_assorter_mean(error_rate_1: Double = 0.0, error_rate_2: Double = 0.0): Double {
        return (1 - error_rate_1 / 2 - error_rate_2) / (2 - this.margin!! / this.assorter.upper_bound)
    }


    fun overstatement_assorter(mvr: Cvr, cvr: Cvr, use_style: Boolean): Double {
        /*
        assorter that corresponds to normalized overstatement error for an assertion

        If `use_style == true`, then if the CVR contains the contest but the MVR does not,
        that is considered to be an overstatement, because the ballot is presumed to contain
        the contest .

        If `use_style == False`, then if the CVR contains the contest but the MVR does not,
        the MVR is considered to be a non -vote in the contest .

        Parameters
        -----------
        mvr: Cvr the manual interpretation of voter intent TODO List or CVR or Contest?
        cvr: Cvr the machine-reported cast vote record .

        Returns
        --------
        over: float
        (1 - o / u) / (2 - v / u), where
        o is the overstatement
        u is the upper bound on the value the assorter assigns to any ballot
        v is the assorter margin
        */
        //         return (1-this.assorter.overstatement(mvr, cvr, use_style)
        //                / this.assorter.upper_bound)/(2-this.margin/this.assorter.upper_bound)
        // TODO possible error in shangrla
        return (1 - this.assorter.overstatement(mvr, cvr, use_style) / this.assorter.upper_bound) /
                (2 - this.margin!! / this.assorter.upper_bound)
    }

    fun set_margin_from_cvrs(audit: Audit, cvr_list: List<Cvr>): Assertion {
        /*
        find assorter margin from cvrs and store it

        Parameters
        ----------
        cvr_list: cvrs from which the sample will be drawn
        use_style: is the sample drawn only from ballots that should contain the contest?

        Side effects
        ------------
        sets assorter.margin
        */

        //         if len(audit.strata) > 1:
        //            raise NotImplementedError('stratified audits not yet supported')
        //        stratum = next(iter(audit.strata.values()))
        //        use_style = stratum.use_style
        //        amean =self.assorter.mean(cvr_list, use_style=use_style)
        //        if amean < 1/2:
        //            warnings.warn(f"assertion {self} not satisfied by CVRs: mean value is {amean}")
        //        self.margin = 2*amean-1
        //        if self.contest.audit_type == Audit.AUDIT_TYPE.POLLING:
        //            self.test.u = self.assorter.upper_bound
        //        elif self.contest.audit_type in [Audit.AUDIT_TYPE.CARD_COMPARISON, Audit.AUDIT_TYPE.ONEAUDIT]:
        //            self.test.u = 2/(2-self.margin/self.assorter.upper_bound)
        //        else:
        //            raise NotImplementedError(f'audit type {self.contest.audit_type} not supported')

        if (audit.strata.size > 1) throw NotImplementedError("stratified audits not yet supported")

        val stratum = audit.strata.values.first()
        val use_style = stratum.use_style
        val amean = this.assorter.mean(cvr_list, use_style = use_style)
        if (amean < .5) {
            println("assertion $this not satisfied by CVRs: mean value is ${amean}")
        }
        this.margin = 2 * amean - 1
        if (this.contest.audit_type == AuditType.POLLING) {
            this.test.u = this.assorter.upper_bound
        } else if (this.contest.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            this.test.u = 2 / (2 - this.margin!! / this.assorter.upper_bound)
        } else {
            throw NotImplementedError("audit type {this.contest.audit_type} not supported")
        }
        return this
    }

    fun find_margin_from_tally(tallyInput: Map<String, Int>? = null) {
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
            throw NotImplementedError("social choice function {this.contest.choice_function} not supported")
        }
    }


    fun make_overstatement(overs: Double, use_style: Boolean = false): Double {
        /*
        Caclulate the numerical value corresponding to an overstatement of `overs` times the assorter upper bound `u`
        **Assumes that the margin has been set.**

        Parameters
        ----------
        overs: the multiple of `u`
        use_style: flag to use style information. Only used if the assorter margin has not been set

        Returns
        -------
        the numerical value corresponding to an overstatement of that multiple
        */

        val result =  (1 - overs / this.assorter.upper_bound) / (2 - this.margin!! / this.assorter.upper_bound)
        println("make_overstatement = $result this.margin = ${this.margin}")
        return result
    }

    fun mvrs_to_data(mvr_sample: List<Cvr>, cvr_sample: List<Cvr>): Pair<DoubleArray, Double> {
        /*
        Process mvrs (and, for comparison audits, cvrs) to create data for the assertion"s test 
        and for sample size simulations.
    
        Creates assorter values for the mvrs, or overstatement assorter values using the mvrs and cvrs,
        according to whether the audit uses ballot polling or card-level comparison
    
        The margin should be set before calling this function.
    
        mvr_sample and cvr_sample should be ordered using CVR.prep_comparison_sample() or
        CVR.prep_polling_sample() before calling this routine
    
        Parameters
        ----------
        mvr_sample: list of CVR objects
        corresponding MVRs
        cvr_sample: list of CVR objects
        sampled CVRs
    
        Returns
        -------
        d: DoubleArray; either assorter values or overstatement assorter values, depending on the audit method
        u: upper bound for the test
        */
        requireNotNull(this.margin)
        val margin = this.margin!!
        val upper_bound = this.assorter.upper_bound
        val con = this.contest
        val use_style = true // con.use_style

        var d: List<Double>
        var u: Double
        if (con.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            require(mvr_sample.size == cvr_sample.size)
            val cvr2: List<Pair<Cvr, Cvr>> = mvr_sample.zip(cvr_sample)

            //       d = np.array(
            //                [this.overstatement_assorter(mvr_sample[i], cvr_sample[i], use_style=use_style)
            //                     for i in range(len(mvr_sample))
            //                         if ((not use_style) or
            //                         (cvr_sample[i].has_contest(con.id) and cvr_sample[i].sample_num <= con.sample_threshold))
            //                ])
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
        // convert to double array
        val fa = DoubleArray(d.size) { d[it]}
        return Pair(fa, u)
    }

    fun find_sample_size(
        data: DoubleArray? = null,
        prefix: Boolean = false,
        rate1: Double? = null,
        rate2: Double? = null,
        reps: Int? = null,
        quantile: Double = 0.5,
        seed: Int = 1234567890
    ): Int {
        /*
        Estimate sample size needed to reject the null hypothesis that the assorter mean is <=1/2,
        for the specified risk function, given:
        - for comparison audits, the assorter margin and assumptions about the rate of overstatement errors
        - for polling audits, either a set of assorter values, or the assumption that the reported tallies
        are correct

        If `data is not None`, uses data to make the estimate. There are three strategies:
        1. if `reps is None`, tile the data to make a list of length N
        2. if `reps is not None and not prefix`, sample from the data with replacement to make `reps` lists of
        length N
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

        This function is for a single assertion.

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
        sample_size: int
        sample size estimated to be sufficient to confirm the outcome if data are generated according to
        the assumptions

        Side effects
        ------------
        sets the sample_size attribute of the assertion

        */
        if (data != null) {
            return this.test.sample_size(
                data, alpha = this.contest.risk_limit, reps = reps,
                prefix = prefix, quantile = quantile, // seed = seed TODO
            )
        }

        require(margin != null) { "Margin as not set"  }
        val amargin = margin!!
        require(amargin > 0.0) { "Margin ${amargin} is nonpositive" }

        /*
        Construct data .
        For POLLING, values are 0, 1/2, and u.
        For CARD_COMPARISON, values are overstatement assorter values corresponding to
        overstatements of 2u(at rate_2), u (at rate_1), or 0.
        */
        //             big = (self.assorter.upper_bound if self.contest.audit_type == Audit.AUDIT_TYPE.POLLING \
        //                  else self.make_overstatement(overs=0))
        //            small = (0 if self.contest.audit_type == Audit.AUDIT_TYPE.POLLING
        //                     else self.make_overstatement(overs=1/2))
        //            rate_1 = rate_1 if rate_1 is not None else (1-self.margin)/2   # rate of small values
        //            x = big*np.ones(self.test.N)
        //            if self.contest.audit_type == Audit.AUDIT_TYPE.POLLING:
        //                if self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.IRV:
        //                    raise NotImplementedError(f'data must be provided to estimate sample sizes for IRV assertions')
        //                else: # get tally
        //                    if self.contest.tally:
        //                        n_0 = self.contest.tally[self.loser]
        //                        n_big = self.contest.tally[self.winner]
        //                        n_half = self.test.N - n_0 - n_big
        //                        x = interleave_values(n_0, n_half, n_big, big=big)
        //                    else:
        //                        raise ValueError(f'contest {self.contest} tally required but not defined')
        //            elif self.contest.audit_type == Audit.AUDIT_TYPE.CARD_COMPARISON: # comparison audit
        //                rate_1_i = np.arange(0, self.test.N, step=int(1/rate_1), dtype=int) if rate_1 else []
        //                rate_2_i = np.arange(0, self.test.N, step=int(1/rate_2), dtype=int) if rate_2 else []
        //                x[rate_1_i] = small
        //                x[rate_2_i] = 0
        //            else:
        //                raise NotImplementedError(f'audit type {self.contest.audit_type} for contest {self.contest} not implemented')
        //            sample_size = self.test.sample_size(x, alpha=self.contest.risk_limit, reps=reps,
        //                                                prefix=prefix, quantile=quantile, seed=seed)
        //        self.sample_size = sample_size
        //        return sample_size
        val big = if (this.contest.audit_type == AuditType.POLLING) this.assorter.upper_bound
                  else this.make_overstatement(overs = 0.0)
        val small = if (this.contest.audit_type == AuditType.POLLING) 0.0
                    else this.make_overstatement(overs = 0.5)
        val rate_1 = rate1 ?: ((1 - amargin) / 2)   // rate of small values
        var x = DoubleArray(this.test.N) { big } // array N floats, all equal to big
        println("  big = $big small = $small")
        // println("  x = ${x.contentToString()}")


        if (this.contest.audit_type == AuditType.POLLING) {
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

        } else if (this.contest.audit_type == AuditType.CARD_COMPARISON) { // # comparison audit
            //     arange([start,] stop[, step,], dtype=None, *, like=None) : Return evenly spaced values within a given interval.
            val rate_1_i = if (rate_1 != null) numpy_arange(0, this.test.N, step = (1.0 / rate_1).toInt()) else IntArray(0)
            val rate_2_i = if (rate2 != null && rate2 != 0.0) numpy_arange(0, this.test.N, step = (1.0 / rate2).toInt())  else IntArray(0)
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

    companion object {
        fun interleave_values(n_small: Int, n_med: Int, n_big: Int,
                              small: Double = 0.0, med: Double = 0.5, big: Double = 1.0): DoubleArray {
            /*
                make an interleaved population of n_s values equal to small, n_m values equal to med, and n_big equal to big
                Start with a small if n_small > 0
            */

            //        N = n_small + n_med + n_big
            //        x = np.zeros(N)
            //        i_small = 0
            //        i_med = 0
            //        i_big = 0
            //        r_small = 1 if n_small else 0
            //        r_med = 1 if n_med else 0
            //        r_big = 1
            //        if r_small:   # start with small
            //            x[0] = small
            //            i_small = 1
            //            r_small = (n_small-i_small)/n_small
            //        elif r_med: # start with 1/2
            //            x[0] = med
            //            i_med = 1
            //            r_med = (n_med-i_med)/n_med
            //        else:
            //            x[0] = big
            //            i_big = 1
            //            r_big = (n_big-i_big)/n_big
            //        for i in range(1, N):
            //            if r_small > r_big:
            //                if r_med > r_small:
            //                    x[i] = med
            //                    i_med += 1
            //                    r_med = (n_med-i_med)/n_med
            //                else:
            //                    x[i] = small
            //                    i_small += 1
            //                    r_small = (n_small-i_small)/n_small
            //            elif r_med > r_big:
            //                x[i] = med
            //                i_med += 1
            //                r_med = (n_med-i_med)/n_med
            //            else:
            //                x[i] = big
            //                i_big += 1
            //                r_big = (n_big-i_big)/n_big
            //        return x
            val N = n_small + n_med + n_big
            val x = DoubleArray(N) // np.zeros(N)
            var i_small = 0
            var i_med = 0
            var i_big = 0
            var r_small = if (n_small != 0) 1.0 else 0.0
            var r_med = if (n_med != 0) 1.0 else 0.0
            var r_big = 1.0
            if (r_small != 0.0) {   //start with small
                x[0] = small
                i_small = 1
                r_small = (n_small - i_small).toDouble() / n_small
            } else if (r_med != 0.0) { // start with 1/2
                x[0] = med
                i_med = 1
                r_med = (n_med - i_med).toDouble() / n_med
            } else {
                x[0] = big
                i_big = 1
                r_big = (n_big - i_big).toDouble() / n_big
            }
            for (i in (1 until N)) {
                if (r_small > r_big) {
                    if (r_med > r_small) {
                        x[i] = med
                        i_med += 1
                        r_med = (n_med - i_med).toDouble() / n_med
                    } else {
                        x[i] = small
                        i_small += 1
                        r_small = (n_small - i_small).toDouble() / n_small
                    }
                } else if (r_med > r_big) {
                    x[i] = med
                    i_med += 1
                    r_med = (n_med - i_med).toDouble() / n_med
                } else {
                    x[i] = big
                    i_big += 1
                    r_big = (n_big - i_big).toDouble() / n_big
                }
            }
            return x
        }

        fun make_plurality_assertions(contest: Contest, winner: List<String>, loser: List<String>,
                                      test: TestFn? = null, estim: EstimatorFn? = null, bet: BetFn? = null): Map<String, Assertion> {
            /*
            Construct assertions that imply the winner(s) got more votes than the loser(s).

            The assertions are that every winner beat every loser: there are
            len(winner)*len(loser) pairwise assertions in all.

            Parameters
            -----------
            contest: instance of Contest; contest to which the assertions are relevant
            winner: list; list of identifiers of winning candidate(s)
            loser: list; list of identifiers of losing candidate(s)

            Returns
            --------
            a dict of Assertions
            */
            //         assertions = {}
            //        test = test if test is not None else contest.test
            //        estim = estim if estim is not None else contest.estim
            //        bet = bet if bet is not None else contest.bet
            //        for winr in winner:
            //            for losr in loser:
            //                wl_pair = winr + ' v ' + losr
            //                _test = NonnegMean(test=test, estim=estim, bet=bet, g=contest.g, u=1, N=contest.cards,
            //                                       t=1/2, random_order=true)
            //                assertions[wl_pair] = Assertion(contest, winner=winr, loser=losr,
            //                                         assorter=Assorter(contest=contest,
            //                                             assort = lambda c, contest_id=contest.id, winr=winr, losr=losr:
            //                                                 (CVR.as_vote(c.get_vote_for(contest.id, winr))
            //                                                 - CVR.as_vote(c.get_vote_for(contest.id, losr))
            //                                                  + 1)/2,
            //                                            upper_bound=1),
            //                                         test=_test)
            //        return assertions
            val assertions = mutableMapOf<String, Assertion>()
            val test = test ?: contest.testFn
            val estim = estim ?: contest.estimFn
            val bet = bet ?: contest.betFn

            for (winr in winner) {
                for (losr in loser) {
                    val wl_pair = winr + " v " + losr
                    val _test = NonnegMean(
                        u = 1.0, N = contest.ncards, t = .5,
                        testFnOverride = test, estimFnOverride = estim, betFnOverride = bet, gOverride = contest.g)

                    assertions[wl_pair] = Assertion(
                        contest,
                        winner = winr,
                        loser = losr,
                        assorter = Assorter(
                            contest = contest,
                            assort =  { cvr: Cvr ->
                                val w = Cvr.as_vote(cvr.get_vote_for(contest.id, winr))
                                val l = Cvr.as_vote(cvr.get_vote_for(contest.id, losr))
                                val calc = (w - l + 1) * 0.5
                                calc},
                            upper_bound = 1.0,
                            ),
                    test = _test)
                }
            }

            return assertions
        }

        fun make_supermajority_assertion(contest: Contest, winner: String, loser: List<String>,
                 test: TestFn? = null, estim: EstimatorFn? = null, bet:BetFn? = null): Map<String, Assertion> {
            /*
            Construct assertion that winner got >= share_to_win \in (0,1) of the valid votes

            **TO DO: This method assumes there was a winner. To audit that there was no winner requires
            flipping things.**

            An equivalent condition is:

            (votes for winner)/(2*share_to_win) + (invalid votes)/2 > 1/2.

            Thus the correctness of a super-majority outcome--where share_to_win >= 1/2--can
            be checked with a single assertion.

            share_to_win < 1/2 might be useful for some social choice functions, including
            primaries where candidates who receive less than some threshold share are
            eliminated.

            A CVR with a mark for more than one candidate in the contest is considered an
            invalid vote.

            Parameters
            -----------
            contest: contest object instance to which the assertion applies
            winner: identifier of winning candidate
            loser: list of identifiers of losing candidate(s)
            share_to_win: float; fraction of the valid votes the winner must get to win
            test: instance of NonnegMean; risk function for the contest
            estim: an estimation method of NonnegMean; estimator the alpha_mart test uses for the alternative
            bet: method to choose the bet for betting_mart risk function

            Returns
            --------
            a dict containing one Assertion
            */
            //         assertions = {}
            //        wl_pair = winner + ' v ' + Contest.CANDIDATES.ALL_OTHERS
            //        cands = loser.copy()
            //        cands.append(winner)
            //        _test = NonnegMean(test=test, estim=estim, bet=bet, u=1/(2*contest.share_to_win),
            //                           N=contest.cards, t=1/2, random_order=True)
            //        assertions[wl_pair] = Assertion(contest, winner=winner, loser=Contest.CANDIDATES.ALL_OTHERS,
            //                                 assorter=Assorter(contest=contest,
            //                                          assort = lambda c, contest_id=contest.id:
            //                                                CVR.as_vote(c.get_vote_for(contest.id, winner))
            //                                                      /(2*contest.share_to_win)
            //                                                if c.has_one_vote(contest.id, cands) else 1/2,
            //                                          upper_bound=1/(2*contest.share_to_win)), test=_test)
            //        return assertions
            val assertions = mutableMapOf<String, Assertion>()
            val wl_pair = winner + " v " + Candidates.ALL_OTHERS
            val cands = mutableListOf<String>()
            cands.addAll(loser)
            cands.add(winner)

            val _test = NonnegMean(
                u = 1 / (2 * contest.share_to_win), N = contest.ncards, t = 0.5,
                testFnOverride = test, estimFnOverride = estim, betFnOverride = bet,
            )
            assertions[wl_pair] = Assertion(
                contest, winner = winner, loser = Candidates.ALL_OTHERS.name,
                assorter = Assorter(
                    contest = contest,
                    assort =  { cvr: Cvr -> if (cvr.has_one_vote(contest.id, cands))
                        (Cvr.as_vote(cvr.get_vote_for(contest.id, winner))
                        / (2 * contest.share_to_win)) else .5 },
                upper_bound = 1 / (2 * contest.share_to_win)),
                test = _test)

            return assertions
        }

        fun make_all_assertions(contests: List<Contest>) {
            /*
            Construct all the assertions to audit the contests and add the assertions to the contest dict

            Side Effects
            ------------
            creates assertions and adds the dict of assertions relevant to each contest to the contest object's `assertions` attribute
            */
            //         for c, con in contests.items():
            //            scf = con.choice_function
            //            winrs = con.winner
            //            losrs = list(set(con.candidates) - set(winrs))
            //            test = con.test
            //            estim = con.estim
            //            bet = con.bet
            //            if scf == Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY:
            //                contests[c].assertions = Assertion.make_plurality_assertions(contest=con, winner=winrs, loser=losrs,
            //                                                                                test=test, estim=estim, bet=bet)
            //            elif scf == Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY:
            //                contests[c].assertions = Assertion.make_supermajority_assertion(contest=con, winner=winrs[0],
            //                                                    loser=losrs, share_to_win=con.share_to_win,
            //                                                    test=test, estim=estim, bet=bet)
            //            elif scf == Contest.SOCIAL_CHOICE_FUNCTION.IRV:
            //                # Assumption: contests[c].assertion_json yields list assertions in JSON format.
            //                contests[c].assertions = Assertion.make_assertions_from_json(contest=con,
            //                                                    candidates=con.candidates,
            //                                                    json_assertions=con.assertion_json,
            //                                                    test=test, estim=estim, bet=bet)
            //            else:
            //                raise NotImplementedError(f'Social choice function {scf} is not implemented.')
            //        return True
            for (contest in contests) {
                val scf = contest.choice_function
                val winrs = contest.reported_winners
                val losrs = contest.candidates.filter { !winrs.contains(it) }
                val test = contest.testFn
                val estim = contest.estimFn
                val bet = contest.betFn
                if (scf == SocialChoiceFunction.PLURALITY) {
                    contest.assertions.putAll(make_plurality_assertions(
                        contest = contest, winner = winrs, loser = losrs,
                        test = test, estim = estim, bet = contest.betFn
                    ))
                } else if (scf == SocialChoiceFunction.SUPERMAJORITY) {
                    contest.assertions.putAll(make_supermajority_assertion(
                        contest = contest, winner = winrs[0],
                        loser = losrs, /* share_to_win = contest.share_to_win, */
                        test = test, estim = estim, bet = bet
                    ))
/*                } else if (scf == SocialChoiceFunction.IRV) {
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

        fun set_all_margins_from_cvrs(audit: Audit, contests: List<Contest>, cvr_list: List<Cvr>): Double {
            /*
            Find all the assorter margins in a set of Assertions. Updates the dict of dicts of assertions
            and the contest dict.

            Appropriate only if cvrs are available. Otherwise, base margins on the reported results.

            This function is primarily about side-effects on the assertions in the contest dict.

            Parameters
            ----------
            audit: Audit; information about the audit
            contests: dict of Contest objects
            cvr_list: collection of CVR objects

            Returns
            -------
            min_margin: float; smallest margin in the audit

            Side effects
            ------------
            sets the margin of every assertion
            sets the assertion.test.u for every assertion, according to whether
            `assertion.contest.audit_type==Audit.AUDIT_TYPE.POLLING`
            or `assertion.contest.audit_type in [Audit.AUDIT_TYPE.CARD_COMPARISON, Audit.AUDIT_TYPE.ONEAUDIT]`
            */
    //         min_margin = np.infty
    //        for c, con in contests.items():
    //            con.margins = {}
    //            for a, asn in con.assertions.items():
    //                asn.set_margin_from_cvrs(audit, cvr_list)
    //                margin = asn.margin
    //                con.margins.update({a: margin})
    //                if con.audit_type==Audit.AUDIT_TYPE.POLLING:
    //                    u = asn.assorter.upper_bound
    //                elif con.audit_type in [Audit.AUDIT_TYPE.CARD_COMPARISON, Audit.AUDIT_TYPE.ONEAUDIT]:
    //                    u = 2/(2-margin/asn.assorter.upper_bound)
    //                else:
    //                    raise NotImplementedError(f'audit type {con.audit_type} not implemented')
    //                asn.test.u = u
    //                min_margin = min(min_margin, margin)
    //        return min_margin
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

        fun set_p_values(contests: List<Contest>, mvr_sample: List<Cvr>, cvr_sample: List<Cvr>): Double {
            /*
            Find the p-value for every assertion and update assertions & contests accordingly

            update p_value, p_history, proved flag, the maximum p-value for each contest.

            Primarily about side-effects.

            Parameters
            ----------
            contests: the contest data structure. outer keys are contest identifiers; inner keys are assertions
            mvr_sample: list of CVR objects; the manually ascertained voter intent from sheets, including entries for phantoms
            cvr_sample: list of CVR objects; the cvrs for the same sheets, for ballot-level comparison audits not needed for polling audits

            Returns
            -------
            p_max: float; largest p-value for any assertion in any contest

            Side-effects
            ------------
            Sets u for every test for every assertion, according to whether the corresponding audit method
            is Audit.AUDIT_TYPE.CARD_COMPARISON, Audit.AUDIT_TYPE.ONEAUDIT, or Audit.AUDIT_TYPE.POLLING.
            Sets contest max_p to be the largest P-value of any assertion for that contest
            Updates p_value, p_history, and proved for every assertion
            */
    //         if cvr_sample is not None:
    //            assert len(mvr_sample) == len(cvr_sample), "unequal numbers of cvrs and mvrs"
    //        p_max = 0
    //        for c, con in contests.items():
    //            con.p_values = {}
    //            con.proved = {}
    //            contest_max_p = 0
    //            for a, asn in con.assertions.items():
    //                d, u = asn.mvrs_to_data(mvr_sample, cvr_sample)
    //                asn.test.u = u       # set upper bound for the test for each assorter
    //                asn.p_value, asn.p_history = asn.test.test(d)
    //                asn.proved = (asn.p_value <= con.risk_limit) or asn.proved
    //                con.p_values.update({a: asn.p_value})
    //                con.proved.update({a: asn.proved})
    //                contest_max_p = np.max([contest_max_p, asn.p_value])
    //            contests[c].max_p = contest_max_p
    //            p_max = np.max([p_max, contests[c].max_p])
    //        return p_max

            if (cvr_sample != null) {
                require(mvr_sample.size == cvr_sample.size) {"unequal numbers of cvrs and mvrs"}
            }
            var p_max = 0.0
            for (con in contests) {
                for ((_, asn) in con.assertions) {
                    val (d, u) = asn.mvrs_to_data(mvr_sample, cvr_sample)
                    asn.test.u = u       // set upper bound for the test for each assorter
                    val (p_value, p_history) = asn.test.test(d)
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