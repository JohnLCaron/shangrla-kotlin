package org.cryptobiotic.shangrla

/*
    Parameters
        ----------
        contest: Contest instance
            contest to which the assorter is relevant
        winner: str
            identifier for the nominal "winner" for this assertion. Can be an element of this.contest.candidates,
            an element of Contest.CANDIDATES, or an arbitrary label.
            Using an element of this.contest.candidates or an element of Contest.CANDIDATES can be useful for
            setting the margin in approval, plurality, and supermajority contests.
        loser: str
            identifier for the nominal "loser" for this assertion. Can be an element of this.contest.candidates,
            an element of Contest.CANDIDATES, or an arbitrary label.
            Using an element of this.contest.candidates or an element of Contest.CANDIDATES can be useful for
            setting the margin in approval, plurality, and supermajority contests.
        assorter: Assorter instance
            the assorter for the assertion
        margin: float
            the assorter margin. Generally this will not be known when the assertion is created, but will be set
            later.
        test: instance of class NonnegMean
            the function to find the p-value of the hypothesis that the assertion is true, i.e., that the
            assorter mean is <=1/2
        p_value: float
            the current p-value for the complementary null hypothesis that the assertion is false
        p_history: list
            the history of p-values, sample by sample. Generally, it is valid only for sequential risk-measuring
            functions.
        proved: boolean
            has the complementary null hypothesis been rejected?
        sample_size: int
            estimated total sample size to complete the audit of this assertion
        tally_pool_means: dict
            dict of reported assorter means for each `tally_pool`, for ONEAudit
 */
class Assertion(
    val contest: Contest,
    val assorter: Assorter,
    val winner: Candidates,
    val loser: Candidates,
    var margin: Float,
    val test: NonnegMean,
    p_value: Float,
    p_history: List<Float>,
    val proved: Boolean,
    var sample_size: Int?,
    tally_pool_means: Map<Any, Any>,
) {
    init {
        require(margin > 0) { "Margin ${margin} is nonpositive" }
    }

    //     find the margin for a list of Cvrs.
    //    By definition, the margin is twice the mean of the assorter, minus 1.}
    fun margin(cvr_list: List<CVR>, use_style: Boolean): Float {
        return 2 * this.assorter.mean(cvr_list, use_style = use_style) - 1
    }

    fun overstatement_assorter_margin(error_rate_1: Float, error_rate_2: Float): Float {
        /*
        find the overstatement assorter margin corresponding to an assumed rate of 1-vote and 2-vote overstatements

        Parameters
        ----------
        error_rate_1: float
        the assumed rate of one-vote overstatement errors in the CVRs
        error_rate_2: float
        the assumed rate of two-vote overstatement errors in the CVRs

        Returns
        -------
        the overstatement assorter margin implied by the reported margin and the assumed rates of overstatements
        */
        return (1 - (error_rate_2 + error_rate_1 / 2) * this.assorter.upper_bound / this.margin) /
                (2 * this.assorter.upper_bound / this.margin - 1)
    }


    fun overstatement_assorter_mean(error_rate_1: Float, error_rate_2: Float): Float {
        return (1 - error_rate_1 / 2 - error_rate_2) / (2 - this.margin / this.assorter.upper_bound)
    }


    // TODO original fun overstatement_assorter(mvr: List<CVR>, cvr: List<CVR>, use_style: Boolean): Float {
    fun overstatement_assorter(mvr: CVR, cvr: CVR, use_style: Boolean): Float {
        /*
        assorter that corresponds to normalized overstatement error for an assertion

        If `use_style == True`, then if the CVR contains the contest but the MVR does not,
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
                (2 - this.margin / this.assorter.upper_bound)
    }

    fun set_margin_from_cvrs(audit: Audit, cvr_list: List<CVR>) {
        /*
        find assorter margin from cvrs and store it

        Parameters
        ----------
        cvr_list: Collection
        cvrs from which the sample will be drawn
        use_style: bool
        is the sample drawn only from ballots that should contain the contest?

        Returns
        -------
        nothing

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
        if (amean < 1 / 2) {
            println("assertion $this not satisfied by CVRs: mean value is ${amean}")
        }
        this.margin = 2 * amean - 1
        if (this.contest.audit_type == AuditType.POLLING) {
            this.test.u = this.assorter.upper_bound
        } else if (this.contest.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            this.test.u = 2 / (2 - this.margin / this.assorter.upper_bound)
        } else {
            throw NotImplementedError("audit type {this.contest.audit_type} not supported")
        }
    }

    fun find_margin_from_tally(tallyInput: Map<Candidates, Int>? = null) {
        /*
        find the assorter margin between implied by a tally.

        Generally useful only for approval, plurality, and supermajority contests.

        Assumes the number of cards containing the contest has been set.

        Parameters
        ----------
        tally: dict
        dict of tallies for the candidates in the contest. Keys are candidates as listed
        in Contest.candidates. If `tally is None` tries to use the contest.tally.

        The margin for a supermajority contest with a winner is (see SHANRGLA section 2.3)
        2(pq/(2f) + (1 − q)/2 - 1/2) = q(p/f-1), where:
        q is the fraction of cards that have valid votes
        p is the fraction of cards that have votes for the winner
        f is the fraction of valid votes required to win.

        Returns
        -------
        nothing

        Side effects
        ------------
        sets this.margin

        */
        //         tally = tally if tally else self.contest.tally
        //        if self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY \
        //             or self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.APPROVAL:
        //            self.margin = (tally[self.winner]-tally[self.loser])/self.contest.cards
        //        elif self.contest.choice_function == Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY:
        //            if self.winner == Contest.CANDIDATES.NO_CANDIDATE or self.loser != Contest.CANDIDATES.ALL_OTHERS:
        //                raise NotImplementedError(f'TO DO: currently only support super-majority with a winner')
        //            else:
        //                q = np.sum([tally[c] for c in self.contest.candidates])/self.contest.cards
        //                p = tally[self.winner]/self.contest.cards
        //                self.margin = q*(p/self.contest.share_to_win - 1)
        //        else:
        //            raise NotImplementedError(f'social choice function {self.contest.choice_function} not supported')
        val tally = tallyInput ?: this.contest.tally
        if (this.contest.choice_function in listOf(SocialChoiceFunction.PLURALITY, SocialChoiceFunction.APPROVAL)) {
            this.margin = (tally[this.winner]!! - tally[this.loser]!!).toFloat() / this.contest.cards // // LOOK check nullable
        } else if (this.contest.choice_function == SocialChoiceFunction.SUPERMAJORITY) {
            if (this.winner == Candidates.NO_CANDIDATE || this.loser != Candidates.ALL_OTHERS) {
                throw NotImplementedError("TO DO: currently only support super-majority with a winner")
            } else {
                // val q = np.sum([tally[c] for c in this.contest.candidates])/this.contest.cards
                val q = this.contest.candidates.map { tally[it]!! }.sum() / this.contest.cards // LOOK check nullable
                val p = tally[this.winner]!! / this.contest.cards // LOOK check nullable
                this.margin = q * (p / this.contest.share_to_win - 1)
            }
        } else {
            throw NotImplementedError("social choice function {this.contest.choice_function} not supported")
        }
    }


    fun make_overstatement(overs: Float, use_style: Boolean = false): Float {
        /*
    return the numerical value corresponding to an overstatement of `overs` times the assorter upper bound `u`

    **Assumes that the margin has been set.**

    Parameters
    ----------
    overs: float
    the multiple of `u`
    use_style: bool
    flag to use style information. Only used if the assorter margin has not been set

    Returns
    -------
    the numerical value corresponding to an overstatement of that multiple

    */
        return (1 - overs / this.assorter.upper_bound) / (2 - this.margin / this.assorter.upper_bound)
    }

    fun mvrs_to_data(mvr_sample: List<CVR>, cvr_sample: List<CVR>): Pair<FloatArray, Float> {
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
        d: np.array; either assorter values or overstatement assorter values, depending on the audit method
        u: upper bound for the test
        */
        val margin = this.margin
        val upper_bound = this.assorter.upper_bound
        val con = this.contest
        val use_style = con.use_style

        var d: List<Float>
        var u: Float
        if (con.audit_type in listOf(AuditType.CARD_COMPARISON, AuditType.ONEAUDIT)) {
            require(mvr_sample.size == cvr_sample.size)
            val cvr2: List<Pair<CVR, CVR>> = mvr_sample.zip(cvr_sample)

            //       d = np.array(
            //                [this.overstatement_assorter(mvr_sample[i], cvr_sample[i], use_style=use_style)
            //                     for i in range(len(mvr_sample))
            //                         if ((not use_style) or
            //                         (cvr_sample[i].has_contest(con.id) and cvr_sample[i].sample_num <= con.sample_threshold))
            //                ])
            d = cvr2.filter { (mvr, cvr) -> !use_style || (cvr.has_contest(con.id) && cvr.sample_num <= con.sample_threshold) }
                    .map { (mvr, cvr) -> this.overstatement_assorter(mvr, cvr, use_style) }
            u = 2 / (2 - margin / upper_bound)
        } else if (con.audit_type == AuditType.POLLING) {  // Assume style information is irrelevant
            // d = np.array([this.assorter.assort(mvr_sample[i]) for i in range(len(mvr_sample))])
            d = mvr_sample.map { this.assorter.assort(it) }
            u = upper_bound
        } else {
            throw NotImplementedError("audit type ${con.audit_type} not implemented")
        }
        val fa = FloatArray(d.size) { d[it]}
        return Pair(fa, u)
    }

    fun find_sample_size(
        data: FloatArray? = null, prefix: Boolean = false, rate1: Float? = null, rate2: Float? = null,
        reps: Int? = null, quantile: Float = 0.5f, seed: Int = 1234567890
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
        data: np.array; observations on which to base the calculation. If `data is not None`, uses them in a bootstrap
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
                  else this.make_overstatement(overs = 0.0f)
        val small = if (this.contest.audit_type == AuditType.POLLING) 0
                    else this.make_overstatement(overs = 0.5f)
        val rate_1 = rate1 ?: ((1 - this.margin) / 2)   // rate of small values
        var x = FloatArray(this.test.N) { big } // array N floats, all equal to big
        if (this.contest.audit_type == AuditType.POLLING) {
            if (this.contest.choice_function == SocialChoiceFunction.IRV) {
                throw NotImplementedError("data must be provided to estimate sample sizes for IRV assertions")
            } else { // get tally
                if (this.contest.tally != null) {
                    val n_0 = this.contest.tally[this.loser]!! // LOOK nulllability
                    val n_big = this.contest.tally[this.winner]!! // LOOK nulllability
                    val n_half = this.test.N - n_0 - n_big
                    //         fun interleave_values(n_small: Int, n_med: Int, n_big: Int, small: Float, med: Float, big: Float): FloatArray {
                    x = interleave_values(n_0, n_half, n_big, big = big)
                } else {
                    throw Exception("contest ${this.contest} tally required but not defined") // LOOK ValueError
                }
            }
        } else if (this.contest.audit_type == AuditType.CARD_COMPARISON) { // # comparison audit
            //val rate_1_i = np.arange(0, this.test.N, step = int(1 / rate_1), dtype = int) if rate_1 else []
            //val rate_2_i = np.arange(0, this.test.N, step = int(1 / rate_2), dtype = int) if rate_2 else []
            //x[rate_1_i] = small
            //x[rate_2_i] = 0
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
                              small: Float = 0.0f, med: Float = 0.5f, big: Float = 1.0f): FloatArray {
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
            val x = FloatArray(N) // np.zeros(N)
            var i_small = 0
            var i_med = 0
            var i_big = 0
            var r_small = if (n_small != 0) 1 else 0 // TODO
            var r_med = if (n_med != 0) 1 else 0 // TODO
            var r_big = 1
            if (r_small != 0) {   //start with small TODO
                x[0] = small
                i_small = 1
                r_small = (n_small - i_small) / n_small
            } else if (r_med != 0) { // start with 1/2 TODO
                x[0] = med
                i_med = 1
                r_med = (n_med - i_med) / n_med
            } else {
                x[0] = big
                i_big = 1
                r_big = (n_big - i_big) / n_big
            }
            for (i in (1..N)) { // TODO wrong, should be  < N
                if (r_small > r_big) {
                    if (r_med > r_small) {
                        x[i] = med
                        i_med += 1
                        r_med = (n_med - i_med) / n_med
                    } else {
                        x[i] = small
                        i_small += 1
                        r_small = (n_small - i_small) / n_small
                    }
                } else if (r_med > r_big) {
                    x[i] = med
                    i_med += 1
                    r_med = (n_med - i_med) / n_med
                } else {
                    x[i] = big
                    i_big += 1
                    r_big = (n_big - i_big) / n_big
                }
            }
            return x
        }
    }
}