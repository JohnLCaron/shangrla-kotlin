package org.cryptobiotic.shangrla.core

import kotlin.math.abs
import kotlin.math.min

typealias TestFn = (x: DoubleArray) -> Pair<Double, DoubleArray>
typealias EstimatorFn = (x: DoubleArray) -> DoubleArray
typealias BetFn = (x: DoubleArray) -> DoubleArray

/*
    Tests of the hypothesis that the mean of a population of values in [0, u] is less than or equal to t.
        Several tests are implemented, all ultimately based on martingales or supermartingales:
            Kaplan-Kolmogorov (with and without replacement)
            Kaplan-Markov (without replacement)
            Kaplan-Wald (without replacement)
            Wald SPRT (with and with replacement)
            ALPHA supermartingale test (with and without replacement)
            Betting martingale tests (with and without replacement)
    Some tests work for all nonnegative populations; others require a finite upper bound `u`.
    Many of the tests have versions for sampling with replacement (`N=np.inf`) and for sampling
    without replacement (`N` finite).
    Betting martingales and ALPHA martingales are different parametrizations of the same tests, but
      lead to different heuristics for selecting the parameters.
 */
class NonnegMean (
    testFn: TestFn?,
    estimFn: EstimatorFn?,
    betFn: BetFn?,
    var u: Double = 1.0,
    val N: Int = Int.MAX_VALUE, // TODO withReplacement
    val t: Double = 0.5,
    val random_order: Boolean = true,
) {
    val uulp = (u * ulpm)
    val test: TestFn = testFn ?: { x -> this.alpha_mart(x) }
    val estim: EstimatorFn = estimFn ?: { x -> this.fixed_alternative_mean(x) }
    val bet: BetFn = betFn ?: { x -> this.fixed_bet(x) }
    val lam = 0.5 // initial fraction of fortune to bet TODO allow to be set

    fun alpha_mart(x: DoubleArray): Pair<Double, DoubleArray> {
        /*
        Finds the ALPHA martingale for the hypothesis that the population
        mean is less than or equal to t using a martingale method,
        for a population of size N, based on a series of draws x.

        **The draws must be in random order**, or the sequence is not a supermartingale under the null

        If N is finite, assumes the sample is drawn without replacement
        If N is infinite, assumes the sample is with replacement

        Parameters
        ----------
        x: list corresponding to the data
        attributes used:
            keyword arguments for estim() and for this function
            u: float > 0 (default 1); upper bound on the population
            eta: float in (t,u] (default u*(1-eps))
                value parametrizing the bet. Use alternative hypothesized population mean for polling audit
                or a value nearer the upper bound for comparison audits


        Returns
        -------
        p: float; sequentially valid p-value of the hypothesis that the population mean is less than or equal to t
        p_history: sample by sample history of p-values. Not meaningful unless the sample is in random order.
        */
        // val atol = kwargs.get("atol", 2 * np.finfo(float).eps)
        // val rtol = kwargs.get("rtol", 10e-6)
        val (_, Stot, _, m) = this.sjm(N, t, x)

        //         with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        //            etaj = self.estim(x)
        //            terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
        val etaj = this.estim(x)
        val xp = DoubleArray(x.size) { (x[it] * etaj[it] / m [it] + (u - x[it]) * (u - etaj[it]) / (u - m[it])) / u }
        val terms = numpy_cumprod(xp)

        // terms[m > u] = 0  # true mean is certainly less than hypothesized
        repeat(terms.size) {if (m[it] > u) terms[it] = 0.0 } // true mean is certainly less than hypothesized
        // terms[np.isclose(0, m, atol=atol)] = 1  # ignore
        repeat(terms.size) {if (numpy_isclose(0.0, m[it])) terms[it] = 1.0 } // ignore
        // terms[np.isclose(0, terms, atol=atol)] = ( 1 } # martingale effectively vanishes; p-value 1
        repeat(terms.size) {if (numpy_isclose(0.0, terms[it])) terms[it] = 1.0 } // martingale effectively vanishes; p-value 1
        //        terms[m < 0] = np.inf  # true mean certainly greater than hypothesized
        repeat(terms.size) {if (m[it] < 0.0) terms[it] = Double.POSITIVE_INFINITY } // true mean is certainly less than hypothesized

        // terms[-1] = ( np.inf if Stot > N * t else terms[-1] )  # final sample makes the total greater than the null
        // -1 is the last element
        if (Stot > N * t) terms[terms.size-1] = N * t // final sample makes the total greater than the null

        // return min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms)
        // np.minimum = element-wise minumum, presumably the smaller of 1 and 1/term
        // np.max = maximum of an array
        // min = min or an iterable

        val npmin = terms.map{ min(1.0, 1.0 / it) }.toDoubleArray()
        val npmax = terms.max()
        return Pair(min(1.0, 1.0 / npmax), npmin)
    }

    data class CumulativeSum(val S: DoubleArray, val Stot: Double, val indices: IntArray, val mean: DoubleArray)

    fun sjm(N: Int, t: Double, x: DoubleArray, withReplacement: Boolean = false): CumulativeSum {
        /*
        This method calculates the cumulative sum of the input array `x`, the total sum of `x`,
        an array of indices, and the mean of the population after each draw if the null hypothesis is true.

        Parameters
        ----------
        N : int or float; The size of the population. If N is np.inf, it means the sampling is with replacement.
        t : float; The hypothesized population mean under the null hypothesis.
        x : np.array; The input data array.

        Returns
        -------
        S : np.array; The cumulative sum of the input array `x`, excluding the last element.
        Stot : float; The total sum of the input array `x`.
        j : np.array; An array of indices from 1 to the length of `x`.
        m : np.array; The mean of the population after each draw if the null hypothesis is true.
        */

//         assert isinstance(N, int) or (math.isinf(N) and N > 0), "Population size is not an integer!"
//        S = np.insert(np.cumsum(x), 0, 0)  # 0, x_1, x_1+x_2, ...,
//        Stot = S[-1]  # sample total
//        S = S[0:-1]  # same length as the data
//        j = np.arange(1, len(x) + 1)  # 1, 2, 3, ..., len(x)
//        assert j[-1] <= N, "Sample size is larger than the population!"
//        m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true (t=eta is the mean)
//        return S, Stot, j, m

        // require( isinstance (N, int) || (math.isinf(N) and N > 0)) { "Population size is not an integer!"}
        //        S = np.insert(np.cumsum(x), 0, 0)  # 0, x_1, x_1+x_2, ...,
        val cum_sum = numpy_cumsum(x)
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Stot = S[S.size-1]  // sample total ""array[-1] means the last element"
        // S = S[0:-1]  // same length as the data. TODO wtf

//        j = np.arange(1, len(x) + 1)  # 1, 2, 3, ..., len(x)
//        assert j[-1] <= N, "Sample size is larger than the population!"
        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        require( x.size <= N) { "Sample size is larger than the population!" }
//        m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true (t=eta is the mean)
        val m = DoubleArray(x.size) { if (!withReplacement) t else (N * t - S[it]) / (N - j[it] + 1)  }
        return CumulativeSum(S, Stot, j, m)
    }

    fun fixed_alternative_mean(x: DoubleArray, etap: Double? = null): DoubleArray {
        /*
        Compute the alternative mean just before the jth draw, for a fixed alternative that the original population mean is eta.
        Throws a warning if the sample implies that the fixed alternative is false (because the population would
        have negative values or values greater than u.

        S_1 := 0
        S_j := \sum_{i=1}^{j-1} x_i, j >= 1
        eta_j := (N*eta-S_j)/(N-j+1) if np.isfinite(N) else t

        Parameters
        ----------
        x: np.array
            input data
        kwargs:
            eta: float in (t, u) (default u*(1-eps))
                alternative hypothethesized value for the population mean
            u: float > 0 (default 1)
                upper bound on the population values
        */
        val eta = etap ?: uulp
        val (_, _, _, m) = this.sjm(N, eta, x)
        val negs = m.filter { it < 0.0 }
        if (negs.sum() > 0) println("Implied population mean is negative in ${negs} of ${x.size} terms")
        val pos = m.filter { it < u }
        if (pos.sum() > 0) println("Implied population mean is greater than ${u} in ${pos} of ${x.size} terms")
        return m
    }

    /*
    fun shrink_trunc(x: DoubleArray): DoubleArray {
        /*
    apply shrinkage/truncation estimator to an array to construct a sequence of "alternative" values

    sample mean is shrunk towards eta, with relative weight d compared to a single observation,
    then that combination is shrunk towards u, with relative weight f/(stdev(x)).

    The result is truncated above at u*(1-eps) and below at m_j+e_j(c,j)

    Shrinking towards eta stabilizes the sample mean as an estimate of the population mean.
    Shrinking towards u takes advantage of low-variance samples to grow the test statistic more rapidly.

    The running standard deviation is calculated using Welford's method.

    S_1 := 0
    S_j := \sum_{i=1}^{j-1} x_i, j >= 1
    m_j := (N*t-S_j)/(N-j+1) if np.isfinite(N) else t
    e_j := c/sqrt(d+j-1)
    sd_1 := sd_2 = 1
    sd_j := sqrt[(\sum_{i=1}^{j-1} (x_i-S_j/(j-1))^2)/(j-2)] \wedge minsd, j>2
    eta_j :=  ( [(d*eta + S_j)/(d+j-1) + f*u/sd_j]/(1+f/sd_j) \vee (m_j+e_j) ) \wedge u*(1-eps)

    Parameters
    ----------
    x: np.array
        input data
    attributes used:
        eta: float in (t, u) (default u*(1-eps))
            initial alternative hypothethesized value for the population mean
        c: positive float
            scale factor for allowing the estimated mean to approach t from above
        d: positive float
            relative weight of eta compared to an observation, in updating the alternative for each term
        f: positive float
            relative weight of the upper bound u (normalized by the sample standard deviation)
        minsd: positive float
            lower threshold for the standard deviation of the sample, to avoid divide-by-zero errors and
            to limit the weight of u

    */
        //         u = self.u
        //        N = self.N
        //        t = self.t
        //        eta = getattr(self, "eta", u * (1 - np.finfo(float).eps))
        //        c = getattr(self, "c", 1 / 2)
        //        d = getattr(self, "d", 100)
        //        f = getattr(self, "f", 0)
        //        minsd = getattr(self, "minsd", 10**-6)
        //        S, _Stot, j, m = self.sjm(N, t, x)
        //        # Welford's algorithm for running mean and running sd
        //        mj = [x[0]]
        //        sdj = [0]
        //        for i, xj in enumerate(x[1:]):
        //            mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
        //            sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
        //        sdj = np.sqrt(sdj / j)
        //        # end of Welford's algorithm.
        //        # threshold the sd, set first two sds to 1
        //        sdj = np.insert(np.maximum(sdj, minsd), 0, 1)[0:-1]
        //        sdj[1] = 1
        //        weighted = ((d * eta + S) / (d + j - 1) + u * f / sdj) / (1 + f / sdj)
        //        return np.minimum(
        //            u * (1 - np.finfo(float).eps),
        //            np.maximum(weighted, m + c / np.sqrt(d + j - 1)),
        //        )

        val eta = getattr(self, "eta", u * (1 - np.finfo(float).eps))
        val c = getattr(self, "c", 1 / 2)
        val d = getattr(self, "d", 100)
        val f = getattr(self, "f", 0)
        val minsd = getattr(self, "minsd", 10**-6)
        val (S, _Stot, j, m) = this.sjm(N, t, x)
        // Welford's algorithm for running mean and running sd
        val mj = [x[0]]
        var sdj = [0]
        for (i, xj in enumerate(x[1:])) {
            mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
            sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
        }
        sdj = np.sqrt(sdj / j)
        // end of Welford's algorithm.
        // threshold the sd, set first two sds to 1
        sdj = np.insert(np.maximum(sdj, minsd), 0, 1)[0:-1]
        sdj[1] = 1
        val weighted = ((d * eta + S) / (d + j - 1) + u * f / sdj) / (1 + f / sdj)
        return np.minimum(u * (1 - np.finfo(float).eps), np.maximum(weighted, m + c / np.sqrt(d + j - 1)))
    }

     */

    fun fixed_bet(x: DoubleArray): DoubleArray
    {
        /*
        Return a fixed value of lambda, the fraction of the current fortune to bet.

        Parameters
        ----------
        x: np.array
            input data

        Assumes the instance variable `lam` has been set.
        */
        return DoubleArray(x.size) { this.lam }
    }

    // TODO guess that x are floats. seperate with reps and without
    fun sample_size(x: DoubleArray, alpha: Double, reps: Int?, prefix: Boolean, quantile: Double): Int {
        /*
        Estimate the sample size to reject the null hypothesis that the population mean of a population of size
        `N` is `<=t` at significance level `alpha`, using pilot data `x`.

        If `reps is None`, tiles copies of `x` to produce a list of length `N`.

        If `reps is not None`, samples at random from `x` to produce `reps` lists of length `N` and reports
        the `quantile` quantile of the resulting sample sizes.

        If `prefix == True`, starts every list with `x` as given, then samples from `x` to produce the
        remaining `N-len(x)` entries.

        Parameters
        ----------
        x: the data for the simulation or calculation
        alpha: float in (0, 1/2); the significance level for the test
        reps: int; the number of replications to use to simulate the sample size. If `reps is None`, estimates deterministically
        prefix: bool; whether to use `x` as the prefix of the data sequence. Not used if `reps is None`
        quantile: float in (0, 1); desired quantile of the sample sizes from simulations. Not used if `reps is None`
        kwargs: keyword args passed to this.test()

        Returns
        -------
        sam_size: int
        estimated sample size
        */

        //         N = self.N
        //         N = self.N
        //        if reps is None:
        //            pop = np.repeat(np.array(x), math.ceil(N/len(x)))[0:N]  # tile data to make the population
        //            p = self.test(pop, **kwargs)[1]
        //            crossed = (p<=alpha)
        //            sam_size = int(N if np.sum(crossed)==0 else (np.argmax(crossed)+1))
        //        else:  # estimate the quantile by a bootstrap-like simulation
        //            seed = kwargs.get('seed',1234567890)
        //            prng = np.random.RandomState(seed)  # use the Mersenne Twister for speed
        //            sams = np.zeros(int(reps))
        //            pfx = np.array(x) if prefix else []
        //            ran_len = (N-len(x)) if prefix else N
        //            for r in range(reps):
        //                pop = np.append(pfx, prng.choice(x, size=ran_len, replace=True))
        //                p = self.test(pop, **kwargs)[1]
        //                crossed = (p<=alpha)
        //                sams[r] = N if np.sum(crossed)==0 else (np.argmax(crossed)+1)
        //            sam_size = int(np.quantile(sams, quantile))
        //        return sam_size

        /* TODO
        val N = this.N
        val sam_size = if (reps == null) {
            val pop = np.repeat(np.array(x), math.ceil(N / x.size))[0:N]  // tile data to make the population
            val p = this.test(pop, kwargs)[1]
            val crossed = (p <= alpha)
            int(N if np.sum(crossed) == 0 else (np.argmax(crossed) + 1))
        } else { // estimate the quantile by a bootstrap-like simulation
            val seed = kwargs.get("seed", 1234567890)
            val prng = np.random.RandomState(seed)  // use the Mersenne Twister for speed
            val sams = np.zeros(int(reps))
            val pfx = if (prefix) np.array(x) else emptyList()
            val ran_len = if (prefix) (N - x.size) else N
            for (r in range(reps)) {
                val pop = np.append(pfx, prng.choice(x, size = ran_len, replace = True))
                val p = this.test(pop, kwargs)[1]
                val crossed = (p <= alpha)
                sams[r] = if (np.sum(crossed) == 0) N else (np.argmax(crossed) + 1)
            }
            int(np.quantile(sams, quantile))
        }
        return sam_size

         */
        return 0
    }

        fun optimal_comparison(x: DoubleArray): Double {
            /*
        The value of eta corresponding to the "bet" that is optimal for ballot-level comparison audits,
        for which overstatement assorters take a small number of possible values and are concentrated
        on a single value when the CVRs have no errors.

        Let p0 be the rate of error-free CVRs, p1=0 the rate of 1-vote overstatements,
        and p2= 1-p0-p1 = 1-p0 the rate of 2-vote overstatements. Then

        eta = (1-u*p0)/(2-2*u) + u*p0 - 1/2, where p0 is the rate of error-free CVRs.

        Translating to p2=1-p0 gives:

        eta = (1-u*(1-p2))/(2-2*u) + u*(1-p2) - 1/2.

        Parameters
        ----------
        x: np.array
            input data
        rate_error_2: float
            hypothesized rate of two-vote overstatements

        Returns
        -------
        eta: float
            estimated alternative mean to use in alpha
        */
            // set the parameters
            // TO DO: double check where rate_error_2 is set
            val p2: Double = 1e-4 // getattr(self, "rate_error_2", 1e-4)  // rate of 2-vote overstatement errors
            return (1 - this.u * (1 - p2)) / (2 - 2 * this.u) + this.u * (1 - p2) - 1 / 2
        }

    companion object {
        val ulp =  java.lang.Math.ulp(1.0)
        val ulpm =  1.0 - ulp

        /*
        fun alpha_mart(x: DoubleArray, N: Int, t: Double = .5, eta: Double, u: Double = 1.0, estimFn: EstimatorFn?):
                Pair<Double, Array<Int>> {

            /*
                Finds the ALPHA martingale for the hypothesis that the population
                mean is less than or equal to mu using a martingale method,
                for a population of size N, based on a series of draws x.

                **The draws must be in random order**, or the sequence is not a supermartingale under the null

                If N is finite, assumes the sample is drawn without replacement
                If N is infinite, assumes the sample is with replacement

                Parameters
                ----------
                x : list corresponding to the data
                N : int; population size for sampling without replacement, or np.infinity for sampling with replacement
                t : float in [0,u); hypothesized fraction of ones in the population
                eta : float in (mu,u]; alternative hypothesized population mean
                estim : function (note: class methods are not of type Callable)
                        estim(x, N, mu, eta, u) -> np.array of length len(x), the sequence of values of eta_j for ALPHA

                Returns
                -------
                p : float; sequentially valid p-value of the hypothesis that the population mean is less than or equal to mu
                p_history : numpy array; sample by sample history of p-values. Not meaningful unless the sample is in random order.
            */
            //         if not estim:
            //            estim = TestNonnegMean.shrink_trunc
            //        S = np.insert(np.cumsum(x),0,0)[0:-1]  # 0, x_1, x_1+x_2, ...,
            //        j = np.arange(1,len(x)+1)              # 1, 2, 3, ..., len(x)
            //        m = (N*t-S)/(N-j+1) if np.isfinite(N) else t   # mean of population after (j-1)st draw, if null is true
            //        etaj = estim(x, N, t, eta, u)
            //        x = np.array(x)
            //        with np.errstate(divide='ignore',invalid='ignore'):
            //            terms = np.cumprod((x*etaj/m + (u-x)*(u-etaj)/(u-m))/u)
            //        terms[m<0] = np.inf                                         # true mean certainly greater than hypothesized
            //        terms[m>u] = 1                                              # true mean certainly less than hypothesized
            //        terms[np.isclose(0, m, atol=2*np.finfo(float).eps)] = 1     # ignore
            //        terms[np.isclose(u, m, atol=10**-8, rtol=10**-6)] = 1       # ignore
            //        terms[np.isclose(0, terms, atol=2*np.finfo(float).eps)] = 1 # martingale effectively vanishes; p-value 1
            //        return min(1, 1/np.max(terms)), np.minimum(1,1/terms)

            val estim: EstimatorFn = estimFn ?: { x -> shrink_trunc(x) }

            val cum_sum = numpy_cumsum(x)
            val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
            val Stot = S[S.size-1]  // sample total ""array[-1] means the last element"
            val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
            // mean of population after (j-1)st draw, if null is true
            val m = DoubleArray(x.size) { if (!withReplacement) t else (N * t - S[it]) / (N - j[it] + 1)  }

            val etaj = estim(x, N, t, eta, u)
            val x = np.array(x)
            with(np.errstate(divide = 'ignore', invalid = 'ignore')) {
                terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            }
            terms[m < 0] =
                np.inf                                         // true mean certainly greater than hypothesized
            terms[m > u] = 1                                              // true mean certainly less than hypothesized
            terms[np.isclose(0, m, atol = 2 * np.finfo(float).eps)] = 1     // ignore
            terms[np.isclose(u, m, atol = 10 * * - 8, rtol = 10 * * - 6)] = 1       // ignore
            terms[np.isclose(0, terms, atol = 2 * np.finfo(float).eps)] =
                1 // martingale effectively vanishes; p - value 1

            return Pair(min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms))
        }

         */
    }
}


// Return the cumulative product of elements
fun numpy_cumprod(a: DoubleArray) : DoubleArray {
    val result = DoubleArray(a.size)
    result[0] = a[0]
    for (i in 1 until a.size) {
        result[i] = result[i-1] * a[i]
    }
    return result
}

// Return the cumulative product of elements
fun numpy_cumsum(a: DoubleArray) : DoubleArray {
    val result = DoubleArray(a.size)
    result[0] = a[0]
    for (i in 1 until a.size) {
        result[i] = result[i-1] + a[i]
    }
    return result
}

// def isclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
fun numpy_isclose(a: Double, b: Double, rtol: Double=1.0e-5, atol:Double=1.0e-8): Boolean {
    //    For finite values, isclose uses the following equation to test whether
    //    two floating point values are equivalent.
    //
    //     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))
    return abs(a - b) <= atol + rtol * abs(b)
}

// Return an array of ones with the same shape and type as a given array.
fun numpy_ones_like(a: DoubleArray): DoubleArray {
    return DoubleArray(a.size) { 1.0 }
}
