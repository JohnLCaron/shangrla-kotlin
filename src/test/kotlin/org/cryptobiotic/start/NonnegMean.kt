package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.*
import org.cryptobiotic.shangrla.core.NonnegMean.CumulativeSum
import kotlin.math.ceil
import kotlin.math.min

// Tests of the hypothesis that the mean of a population of values in [0, u] is less than or equal to t
class NonnegMean (
        val N: Int,
        var u: Double = 1.0,        // TODO mutable
        val t: Double = 0.5,        // TODO is it ever anything other than .5 ??
        val withReplacement: Boolean = false,

        testFnType: TestFnType = TestFnType.ALPHA_MART,
        estimFnType: EstimFnType = EstimFnType.OPTIMAL,
        gOverride: Double = 0.1, // only used by kaplan_wald
        random_order: Boolean = true, // only used by kaplan_wald
    ) {
    val eta = (t + (u - t) / 2)
    val estimFn: EstimatorFn
    val testFn: TestFn

    init {
        val testAlpha: TestFn = { x -> this.alpha_mart(x) }
        val testWald: TestFn = { x -> this.kaplan_wald(x, gOverride, random_order) }
        testFn  = if (testFnType == TestFnType.KAPLAN_WALD) testWald else testAlpha

        // TODO only place this is used is in alpha_mart(), betFn in net_mart, so could be added to the testFn and passed in.
        val estimFixed: EstimatorFn = { x -> this.fixed_alternative_mean(x) }
        val estimOptimal: EstimatorFn = { x -> this.optimal_comparison() } // TODO pass error_2_rate ?
        estimFn = if (estimFnType == EstimFnType.OPTIMAL) estimOptimal else estimFixed
    }

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
        fun sample_size(x: DoubleArray, alpha: Double, reps: Int? = null, prefix: Boolean, quantile: Double): Int {

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

        require (quantile in 0.0..1.0)

        val sam_size = if (reps == null) {
            val repeats = ceil(N.toDouble() / x.size).toInt()
            val pop = numpy_repeat(x, repeats) // tile data to make the population
            require( pop.size == N)
            val (p, p_history) = this.testFn(pop)
//            crossed = p <= alpha : Array of true or false
//            sum = np.sum(crossed) : count of true
//            argmax = np.argmax(crossed): ?? the first index that has true
//            sam_size = int(N if sum == 0 else (argmax + 1))
            // int(N if np.sum(crossed) == 0 else (np.argmax(crossed) + 1))
            // i guess you could just find the first index thats true.
            val firstIndex = p_history.indexOfFirst{ it <= alpha }
            if (firstIndex < 0) N else firstIndex + 1

        } else { // estimate the quantile by a bootstrap-like simulation
            //            seed = kwargs.get("seed", 1234567890)
            //            prng = np.random.RandomState(seed)  # use the Mersenne Twister for speed
            //            sams = np.zeros(int(reps))
            //            pfx = np.array(x) if prefix else []
            //            ran_len = (N - len(x)) if prefix else N
            //            for r in range(reps):
            //                pop = np.append(pfx, prng.choice(x, size=ran_len, replace=True))
            //                p = self.test(pop, **kwargs)[1]
            //                crossed = p <= alpha
            //                sams[r] = N if np.sum(crossed) == 0 else (np.argmax(crossed) + 1)
            //            sam_size = int(np.quantile(sams, quantile))
            //        return sam_size
            val sams = IntArray(reps)
            val pfx = if (prefix) x else DoubleArray(0)
            val ran_len = if (prefix) (N - x.size) else N
            repeat(reps) {
                val choices = python_choice(x, size = ran_len)
                val pop = numpy_append(pfx, choices) // tile data to make the population
                val (p, p_history) = this.testFn(pop)
                val crossed = p_history.map{ it <= alpha }
                val crossedCount = crossed.filter { it }.count()
                sams[it] = if (crossedCount == 0) N else (indexFirstTrue(crossed) + 1)
            }
            sams.sort() // sort in place
            numpy_quantile(sams, quantile)
        }
        return sam_size
    }

    // testFn
    // return (p, p_history)
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
        val etaj = this.estimFn(x)
        val xp = if (etaj.size == x.size) DoubleArray(x.size) { (x[it] * etaj[it] / m [it] + (u - x[it]) * (u - etaj[it]) / (u - m[it])) / u }
        else if (etaj.size == 1) DoubleArray(x.size) { (x[it] * etaj[0] / m [it] + (u - x[it]) * (u - etaj[0]) / (u - m[it])) / u }
        else throw RuntimeException("NonnegMean alpha_mart")

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
        if (Stot > N * t) terms[terms.size-1] = Double.POSITIVE_INFINITY // final sample makes the total greater than the null

        // return min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms)
        // np.minimum = element-wise minumum, presumably the smaller of 1 and 1/term
        // np.max = maximum of an array
        // min = min or an iterable

        val npmin = terms.map{ min(1.0, 1.0 / it) }.toDoubleArray()
        val npmax = terms.max()
        val first = min(1.0, 1.0 / npmax)
        return Pair(first, npmin)
    }

    fun sjm(N: Int, t: Double, x: DoubleArray): CumulativeSum {
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
        S    : np.array; The cumulative sum of the input array `x`, excluding the last element.
        Stot : float; The total sum of the input array `x`.
        j    : np.array; An array of indices from 1 to the length of `x`.
        m    : np.array; The mean of the population after each draw if the null hypothesis is true.
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
        val Stot = S.last()  // sample total ""array[-1] means the last element"
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

//        j = np.arange(1, len(x) + 1)  # 1, 2, 3, ..., len(x)
//        assert j[-1] <= N, "Sample size is larger than the population!"
        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        require( withReplacement || x.size <= N) { "Sample size is larger than the population!" }
//        m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true (t=eta is the mean)
        // val m = if (withReplacement) doubleArrayOf(t) else DoubleArray(x.size) { (N * t - Sp[it]) / (N - j[it] + 1)  }
        val m = DoubleArray(x.size) { if (withReplacement) t else (N * t - Sp[it]) / (N - j[it] + 1)  }
        return CumulativeSum(Sp, Stot, j, m)
    }

    // estimFn
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
        // TODO eta based on ulp
        //        ulp = np.finfo(float).eps
        //        uulp = (1 - ulp)
        //        eta = getattr(self, "eta", u * (1 - np.finfo(float).eps))
        //        _S, _Stot, _j, m = self.sjm(N, eta, x)
        //        if (negs := np.sum(m < 0)) > 0:
        //            warnings.warn(
        //                f"Implied population mean is negative in {negs} of {len(x)} terms"
        //            )
        //        if (pos := np.sum(m > u)) > 0:
        //            warnings.warn(
        //                f"Implied population mean is greater than {u=} in {pos} of {len(x)} terms"
        //            )
        //        return m

        val (_, _, _, m) = this.sjm(N, eta, x)
        val negs = m.filter { it < 0.0 }
        if (negs.count() > 0) {
            println("Implied population mean is negative in ${negs.size} of ${x.size} terms")
        }
        val pos = m.filter { it > u }
        if (pos.count() > 0) {
            println("Implied population mean is greater than ${u} in ${pos.size} of ${x.size} terms")
        }
        return m
    }

    fun optimal_comparison(rate_error_2: Double = 1e-4): DoubleArray {
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
        // TODO: where is rate_error_2 set?
        if (this.u == 1.0)
            throw RuntimeException("optimal_comparison: u ${this.u} must be < 1")
        val p2 = rate_error_2 // getattr(self, "rate_error_2", 1e-4)  // rate of 2-vote overstatement errors
        val result = (1 - this.u * (1 - p2)) / (2 - 2 * this.u) + this.u * (1 - p2) - .5
        return doubleArrayOf(result)
    }

    fun kaplan_wald(x: DoubleArray, g : Double = 0.1, random_order: Boolean = true): Pair<Double, DoubleArray> {
        /*
        Kaplan-Wald p-value for the hypothesis that the sample x is drawn IID from a population
        with mean t against the alternative that the mean is less than t.

        If there is a possibility that x has elements equal to zero, set g \in (0, 1);
        otherwise, the p-value will be 1.

        If the order of the values in the sample is random, you can set random_order = True to use
        optional stopping to increase the power. If the values are not in random order or if you want
        to use all the data, set random_order = False

        Parameters:
        -----------
        x: array-like the sample
        g: "padding" in case there any values in the sample are zero
        random_order: Boolean
                if the sample is in random order, it is legitimate to stop early, which
                can yield a more powerful test. See above.

        Returns:
        --------
        p: float p-value
        p_history: sample by sample history of p-values. Not meaningful unless the sample is in random order.

        */
        if (g < 0.0 || g > 1.0) {
            throw Exception ("{$g}, but g must be between 0 and 1. ")
        }
        x.forEach { if( it < 0.0) throw Exception ("Negative value in sample from a nonnegative population.") }
        val xp = x.map { (1 - g) * it / t + g }.toDoubleArray()
        val p_history = numpy_cumprod(xp)

        // return np.min([1, 1 / np.max(p_history) if random_order else 1 / p_history[-1]]), np.minimum(1 / p_history, 1)
        val npmin = p_history.map{ min(1.0, 1.0 / it) }.toDoubleArray()
        val term = if (random_order) p_history.max() else (1.0 / p_history.last())
        val first = min (1.0, term)
        return Pair(first, npmin)

    }

    companion object {
        val PLURALITY = NonnegMean()
    }
}