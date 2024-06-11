package org.cryptobiotic.shangrla

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
    val test: callable?,
    val estim: callable?,
    val bet: callable?,
    var u: Float = 1.0f,
    val N: Int,
    val t: Float = 0.5f,
    val random_order: Boolean = true,
    val kwargs: Map<String, Any>,
) {

    // TODO guess that x are floats. seperate with reps and without
    fun sample_size(x: FloatArray, alpha: Float, reps: Int?, prefix: Boolean, quantile: Float): Int {
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
        x: list or np.array; the data for the simulation or calculation
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

    companion object {

        // TODO eta: Float=1-np.finfo(float).eps
        fun alpha_mart(x: Array<Int>, N: Int, t: Float=.5f, eta: Float, u: Float=1.0f, estim: callable=None):
                Pair<Float, Array<Int>> {

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

            if (!estim) {
                estim = TestNonnegMean.shrink_trunc
            }
            val S = np.insert(np.cumsum(x), 0, 0)[0:-1]  # 0, x_1, x_1+x_2, ...,
            val j = np.arange(1, len(x) + 1)              # 1, 2, 3, ..., len(x)
            val m =
                (N * t - S) / (N - j + 1) if np.isfinite(N) else t   # mean of population after (j-1)st draw, if null is true
            val etaj = estim(x, N, t, eta, u)
            val x = np.array(x)
            with (np . errstate (divide = 'ignore', invalid = 'ignore')) {
                terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            }
            terms[m < 0] =
                np.inf                                         # true mean certainly greater than hypothesized
            terms[m > u] = 1                                              # true mean certainly less than hypothesized
            terms[np.isclose(0, m, atol = 2 * np.finfo(float).eps)] = 1     # ignore
            terms[np.isclose(u, m, atol = 10 * * - 8, rtol = 10 * * - 6)] = 1       # ignore
            terms[np.isclose(0, terms, atol = 2 * np.finfo(float).eps)] =
                1 # martingale effectively vanishes; p - value 1

            return Pair(min(1, 1 / np.max(terms)), np.minimum(1, 1/terms))
        }
    }
}
