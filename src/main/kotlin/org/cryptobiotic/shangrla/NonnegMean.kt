package org.cryptobiotic.shangrla

class NonnegMean (
    //test: callable=None,
    //estim: callable=None,
    //bet: callable=None,
    var u: Float = 1.0f,
    val N: Int,
    t: Float = 0.5f,
    random_order: Boolean = true,
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
}
