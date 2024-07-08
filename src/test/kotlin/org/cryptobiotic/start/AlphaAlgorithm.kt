package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.numpy_cumsum
import org.cryptobiotic.start.AlphaMart.CumulativeSum
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

//         return AlgoValues(etajsa, populationMeanValues, tjs, tjs
data class AlgoValues(val etaj: List<Double>, val populationMeanValues: List<Double>, val tjs: List<Double>, val tstat: List<Double>)


// ALPHA paper, section 3
class AlphaAlgorithm(
    val risk_limit: Double = 0.05, // α ∈ (0, 1)
    val withoutReplacement: Boolean = false,
    val upperBound: Double = 1.0,  // aka u
    val N: Int, // number of ballot cards in the population of cards from which the sample is drawn
                // If N is np.inf, it means the sampling is with replacement

    val isPolling: Boolean = false,
    val d: Double,
    val t: Double = 0.5,        // the hypothesized mean "under the null".
) {
    var eta0 = setEta0()
    val c = (eta0 - upperBound / 2)

    init {
        //• Set audit parameters:
        // – Select the risk limit α ∈ (0, 1)
        // - decide whether to sample with or without replacement.
        // – Set u as appropriate for the assertion under audit.
        // – Set N to the number of ballot cards in the population of cards from which the sample is drawn.
        // – Set η0
        //    For polling audits, η0 could be the reported mean value of the assorter.
        //	    For instance, for the assertion corresponding to checking whether w got more votes than ℓ,
        //	      η0 = (Nw + Nc /2)/N , where Nw is the number of votes reported for w , Nℓ is the
        //	   number of votes reported for ℓ, and Nc = N − Nw − Nℓ is the number of ballot cards
        //	   reported to have a vote for some other candidate or no valid vote in the contest.
        //    For comparison audits, η0 can be based on assumed or historical rates of overstatement errors.
        //
        // – Define the function to update η based on the sample,
        //	  e.g, η(i, X i−1 ) = ((d * η0 + S)/(d + i − 1) ∨ (eps(i) + µi )) ∧ u,    (2.5.2, eq 14, "truncated shrinkage")
        //	    where S = Sum i−1 k=1 (Xk) is the sample sum of the first i − 1 draws
        //	    and eps(i) = c/ sqrt(d + i − 1)
        //	  set any free parameters in the function (e.g., d and c in this example). The only requirement is that
        //	 η(i, X i−1 ) ∈ (µi , u), where µi := E(Xi |X i−1 ) is computed under the null.
    }

    fun run2(maxSample: Int, drawSample : () -> Double) : AlgoValues {

        // Initialize variables
        // – j ← 0: sample number
        // – T ← 1: test statistic
        // – S ← 0: sample sum
        // – m = µ_j = 1/2: population mean under the null

        var sampleNumber = 0
        var testStatistic = 1.0
        var sampleSum = 0.0
        var populationMean = .5

        // Let
        //  theta = 1/N * Sum(xj)
        //  xi = 1/2 otherwise. The reported winner really beat the reported loser if
        //  xi in (0,upperBound) // the sampled assort values
        //  X^j := (x1, x2, ... xj)
        //  µ_j := E(X_j | X^j−1 ) assuming the null hypothesis
        //  eta_j := eta_j(X^j−1) estimate of the population mean, (depends on X^j−1, but not on x_k for k ≥ j)

        val sampleAssortValues = mutableListOf<Double>()
        val populationMeanValues = mutableListOf<Double>()
        val etajsa = mutableListOf<Double>()
        val tjs = mutableListOf<Double>()
        val tstat = mutableListOf<Double>()

        //
        // While T < 1/α and not all ballot cards have been audited:
        // while (testStatistic < 1/risk_limit && sampleNumber < maxSample) {
        while (sampleNumber < maxSample) {

            // – Draw a ballot at random
            // – Determine Xj by applying the assorter to the selected ballot card (and the CVR, for comparison audits)
            val xj: Double = drawSample()
            sampleAssortValues.add(xj)
            val sampleMean = sampleAssortValues.average()
            sampleNumber++
            sampleSum += xj // – S ← S + Xj

            val m = this.sjm2(N, t, sampleAssortValues.toDoubleArray(), sampleNumber)
            populationMeanValues.add(m)

            val eta = (t + (upperBound - t) / 2) // initial estimate of the population mean
            val etajs = this.fixedAlternativeMean(sampleAssortValues.toDoubleArray(), eta) // estimFixed returns single "fixed" value; estimOptimal an array of size x
            val etaj = etajs.last()
            etajsa.add(etaj)

            // estimate the populationMean
            //     fun updateEta(c: Double, d: Double, eta0: Double, sumkm1: Double, sampleNum: Int, populationMean: Double, upperBound: Double) : Double {
            // val etaj : Double = truncShrinkage(c=c, d=d, eta0=eta0, sumkm1=sampleSum, sampleNum=sampleNumber) // , populationMean=populationMean, upperBound=upperBound)

            // – If m < 0, T ← ∞. Otherwise, T ← T / u * ( Xj * η(j,S)/m + (u - Xj) * (u−η(j,S))/(u-m))
            // TODO This is eq 4 of ALPHA, p.5 :
            //      T_j = T_j-1 * (X_j * eta_j / µ_j + (u - X_j) * (u - eta_j) / ( u - µ_j)) / u
            // python:  terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            // python: m = µ_j = The mean of the population after each draw if the null hypothesis is true. TODO mean of the sample?
            //          m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            //        draw, if null is true (t=eta is the mean)
            // u = upperBound
            // m = The mean of the population after each draw if the null hypothesis is true.

            //val m2 = (N * t - sampleSum)
            //val m3 = (N - sampleNumber + 1)
            //val m4 = m2 / m3

            //val epsi = c / sqrt(d + sampleNumber - 1)
            //val etajBounded = min( max(epsi + populationMean, etaj), upperBound)

            val tj = if (m < 0.0) Double.POSITIVE_INFINITY else {
                (xj * etaj / m + (upperBound - xj) * (upperBound - etaj) / (upperBound - m))/ upperBound
            }
            tjs.add(tj)
            testStatistic *= tj
            println(" $sampleNumber = $xj etaj = $etaj tj=$tj, T = $testStatistic")
            tstat.add(testStatistic)

            // update the population mean
            // – If the sample is drawn without replacement, m ← (N/2 − S)/(N − j + 1)
            //   m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            if (withoutReplacement) {
                populationMean = (.5*N - sampleSum)/(N - sampleNumber + 1)
                if (populationMean < 0.0)
                    println("wtf")
            }

            // – If desired, break and conduct a full hand count instead of continuing to audit.
        }
        println(" etajs ${etajsa}")
        println(" populationMeanValues ${populationMeanValues}")
        println(" tjs ${tjs}")
        println(" tstat ${tstat}")

        //
        //• If a full hand count is conducted, its results replace the reported results if they differ.
        return AlgoValues(etajsa, populationMeanValues, tjs, tstat)
    }

    // sampleNum starts at 1
    fun sjm2(N: Int, t: Double, x: DoubleArray, sampleNum: Int): Double {
        val cum_sum = numpy_cumsum(x)
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

        return if (!withoutReplacement) t else {
            (N * t - Sp[sampleNum]) / (N - sampleNum)
        }
    }

    fun run(maxSample: Int, drawSample : () -> Double) : AlgoValues {

        // Initialize variables
        // – j ← 0: sample number
        // – T ← 1: test statistic
        // – S ← 0: sample sum
        // – m = µ_j = 1/2: population mean under the null

        var sampleNumber = 0
        var testStatistic = 1.0
        var sampleSum = 0.0
        var populationMean = .5

        // Let
        //  theta = 1/N * Sum(xj)
        //  xi = 1/2 otherwise. The reported winner really beat the reported loser if
        //  xi in (0,upperBound) // the sampled assort values
        //  X^j := (x1, x2, ... xj)
        //  µ_j := E(X_j | X^j−1 ) assuming the null hypothesis
        //  eta_j := eta_j(X^j−1) estimate of the population mean, (depends on X^j−1, but not on x_k for k ≥ j)

        val sampleAssortValues = mutableListOf<Double>()
        val populationMeanValues = mutableListOf<Double>()
        val etajsa = mutableListOf<Double>()
        val tjs = mutableListOf<Double>()
        val tstat = mutableListOf<Double>()

        //
        // While T < 1/α and not all ballot cards have been audited:
        // while (testStatistic < 1/risk_limit && sampleNumber < maxSample) {
        while (sampleNumber < maxSample) {

            // – Draw a ballot at random
            // – Determine Xj by applying the assorter to the selected ballot card (and the CVR, for comparison audits)
            val xj: Double = drawSample()
            sampleAssortValues.add(xj)
            val sampleMean = sampleAssortValues.average()
            sampleNumber++
            sampleSum += xj // – S ← S + Xj

            val (_, Stot, _, ms) = this.sjm(N, t, sampleAssortValues.toDoubleArray())
            val m = ms.last()
            populationMeanValues.add(m)

            val eta = (t + (upperBound - t) / 2) // initial estimate of the population mean
            val etajs = this.fixedAlternativeMean(sampleAssortValues.toDoubleArray(), eta) // estimFixed returns single "fixed" value; estimOptimal an array of size x
            val etaj = etajs.last()
            etajsa.add(etaj)

            // estimate the populationMean
            //     fun updateEta(c: Double, d: Double, eta0: Double, sumkm1: Double, sampleNum: Int, populationMean: Double, upperBound: Double) : Double {
            // val etaj : Double = truncShrinkage(c=c, d=d, eta0=eta0, sumkm1=sampleSum, sampleNum=sampleNumber) // , populationMean=populationMean, upperBound=upperBound)

            // – If m < 0, T ← ∞. Otherwise, T ← T / u * ( Xj * η(j,S)/m + (u - Xj) * (u−η(j,S))/(u-m))
            // TODO This is eq 4 of ALPHA, p.5 :
            //      T_j = T_j-1 * (X_j * eta_j / µ_j + (u - X_j) * (u - eta_j) / ( u - µ_j)) / u
            // python:  terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            // python: m = µ_j = The mean of the population after each draw if the null hypothesis is true. TODO mean of the sample?
            //          m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            //        draw, if null is true (t=eta is the mean)
            // u = upperBound
            // m = The mean of the population after each draw if the null hypothesis is true.

            //val m2 = (N * t - sampleSum)
            //val m3 = (N - sampleNumber + 1)
            //val m4 = m2 / m3

            //val epsi = c / sqrt(d + sampleNumber - 1)
            //val etajBounded = min( max(epsi + populationMean, etaj), upperBound)

            val tj = if (m < 0.0) Double.POSITIVE_INFINITY else {
                (xj * etaj / m + (upperBound - xj) * (upperBound - etaj) / (upperBound - m))/ upperBound
            }
            tjs.add(tj)
            testStatistic *= tj
            println(" $sampleNumber = $xj etaj = $etaj tj=$tj, T = $testStatistic")
            tstat.add(testStatistic)

            // update the population mean
            // – If the sample is drawn without replacement, m ← (N/2 − S)/(N − j + 1)
            //   m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            if (withoutReplacement) {
                populationMean = (.5*N - sampleSum)/(N - sampleNumber + 1)
                if (populationMean < 0.0)
                    println("wtf")
            }

            // – If desired, break and conduct a full hand count instead of continuing to audit.
        }
        println(" etajs ${etajsa}")
        println(" populationMeanValues ${populationMeanValues}")
        println(" tjs ${tjs}")
        println(" tstat ${tstat}")

        //
        //• If a full hand count is conducted, its results replace the reported results if they differ.
        return AlgoValues(etajsa, populationMeanValues, tjs, tstat)
    }

    fun setEta0() : Double {
        //    For polling audits, eta0 could be the reported mean value of the assorter.
        //	    For instance, for the assertion corresponding to checking whether w got more votes than ℓ,
        //	      η0 = (Nw + Nc/2)/N , where Nw is the number of votes reported for w , Nℓ is the
        //	   number of votes reported for ℓ, and Nc = N − Nw − Nℓ is the number of ballot cards
        //	   reported to have a vote for some other candidate or no valid vote in the contest.

        //    For comparison audits, eta0 can be based on assumed or historical rates of overstatement errors.
        val eps = 0.0001  // Generic small value
        val eta0 = (eps + (upperBound - eps) / 2) // initial estimate of the population mean
        return eta0
    }

    fun truncShrinkage(c: Double, d: Double, eta0: Double, sumkm1: Double, sampleNum: Int) : Double {
        // – Define the function to update eta based on the sample,
        //	  e.g, eta(i, X^i−1 ) = ((d * η0 + S)/(d + i − 1) ∨ (eps(i) + µi )) ∧ u,    (2.5.2, eq 14, "truncated shrinkage")
        //	    where S = Sum(Xk) k=1..i-1  is the sample sum of the first i-1 draws
        //	    and eps(i) = c / sqrt(d + i − 1)

        // 2.5.2 Choosing ǫi . To allow the estimated winner’s share ηi to approach
        //√ µi as the sample grows
        //(if the sample mean approaches µi or less), we shall take ǫi := c/ d + i − 1 for a nonnega-
        //tive constant c, for instance c = (η0 − µ)/2. The estimate ηi is thus the sample mean, shrunk

        var result = ((d * eta0 + sumkm1) / (d + sampleNum - 1)) //  ∨ (eps(i) + µi )) ∧ u

        // The only requirement is that eta(i, X^i−1 ) ∈ (µi .. u), where µi := E(Xi |X i−1 ) is computed under the null.
        // result = min( max(epsi + populationMean, result), upperBound)
        return result
    }

    fun sjm(N: Int, t: Double, x: DoubleArray): CumulativeSum {
        val cum_sum = numpy_cumsum(x)
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Stot = S.last()  // sample total ""array[-1] means the last element"
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

//        j = np.arange(1, len(x) + 1)  # 1, 2, 3, ..., len(x)
//        assert j[-1] <= N, "Sample size is larger than the population!"
        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        require( !withoutReplacement || x.size <= N) { "Sample size is larger than the population!" }
//        m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true (t=eta is the mean)
        // val m = if (withReplacement) doubleArrayOf(t) else DoubleArray(x.size) { (N * t - Sp[it]) / (N - j[it] + 1)  }
        val m = DoubleArray(x.size) {
            val m1 = (N * t - Sp[it])
            val m2 = (N - j[it] + 1)
            m1 / m2
            // if (!withoutReplacement) t else (N * t - Sp[it]) / (N - j[it] + 1)
        }
        return CumulativeSum(Sp, Stot, j, m)
    }

    fun sjmFromFixed(N: Int, t: Double, x: DoubleArray): CumulativeSum {
        val cum_sum = numpy_cumsum(x)
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Stot = S.last()  // sample total ""array[-1] means the last element"
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

//        j = np.arange(1, len(x) + 1)  # 1, 2, 3, ..., len(x)
//        assert j[-1] <= N, "Sample size is larger than the population!"
        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        require( !withoutReplacement || x.size <= N) { "Sample size is larger than the population!" }
//        m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true (t=eta is the mean)
        // val m = if (withReplacement) doubleArrayOf(t) else DoubleArray(x.size) { (N * t - Sp[it]) / (N - j[it] + 1)  }
        val m = DoubleArray(x.size) {
            val m1 = (N * t - Sp[it])
            val m2 = (N - j[it] + 1)
            val m3 = m1 / m2
            m3
            // if (!withoutReplacement) t else (N * t - Sp[it]) / (N - j[it] + 1)
        }
        return CumulativeSum(Sp, Stot, j, m)
    }

    fun fixedAlternativeMean(x: DoubleArray, eta: Double): DoubleArray {
        val (_, _, _, m) = this.sjmFromFixed(N, eta, x)
        val negs = m.filter { it < 0.0 }
        if (negs.count() > 0) {
            println("Implied population mean is negative in ${negs.size} of ${x.size} terms")
        }
        val pos = m.filter { it > upperBound}
        if (pos.count() > 0) {
            println("Implied population mean is greater than ${upperBound} in ${pos.size} of ${x.size} terms")
        }
        return m
    }
}
