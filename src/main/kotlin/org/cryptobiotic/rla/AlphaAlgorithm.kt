package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.numpy_cumsum
import org.cryptobiotic.shangrla.core.numpy_isclose
import kotlin.math.min

//         return AlgoValues(etajsa, populationMeanValues, tjs, tjs
data class AlgoValues(val etaj: List<Double>, val populationMeanValues: List<Double>, val tjs: List<Double>,
                      val tstat: List<Double>, val phistory: List<Double>)

interface EstimFn {
    // return eta estimate
    fun eta(j: Int, sampleSum: Double): Double
}

// ALPHA paper, section 3, p.9
class AlphaAlgorithm(
    val estimFn : EstimFn,
    val N: Int, // number of ballot cards in the population of cards from which the sample is drawn
                // (old) "If N is np.inf, it means the sampling is with replacement"
    val withoutReplacement: Boolean = true,
    val risk_limit: Double = 0.05, // α ∈ (0, 1)
    val upperBound: Double = 1.0,  // aka u
    val t: Double = 0.5,        // the hypothesized mean "under the null". TODO is it ever not 1/2 ?
    val isPolling: Boolean = false,
) {
    //var eta0 = setEta0()
    //val c = (eta0 - upperBound / 2)

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
        val sampleMeanValues = mutableListOf<Double>()
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
            // val sampleMean = sampleAssortValues.average()
            sampleNumber++

            //// use previous sample sum
            // population mean if theta = t
            populationMean = this.populationMean(sampleNumber, sampleSum)
            sampleMeanValues.add(populationMean)

            // fixed alternative that the original population mean is eta
            val eta = (t + (upperBound - t) / 2)
            val etaj = this.alternativeFixedMean(eta, sampleNumber, sampleSum)
            etajsa.add(etaj)

            sampleSum += xj // – S ← S + Xj

            // estimate of the populationMean
            // – If m < 0, T ← ∞. Otherwise, T ← T / u * ( Xj * η(j,S)/m + (u - Xj) * (u−η(j,S))/(u-m))
            // TODO This is eq 4 of ALPHA, p.5 :
            //      T_j = T_j-1 * (X_j * eta_j / µ_j + (u - X_j) * (u - eta_j) / ( u - µ_j)) / u
            // python:  terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            // python: m = µ_j = The mean of the population after each draw if the null hypothesis is true.
            //          m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            //        draw, if null is true
            // u = upperBound
            // m = The mean of the population after each draw if the null hypothesis is true.

            val tj = if (populationMean < 0.0) Double.POSITIVE_INFINITY else {
                (xj * etaj / populationMean + (upperBound - xj) * (upperBound - etaj) / (upperBound - populationMean))/ upperBound
            }
            tjs.add(tj)
            testStatistic *= tj
            println(" $sampleNumber = $xj etaj = $etaj tj=$tj, T = $testStatistic")
            tstat.add(testStatistic)

            /* update the population mean
            // – If the sample is drawn without replacement, m ← (N/2 − S)/(N − j + 1)
            //   m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st
            if (withoutReplacement) {
                populationMean = (.5*N - sampleSum)/(N - sampleNumber + 1)
                if (populationMean < 0.0)
                    println("wtf")
            }

             */

            // – If desired, break and conduct a full hand count instead of continuing to audit.
        }
        println(" etajs ${etajsa}")
        println(" populationMeanValues ${sampleMeanValues}")
        println(" tjs ${tjs}")

        //// exceptional conditions
        // terms[m > u] = 0  # true mean is certainly less than hypothesized
        repeat(tstat.size) { if (sampleMeanValues[it] > upperBound) tstat[it] = 0.0 } // true mean is certainly less than hypothesized

        // terms[np.isclose(0, m, atol=atol)] = 1  # ignore
        repeat(tstat.size) { if (numpy_isclose(0.0, sampleMeanValues[it])) tstat[it] = 1.0 } // ignore

        // terms[np.isclose(0, terms, atol=atol)] = ( 1 } # martingale effectively vanishes; p-value 1
        repeat(tstat.size) { if (numpy_isclose(0.0, tstat[it])) tstat[it] = 1.0 } // martingale effectively vanishes; p-value 1

        // terms[m < 0] = np.inf  # true mean certainly greater than hypothesized
        repeat(tstat.size) { if (sampleMeanValues[it] < 0.0) tstat[it] = Double.POSITIVE_INFINITY } // true mean is certainly less than hypothesized

        // terms[-1] = ( np.inf if Stot > N * t else terms[-1] )  # final sample makes the total greater than the null
        if (sampleSum > N * t) tstat[tstat.size - 1] = Double.POSITIVE_INFINITY // final sample makes the total greater than the null

        println(" tstat ${tstat}")

        // return min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms)
        val phistory = tstat.map { min(1.0, 1.0 / it) }
        // val p = min(1.0, 1.0 / tstat.max()) // seems wrong
        println(" phistory ${phistory}")

        //
        //• If a full hand count is conducted, its results replace the reported results if they differ.
        return AlgoValues(etajsa, sampleMeanValues, tjs, tstat, phistory)
    }

    // Compute the alternative mean just before the jth draw, for a fixed alternative that the original population mean is eta. ??
    // sampleNum starts at 1
    fun alternativeFixedMean(t: Double, sampleNum: Int, sampleSumM1: Double): Double {
        return (N * t - sampleSumM1) / (N - sampleNum + 1)
    }

    // population mean under the null hypothesis that “the average of this list is not greater than 1/2”
    // TODO seems to be the "mean of the remaining sample", which is why it depends on whether you are replacing or not
    // TODO detect if it goes negetive
    fun populationMean(sampleNum: Int, sampleSumM1: Double): Double {
        return if (withoutReplacement) (N * t - sampleSumM1) / (N - sampleNum + 1) else t
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

    data class CumulativeSum(val S: DoubleArray, val Stot: Double, val indices: IntArray, val mean: DoubleArray)

}

// Compute the alternative mean just before the jth draw, for a fixed alternative that the original population mean is eta.
// sampleNum starts at 1
class FixedAlternativeMean(val N: Int, val eta0:Double): EstimFn {
    override fun eta(j: Int, sampleSum: Double): Double {
        return (N * eta0 - sampleSum) / (N - j + 1)
    }
}
