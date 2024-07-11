package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.numpy_cumsum
import org.cryptobiotic.shangrla.core.numpy_isclose
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

//         return AlgoValues(etajsa, populationMeanValues, tjs, tjs
data class AlgoValues(val etaj: List<Double>, val populationMeanValues: List<Double>, val tjs: List<Double>,
                      val tstat: List<Double>, val phistory: List<Double>)

// estimate the population mean for the jth sample from the previous j-1 samples
interface EstimFn {
    // return eta estimate
    // dont pass the empty list, start with j = 2 with prevSamples.size = 1
    fun eta(prevSamples: List<Double>): Double
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
            // val sampleMean = sampleAssortValues.average()
            sampleNumber++

            //// use previous sample sum
            populationMean = this.populationMean(sampleNumber, sampleSum)
            sampleMeanValues.add(populationMean)

            //// use previous sample list
            val etaj = estimFn.eta(sampleAssortValues)
            etajsa.add(etaj)

            sampleAssortValues.add(xj)
            sampleSum += xj // – S ← S + Xj

            // T is the "ALPHA supermartingale"
            // – If m < 0, T ← ∞. Otherwise, T ← T / u * ( Xj * η(j,S)/m + (u - Xj) * (u−η(j,S))/(u-m))
            // TODO This is eq 4 of ALPHA, p.5 :
            //      T_j = T_j-1 * (X_j * eta_j / µ_j + (u - X_j) * (u - eta_j) / ( u - µ_j)) / u
            // python:  terms = np.cumprod((x * etaj / m + (u - x) * (u - etaj) / (u - m)) / u)
            // python: m = µ_j = The mean of the population after each draw if the null hypothesis is true.
            //          m = ( (N * t - S) / (N - j + 1) if np.isfinite(N) else t )  # mean of population after (j-1)st draw, if null is true
            // X_j = jthe sample
            // eta_j = estimate of population mean = eta(X^j-a)
            // u = upperBound
            // µ_j = m = The mean of the population after each draw if the null hypothesis is true.

            val tj = if (populationMean < 0.0) Double.POSITIVE_INFINITY else {
                (xj * etaj / populationMean + (upperBound - xj) * (upperBound - etaj) / (upperBound - populationMean))/ upperBound
            }
            tjs.add(tj)
            testStatistic *= tj
            println("   $sampleNumber = $xj etaj = $etaj tj=$tj, T = $testStatistic")
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
        println("   etajs ${etajsa}")
        println("   populationMeanValues ${sampleMeanValues}")
        println("   tjs ${tjs}")

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

        println("   tstat ${tstat}")

        // return min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms)
        val phistory = tstat.map { min(1.0, 1.0 / it) }
        // val p = min(1.0, 1.0 / tstat.max()) // seems wrong
        println("   phistory ${phistory}")

        //
        //• If a full hand count is conducted, its results replace the reported results if they differ.
        return AlgoValues(etajsa, sampleMeanValues, tjs, tstat, phistory)
    }

    // population mean under the null hypothesis that “the average of this list is not greater than 1/2”
    // TODO seems to be the "mean of the remaining sample", which is why it depends on whether you are replacing or not
    // TODO detect if it goes negetive
    fun populationMean(sampleNum: Int, sampleSumM1: Double): Double {
        return if (withoutReplacement) (N * t - sampleSumM1) / (N - sampleNum + 1) else t
    }
}

data class CumulativeSum(val S: DoubleArray, val Stot: Double, val indices: IntArray, val mean: DoubleArray)


// Compute the alternative mean just before the jth draw, for a fixed alternative that the original population mean is eta.
class FixedAlternativeMean(val N: Int, val eta0:Double): EstimFn {

    override fun eta(prevSamples: List<Double>): Double {
        val j = prevSamples.size + 1
        val sampleSum = prevSamples.sum()
        return (N * eta0 - sampleSum) / (N - j + 1)
    }

}

class TruncShrinkage(val N: Int, val u: Double, val t: Double, val minsd : Double, val d: Int, val eta0: Double, val f: Double, val c: Double, val eps: Double): EstimFn {
    val welford = Welford()

    override fun eta(prevSamples: List<Double>): Double {
        if (prevSamples.size == 0) return eta0
        welford.update(prevSamples.last())
        // if (prevSamples.size == 1) return eta0

        val lastj = prevSamples.size
        val sampleSum = prevSamples.subList(0, lastj - 1).sum()

        val (_, _, std) = welford.result()
        val sdj3 = if (lastj < 2) 1.0 else max(std, minsd)

        val mean = mean(N, t, prevSamples)
        val mean2 = mean2(N, t, prevSamples)
        require(mean == mean2)

        val weighted = ((d * eta0 + sampleSum) / (d + lastj - 1) + u * f / sdj3) / (1 + f / sdj3) // (2.5.2, eq 14, "truncated shrinkage")
        // Choosing ǫi . To allow the estimated winner’s share ηi to approach √ µi as the sample grows
        // (if the sample mean approaches µi or less), we shall take ǫi := c/ d + i − 1 for a nonnegative constant c, for instance c = (η0 − µ)/2.
        // The estimate ηi is thus the sample mean, shrunk towards η0 and truncated to the interval [µi + ǫi , 1), where ǫi → 0 as the sample size grows.
        val npmax = max( weighted, mean2 + c / sqrt((d + lastj - 1).toDouble()))  // 2.5.2 "choosing ǫi"
        val gold = min(u * (1 - eps), npmax)
        return gold
    }

    fun mean(N: Int, t: Double, x: List<Double>): Double {
        val cum_sum = numpy_cumsum(x.toDoubleArray())
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        val m = DoubleArray(x.size) {
            val m1 = (N * t - Sp[it])
            val m2 = (N - j[it] + 1)
            val m3 = m1 / m2
            m3
        }
        println("   j = ${j.last()} m = ${m.last()}")
        return m.last()
    }

    fun mean2(N: Int, t: Double, x: List<Double>): Double {
        val sum = x.subList(0, x.size-1).sum()
        val m1 = (N * t - sum)
        val m2 = (N - x.size + 1)
        val m3 = m1 / m2
        return m3
    }

    fun etaOld(prevSamples: List<Double>): Double {
        // if (prevSamples.size == 0) return eta0

        val mj = mutableListOf<Double>()
        var sdj = mutableListOf<Double>()

        mj.add(prevSamples[0])
        sdj.add(0.0)

        // Welford's algorithm for running mean and running sd
        val (S, _, j, m) = this.sjm(N, t, prevSamples)
        // Welford's algorithm for running mean and running sd
        //        for i, xj in enumerate(x[1:]):
        //            mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
        //            sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
        // enumerate returns Pair(index, element)
        for (idx in 0 until prevSamples.size-1) {
            val xj = prevSamples[idx+1]
            // mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
            mj.add(mj.last() + (xj - mj.last()) / (idx + 2)) // fixed by PR # 89
            // sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
            sdj.add(sdj.last() + (xj - mj[mj.size - 2]) * (xj - mj.last()))
        }
        // sdj = np.sqrt(sdj / j)
        val sdj2 = sdj.mapIndexed { idx, it -> sqrt(it / j[idx]) }
        // end of Welford's algorithm.

        // threshold the sd, set first two sds to 1
        // sdj = np.insert(np.maximum(sdj, minsd), 0, 1)[0:-1]
        val sdj3 = DoubleArray(sdj2.size) { if (it < 2) 1.0 else max(sdj2[it-1], minsd) }

        // weighted = ((d * eta + S) / (d + j - 1) + u * f / sdj) / (1 + f / sdj)
        val weighted = sdj3.mapIndexed { idx, it -> ((d * eta0 + S[idx]) / (d + j[idx] - 1) + u * f / it) / (1 + f / it) }

        // return min(1, 1 / np.max(terms)), np.minimum(1, 1 / terms)
        // np.minimum = element-wise minumum, presumably the smaller of 1 and 1/term
        // np.max = maximum of an array
        // min = min of an iterable
        //val npmin = terms.map { min(1.0, 1.0 / it) }.toDoubleArray()
        //val p = min(1.0, 1.0 / terms.max()) // seems wrong

        // return np.minimum( u * (1 - np.finfo(float).eps), np.maximum(weighted, m + c / np.sqrt(d + j - 1)), )
        val npmax = weighted.mapIndexed { idx, it ->  max( it, m[idx] + c / sqrt((d + j[idx] - 1).toDouble())) }
        return npmax.map { min(u * (1 - eps), it) }.last()
    }

    fun sjm(N: Int, t: Double, x: List<Double>): CumulativeSum {
        val cum_sum = numpy_cumsum(x.toDoubleArray())
        val S = DoubleArray(x.size+1) { if (it == 0) 0.0 else cum_sum[it-1] }   // 0, x_1, x_1+x_2, ...,
        val Stot = S.last()  // sample total ""array[-1] means the last element"
        val Sp = DoubleArray(x.size) { S[it] } // same length as the data.

        val j = IntArray(x.size) { it+1 } // 1, 2, 3, ..., len(x)
        val m = DoubleArray(x.size) {
            val m1 = (N * t - Sp[it])
            val m2 = (N - j[it] + 1)
            val m3 = m1 / m2
            m3
        }
        println("   j = ${j.last()} m = ${m.last()}")
        return CumulativeSum(Sp, Stot, j, m)
    }
}

class Welford() {
    var count = 0
    var mean = 0.0 // mean accumulates the mean of the entire dataset
    var M2 = 0.0 // M2 aggregates the squared distance from the mean

    // For a new value new_value, compute the new count, new mean, the new M2.
    fun update(new_value: Double) {
        count++
        val delta = new_value - mean
        mean += delta / count
        val delta2 = new_value - mean
        M2 += delta * delta2
    }

    // Retrieve the mean, variance and sample variance from an aggregate
    fun result() : Triple<Double, Double, Double> {
        if (count < 2) return Triple(mean, Double.NaN, Double.NaN)
        val variance = M2 / count
        val sample_variance = M2 / (count - 1)
        val std_dev = Math.sqrt(variance)
        return Triple(mean, variance, std_dev)
    }
}
