package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.Bernoulli
import org.cryptobiotic.start.AlphaMart
import org.cryptobiotic.start.EstimFnType
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt
import kotlin.test.Test

fun makeEstimator() {
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
}

// Compare AlphaAlgorithm with output from start/TestNonnegMean
class TestAlphaShrinkTrunc {
    val eps = 0.0001  // Generic small value

    @Test
    fun testShrinkTrunc() {
        val bernoulliDist = Bernoulli(.5)
        val bernoulliList = DoubleArray(20) { bernoulliDist.get() }.toList()

        val v = listOf(
            listOf(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            listOf(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            bernoulliList
        )

        // x = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        //j = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        //running mean = [0.0, 0.0, 0.0, 0.0, 0.0, 0.16666666666666666, 0.2857142857142857, 0.375, 0.4444444444444444, 0.5]
        //running std = [0.0, 0.0, 0.0, 0.0, 0.0, 0.37267799624996495, 0.45175395145262565, 0.4841229182759271, 0.4969039949999533, 0.5]
        //xfin = [0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999]
        //yfin = [0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999]

        testShrinkTrunc(
            .51,
            v[0],
            listOf(0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999)
        )
        val etas = listOf(.51, .55, .6) // alternative means
        for (eta in etas) {
            for (x: List<Double> in v) {
                testShrinkTrunc(
                    eta,
                    x,
                    listOf(0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999)
                )
            }
        }
    }

    @Test
    fun testTruncShrinkage() {
        val eta0 = .51
        val x = listOf(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        val expected = listOf(0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999)

        val t = .5
        val u = 1.0
        val d = 10
        val f = 0.0
        val minsd = 1.0e-6
        val c = (eta0 - t) / 2

        val N = x.size

        val estimFn = TruncShrinkage(N = N, u=u, t=t, minsd=minsd, d=d, eta0=eta0, f=f, c=c, eps=eps)
        val result = mutableListOf<Double>()
        repeat(x.size) { it ->
            val sublist = x.subList(0, it+1)
            val estim = estimFn.eta(sublist)
            println("estimOld = $estim")
            result.add(estim)
        }
        doublesAreEqual(expected, result)
    }

    @Test
    fun testTruncShrinkageOld() {
        val eta0 = .51
        val x = listOf(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        val expected = listOf(0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999)

        val t = .5
        val u = 1.0
        val d = 10
        val f = 0.0
        val minsd = 1.0e-6
        val c = (eta0 - t) / 2

        val N = x.size
        val estimFn = TruncShrinkage(N = N, u=u, t=t, minsd=minsd, d=d, eta0=eta0, f=f, c=c, eps=eps)

        val result = mutableListOf<Double>()
        repeat(x.size) { it ->
            val sublist = x.subList(0, it+1)
            val estim = estimFn.etaOld(sublist)
            println("estimOld = $estim")
            result.add(estim)
        }
        doublesAreEqual(expected, result)
    }


    fun testShrinkTrunc(eta0: Double, x: List<Double>, expected: List<Double>) {
        val t = .5
        val u = 1.0
        val d = 10
        val f = 0.0
        val minsd = 1.0e-6
        val c = (eta0 - t) / 2

        val N = x.size


        val estimFn = TruncShrinkage(N = N, u=u, t=t, minsd=minsd, d=d, eta0=eta0, f=f, c=c, eps=eps)
        val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, upperBound=u, t=t)

        val sampler = SampleFromList(x.toDoubleArray())
        val algoValues = alpha.run(x.size) { sampler.sample() }

        doublesAreEqual(expected, algoValues.etaj)
    }

    @Test
    fun test_shrink_trunc() {
        val t = .5
        val u = 1.0
        val d = 10
        val f = 0.0
        val minsd = 1.0e-6

        val bernoulliDist = Bernoulli(.5)
        val bernoulliList = DoubleArray(20) { bernoulliDist.get() }.toList()

        val v = listOf(
            listOf(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            listOf(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            bernoulliList
        )

        val etas = listOf(.51, .55, .6) // alternative means
        for (eta in etas) {
            val c = (eta - t) / 2

            for (x: List<Double> in v) {
                val N = x.size
                // class AlphaMart(val N: Int, val withReplacement: Boolean, val t: Double, val u: Double, val estimFnType: EstimFnType) : TestFn {
                val alphamart = AlphaMart(N = N, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)

                // fun shrink_trunc(x: DoubleArray, minsd : Double, d: Int, eta: Double, f: Double, c: Double, eps: Double): DoubleArray {
                val xfin = alphamart.shrink_trunc(x.toDoubleArray(), minsd=minsd, d = d, eta=eta, f=f, c=c, eps=eps)
                val yfin = DoubleArray(N)
//                 for j in range(1,N+1):
//                    est = (d*eta + Sj(x,j))/(d+j-1)
//                    most = u*(1-np.finfo(float).eps)
//                    yinf[j-1] = np.minimum(np.maximum(t+epsj(c,d,j), est), most)
//                    yfin[j-1] = np.minimum(np.maximum(tj(N,t,x,j)+epsj(c,d,j), est), most)

                repeat(N) {
                    val j = it + 1
                    val sj = Sj(x, j)
                    val est = (d * eta + Sj(x, j)) / (d + j - 1)
                    val most = u * (1 - eps)
                    val tjv = tj(N, t, x, j)
                    val epsjv = epsj(c, d, j)
                    yfin[it] = min( max(tj(N, t, x, j) + epsj(c, d, j), est), most)
                }
                println("xfin = ${xfin.contentToString()}")
                println("yfin = ${yfin.contentToString()}")
                doublesAreClose(xfin, yfin)
                println()
            }
        }
    }

    fun epsj(c: Double, d: Int, j:Int): Double =  c/ sqrt(d+j-1.0)
    fun Sj(x: List<Double>, j:Int): Double = if (j == 1) 0.0 else x.subList(0,j-1).sum()
    fun tj(N:Int, t: Double, x: List<Double>, j:Int) =  (N*t-Sj(x, j))/(N-j+1)
}