package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.Bernoulli
import org.cryptobiotic.start.AlphaMart
import org.cryptobiotic.start.EstimFnType
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt
import kotlin.test.Test

// Compare AlphaAlgorithm with output from start/TestNonnegMean
class TestAlphaShrinkTrunc {
    val eps = 0.0001  // Generic small value

    /*
    @Test
    fun testShrinkTrunc() {
        val bernoulliDist = Bernoulli(.5)
        val bernoulliList = DoubleArray(20) { bernoulliDist.get() }.toList()

        val v = listOf(
            listOf(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            listOf(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            bernoulliList
        )

        testShrinkTrunc(
            .51,
            v[0],
            listOf(0.51, 0.5570631122784444, 0.626443375672974, 0.7156724647762773, 0.8346696395428955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999)
        )
        /*
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
         */
    }


    fun testShrinkTrunc(eta: Double, x: List<Double>, expected: List<Double>) {
        val t = .5
        val u = 1.0
        val d = 10
        val f = 0.0
        val minsd = 1.0e-6

        val N = x.size

        val alpha = AlphaAlgorithm(N = N, d = d, withoutReplacement = false, t = .0001, upperBound = 2.0)
        val sampler = SampleFromList(x)
        val algoValues = alpha.run(x.size) { sampler.sample() }

        val algoValues = alpha.shrink_trunc(x.toDoubleArray(), minsd=minsd, d = d, eta=eta, f=f, c=c, eps=eps)

        doublesAreEqual(expected, algoValues.etaj)
    }

     */

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