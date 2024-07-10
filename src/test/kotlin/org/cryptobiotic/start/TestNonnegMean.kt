package org.cryptobiotic.start

import org.cryptobiotic.rla.doublesAreClose
import org.cryptobiotic.shangrla.Bernoulli
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull

class TestNonnegMean {
    val eps = 0.0001  // Generic small value

    //     var u: Double = 1.0,        // TODO mutable
    //    val N: Int = Int.MAX_VALUE, // TODO withReplacement: If N is np.inf, it means the sampling is with replacement
    //    val t: Double = 0.5, // TODO is it ever anything other than .5 ??
    //    val random_order: Boolean = true,
    //    val gOverride: Double? = null, // used only in the 3 KAPLAN_* tests.
    //    val withReplacement: Boolean = false,
    //    testFnType : TestFnType? = TestFnType.ALPHA_MART,

    @Test
    fun testAlphaMartAllHalf() {
        //        # When all the items are 1/2, estimated p for a mean of 1/2 should be 1.
        //        s = np.ones(5)/2
        //        test = NonnegMean(N=int(10**6))
        //        np.testing.assert_almost_equal(test.alpha_mart(s)[0],1.0)
        val eta = .5
        val s1 = DoubleArray(5) { eta }
        val t = .5
        val u = 1.0
        val d = 10.0

        // class AlphaMart(val N: Int, val withReplacement: Boolean, val t: Double, val u: Double, val estimFnType: EstimFnType) : TestFn {
        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)
        val (p, p_history) = alphamart.test(s1)
        assertEquals(1.0, p)
        println("testAlphaMartAllHalf p_history=${p_history.contentToString()}")
    }

    @Test
    fun testAlphaMartEps() {
        // test.t = eps
        // np.testing.assert_array_less(test.alpha_mart(s)[1][1:],[eps]*(len(s)-1))
        val eta = .5
        val s1 = DoubleArray(5) { eta }
        val u = 1.0

        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = eps, u = u, estimFnType = EstimFnType.FIXED)
        val (_, p_history) = alphamart.test(s1)
        p_history.forEach { it < eps * (s1.size - 1) }
        println("testAlphaMartEps p_history=${p_history.contentToString()}")

        // 3.99920016e-04 1.59136578e-07 6.30056727e-11 2.48193829e-14 9.72731098e-18
        doublesAreClose(listOf(3.99920016e-04, 1.59136578e-07, 6.30056727e-11, 2.48193829e-14, 9.72731098e-18), p_history.toList())
    }

    @Test
    fun testAlphaMartU2() {
        val u = 2.0
        val t = eps

        //        s = [0.6,0.8,1.0,1.2,1.4]
        //        test.u=2 # TODO eta not recalculated
        //        np.testing.assert_array_less(test.alpha_mart(s)[1][1:],[eps]*(len(s)-1))
        val s2 = doubleArrayOf(0.6, 0.8, 1.0, 1.2, 1.4)
        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)
        val (_, p_history) = alphamart.test(s2)
        p_history.forEach { it < eps * (s2.size - 1) }
        println("testAlphaMartU2 p_history=${p_history.contentToString()}")

        // 3.33277787e-04 8.28092659e-08 1.63283887e-11 2.65587174e-15 3.65726954e-19
        doublesAreClose(listOf(3.33277787e-04, 8.28092659e-08, 1.63283887e-11, 2.65587174e-15, 3.65726954e-19), p_history.toList())
    }

    @Test
    fun test_lastObs() {
        val u = 1.0
        val t = 3.0 / 7.0

        val alphamart = AlphaMart(N = 7, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)

        val s1 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0)
        val (_, alpha_mart1) = alphamart.test(s1)
        // p - values should be big until the last, which should be 0
        println("test_lastObs alpha_mart1=${alpha_mart1.contentToString()}")
        assertNull(alpha_mart1.find { it == Double.POSITIVE_INFINITY})
        assertEquals(0.0, alpha_mart1.last())

         // alpha_mart1=array([0.6, 1. , 0.6, 0.2, 1. , 1. , 0. ])
        // test_lastObs alpha_mart1=[0.6000000000000001, 1.0, 0.6000000000000001, 0.20000000000000007, 1.0, 1.0, 0.0]
        doublesAreClose(listOf(0.6, 1.0 , 0.6, 0.2, 1.0 , 1.0 , 0.0), alpha_mart1.toList())

        val s2 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        val (_, alpha_mart2) = alphamart.test(s2)

        println("test_lastObs alpha_mart1=${alpha_mart2.contentToString()}")
        assertNull(alpha_mart2.find { it == Double.POSITIVE_INFINITY})
        assertEquals(1.0, alpha_mart2.last())

        // alpha_mart2=array([0.6, 1. , 0.6, 0.2, 1. , 1. , 1. ])
        doublesAreClose(listOf(0.6, 1.0 , 0.6, 0.2, 1.0 , 1.0 , 1.0), alpha_mart2.toList())
    }

    //     def test_shrink_trunc(self):
    //        epsj = lambda c, d, j: c/math.sqrt(d+j-1)
    //        Sj = lambda x, j: 0 if j==1 else np.sum(x[0:j-1])
    //        tj = lambda N, t, x, j: (N*t - Sj(x, j))/(N-j+1) if np.isfinite(N) else t
    //        etas = [.51, .55, .6]  # alternative means
    //        t = 1/2
    //        u = 1
    //        d = 10
    //        f = 0
    //        vrand =  sp.stats.bernoulli.rvs(1/2, size=20)
    //        v = [
    //            np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
    //            np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0]),
    //            vrand
    //        ]
    //        test_inf = NonnegMean(N=np.inf, t=t, u=u, d=d, f=f)
    //        test_fin = NonnegMean(          t=t, u=u, d=d, f=f)
    //        for eta in etas:
    //            c = (eta-t)/2
    //            test_inf.c = c
    //            test_inf.eta = eta
    //            test_fin.c=c
    //            test_fin.eta=eta
    //            for x in v:
    //                N = len(x)
    //                test_fin.N = N
    //                xinf = test_inf.shrink_trunc(x)
    //                xfin = test_fin.shrink_trunc(x)
    //                yinf = np.zeros(N)
    //                yfin = np.zeros(N)
    //                for j in range(1,N+1):
    //                    est = (d*eta + Sj(x,j))/(d+j-1)
    //                    most = u*(1-np.finfo(float).eps)
    //                    yinf[j-1] = np.minimum(np.maximum(t+epsj(c,d,j), est), most)
    //                    yfin[j-1] = np.minimum(np.maximum(tj(N,t,x,j)+epsj(c,d,j), est), most)
    //                np.testing.assert_allclose(xinf, yinf)
    //                np.testing.assert_allclose(xfin, yfin)
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

    fun epsj(c: Double, d: Int, j:Int): Double =  c/sqrt(d+j-1.0)
    fun Sj(x: List<Double>, j:Int): Double = if (j == 1) 0.0 else x.subList(0,j-1).sum()
    fun tj(N:Int, t: Double, x: List<Double>, j:Int) =  (N*t-Sj(x, j))/(N-j+1)

}

