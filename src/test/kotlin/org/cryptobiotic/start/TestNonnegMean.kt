package org.cryptobiotic.start

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class TestNonnegMean {

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
        val c = (eta - t) / 2

        // class AlphaMart(val N: Int, val withReplacement: Boolean, val t: Double, val u: Double, val estimFnType: EstimFnType) : TestFn {
        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)
        val (p, p_history) = alphamart.test(s1)
        assertEquals(1.0, p)
        println("testAlphaMartAllHalf p_history=${p_history.contentToString()}")

        println("\nAlphaAlgorithm")
        val alpha = AlphaAlgorithm(N = 100, d = d, withoutReplacement = false)
        val sampler = SampleFromList(s1)
        alpha.run(s1.size) { sampler.sample() }
    }

    @Test
    fun testAlphaMartEps() {
        // test.t = eps
        // np.testing.assert_array_less(test.alpha_mart(s)[1][1:],[eps]*(len(s)-1))
        val eta = .5
        val s1 = DoubleArray(5) { eta }
        val t = .5
        val u = 1.0
        val d = 10.0
        val c = (eta - t) / 2

        val eps = 0.0001  // Generic small value
        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = eps, u = u, estimFnType = EstimFnType.FIXED)
        val (_, p_history) = alphamart.test(s1)
        p_history.forEach { it < eps * (s1.size - 1) }
        println("testAlphaMartEps p_history=${p_history.contentToString()}")
    }

    @Test
    fun testAlphaMartU2() {
        val eps = 0.0001  // Generic small value
        val eta = .5
        val u = 2.0
        val d = 10.0
        val t = eps
        val c = (eta - t) / 2

        //        s = [0.6,0.8,1.0,1.2,1.4]
        //        test.u=2 # TODO eta not recalculated
        //        np.testing.assert_array_less(test.alpha_mart(s)[1][1:],[eps]*(len(s)-1))
        val s2 = doubleArrayOf(0.6, 0.8, 1.0, 1.2, 1.4)
        val alphamart = AlphaMart(N = 1_000_000, withReplacement = false, t = t, u = u, estimFnType = EstimFnType.FIXED)
        val (_, p_history) = alphamart.test(s2)
        p_history.forEach { it < eps * (s2.size - 1) }
        println("testAlphaMartU2 p_history=${p_history.contentToString()}")

    }

        /*
            s1 = [1, 0, 1, 1, 0, 0, 1]
    //        test.u=1
    //        test.N = 7
    //        test.t = 3/7
    //        alpha_mart1 = test.alpha_mart(s1)[1]
    //        # p-values should be big until the last, which should be 0
    //        print(f'{alpha_mart1=}')
    //        assert(not any(np.isnan(alpha_mart1)))
    //        assert(alpha_mart1[-1] == 0)
    val s3 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0)
    val test3 = NonnegMean(N = 7, u = 1.0, t = 3.0 / 7)
    val (_, p3_history)  = test3.alpha_mart(s3)
    // p-values should be big until the last, which should be 0
    // assert(not any (np.isnan(alpha_mart1)))
    p3_history.forEach{ !it.isNaN() }
    assert(p3_history.last() == 0.0)

    //         s2 = [1, 0, 1, 1, 0, 0, 0]
    //        alpha_mart2 = test.alpha_mart(s2)[1]
    //        # Since s1 and s2 only differ in the last observation,
    //        # the resulting martingales should be identical up to the next-to-last.
    //        # Final entry in alpha_mart2 should be 1
    //        assert(all(np.equal(alpha_mart2[0:(len(alpha_mart2)-1)],
    //                            alpha_mart1[0:(len(alpha_mart1)-1)])))
    //        print(f'{alpha_mart2=}')
    val s4 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)
    val (_, p4_history) = test3.alpha_mart(s4)
    // Since s3 and s4 only differ in the last observation,
    // the resulting martingales should be identical up to the next-to-last.
    // Final entry in alpha_mart2 should be 1
    assertEquals(p3_history.size, p4_history.size)
    repeat(p3_history.size-1) { assertEquals(p3_history[it], p4_history[it])}
    assertEquals(1.0, p4_history.last())
    println("p3_history=${p3_history.contentToString()}")
    println("p4_history=${p4_history.contentToString()}")

     */

    class SampleFromList(val list: DoubleArray) {
        var index = 0
        fun sample() = list[index++]
    }
}