package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.NonnegMean
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals

class TestAlpha {

    @Test
    fun testWithReplacement() {
        val etas = listOf(.66, .6, .55, .51, .501)  // alternative means
        val t = .5
        val u = 1
        val d = 10.0

        for (eta in etas) {
            val c = (eta - t) / 2
            val alpha = AlphaAlgorithm(N = 100, d = d, withoutReplacement = false)
            testAlpha(alpha, eta)
        }
    }

    fun testAlpha(alpha: AlphaAlgorithm, theta: Double) {
        println("Test theta = $theta")
        alpha.run( 100) { sample(theta) }
    }

    fun sample(theta: Double): Double {
        val random = Random.nextDouble(1.0)
        val vote =  if (random < theta) 1 else 0
        return assort(vote)
    }

    fun assort(vote: Int): Double {
        val w = if (vote == 1) 1.0 else 0.0
        val l = if (vote == 0) 1.0 else 0.0
        val a =  (w - l + 1) * 0.5 // eq 1.
        return a
    }

    @Test
    fun test_alpha_mart() {
        //        # When all the items are 1/2, estimated p for a mean of 1/2 should be 1.
        //        s = np.ones(5)/2
        //        test = NonnegMean(N=int(10**6))
        //        np.testing.assert_almost_equal(test.alpha_mart(s)[0],1.0)
        val s1 = DoubleArray(5) { .5 }
        val test = NonnegMean(N = 1_000_000)
        val (p, _) = test.alpha_mart(s1)
        assertEquals(1.0, p)
    }

}