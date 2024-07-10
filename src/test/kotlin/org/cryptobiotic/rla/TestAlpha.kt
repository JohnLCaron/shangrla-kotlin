package org.cryptobiotic.rla

import kotlin.random.Random
import kotlin.test.Test

class TestAlpha {

    @Test
    fun testWithReplacement() {
        val etas = listOf(.66, .6, .55, .51, .501)  // alternative means
        val t = .5
        val u = 1.0
        val N = 100

        for (eta in etas) {
            val c = (eta - t) / 2
            val estimFn = FixedAlternativeMean(N, setEta0(u))
            val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, upperBound=u)
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

}