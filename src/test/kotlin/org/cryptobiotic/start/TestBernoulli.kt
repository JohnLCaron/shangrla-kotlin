package org.cryptobiotic.start

import org.cryptobiotic.shangrla.Bernoulli
import org.junit.jupiter.api.Test

class TestBernoulli {

    @Test
    fun testDist() {
        val bernoulliDist = Bernoulli(.5)
        var sum = 0.0
        repeat(100) {
            var sumInner = 0.0
            repeat(100) {
                val s = bernoulliDist.get()
                // print(" $s,")
                sumInner += s
            }
            print("  ${sumInner/100}")
            sum += (sumInner/100)
        }
        println("\nsum = ${sum/100}")
    }
}