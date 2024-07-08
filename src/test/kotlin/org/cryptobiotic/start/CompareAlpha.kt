package org.cryptobiotic.start

import org.cryptobiotic.shangrla.core.numpy_isclose
import org.cryptobiotic.start.TestNonnegMean.SampleFromList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CompareAlpha {
    val epsilon = .0000001

    @Test
    fun testAlphaMartEps() {
        val x = DoubleArray(5) { .5 }
        val d = 10.0

        println("\nAlphaAlgorithm")
        val alpha = AlphaAlgorithm(N = 1_000_000, d = d, withoutReplacement = false, t = .0001)
        val sampler = SampleFromList(x)
        val algoValues = alpha.run(x.size) { sampler.sample() }

        // data class AlgoValues(val etaj: List<Double>, val populationMeanValues: List<Double>, val tjs: List<Double>, val tstat: List<Double>)
        doublesAreClose(listOf(0.50005, 0.50005, 0.50005, 0.50005, 0.50005), algoValues.etaj)
        doublesAreClose(listOf(1.0E-4, 9.95000995000995E-5, 9.900019800039601E-5, 9.85002955008865E-5, 9.8000392001568E-5), algoValues.populationMeanValues)
        doublesAreClose(listOf(2500.5, 2513.0615574127005, 2525.7499992449743, 2538.5672577468645, 2551.5153046020896), algoValues.tjs)
        doublesAreClose(listOf(2500.5, 6283910.424310458, 1.5871586749457624E10, 4.029109045066211E13, 1.0280333392397147E17), algoValues.tstat)
    }

    @Test
    fun testAlphaMartEps2() {
        val x = DoubleArray(5) { .5 }
        val d = 10.0

        println("\nAlphaAlgorithm2")
        val alpha = AlphaAlgorithm(N = 1_000_000, d = d, withoutReplacement = false, t = .0001)
        val sampler = SampleFromList(x)
        val algoValues = alpha.run2(x.size) { sampler.sample() }

        // data class AlgoValues(val etaj: List<Double>, val populationMeanValues: List<Double>, val tjs: List<Double>, val tstat: List<Double>)
        doublesAreClose(listOf(0.50005, 0.50005, 0.50005, 0.50005, 0.50005), algoValues.etaj)
        doublesAreClose(listOf(1.0E-4, 9.95000995000995E-5, 9.900019800039601E-5, 9.85002955008865E-5, 9.8000392001568E-5), algoValues.populationMeanValues)
        doublesAreClose(listOf(2500.5, 2513.0615574127005, 2525.7499992449743, 2538.5672577468645, 2551.5153046020896), algoValues.tjs)
        doublesAreClose(listOf(2500.5, 6283910.424310458, 1.5871586749457624E10, 4.029109045066211E13, 1.0280333392397147E17), algoValues.tstat)
    }
}

//  etajs [0.50005, 0.50005, 0.50005, 0.50005, 0.50005]
// populationMeanValues [1.0E-4, 9.95000995000995E-5, 9.900019800039601E-5, 9.85002955008865E-5, 9.8000392001568E-5]
// tjs [2500.5, 2513.0615574127005, 2525.7499992449743, 2538.5672577468645, 2551.5153046020896]
// tstat [2500.5, 6283910.424310458, 1.5871586749457624E10, 4.029109045066211E13, 1.0280333392397147E17]

fun doublesAreClose(a: List<Double>, b: List<Double>, rtol: Double=1.0e-5, atol:Double=1.0e-8) {
    //    For finite values, isclose uses the following equation to test whether
    //    two floating point values are equivalent.
    //
    //     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    assertEquals(a.size, b.size, "size differs")
    repeat(a.size) { assertTrue(numpy_isclose(a[it], b[it])) }
}