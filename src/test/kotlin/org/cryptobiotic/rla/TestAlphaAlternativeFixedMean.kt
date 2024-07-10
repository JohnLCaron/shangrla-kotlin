package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.core.numpy_isclose
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

fun setEta0(upperBound: Double) : Double {
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

// Compare AlphaAlgorithm with output from start/TestNonnegMean
class TestAlphaAlternativeFixedMean  {

    @Test
    fun testAlphaMartAllHalf2() {
        //        # When all the items are 1/2, estimated p for a mean of 1/2 should be 1.
        //        s = np.ones(5)/2
        //        test = NonnegMean(N=int(10**6))
        //        np.testing.assert_almost_equal(test.alpha_mart(s)[0],1.0)
        val s1 = DoubleArray(5) { .5 }
        val u = 1.0
        val N = 100

        val estimFn = FixedAlternativeMean(N, setEta0(u))
        val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, upperBound=u)
        val sampler = SampleFromList(s1)
        val algoValues = alpha.run(s1.size) { sampler.sample() }

        // run 0.75, 0.75000025000025, 0.750000500001, 0.75000075000225, 0.750001000004
        // run2 0.75, 0.7525252525252525, 0.7551020408163265, 0.7577319587628866, 0.7604166666666666
        doublesAreClose(listOf(0.75, 0.7525252525252525, 0.7551020408163265, 0.7577319587628866, 0.7604166666666666), algoValues.etaj)
        algoValues.phistory.forEach { assertEquals(1.0, it)}
    }

    @Test
    fun testAlphaMartEps() {
        val N = 1_000_000
        val t = .0001
        val x = DoubleArray(5) { .5 }

        val estimFn = FixedAlternativeMean(N, setEta0(1.0))
        val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, t=t)
        val sampler = SampleFromList(x)
        val algoValues = alpha.run(x.size) { sampler.sample() }

        doublesAreEqual(listOf(0.50005, 0.50005000005, 0.5000500001000002, 0.5000500001500005, 0.5000500002000008), algoValues.etaj)
        doublesAreEqual(listOf(1.0E-4, 9.95000995000995E-5, 9.900019800039601E-5, 9.85002955008865E-5, 9.8000392001568E-5), algoValues.populationMeanValues)
        doublesAreEqual(listOf(2500.5, 2513.0615576639316, 2525.749999749975, 2538.5672585082107, 2551.515305622398), algoValues.tjs)
        doublesAreEqual(listOf(2500.5, 6283910.424938661, 1.5871586754217688E10, 4.0291090474829625E13, 1.028033340267446E17), algoValues.tstat)
        doublesAreEqual(listOf(3.999200159968006E-4, 1.591365777639584E-7, 6.300567268324711E-11, 2.4819382851519325E-14, 9.727310981371936E-18), algoValues.phistory)
    }

    // etaj = [1.00005, 1.0000504000504002, 1.0000506001012004, 1.0000506001518006, 1.0000504002016009]
    //m = [1.0E-4, 9.94000994000994E-5, 9.86001972003944E-5, 9.760029280087839E-5, 9.640038560154241E-5]
    //tj = [3000.5, 4024.643661761821, 5071.490364786546, 6148.03360619736, 7261.897717512307]
    //T = [3000.5, 1.2075943307116345E7, 6.1243030127749115E10, 3.76524207370759E14, 2.7342802820938455E18]
    @Test
    fun testAlphaMartU2() {
        val N = 1_000_000
        val u = 2.0
        val t = .0001
        val x = doubleArrayOf(0.6, 0.8, 1.0, 1.2, 1.4)

        val estimFn = FixedAlternativeMean(N, setEta0(u))
        val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, upperBound=u, t=t)

        val sampler = SampleFromList(x)
        val algoValues = alpha.run(x.size) { sampler.sample() }

        doublesAreEqual(listOf(1.00005, 1.0000504000504002, 1.0000506001012004, 1.0000506001518006, 1.0000504002016009), algoValues.etaj)
        doublesAreEqual(listOf(1.0E-4, 9.94000994000994E-5, 9.86001972003944E-5, 9.760029280087839E-5, 9.640038560154241E-5), algoValues.populationMeanValues)
        doublesAreEqual(listOf(3000.5, 4024.643661761821, 5071.490364786546, 6148.03360619736, 7261.897717512307), algoValues.tjs)
        doublesAreEqual(listOf(3000.5, 1.2075943307116345E7, 6.1243030127749115E10, 3.76524207370759E14, 2.7342802820938455E18), algoValues.tstat)
        doublesAreEqual(listOf(3.332777870354941E-4, 8.280926587413679E-8, 1.6328388682174326E-11, 2.655871735267506E-15, 3.6572695438311987E-19), algoValues.phistory)
    }

    @Test
    fun testLastObs() {
        //         s1 = [1, 0, 1, 1, 0, 0, 1]
        //        test = NonnegMean(N=7, u=1, t = 3/7)
        val N = 7
        val u = 1.0
        val t = 3.0 / 7.0
        val x1 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0)

        val estimFn = FixedAlternativeMean(N, setEta0(u))
        val alpha = AlphaAlgorithm(estimFn=estimFn, N=N, upperBound=u, t=t)

        val sampler = SampleFromList(x1)
        val algoValues = alpha.run(x1.size) { sampler.sample() }

        // test_lastObs()
        // etaj = [0.7142857142857142, 0.6666666666666665, 0.7999999999999998, 0.7499999999999998, 0.6666666666666664, 0.9999999999999996, 1.9999999999999991]
        //m = [0.42857142857142855, 0.3333333333333333, 0.4, 0.25, 0.0, 0.0, 0.0]
        //tj = [1.6666666666666665, 0.5000000000000001, 1.9999999999999996, 2.999999999999999, NaN, NaN, Infinity]
        //T = [1.6666666666666665, 0.8333333333333335, 1.6666666666666665, 4.999999999999998, NaN, NaN, NaN]
        //test_lastObs alpha_mart1=[0.6000000000000001, 1.0, 0.6000000000000001, 0.20000000000000007, 1.0, 1.0, 0.0]

        // actual
        //  etajs [0.7142857142857142, 0.6666666666666665, 0.7999999999999998, 0.7499999999999998, 0.6666666666666664, 0.9999999999999996, 1.9999999999999991]
        // populationMeanValues [0.5, 0.4, 0.5, 0.3333333333333333, 0.0, 0.0, NaN]
        // tjs [1.4285714285714284, 0.5555555555555558, 1.5999999999999996, 2.2499999999999996, NaN, NaN, NaN]
        // tstat [1.4285714285714284, 0.7936507936507939, 1.26984126984127, 2.8571428571428568, NaN, NaN, NaN]
        // phistory [0.7000000000000001, 1.0, 0.7874999999999999, 0.35000000000000003, NaN, NaN, NaN]

        doublesAreEqual(listOf(0.7142857142857142, 0.6666666666666665, 0.7999999999999998, 0.7499999999999998, 0.6666666666666664, 0.9999999999999996, 1.9999999999999991), algoValues.etaj)
        doublesAreEqual(listOf(0.42857142857142855, 0.3333333333333333, 0.4, 0.25, 0.0, 0.0, 0.0), algoValues.populationMeanValues)
        doublesAreEqual(listOf(1.6666666666666665, 0.5000000000000001, 1.9999999999999996, 2.999999999999999, Double.NaN, Double.NaN, Double.POSITIVE_INFINITY), algoValues.tjs)
        doublesAreEqual(listOf(1.6666666666666665, 0.8333333333333335, 1.6666666666666665, 4.999999999999998, 1.0, 1.0, Double.POSITIVE_INFINITY), algoValues.tstat)
        doublesAreEqual(listOf(0.6000000000000001, 1.0, 0.6000000000000001, 0.20000000000000007, 1.0, 1.0, 0.0), algoValues.phistory)

        val x2 = doubleArrayOf(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        val sampler2 = SampleFromList(x2)
        val algoValues2 = alpha.run(x2.size) { sampler2.sample() }

        doublesAreEqual(listOf(0.7142857142857142, 0.6666666666666665, 0.7999999999999998, 0.7499999999999998, 0.6666666666666664, 0.9999999999999996, 1.9999999999999991), algoValues2.etaj)
        doublesAreEqual(listOf(0.42857142857142855, 0.3333333333333333, 0.4, 0.25, 0.0, 0.0, 0.0), algoValues2.populationMeanValues)
        doublesAreEqual(listOf(1.6666666666666665, 0.5000000000000001, 1.9999999999999996, 2.999999999999999, Double.NaN, Double.NaN, Double.NaN), algoValues2.tjs)
        doublesAreEqual(listOf(1.6666666666666665, 0.8333333333333335, 1.6666666666666665, 4.999999999999998, 1.0, 1.0, 1.0), algoValues2.tstat)
        doublesAreEqual(listOf(0.6000000000000001, 1.0, 0.6000000000000001, 0.20000000000000007, 1.0, 1.0, 1.0), algoValues2.phistory)
    }
}

class SampleFromList(val list: DoubleArray) {
    var index = 0
    fun sample() = list[index++]
}

fun doublesAreClose(a: List<Double>, b: List<Double>, rtol: Double=1.0e-5, atol:Double=1.0e-8) {
    //    For finite values, isclose uses the following equation to test whether
    //    two floating point values are equivalent.
    //
    //     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    assertEquals(a.size, b.size, "size differs")
    repeat(a.size) { assertTrue(numpy_isclose(a[it], b[it], rtol, atol), "$it: ${a[it]} !~ ${b[it]}") }
}

fun doublesAreClose(a: DoubleArray, b: DoubleArray, rtol: Double=1.0e-5, atol:Double=1.0e-8) {
    //    For finite values, isclose uses the following equation to test whether
    //    two floating point values are equivalent.
    //
    //     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    assertEquals(a.size, b.size, "size differs")
    repeat(a.size) { assertTrue(numpy_isclose(a[it], b[it], rtol, atol), "$it: ${a[it]} !~ ${b[it]}") }
}

fun doublesAreEqual(a: List<Double>, b: List<Double>) {
    //    For finite values, isclose uses the following equation to test whether
    //    two floating point values are equivalent.
    //
    //     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    assertEquals(a.size, b.size, "size differs")
    repeat(a.size) { assertEquals(a[it], b[it], "$it: ${a[it]} != ${b[it]}") }
}