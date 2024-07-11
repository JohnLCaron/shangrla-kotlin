package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.Bernoulli
import kotlin.math.pow
import kotlin.math.sqrt
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

    // oops python caught this at PR #89
    @Test
    fun testWelford() {
        val sample = DoubleArray(10) { if (it % 2 == 0) 1.0 else 2.0 }
        w2(sample)

        println("mean = ${sample.average()}")
        println("variance = ${sample.variance()}")

        println("Welford")
        val welford = Welford()
        sample.forEachIndexed { idx, it ->
            welford.update(it)
            println(" $idx = ${welford.result()}")
        }
    }

    @Test
    fun testWelfordRandom() {
        val b = Bernoulli( .45)
        val sample = DoubleArray(10) { b.get() }
        w2(sample)
    }

}

fun w2(x: DoubleArray) {
    // Welford's algorithm for running mean and running sd
    val mj = mutableListOf<Double>()
    mj.add(x[0])
    var sdj = mutableListOf<Double>()
    sdj.add(0.0)

    //        for i, xj in enumerate(x[1:]):
    //            mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
    //            sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
    // enumerate returns Pair(index, element)
    for (idx in 0 until x.size-1) {
        val xj = x[idx+1]
        // mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
        mj.add(mj.last() + (xj - mj.last()) / (idx + 2))
        // sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
        sdj.add(sdj.last() + (xj - mj[mj.size - 2]) * (xj - mj.last()))
    }
    // sdj = np.sqrt(sdj / j)
    val sdj2 = sdj.mapIndexed { idx, it -> sqrt(it / (idx + 1)) }
    // end of Welford's algorithm.

    println("running mean = ${mj}")
    println("running variance = ${sdj2}")
}

/* Welford's algorithm for running mean and running sd
fun welford(x: DoubleArray): Pair<List<Double>, List<Double>> {
    val mj = mutableListOf<Double>()
    mj.add(x[0])
    var sdj = mutableListOf<Double>()
    sdj.add(0.0)

    for (idx in 0 until x.size - 1) {
        val xj = x[idx + 1]
        // mj.append(mj[-1] + (xj - mj[-1]) / (i + 1))
        mj.add(mj.last() + (xj - mj.last()) / (idx + 1))
        // sdj.append(sdj[-1] + (xj - mj[-2]) * (xj - mj[-1]))
        sdj.add(sdj.last() + (xj - mj[mj.size - 2]) * (xj - mj.last()))
    }
    val sdj2 = sdj.mapIndexed { idx, it -> it / (idx + 1) }

    println("running mean = ${mj}")
    println("running std = ${sdj2}")
    return Pair(mj, sdj2)
}

 */

fun DoubleArray.variance(): Double {
    val mean = this.average()
    val variance = this.map { (it - mean).pow(2) }.average()
    return Math.sqrt(variance)
}



