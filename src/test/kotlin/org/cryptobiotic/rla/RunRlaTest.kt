package org.cryptobiotic.rla

import io.github.oshai.kotlinlogging.KotlinLogging
import kotlinx.cli.ArgParser
import kotlinx.cli.ArgType
import kotlinx.cli.required
import kotlin.test.Test

class RunRlaTest {

    companion object {
        val logger = KotlinLogging.logger("RunRlaTest")

        @JvmStatic
        fun main(args: Array<String>) {
            val parser = ArgParser("RunRlaTest")
            val nballots by parser.option(
                ArgType.Int,
                shortName = "nb",
                description = "number of ballots"
            ).required()
            val ntrials by parser.option(
                ArgType.Int,
                shortName = "nt",
                description = "number of trials"
            ).required()

            parser.parse(args)

            val info = "starting RunRlaTest ntrials= $ntrials nballots= $nballots"
            logger.info { info }

            runTest(ntrials, nballots)
        }

        fun runTest(ntrials: Int, nballots: Int, show: Boolean = false, plot: Boolean = false): List<Double> {
            val hist = mutableListOf<Double>()
            repeat(ntrials) {
                val cvrs = makeSimpleCvrs(nballots)
                val votes = cvrs.map { if (it.votedForA()) 1.0 else 0.0 }.sum()
                val pct = votes/nballots
                hist.add(pct)
                if (show) println("trial $it count ${votes/nballots}")
            }
            return hist
        }
    }

    @Test
    fun testRunTest() {
        RunRlaTest.runTest(10, 1000)
    }
}