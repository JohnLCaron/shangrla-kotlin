package org.cyrptobiotic.shangrla.reader

import org.cryptobiotic.shangrla.reader.readRaireBallots
import org.cryptobiotic.shangrla.reader.makeCvrsFromRaireBallots
import kotlin.test.Test

class TestRaireReader {
    val dataDir = "src/test/data/"

    @Test
    fun test_read_raire_SFDA() {
        val filename = dataDir + "rla/SFDA2019_PrelimReport12VBMJustDASheets.raire"
        testReadRaireBallots(filename, 11)
        testMakeCvrsFromRaireBallots(filename, 11)
    }

    @Test
    fun test_read_raire_AspenMayor() {
        testReadRaireBallots(dataDir + "raire/Aspen_2009_Mayor.raire")
    }

    fun testMakeCvrsFromRaireBallots(filename: String, limit: Int = Integer.MAX_VALUE) {
        val (raireContests, ballots) = readRaireBallots(filename)
        val cvrs = makeCvrsFromRaireBallots(ballots, 11)
        cvrs.forEach { println(it) }
    }

    fun testReadRaireBallots(filename: String, limit: Int = Integer.MAX_VALUE) {
        val (raireContests, ballots) = readRaireBallots(filename)
        println("RaireContests ${raireContests}")
        println("Cvrs Records")

        // Map<String, MutableMap<String, MutableMap<String, Int>>> ballotId: contestId : candidateId : count
        var ballotIter = ballots.entries.iterator()
        for (idx in 0..limit) {
            if (!ballotIter.hasNext()) break
            val (ballotId, cvr) = ballotIter.next()

            println(buildString {
                append(" Ballot '$ballotId'= ")
                for ((contId, mm2) in cvr) {
                    append("Contest '$contId': ")
                    for ((candId, pref) in mm2) {
                        append("'$candId':$pref, ")
                    }
                }
            })
        }
        if (ballotIter.hasNext()) println(" ...")
    }

}