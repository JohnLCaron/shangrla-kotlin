package org.cyrptobiotic.shangrla.reader

import org.cryptobiotic.shangrla.reader.load_contests_from_raire
import kotlin.test.Test

class TestRaireReader {
    val dataDir = "src/test/data/raire/"

    @Test
    fun test_read_raire_file() {
        val wtf = load_contests_from_raire(dataDir + "Aspen_2009_Mayor.raire")
        println("RaireContests ${wtf.first}")
        println()

        // Map<String, MutableMap<String, MutableMap<String, Int>>>
        println( buildString {
            appendLine("Cvrs Records")
            for ((bid, mm1) in wtf.second) {
                appendLine(" Ballot Id $bid")
                for ((contId, mm2) in mm1) {
                    append("  ContestId $contId :")
                    for ((candId, pref) in mm2) {
                        append("    $candId:$pref, ")
                    }
                    appendLine()
                }
            }
        })
    }
}