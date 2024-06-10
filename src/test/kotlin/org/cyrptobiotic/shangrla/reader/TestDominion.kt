package org.cyrptobiotic.shangrla.reader

import org.cryptobiotic.shangrla.reader.readDominionJsonFromFile
import kotlin.test.*

class TestDominion {
    val dataDir = "src/test/data/Dominion/json/"

    @Test
    fun testReadDominionJsonFromFile() {
        readDominionJsonFromFile(dataDir + "CvrExport_0.json")
        readDominionJsonFromFile(dataDir + "CvrExport_1.json")
    }
}