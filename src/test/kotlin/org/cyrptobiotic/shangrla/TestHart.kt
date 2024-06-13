package org.cyrptobiotic.shangrla

import org.cryptobiotic.shangrla.reader.readHartXmlFromFile
import kotlin.test.*

class TestHart {
    val dataDir = "src/test/data/Hart/xml/"

    @Test
    fun testReadHartXml() {
        readHartXmlFromFile(dataDir + "test_Hart_CVR_1.xml")
        readHartXmlFromFile(dataDir + "test_Hart_CVR_2.xml")
    }
}