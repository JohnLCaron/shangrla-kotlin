package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.reader.readRaireBallots
import org.cryptobiotic.shangrla.reader.showRaireBallots
import org.cryptobiotic.shangrla.reader.makeCvrsFromRaireBallots

import kotlin.test.Test

class TestWorkflow {

    @Test
    fun testWorkflow() {
        //// Read overall audit information (including the seed) and contest information
        val audit: AuditSimple = AuditSimple(
            cvr_file=      "/home/stormy/dev/github/rla/shangrla-kotlin/src/test/data/rla/SFDA2019_PrelimReport12VBMJustDASheets.raire",
            manifest_file= "/home/stormy/dev/github/rla/shangrla-kotlin/src/test/data/rla/N19.ballotmanifest.VBM.11-14.xlsx",
            max_cards=293555)

        //// Read ballot manifest, get total ballots
        //val manifest = DataFrame.read(audit.manifest_file)
        //manifest.schema()

        //val totalBallots = manifest["Total Ballots"] as ValueColumn<Double>
        //totalBallots.toList().sum()

        //// Read cvrs.
        val raireBallots = readRaireBallots(audit.cvr_file)
        showRaireBallots(raireBallots, 11)
        val cvrs : List<CvrSimple> = makeCvrsFromRaireBallots(raireBallots.second, 11)
        cvrs.forEach { println(it) }

        // make contests
        // // TODO: 293555.0 vs 146662; python has 293555
        val votes: Map<String, Map<String, Int>> = CvrSimple.tabulate_votes(cvrs)
        val styles: Map<Set<String>, Int> = CvrSimple.tabulate_styles(cvrs)
        val cards: Map<String, Int> = CvrSimple.tabulate_cards_contests(cvrs)

        // TODO just use RaireContests?
        val contests: List<ContestSimple> = ContestSimple.fromVotes(audit, votes, cards)
        println("Contests")
        contests.forEach { println("  ${it}") }

        //// Create Assertions for every Contest, including an Assorter and NonnegMean for every Assertion

        make_all_assertions(contests)
        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }
    }
}