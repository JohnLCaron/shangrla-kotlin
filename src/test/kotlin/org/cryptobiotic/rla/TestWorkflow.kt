package org.cryptobiotic.rla

import org.cryptobiotic.shangrla.reader.readRaireBallots
import org.cryptobiotic.shangrla.reader.showRaireBallots
import org.cryptobiotic.shangrla.reader.makeCvrsFromRaireBallots
import org.cryptobiotic.shangrla.reader.makeCvrsFromRaireBallotsPlurality
import org.cyrptobiotic.shangrla.core.AssertionUtils.Companion.check_audit_parameters

import org.cryptobiotic.shangrla.core.*


import kotlin.test.Test

class TestWorkflow {

    @Test
    fun testWorkflow() {
        //// Read overall audit information (including the seed) and contest information
        val audit: AuditSimple = AuditSimple(
            cvr_file=      "/home/stormy/dev/github/rla/shangrla-kotlin/src/test/data/rla/SFDA2019_PrelimReport12VBMJustDASheets.raire",
            manifest_file= "/home/stormy/dev/github/rla/shangrla-kotlin/src/test/data/rla/N19.ballotmanifest.VBM.11-14.xlsx",
            error_rate_1=  0.001,
            max_cards=146662)

        //// Read ballot manifest, get total ballots
        //val manifest = DataFrame.read(audit.manifest_file)
        //manifest.schema()

        //val totalBallots = manifest["Total Ballots"] as ValueColumn<Double>
        //totalBallots.toList().sum()

        //// Read cvrs.
        val raireBallots = readRaireBallots(audit.cvr_file)
        showRaireBallots(raireBallots, 11)
        val cvrs : List<CvrSimple> = makeCvrsFromRaireBallotsPlurality(raireBallots.second)
        println("ncvrs = ${cvrs.size}")

        // cvrs.forEach { println(it) }

        // make contests
        // // TODO: 293555.0 vs 146662; python has 293555
        val votes: Map<String, Map<String, Int>> = CvrSimple.tabulate_votes(cvrs)
        val styles: Map<Set<String>, Int> = CvrSimple.tabulate_styles(cvrs)
        val cards: Map<String, Int> = CvrSimple.tabulate_cards_contests(cvrs)

        println("Cards")
        cards.forEach { println("  ${it}") }

        println("styles")
        styles.forEach { println("  ${it}") }

        // TODO just use RaireContests?
        val contests: List<ContestSimple> = ContestSimple.fromVotes(audit, votes, cards)
        contests.forEach { it.estimOverride = EstimFnType.OPTIMAL }
        println("Contests")
        contests.forEach { println("  ${it}") }

        //// Create Assertions for every Contest, including an Assorter and NonnegMean for every Assertion

        make_all_assertions(contests)
        audit.check_audit_parameters(contests)

        // Create CVRs for phantom cards
        // skip for now, no phantoms

        // sets margins on the assertions
        val min_margin = set_all_margins_from_cvrs(audit=audit, contests=contests, cvr_list=cvrs)
        println("minimum assorter margin = ${min_margin}")

        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }

        //// Set up for sampling
        val sample_size = audit.find_sample_size(contests, cvrs=cvrs)
        println("sample_size = ${sample_size}")

    }
}