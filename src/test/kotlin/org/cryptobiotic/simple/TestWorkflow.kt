package org.cryptobiotic.simple

import org.cryptobiotic.shangrla.reader.readRaireBallots
import org.cryptobiotic.shangrla.reader.showRaireBallots

import org.cryptobiotic.shangrla.core.*


import kotlin.test.Test

class TestWorkflow {

    @Test
    fun testWorkflow() {
        //// overall audit information (including the seed) and contest information
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

        //// Read contests, cvrs.
        val raireBallots = readRaireBallots(audit.cvr_file)
        showRaireBallots(raireBallots, 11)
        // // TODO: 293555.0 vs 146662; python has 293555
        val cvrs : List<CvrSimple> = makeCvrsFromRaireBallotsPlurality(raireBallots.second)
        println("ncvrs = ${cvrs.size}")

        val votes: Map<String, Map<String, Int>> = CvrSimple.tabulate_votes(cvrs) // contest -> candidate -> count
        val styles: Map<Set<String>, Int> = CvrSimple.tabulate_styles(cvrs) // style -> cvr count
        val cards: Map<String, Int> = CvrSimple.tabulate_cards_contests(cvrs) // contestId -> ncards aka cvrs

        votes.forEach { key, cands ->
            println("contest ${key} ")
            cands.forEach { println("  ${it}") }
        }

        println("styles")
        styles.forEach { println("  ${it}") }

        println("Cards")
        cards.forEach { println("  ${it}") }

        // make contests from cvrs
        val contests: List<ContestSimple> = ContestSimple.fromVotes(audit, votes, cards)
        contests.forEach { // overide for testing against python
            it.candidates = listOf("16", "17", "18")
            it.estimOverride = EstimFnType.OPTIMAL
        }
        println("Contests")
        contests.forEach { println("  ${it}") }

        //// Create Assertions for every Contest, including an Assorter and NonnegMean for every Assertion
        make_all_assertions(contests)
        audit.check_audit_parameters(contests)

        // Create CVRs for phantom cards
        // skip for now, no phantoms

        // sets margins on the assertions
        val min_margin = set_all_margins_from_cvrs(audit=audit, contests=contests, cvr_list=cvrs)
        // println("minimum assorter margin = ${min_margin}")

        contests.map { contest ->
            println("Assertions for Contest ${contest.id}")
            contest.assertions.forEach { println("  ${it}") }
        }

        //// Set up for sampling
        val sample_size = audit.find_sample_size(contests, cvrs=cvrs)
        println("sample_size = ${sample_size}")

        CvrSimple.assign_sample_nums(cvrs)
        //val sampled_cvr_indices = CvrSimple.consistent_sampling(cvrs, contests=contests)

    }
}