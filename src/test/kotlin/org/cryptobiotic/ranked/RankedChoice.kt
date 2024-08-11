package org.cryptobiotic.ranked

import org.junit.jupiter.api.Test

class RankedChoice {
    val john = listOf("booker", "buttigieg", "whitmer", "klobuchar", "bernie", "harris", "polis",)
    val jerry = listOf("buttigieg", "whitmer", "harris", "newsome", "kelly",)
    val jeff = listOf("whitmer", "harris", "raimondo", "buttigieg","newsome", "shapiro")
    val jean = listOf("whitmer", "harris", "shapiro", "booker", "buttigieg",)
    val jim = listOf("whitmer", "warren", "klobuchar", "harris", "inslee", "newsom","booker","brown","coons",)
    val joe = listOf("bernie", "kusinich", "gabbard", "inslee")
    val mark = listOf("mark", "bernie", "harris", "buttigieg", "whitmer", "crockett")

    val names = listOf("john", "jerry", "jeff", "jean", "jim", "joe", "mark")

    val ranks = listOf(john, jerry, jeff, jean, jim, joe, mark)
    val voters = ranks.mapIndexed { idx, it -> Voter(names[idx], it) }

    @Test
    fun runTest() {
        println("Voters")
        voters.forEach {
            println("  ${it.name}: ${it.votes}")
        }
        println()

        val candidates = mutableSetOf<String>()
        voters.forEach { candidates.addAll( it.votes )}
        val scores = mutableMapOf<String, Double>()
        candidates.forEach { candidate ->
            scores[candidate] = voters.map { it.getCandidateScore(candidate) }.sum()
        }
        val sorted : List<Pair<String, Double>> = scores.toList().sortedBy { it.second }
        println("Candidates Global Score")
        sorted.forEach {
            println("  ${it.first}: ${it.second}")
        }
        println()

        runRankedChoice()
    }

    fun runRankedChoice() {
        var winner: String? = null
        var round = 1
        while (winner == null) {
            val sorted : List<Pair<String, Int>> = runRound(round)
            winner = checkForWinner(sorted)
            if (winner == null) {
                val loser = chooseLoser(sorted)
                println("Eliminate $loser")
                voters.forEach { it.eliminate(loser) }
                println()
            }
            round++
        }

        println("*** The winner is $winner")
    }

    fun runRound(round: Int): List<Pair<String, Int>>  {
        val votes = voters.filter { it.getCurrentVote() != null }.map { it.getCurrentVote()!! }
        val voteCount = mutableMapOf<String, Int>()
        votes.forEach { vote ->
            val count = voteCount.getOrPut(vote) { 0 }
            voteCount[vote] = count + 1
        }
        val sorted : List<Pair<String, Int>> = voteCount.toList().sortedBy { -it.second }
        println("Round $round")
        sorted.forEach { println(" ${it}") }
        return sorted
    }

    fun checkForWinner(sorted : List<Pair<String, Int>>) : String? {
        val total = voters.size
        val best = sorted[0].second
        if (best > (total + 1) / 2) {
            return sorted[0].first
        }
        return null
    }

    fun chooseLoser(sorted : List<Pair<String, Int>>) : String {
        val loserScore = mutableMapOf<String, Double>()
        sorted.forEach{ (candidate, _) ->
            loserScore[candidate] = voters.map { it.getCandidateScore(candidate) }.sum()
        }
        println(" Losers")
        val sorted : List<Pair<String, Double>> = loserScore.toList().sortedBy { it.second }
        sorted.forEach { println("  ${it}") }
        return sorted[0].first
    }

    class Voter( val name: String, val votes: List<String>) {
        var currentVote = 0
        fun getCurrentVote() : String? = if (currentVote < votes.size) votes.get(currentVote) else {
            null
        }

        // TODO this is fixed, could be computed once
        fun getCandidateScore(candidate: String): Double {
            val idx: Int = votes.indexOf(candidate)
            val score = if (idx < 0) 0.0 else 1.0 / (idx + 1)
            // println("  voter $name candidate $candidate score $score ")
            return score
        }

        fun eliminate(loser: String) {
            if ((currentVote < votes.size) && (votes[currentVote] == loser)) {
                currentVote++
                if (currentVote == votes.size) println("Voter $name is dead")
            }
        }
    }
}