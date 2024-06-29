package org.cryptobiotic.rla

// An assorter A assigns a nonnegative value to each ballot card, depending on the marks
// the voter made on that ballot card.

// value “1” to a ballot card if it has a mark for Alice but not for Bob; assign the value “0”
//if the card has a mark for Bob but not for Alice; assign the value 1/2, otherwise (e.g., if
//the card has an overvote or an undervote in this contest or does not contain the contest).
//Then Alice beat Bob iff the average value of the assorter for the full set of cast ballot
//cards is greater than 1/2: then Alice got more than 50% of the valid votes.
data class AssorterFn(
    val contest: ContestSimple,
    val upper_bound: Double, // a priori upper bound on the value the assorter can take
    val assort: (CvrSimple) -> Double, // maps votes into [0, upper_bound]
    // tally_pool_means: mean value of the assorter over each tally_pool of CVRs
) {

    //         Compute the arithmetic mean along the specified axis.
    fun mean(cvr_list : List<CvrSimple>, use_style: Boolean = true): Double {
          val result = cvr_list.filter { cvr -> if (use_style) cvr.has_contest(this.contest.id) else true }
            .map { this.assort( it) }
            .average()
        return result
    }

    fun overstatement(mvr : CvrSimple, cvr: CvrSimple, use_style: Boolean = true): Double {
        // sanity check
        if (use_style && !cvr.has_contest(this.contest.id)) {
            throw Exception("use_style==True but Cvr '${cvr.id}' does not contain contest '${this.contest.id}'")
        }
        // assort the MVR
        val mvr_assort = if (mvr.phantom || (use_style && !mvr.has_contest (this.contest.id))) 0.0
        else this.assort(mvr)

        // assort the CVR
        //        cvr_assort = (
        //            self.tally_pool_means[cvr.tally_pool]
        //            if cvr.pool and self.tally_pool_means is not None
        //            else int(cvr.phantom) / 2 + (1 - int(cvr.phantom)) * self.assort(cvr)
        //        )
        // val cvr_assort: Double = if (cvr.pool && this.tally_pool_means != null) this.tally_pool_means!![cvr.tally_pool]!!
        // else phantomValue / 2 + (1 - phantomValue) * this.assort(cvr)
        val phantomValue = if (cvr.phantom) 1.0 else 0.0 // TODO really ? int(cvr.phantom)
        val temp = phantomValue / 2 + (1 - phantomValue) * this.assort(cvr)
        val cvr_assort = if (cvr.phantom) .5 else this.assort(cvr)
        require(temp == cvr_assort)

        return cvr_assort - mvr_assort
    }

}

/**
    Construct assertions that imply the winner(s) got more votes than the loser(s).
    The assertions are that every winner beat every loser.
    There are len(winner)*len(loser) pairwise assertions in all.
    Return Map<String, AssorterFn>, where key = winr + ' v ' + losr
*/

