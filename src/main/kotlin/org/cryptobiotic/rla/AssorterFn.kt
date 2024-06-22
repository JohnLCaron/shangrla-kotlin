package org.cryptobiotic.rla

data class AssorterFn(
    val contest: ContestSimple,
    val upper_bound: Double, // a priori upper bound on the value the assorter can take
    val assort: (CvrSimple) -> Double, // maps votes into [0, upper_bound]
    // tally_pool_means: mean value of the assorter over each tally_pool of CVRs
) {

    //         if use_style and not cvr.has_contest(self.contest.id):
    //            raise ValueError(
    //                f"use_style==True but {cvr=} does not contain contest {self.contest.id}"
    //            )
    //        # assort the MVR
    //        mvr_assort = (
    //            0
    //            if mvr.phantom or (use_style and not mvr.has_contest(self.contest.id))
    //            else self.assort(mvr)
    //        )
    //        # assort the CVR
    //        cvr_assort = (
    //            self.tally_pool_means[cvr.tally_pool]
    //            if cvr.pool and self.tally_pool_means is not None
    //            else int(cvr.phantom) / 2 + (1 - int(cvr.phantom)) * self.assort(cvr)
    //        )
    //        return cvr_assort - mvr_assort

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

