package org.cryptobiotic.start

class Assorter(
    val contest: Contest,
    val upper_bound: Double, // a priori upper bound on the value the assorter can take
    val assort: (Cvr) -> Double, // maps votes into [0, upper_bound]
) {
    // Compute the arithmetic mean of the assort value over the cvrs that have this contest, // eq 2
    fun mean(cvrs : List<Cvr>): Double {
        return cvrs.filter { cvr -> cvr.has_contest(this.contest.id) }
            .map { this.assort( it) }
            .average()
    }

    fun overstatement(mvr : Cvr, cvr: Cvr): Double {
        // TODO paper reference
        val mvr_assort = if (!mvr.has_contest (this.contest.id)) 0.0 else this.assort(mvr)
        return this.assort(cvr) - mvr_assort

        // assort the CVR
        //        cvr_assort = (
        //            self.tally_pool_means[cvr.tally_pool]
        //            if cvr.pool and self.tally_pool_means is not None
        //            else int(cvr.phantom) / 2 + (1 - int(cvr.phantom)) * self.assort(cvr)
        //        )
        // val cvr_assort: Double = if (cvr.pool && this.tally_pool_means != null) this.tally_pool_means!![cvr.tally_pool]!!
        // else phantomValue / 2 + (1 - phantomValue) * this.assort(cvr)
        // val phantomValue = if (cvr.phantom) 1.0 else 0.0 // TODO really ? int(cvr.phantom)
        //val temp = phantomValue / 2 + (1 - phantomValue) * this.assort(cvr)
        //val cvr_assort = if (cvr.phantom) .5 else this.assort(cvr)

    }
}