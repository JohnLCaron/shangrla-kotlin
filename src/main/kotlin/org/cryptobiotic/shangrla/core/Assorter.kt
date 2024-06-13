package org.cryptobiotic.shangrla.core

class Assorter(
    val contest: Contest,
    val assort: (CVR) -> Double, // maps a dict of votes into [0, upper_bound]
    val upper_bound: Double, // a priori upper bound on the value the assorter can take
) {
    var winner:  ((CVR) -> Int)? = null
    var loser: ((CVR) -> Int)? = null
    var tally_pool_means: MutableMap<String, Double>? = null

    fun sum(cvr_list: List<CVR>, use_style: Boolean): Double {
        /*
        find the mean of the assorter applied to a list of CVRs

        Parameters
        ----------
        cvr_list: a collection of cast-vote records
        use_style: does the audit use card style information? If so, apply the assorter only to CVRs
        that contain the contest in question.

        Returns
        -------
        the mean value of the assorter over the collection of cvrs. If use_style, ignores CVRs that
        do not contain the contest.
        */

        //     if use_style:
        //    filtr = lambda c: c.has_contest(self.contest.id)
        //    else:
        //    filtr = lambda c: True
        //    return np.sum([self.assort(c) for c in cvr_list if filtr(c)])

        return cvr_list.filter{cvr -> if (use_style) { cvr.has_contest(this.contest.id) } else true }
            .map { this.assort( it) }
            .sum()
    }

    fun overstatement(mvr : CVR, cvr: CVR, use_style: Boolean = true): Double {
        /*
        overstatement error for a CVR compared to the human reading of the ballot

        If use_style, then if the CVR contains the contest but the MVR does
        not, treat the MVR as having a vote for the loser (assort()=0)

        If not use_style, then if the CVR contains the contest but the MVR does not,
        the MVR is considered to be a non-vote in the contest (assort()=1/2).

        Phantom CVRs and MVRs are treated specially:
        A phantom CVR is considered a non-vote in every contest (assort()=1/2).
        A phantom MVR is considered a vote for the loser (i.e., assort()=0) in every contest.

        Parameters
        ----------
        mvr: Cvr the manual interpretation of voter intent
        cvr: Cvr the machine-reported cast vote record

        Returns
        -------
        overstatement: float
        the overstatement error
        */
        //         if use_style and not cvr.has_contest(self.contest.id):
        //            raise ValueError(f'use_style==True but {cvr=} does not contain contest {self.contest.id}')
        //        # assort the MVR
        //        mvr_assort = (0 if mvr.phantom or (use_style and not mvr.has_contest(self.contest.id))
        //                      else self.assort(mvr)
        //                     )
        //         # assort the CVR
        //        cvr_assort = (self.tally_pool_means[cvr.tally_pool] if cvr.pool and self.tally_pool_means is not None
        //                      else int(cvr.phantom)/2 + (1-int(cvr.phantom))*self.assort(cvr)
//                     )
        //        return cvr_assort - mvr_assort

        // sanity check
        if (use_style && !cvr.has_contest(this.contest.id)) {
            throw Exception("use_style==True but ${cvr} does not contain contest ${this.contest.id}")
        }
        // assort the MVR
        val mvr_assort = if (mvr.phantom || (use_style && !mvr.has_contest (this.contest.id))) 0.0
                         else this.assort(mvr)

        // assort the CVR
        val phantomValue = if (cvr.phantom) 1 else 0 // TODO really ? int(cvr.phantom)
        val cvr_assort: Double = if (cvr.pool && this.tally_pool_means != null) this.tally_pool_means!![cvr.tally_pool]!!
                         else phantomValue / 2 + (1 - phantomValue) * this.assort(cvr)

        return cvr_assort - mvr_assort
    }

    fun mean(cvr_list : List<CVR>, use_style: Boolean = true): Double {
        /*
        find the mean of the assorter applied to a list of CVRs

                Parameters:
        ----------
        cvr_list : list
        a list of cast -vote records
                use_style : Boolean
        does the audit use card style information? If so, apply the assorter only to CVRs
        that contain the contest in question .

        Returns:
        ----------
        mean : float
        the mean value of the assorter over the list of cvrs.If use_style, ignores CVRs that
        do not contain the contest .
        */
//        if use_style:
//            filtr = lambda c: c.has_contest(self.contest)
//        else:
//            filtr = lambda c: True
//        return np.mean([self.assorter.assort(c) for c in cvr_list if filtr(c)])

        return cvr_list.filter { cvr -> if (use_style) cvr.has_contest(this.contest.id) else true }
            .map { this.assort( it) }
            .average()
    }

    fun set_tally_pool_means(cvr_list: List<CVR>, tally_pool: Set<String>?, use_style: Boolean) {
        /*
        create dict of pool means for the assorter from a set of CVRs

        Parameters
        ----------
        cvr_list: cvrs from which the sample will be drawn
        tally_pools: Collection [optional]; the labels of the tally groups

        Returns
        -------
        nothing

        Side effects
        ------------
        sets self.tally_pool_means
        */
        //         if not tally_pool:
        //            tally_pool = set(c.tally_pool for c in cvr_list)
        //        tally_pool_dict = {}
        //        for p in tally_pool:
        //            tally_pool_dict[p] = {}
        //            tally_pool_dict[p]['n'] = 0
        //            tally_pool_dict[p]['tot'] = 0
        //        if use_style:
        //            filtr = lambda c: c.has_contest(self.contest.id)
        //        else:
        //            filtr = lambda c: True
        //        for c in [cvr for cvr in cvr_list if filtr(cvr)]:
        //            tally_pool_dict[c.tally_pool]['n'] += 1
        //            tally_pool_dict[c.tally_pool]['tot'] += self.assort(c)
        //        self.tally_pool_means = {}
        //        for p in tally_pool:
        //            self.tally_pool_means[p] = (np.nan if tally_pool_dict[p]['n'] == 0
        //                                  else tally_pool_dict[p]['tot']/tally_pool_dict[p]['n']
        //                                 )

        val tally_set = tally_pool ?: cvr_list.map { it.tally_pool }.toSet()
        val tally_pool_dict = mutableMapOf<String, Pair<Int, Double>>()
        for (p in tally_set) {
            tally_pool_dict[p] = Pair(0, 0.0)
        }

        cvr_list.filter { cvr -> if (use_style) cvr.has_contest(this.contest.id) else true }
            .forEach { cvr ->
                val (n, tot) = tally_pool_dict[cvr.tally_pool]!!
                tally_pool_dict[cvr.tally_pool] = Pair(n+1, tot + this.assort(cvr))
            }

        val pool_means = mutableMapOf<String, Double>()
        for (p in tally_set) {
            pool_means[p] = if (tally_pool_dict[p]!!.first == 0) Double.NaN
                else tally_pool_dict[p]!!.first  / tally_pool_dict[p]!!.second
        }
        this.tally_pool_means = pool_means
    }
}