package org.cyrptobiotic.shangrla.core

import org.cryptobiotic.shangrla.core.Cvr
import org.cryptobiotic.shangrla.core.Contest
import org.cryptobiotic.shangrla.core.SocialChoiceFunction
import kotlin.math.max
import kotlin.math.min

class AssertionUtils {

    companion object {
        fun find_sample_size(contests: List<Contest>, sample_size_function: (Double, Double) -> Int): Int {
            /*
            Find initial sample size : maximum across assertions for all contests .

            Parameters:
            -----------
            contests : dict of dicts
            assertion : dict of dicts
            sample_size_function : callable; takes two parameters, the margin and the risk limit; returns a sample size

                    Returns:
            --------
            size : int; sample size expected to be adequate to confirm all assertions
            */

            var sample_size = 0
            for (c in contests) {
                val risk = c.risk_limit
                for (a in c.assertions.values) {
                    val margin = a.margin
                    val n = sample_size_function(margin!!, risk)
                    sample_size = max(sample_size, n)
                }
            }
            return sample_size
        }

        fun check_audit_parameters(risk_function: String, g: Double, error_rate: Double, contests: Map<String, Contest>) {
            /*
            Check whether the audit parameters are valid; complain if not.
    
            Parameters:
            ---------
            risk_function : string; the name of the risk-measuring function for the audit
            g : float in [0, 1); padding for Kaplan-Markov or Kaplan-Wald
            error_rate : float; expected rate of 1-vote overstatements
            contests : dict of dicts; contest-specific information for the audit
    
            Returns:
            --------
        */
            if (risk_function in listOf("kaplan_markov", "kaplan_wald")) {
                assert (g >=0) { "g must be at least 0" }
                assert (g < 1) { "g must be less than 1" }
            }
            assert (error_rate >= 0) {"expected error rate must be nonnegative"}
            for (contestId in contests.keys) {
                val contest = contests[contestId]!!
                assert (contest.risk_limit > 0) {"risk limit must be nonnegative in $contest" }
                assert (contest.risk_limit < 1) { "risk limit must be less than 1 in $contest" }
                assert (contest.choice_function in listOf(
                    SocialChoiceFunction.IRV, SocialChoiceFunction.PLURALITY, SocialChoiceFunction.SUPERMAJORITY))
                    { "unsupported choice function ${ contest.choice_function} in $contest" }
                assert (contest.n_winners <= contest.candidates.size) { "fewer candidates than winners in $contest" }
                assert (contest.reported_winners.size == contest.n_winners)
                    { "number of reported winners does not equal n_winners in $contest" }
                for (w in contest.reported_winners) {
                    assert (w in contest.candidates) { "reported winner $w is not a candidate in $contest" }
                }
                if (contest.choice_function in listOf(SocialChoiceFunction.IRV, SocialChoiceFunction.SUPERMAJORITY)) {
                    assert(contest.n_winners == 1) { "${contest.choice_function} can have only 1 winner in $contest" }
                }
                /* if (contest.choice_function == SocialChoiceFunction.IRV) {
                    assert(contest.assertion_file != null) { "IRV contest $contest requires an assertion file" }
                }
                // TODO
                 */
                if (contest.choice_function == SocialChoiceFunction.SUPERMAJORITY) {
                    assert(contest.share_to_win >= 0.5) { "super-majority contest requires winning at least 50% of votes in $contest" }
                }
            }
        }

        fun find_margins(contests: List<Contest>, cvr_list : List<Cvr>, use_style : Boolean): Double {
            /*
            Find all the assorter margins in a set of Assertions. Updates the dict of dicts of assertions,
            and the contest dict.

            Appropriate only if cvrs are available. Otherwise, base margins on the reported results.

            This function is primarily about side-effects.

            Parameters:
            -----------
            contests : dict of contest data, including assertions
            cvr_list : list list of cvr objects
            use_style : bool flag indicating the sample will use style information to target the contest

            Returns:
            --------
            min_margin : float smallest margin in the audit
            */
            //     in_margin = np.infty
            //    for c in contests:
            //        contests[c]['margins'] = {} // TODO adding a field?
            //        for a in contests[c]['assertions']:
            //            # find mean of the assertion for the CVRs
            //            amean = contests[c]['assertions'][a].assorter_mean(cvr_list, use_style=use_style)
            //            if amean < 1/2:
            //                warn(f"assertion {a} not satisfied by CVRs: mean value is {amean}")
            //            margin = 2*amean-1
            //            contests[c]['assertions'][a].margin = margin // TODO wtf?
            //            contests[c]['margins'].update({a: margin})
            //            min_margin = np.min([min_margin, margin])
            //    return min_margin

            var min_margin = Double.POSITIVE_INFINITY
            for (contest in contests) {
                //contest.margins = {} TODO
                for (a in contest.assertions.values) {
                    // find mean of the assertion for the CVRs
                    val amean = a.assorter.mean(cvr_list, use_style = use_style)
                    if (amean < .5) {
                        println("assertion ${a} not satisfied by CVRs: mean value is ${amean}")
                    }
                    val margin = 2 * amean - 1
                    a.margin = margin
                    // contest.margins.update({ a: margin }) TODO
                    min_margin = min(min_margin, margin)
                }
            }
            return min_margin
        }

    }
}