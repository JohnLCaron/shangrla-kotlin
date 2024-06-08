package org.cryptobiotic.shangrla

import kotlin.math.max

enum class SocialChoiceFunction{ APPROVAL, PLURALITY, SUPERMAJORITY, IRV }

enum class Candidates { ALL, ALL_OTHERS, WRITE_IN, NO_CANDIDATE }

class Contest(
    val id: String,
    name: String,
    val risk_limit: Float,
    val cards: Int,
    val choice_function: SocialChoiceFunction,
    n_winners: Int,
    val share_to_win: Float,
    val candidates: List<Candidates>,
    winner: List<Any>,
    assertion_file: String,
    val audit_type: AuditType,
    //test: callable=None,
    g: Float,
    //estim: callable=None,
    //bet: callable=None,
    val use_style: Boolean,
    val assertions: Map<String, Assertion>, // TODO why is this a map?
    val tally: Map<Candidates, Int>,
    var sample_size: Int,
    val sample_threshold: Float,
) {

    fun find_margins_from_tally() {
        /*
    Use the `Contest.tally` attribute to set the margins of the contest's assorters.

    Appropriate only for the social choice functions
    Contest.SOCIAL_CHOICE_FUNCTION.PLURALITY,
    Contest.SOCIAL_CHOICE_FUNCTION.SUPERMAJORITY,
    Contest.SOCIAL_CHOICE_FUNCTION.APPROVAL


    Parameters
    ----------
    None

    Returns
    -------
    None

    Side effects
    ------------
    sets Assertion.margin for all Assertions in the Contest
    */
        for ((_, assn) in this.assertions) {
            assn.find_margin_from_tally()
        }
    }

    fun find_sample_size(audit: Audit, mvr_sample: List<CVR>, cvr_sample: List<CVR>): Int {
        /*
        Estimate the sample size required to confirm the contest at its risk limit.

        This function can be used with or without data, for Audit.AUDIT_TYPE.POLLING,
        Audit.AUDIT_TYPE.CARD_COMPARISON, and Audit.AUDIT_TYPE.ONEAUDIT audits.

        The simulations in this implementation are inefficient because the randomization happens separately
        for every assorter, rather than in parallel.

        Parameters
        ----------
        cvr_sample: list of CVRs; data (or simulated data) to base the sample size estimates on
        mvr_sample: list of MVRs (CVR objects); manually read votes to base the sample size estimates on, if data are available.

        Returns
        -------
        estimated sample size

        Side effects
        ------------
        sets this.sample_size to the estimated sample size
        */

        this.sample_size = 0
        for (a in this.assertions.values) {
            val data = null
            if (mvr_sample != null) {  // process the MVRs/CVRs to get data appropriate to each assertion
                val (data, u) = a.mvrs_to_data(mvr_sample, cvr_sample)
            }
            this.sample_size = max(
                this.sample_size,
                a.find_sample_size(
                    data = data, rate1 = audit.error_rate_1, rate2 = audit.error_rate_2,
                    reps = audit.reps, quantile = audit.quantile, seed = audit.sim_seed
                )
            )
        }
        return this.sample_size
    }
}