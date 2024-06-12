package org.cryptobiotic.shangrla

import org.cryptobiotic.shangrla.core.SocialChoiceFunction

class ContestSimple(
    val id: String,
    val risk_limit: Double,
    val cards: Int?,
    val choice_function: SocialChoiceFunction = SocialChoiceFunction.PLURALITY,
    val n_winners: Int,
    val candidates: List<String>,
    val reported_winners: List<String>, // TODO was called winners?
    val share_to_win: Double? = null,
)