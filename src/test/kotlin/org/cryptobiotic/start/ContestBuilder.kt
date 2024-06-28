package org.cryptobiotic.start

class ContestBuilder(
    val id: String,
    val risk_limit: Double,
    var ncvrs: Int? = 0,
    var cards: Int?,
    val n_winners: Int,
    val candidates: List<String>,
    val winners: List<String>,
    val share_to_win: Double? = null,
)