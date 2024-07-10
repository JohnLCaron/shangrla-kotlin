package org.cryptobiotic.simple

import org.cryptobiotic.shangrla.core.*

class ContestCalc(
    val contest: ContestSimple,
    val betFn: BetFn? = null,
    val estimFn: EstimatorFn? = null,
    val testFn: TestFn? = null,
    val g: Double = 0.1,
    val winners: List<String>,
    val risk_limit: Double = 0.05,
    val share_to_win: Double = 0.5,
    val use_style: Boolean = true,

    var assertions: MutableMap<String, AssertionSimple> = mutableMapOf(), // key = winr + " v " + losr
    var ncvrs: Int = 0,
    var ncards: Int = 0,
    var sample_size: Int? = null,
    var sample_threshold: Double? = null,
    val tally: MutableMap<String, Int> = mutableMapOf(), // candidate name -> vote count
)