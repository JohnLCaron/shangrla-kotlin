# shangrla-kotlin

A port of [SHANGRLA](https://github.com/pbstark/SHANGRLA) to kotlin, for the purpose of learning RLA.

## Code
    shangrla/core is more or less a port from the SHANGRLA python
    test/start starts to refactors the core classes
        AlphaMart, ShrinkTrunc
        KaplanWald (not used)
    rla directly implements the algorithm described in ALPHA paper

    so, compare shangrla/core or test/start against python
    compare rla against test/start. test/rla/Test*

## Papers

    SHANGRLA	Sets of Half-Average Nulls Generate Risk-Limiting Audits: SHANGRLA	Stark; 24 Mar 2020

    ALPHA	    ALPHA: Audit that Learns from Previously Hand-Audited Ballots	Stark; Jan 7, 2022

    Sweeter     Sweeter than SUITE: Supermartingale Stratified Union-Intersection Tests of Elections Jacob V. Spertus and Philip B. Stark; 25 Jul 2022
            combine the union-intersection tests in SUITE and the nonnegative supermartingale (NNSM) tests in ALPHA.
            Seems to be about stratified audits
            https://arxiv.org/abs/2207.03379

    ONEAudit    ONEAudit: Overstatement-Net-Equivalent Risk-Limiting Audit	Stark, P.B., 6 Mar 2023.
            https://github.com/pbstark/ONEAudit

## Notes


SHANGRLA is the framework for RLAs
    assorter: assign a numeric value to a Cvr: assort: (Cvr) -> Double in (0, upper)
    estimator: estimate the population mean of the assorter values: estimate: (x: DoubleArray) -> DoubleArray // LOOK why array?
    assertions : assert that the mean is greater than 1/2: "half-average assertions".
    each audit has a list of half-average assertions
    audit (for a contest): iterative process of picking ballots and checking if all the assertions are true.
    the "test function" is the statistical method to test if the assertion is true. aka "risk function".

enum class TestFnType {
    ALPHA_MART,         // test=NonnegMean.alpha_mart
    BETTING_MART,       // bet=NonnegMean.agrapa, bet=NonnegMean.fixed_bet
    KAPLAN_KOLMOGOROV,  // NonnegMean.kaplan_kolmogorov
    KAPLAN_MARKOV,      // NonnegMean.kaplan_markov
    KAPLAN_WALD,        // NonnegMean.kaplan_wald
    WALD_SPRT,          // NonnegMean.wald_sprt
}

ALPHA is a test function, can use different assorter functions?, and use different estimators ??
I think theyre all tied in together, you dont just mix and match.

ALPHA takes EstimFn that must only depend on the previous j-1 samples
    // estimate the population mean for the jth sample from the previous j-1 samples

estimators
    TruncShrinkage truncated shrinkage // estim=NonnegMean.shrink_trunc
    FixedAlternativeMean: fixed alternative that the original population mean is eta // estim=NonnegMean.fixed_alternative_mean

enum class AuditType { POLLING, CARD_COMPARISON, ONEAUDIT }
enum class SocialChoiceFunction{ APPROVAL, PLURALITY, SUPERMAJORITY, IRV }

