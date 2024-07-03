package org.cryptobiotic.start

class Algorithm(
    val risk_limit: Double,
    val u: Double, // max value of assorter
    val N: Int,

) {

    fun run() {
        //• Set audit parameters:
        // – Select the risk limit α ∈ (0, 1); decide whether to sample with or without replacement.
        // – Set u as appropriate for the assertion under audit.
        // – Set N to the number of ballot cards in the population of cards from which the sample is drawn.
        // – Set η0 .
        //    For polling audits, η0 could be the reported mean value of the assorter.
        //	For instance, for the assertion corresponding to checking whether w got more votes than ℓ,
        //	 η0 = (Nw + Nc /2)/N , where Nw is the number of votes reported for w , Nℓ is the
        //	number of votes reported for ℓ, and Nc = N − Nw − Nℓ is the number of ballot cards
        //	reported to have a vote for some other candidate or no valid vote in the contest.
        //    For comparison audits, η0 can be based on assumed or historical rates of overstatement errors.
        // – Define the function to update η based on the sample,
        //	e.g, η(i, X i−1 ) = ((d * η0 + S)/(d + i − 1) ∨ (eps(i) + µi )) ∧ u,    (2.5.2, eq 14, "truncated shrinkage")
        //	where S = Sum i−1 k=1 (Xk) is the sample sum of the first i − 1 draws
        //	and eps(i) = c/ sqrt(d + i − 1)
        //	set any free parameters in the function (e.g., d and c in this example). The only requirement is that
        //	η(i, X i−1 ) ∈ (µi , u), where µi := E(Xi |X i−1 ) is computed under the null.
        //
        //• Initialize variables
        // – j ← 0: sample number
        // – T ← 1: test statistic
        // – S ← 0: sample sum
        // – m = 1/2: population mean under the null
        //
        //• While T < 1/α and not all ballot cards have been audited:
        // – Draw a ballot at random
        // – j ←j +1
        // – Determine Xj by applying the assorter to the selected ballot card (and the CVR, for comparison audits)
        // – If m < 0, T ← ∞. Otherwise, T ← T / u * ( Xj * η(j,S)/m + (u - Xj) * (u−η(j,S))/(u-m))
        // – S ← S + Xj
        // – If the sample is drawn without replacement, m ← (N/2 − S)/(N − j + 1)
        // – If desired, break and conduct a full hand count instead of continuing to audit.
        //
        //• If a full hand count is conducted, its results replace the reported results if they differ.
    }
}