# stratified_athena

This repository is my exploratory work on stratified audits, under Professor Vora. Specifically I'm working on 2-strata audits that make use of the Athena class of ballot polling audits for the ballot polling stratum in order to increase the efficiency (and therefore feasability) of stratified audits.

A lot of this code is pulled straight from [the SUITE work](https://github.com/pbstark/CORLA18/tree/master/code) and some is pulled from [the minerva branch of r2b2](https://github.com/gwexploratoryaudits/r2b2/tree/minerva).

This code is messy and ugly. For now I'm more focused on getting things working and producing meaningful results, less focused on code ready to be integrated into r2b2....

Please ask me any questions that come up. obroadrick@gwu.edu

# Examples

## Stopping probability for round sizes

Below is a plot of stopping probability as a function of ballot polling stratum round size for a 2-strata audit using Minerva in the ballot polling stratum. This contest has a margin of 2% for a total of 104000 relevant votes with only 5% of those votes in the polling stratum; this allows a better view of the jagged curve on the plot.

![Stopping Probability](imgs/stopping_probability.png)

## Finding round sizes

Here is a plot illustrating the method I use to find the polling stratum round size that achieves a desired probability of stopping. Plotted over possible values for k (winner votes drawn in the round) are both the combined pvalue produced by that k (using 2-strata SUITE with Minerva for the polling stratum and a fixed 750 comparisons) as well as the probability of that k under the alternative hypothesis (a binomial distribution centered around the announced tally for the winner. All the numbers are based on the first example SUITE uses, a relatively tight margin of .0192.

![Round Size Visualization Plot](imgs/round_size_visualization.png)

From the plot, you can tell that this round size of 84 achieves a 90% probability of stopping because the value kmax such that Pr[k>=kmax | H_a] >= .9 and the pvalue corresponding to kmax <= .1 = alpha, and as k increases from kmax, the pvalue only decreases. Another interesting way to look at this is a probability distribution over possible combined pvalues for the given round sizes as seen in the not-so-well-labelled graph below.

![Distribution over Pvalues](imgs/dist_over_pvalues.png)

## Round sizes for various polling stratum sizes

Here is a plot of the ballot polling stratum minimum first round sizes that achieve 90% stopping probability for a 2-strata audit with ~2% margin, a fixed 750 ballots sampled in the comparison stratum, and various stratum sizes (on x-axis as percantage of total relevant ballots). Note that the 2-strata audit can reduce the number of ballots polled if a contest-wide polling audit were performed for significant polling stratum sizes (over half in this example). A tighter margin reduces the polling stratum size that still yields gains in efficiency, but even for this relatively tight race (~2% margin) more than half the relevant ballots can be in the polling stratum and still reduce the polling round size.

![Round Size Plot](imgs/stratified_plot.png)

## Variables that affect round size
Below are the results of tests in which I used the same numbers as in the previous example, but I have isolated and changed one variable to see how the round size is affected.

### Affect of margin on round size
Below is plot of Minerva stratum round sizes for several different margins. 

![Many Margins Minerva](imgs/many_margins_minerva.png)

And here is the same plot for R2 Bravo; it has similar traits, just with higher round sizes.

![Many Margins R2 Bravo](imgs/many_margins_r2bravo.png)

### AFfect of more comparisons on polling round size
While it is reasonable to expect that increasing the number of comparisons would decrease the required polling sample size, this is often not the case. For relatively high perctentage polling stratum sizes, the SUITE overstatement error allocation, lambda, is 0. Therefore, increasing the sample size in the comparison stratum doesn't actually reduce the comparison stratum pvalue, rather it remains near 1, while the pressure of handling the whole of the error is left to the less efficient polling stratum. The following plot of the same contest round sizes for both the original 750 ballot comparison sample as well as a 1500 ballot comparison sample illustrates the point. While a reduction in round size is seen for small polling stratum sizes, the rest of the round sizes are the same for both 750 comparisons and 1500 comparisons.

![Double Comparison Sample Size Plot](imgs/double_comparisons.png)

### Plot of polling round sizes for different comparison sample sizes
As seen in the previous graph, for a sufficient number of comparisons, increasing the number of comparisons no longer decreases the polling stratum first round size. For such audits, lambda = 0 and all overstatement error is allocated to the polling stratum. Below is a plot of polling stratum first round size for various comparison sample sizes. Such a plot can be used to determine roughly the minimum comparison sample size required to achieve lambda = 0, thus minimizing the polling stratum sample size with as few comparisons as possible.

![Various Comparison Sample Sizes](imgs/comparisons.png)

