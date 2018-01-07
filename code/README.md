# Code provided for the workshop

## The TOY example

The toy example is 100% fictional, and shows a (relatively) simple way to do some basic propensity score analyses, using a quantitative outcome, a binary outcome, and a time-to-event outcome. In no way would I claim that the approach shown here is in any way optimal. It's just a first demo.

- [toy.csv](https://raw.githubusercontent.com/THOMASELOVE/ichps2018/master/code/toy.csv) data file
- [toy_ICHPS.pdf](https://github.com/THOMASELOVE/ichps2018/blob/master/code/toy_ICHPS.pdf) documentation (result of running the R Markdown file below)
- [toy_ICHPS.Rmd](https://raw.githubusercontent.com/THOMASELOVE/ichps2018/master/code/toy_ICHPS.Rmd) R Markdown file

The toy example addresses the following thirteen tasks:

1. Ignoring the covariate information, what is the unadjusted point estimate (and 95\% confidence interval) for the effect of the treatment on each of the three outcomes (`out1.cost`, `out2.event`, and `out3.time`)?
2. Assume that theory suggests that the square of `covA`, as well as the interactions of `covB` with `covC` and `covB` with `covD` should be related to treatment assignment. Fit a propensity score model to the data, using the six covariates (A-F) and the three transformations (A^2^, and the B-C and B-D interactions.) Plot the resulting propensity scores, by treatment group, in an attractive and useful way.
3. Use Rubin's Rules to assess the overlap of the propensity scores and the individual covariates prior to the use of any propensity score adjustments.
4. Use 1:1 greedy matching to match all 70 treated subjects  to control subjects without replacement on the basis of the linear propensity for treatment. Evaluate the degree of covariate imbalance before and after propensity matching for each of the six covariates, and present the pre- and post-match standardized differences and variance ratios for the covariates, as well as the square term and interactions, as well as both the raw and linear propensity score in appropriate plots. Now, build a new data frame containing the propensity-matched sample, and use it to first check Rubin's Rules after matching.
5. Now, use the matched sample data set to evaluate the treatment's average causal effect on each of the three outcomes. In each case, specify a point estimate (and associated 95\% confidence interval) for the effect of being treated (as compared to being a control subject) on the outcome. Compare your results to the automatic versions reported by the Matching package when you include the outcome in the matching process.
6. Now, instead of matching, instead subclassify the subjects into quintiles by the raw propensity score. Display the balance in terms of standardized differences by quintile for the covariates, their transformations, and the propensity score in an appropriate table or plot(s). Are you satisfied? 
7. Regardless of your answer to the previous question, use the propensity score quintile subclassification approach to find a point estimate (and 95% confidence interval) for the effect of the treatment on each outcome. 
8. Now using a reasonable propensity score weighting strategy, assess the balance of each covariate, the transformations and the linear propensity score prior to and after propensity weighting. Is the balance after weighting satisfactory?
9. Using propensity score weighting to evaluate the treatment's effect, developing a point estimate and 95% CI for the average causal effect of treatment on each outcome.
10. Finally, use direct adjustment for the linear propensity score on the entire sample to evaluate the treatment's effect, developing a point estimate and 95\% CI for each outcome.
11. Now, try a double robust approach. Weight, then adjust for linear propensity score.
12. Compare your conclusions about the average causal effect obtained in the following six ways to each other. What happens and why? Which of these methods seems most appropriate given the available information?
    + without propensity adjustment, 
    + after propensity matching, 
    + after propensity score subclassification, 
    + after propensity score weighting, 
    + after adjusting for the propensity score directly, and 
    + after weighting then adjusting for the PS.  
13. Perform a sensitivity analysis for your matched samples analysis and the first outcome (`out1.cost`) if it turns out to show a statistically significant treatment effect.

## Sensitivity Spreadsheet: Beware - this was built in 2008

- [Spreadsheet](https://github.com/THOMASELOVE/ichps2018/blob/master/code/sensitivityspreadsheet_2008.xls) 
- [Documentation](https://github.com/THOMASELOVE/ichps2018/blob/master/code/sensitivitydocumentation_2008.pdf) (10 pages)

# Work from Other People

## Elizabeth Stuart's Propensity Score Software Page

- Elizabeth Stuart's [Propensity Score software page](http://www.biostat.jhsph.edu/~estuart/propensityscoresoftware.html) remains incredibly useful.

## Videos To Learn More About Propensity Scores (and if you visit any of these, you'll find many more)

- Jason Roy (Penn) has a Crash Course in Causality [at Coursera](https://www.coursera.org/learn/crash-course-in-causality).
- Kristin Rankin has a session on [Propensity Score Matching](https://www.youtube.com/watch?v=6ZVIfW6_v1U)
- Ben Lambert has a series of YouTube videos on [Propensity Score Matching](https://www.youtube.com/watch?v=h0UU6trKR0E)
- Doug McKee has posted [An Intuitive Introduction to Propensity Score Matching](https://www.youtube.com/watch?v=ACVyPp1Fy6Y)
- Gary King has a lecture entitled [Why Propensity Scores Should Not Be Used for Matching](https://www.youtube.com/watch?v=rBv39pK1iEs)

## Software

- Noah Griefer's [cobalt vignette](https://cran.r-project.org/web/packages/cobalt/vignettes/cobalt_basic_use.html)
- RAND corporation's [TWANG package materials](https://cran.r-project.org/web/packages/twang/index.html)
- Luke Keele's [Rbounds package](http://www.personal.psu.edu/ljk20/rbounds.html)
- A student of mine, David Ngendahimana, is working on the problem of improving the link between covariate balance and improved sensitivity in selecting a matching procedure. We hope to have something arXived and available for testing on github in the Spring, but David is finishing the other two-thirds of his dissertation first. David is at `dkn18 at case dot edu`.
