{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "ebdid##syntax"}{...}
{viewerjumpto "Description" "ebdid##description"}{...}
{viewerjumpto "Options" "ebdid##options"}{...}
{viewerjumpto "Stored results" "ebdid##stored"}{...}
{viewerjumpto "Examples" "ebdid##examples"}{...}

{title:Title}

{phang}
{bf:ebdid} {hline 2} Empirical Bayes shrinkage for difference-in-differences


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:ebdid}
{depvar}
{ifin}{cmd:,}
{cmdab:g:roup(}{varname}{cmd:)}
{cmdab:t:ime(}{varname}{cmd:)}
{cmdab:tr:eated(}{varname}{cmd:)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent :* {opth g:roup(varname)}}panel group identifier{p_end}
{p2coldent :* {opth t:ime(varname)}}time period variable{p_end}
{p2coldent :* {opth tr:eated(varname)}}binary treatment indicator (0/1){p_end}
{synopt:{opt br:eps(#)}}number of bootstrap replications; default is {cmd:breps(200)}{p_end}
{synopt:{opt vc:alib(string)}}noise calibration method; default is {cmd:vcalib(untreated)}{p_end}
{synopt:{opt l:evel(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{p2colreset}{...}
{pstd}* {opt group()}, {opt time()}, and {opt treated()} are required.


{marker description}{...}
{title:Description}

{pstd}
{cmd:ebdid} estimates group-level treatment effects in a
difference-in-differences design using James-Stein / empirical Bayes
shrinkage.  Each treated group's DiD estimate is shrunk toward the
grand mean, with the degree of shrinkage determined by the
signal-to-noise ratio.

{pstd}
The estimator reduces mean squared error relative to OLS when there are
4 or more treated groups (the SURE-optimal James-Stein constant for
shrinkage toward an estimated grand mean requires p >= 4).  It is
adaptive: when treatment effects are homogeneous, estimates are pooled
aggressively; when effects are heterogeneous, the estimator recovers
OLS.

{pstd}
Inference is based on a parametric bootstrap that accounts for
uncertainty in the shrinkage factor, the grand mean, and the individual
group estimates.

{pstd}
Staggered treatment adoption is supported: {cmd:treated()} should be 1
in periods when a group is treated and 0 otherwise.  The first period
with {cmd:treated}==1 defines treatment onset for each group.  Groups
that are never treated serve as controls.


{marker options}{...}
{title:Options}

{phang}
{opth group(varname)} specifies the panel group identifier (e.g., state,
country, district).  Required; must be numeric.

{phang}
{opth time(varname)} specifies the time period variable.  Required; must
be numeric.

{phang}
{opth treated(varname)} specifies a binary (0/1) indicator equal to 1
when and where treatment is active.  For staggered designs, this
variable is 0 in pre-treatment periods and 1 from the treatment start
onward for each treated group.  Groups where {cmd:treated} is always 0
are never-treated controls.  Required; must be numeric.

{phang}
{opt breps(#)} sets the number of parametric bootstrap replications
used to compute standard errors.  Default is 200.

{phang}
{opt vcalib(string)} selects the noise variance calibration method:

{p 12 16 2}
{cmd:untreated} (default) uses leave-one-out placebos from never-treated
groups across all time periods.  This is the most data-rich estimator
and captures post-period noise conditions.

{p 12 16 2}
{cmd:treated} uses placebos from treated groups in pre-treatment periods
only.

{p 12 16 2}
{cmd:both} pools both sources.

{phang}
{opt level(#)} sets the confidence level for confidence intervals.
Default is {cmd:c(level)}, typically 95.


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:ebdid} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(B_hat)}}shrinkage factor (0 = full pooling, 1 = OLS){p_end}
{synopt:{cmd:e(V_hat)}}estimated noise variance{p_end}
{synopt:{cmd:e(tau_bar)}}grand mean of group DiD estimates{p_end}
{synopt:{cmd:e(S2)}}sum of squared deviations of group DiDs{p_end}
{synopt:{cmd:e(G_treat)}}number of treated groups{p_end}
{synopt:{cmd:e(G_ctrl)}}number of never-treated control groups{p_end}
{synopt:{cmd:e(breps)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ebdid}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(vcalib)}}noise calibration method used{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x G row vector of EB-shrunk estimates (tau_eb){p_end}
{synopt:{cmd:e(V)}}G x G diagonal matrix of bootstrap variances{p_end}
{synopt:{cmd:e(tau_ols)}}G x 1 vector of unshrunk OLS DiD estimates{p_end}
{synopt:{cmd:e(group_id)}}G x 1 vector of group identifiers{p_end}
{synopt:{cmd:e(se_boot)}}G x 1 vector of bootstrap standard errors{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup: panel of 20 groups, 10 time periods, first 10 groups treated
from period 7 onward{p_end}

{phang2}{cmd:. ebdid y, group(state) time(year) treated(D)}{p_end}

{pstd}With more bootstrap replications:{p_end}

{phang2}{cmd:. ebdid y, group(state) time(year) treated(D) breps(500)}{p_end}

{pstd}Using treated pre-period placebos for noise calibration:{p_end}

{phang2}{cmd:. ebdid y, group(state) time(year) treated(D) vcalib(treated)}{p_end}

{pstd}Access stored results after estimation:{p_end}

{phang2}{cmd:. matrix list e(b)}{p_end}
{phang2}{cmd:. di "Shrinkage factor: " e(B_hat)}{p_end}
{phang2}{cmd:. matrix list e(tau_ols)}{p_end}


{title:References}

{pstd}
James, W. and C. Stein. 1961. Estimation with quadratic loss.
{it:Proceedings of the Fourth Berkeley Symposium on Mathematical}
{it:Statistics and Probability} 1: 361-379.

{pstd}
Morris, C. N. 1983. Parametric empirical Bayes inference: Theory and
applications. {it:Journal of the American Statistical Association}
78(381): 47-55.


{title:Author}

{pstd}
Michael V. Elwell  
University of North Texas  
Email: MichaelElwell@my.unt.edu
{p_end}

