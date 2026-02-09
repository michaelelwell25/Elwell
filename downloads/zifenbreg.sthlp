{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "zifenbreg##syntax"}{...}
{viewerjumpto "Description" "zifenbreg##description"}{...}
{viewerjumpto "Options" "zifenbreg##options"}{...}
{viewerjumpto "Remarks" "zifenbreg##remarks"}{...}
{viewerjumpto "Examples" "zifenbreg##examples"}{...}
{viewerjumpto "Stored results" "zifenbreg##results"}{...}
{viewerjumpto "References" "zifenbreg##references"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{bf:zifenbreg} {hline 2}}Zero-inflated fixed-effects negative binomial regression{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:zifenbreg} {depvar} {indepvars} {ifin}
{cmd:,} {opt inf:late(varlist)} [{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opth inf:late(varlist)}}variables for zero-inflation equation; required{p_end}
{synopt:{opt nb1}}use NB1 parameterization; default is NB2{p_end}
{synopt:{opth exp:osure(varname)}}include ln({it:varname}) as offset{p_end}
{synopt:{opth off:set(varname)}}include {it:varname} as offset{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(robust)}}sandwich SEs clustered at panel level{p_end}
{synopt:{opt vce(cluster} {it:clustvar}{opt )}}sandwich SEs clustered by {it:clustvar}{p_end}

{syntab:Reporting}
{synopt:{opt irr}}report incidence-rate ratios for count equation{p_end}
{synopt:{opt vuong}}report Vuong test comparing ZINB vs NB{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt nolog}}suppress iteration log{p_end}

{syntab:Maximization}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is {cmd:tolerance(1e-8)}{p_end}
{synopt:{opt iter:ate(#)}}maximum iterations; default is {cmd:iterate(1000)}{p_end}
{synoptline}
{p 4 6 2}
A panel variable and a time variable must be set using {helpb xtset}.{p_end}
{p 4 6 2}
Factor variables and time-series operators are allowed; see {help fvvarlist}.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:zifenbreg} fits a zero-inflated fixed-effects negative binomial
regression model. The count equation uses a concentrated (profile)
likelihood to absorb unit-level fixed effects, as in {helpb fenbreg}.
The zero-inflation (inflate) equation is a logit model without unit
fixed effects.

{pstd}
The model is:

{pmore}
P(y_it = 0) = pi_it + (1 - pi_it) * NB(0 | mu_it, alpha){p_end}
{pmore}
P(y_it = j) = (1 - pi_it) * NB(j | mu_it, alpha), for j > 0{p_end}

{pmore}
Count equation:    ln(mu_it) = a_i + x_it * beta + offset{p_end}
{pmore}
Inflate equation:  logit(pi_it) = z_it * gamma{p_end}

{pstd}
where a_i are unit fixed effects concentrated out of the likelihood
(as in {cmd:fenbreg}) and gamma are estimated by maximum likelihood
jointly with beta and alpha.

{pstd}
The inflate equation does {bf:not} include unit fixed effects because
logit FE suffers from the incidental parameters problem (O(1/T) bias)
that is not resolved by the concentration approach. This is the
standard approach in applied work for short panels. Users can include
unit-level covariates (e.g., population, ideology scores, region
dummies) in {opt inflate()} to proxy for unobserved heterogeneity.

{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth inflate(varlist)} specifies the variables that determine whether
an observation has a zero from the zero-inflation process. A constant
is automatically included. This option is required.

{phang}
{opt nb1} specifies the NB1 (linear variance) parameterization where
Var(y|x,c) = mu*(1+alpha). The default is NB2 where
Var(y|x,c) = mu*(1+alpha*mu).

{phang}
{opth exposure(varname)} specifies a variable whose log is included as
an offset in the count equation. Cannot be combined with {cmd:offset()}.

{phang}
{opth offset(varname)} specifies a variable to be included directly as
an offset in the count equation. Cannot be combined with {cmd:exposure()}.

{dlgtab:SE/Robust}

{phang}
{opt vce(robust)} computes sandwich standard errors clustered at the
panel level. This is the default.

{phang}
{opt vce(cluster} {it:clustvar}{opt )} clusters standard errors by
{it:clustvar}.

{dlgtab:Reporting}

{phang}
{opt irr} reports exponentiated coefficients (incidence-rate ratios)
for the count equation.

{phang}
{opt vuong} reports a Vuong (1989) test comparing the zero-inflated
model against a standard negative binomial. A positive test statistic
favors the ZINB specification.

{phang}
{opt level(#)} specifies the confidence level for confidence intervals.

{phang}
{opt nolog} suppresses the iteration log.

{dlgtab:Maximization}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for both the
inner Newton-Raphson loop and the outer optimization. Default is 1e-8.

{phang}
{opt iterate(#)} specifies the maximum number of outer iterations.
Default is 1000.

{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Why no fixed effects in the inflate equation?}

{pstd}
The concentrated likelihood approach solves the incidental parameters
problem (IPP) for the count equation's NB/Poisson component. However,
this concentration does not resolve the IPP for the logit inflate
equation. Including unit FE in a logit model with short T produces
O(1/T) bias (Neyman and Scott, 1948; Lancaster, 2000). For a typical
panel with T around 10, this bias can be substantial.

{pstd}
The practical solution (standard in applied work and matching Stata's
{cmd:zinb} approach) is to exclude unit FE from the inflate equation
and instead include observable unit-level covariates to capture
cross-sectional variation in the probability of structural zeros.

{pstd}
{bf:Algorithm}

{pstd}
The estimator uses nested optimization identical in philosophy to
{cmd:fenbreg}:

{phang2}1. {it:Inner loop}: For given (beta, gamma, alpha), solve for
each unit's log-effect a_i via Newton-Raphson. The FOC is the same as
{cmd:fenbreg}'s but with posterior ZI weights on y=0 observations.

{phang2}2. {it:Outer loop}: Maximize the concentrated ZI NB log-likelihood
over (beta, gamma, ln(alpha)) using NR/BFGS with analytical gradients.

{phang2}3. {it:Variance}: Sandwich VCE with panel-clustered scores of
dimension k_count + k_inflate + 1.

{pstd}
Starting values are obtained from a Poisson FE regression (for beta),
a pooled logit of I(y=0) on inflate variables (for gamma), and
lnalpha=0.

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang}{cmd:. webuse ships}{p_end}
{phang}{cmd:. xtset ship year}{p_end}

{pstd}Basic ZINB FE model with zero-inflation depending on covariates{p_end}
{phang}{cmd:. zifenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, inflate(op_75_79 co_75_79)}{p_end}

{pstd}With exposure and IRR{p_end}
{phang}{cmd:. zifenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, inflate(op_75_79) exposure(service) irr}{p_end}

{pstd}NB1 parameterization{p_end}
{phang}{cmd:. zifenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, inflate(op_75_79) nb1}{p_end}

{pstd}With Vuong test{p_end}
{phang}{cmd:. zifenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, inflate(op_75_79) vuong}{p_end}

{pstd}Post-estimation predictions{p_end}
{phang}{cmd:. predict yhat, yhat}{p_end}
{phang}{cmd:. predict mu_hat, mu}{p_end}
{phang}{cmd:. predict pi_hat, pr}{p_end}
{phang}{cmd:. predict xb_hat, xb}{p_end}
{phang}{cmd:. predict zg_hat, zg}{p_end}
{phang}{cmd:. predict fe_hat, fe}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:zifenbreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of groups{p_end}
{synopt:{cmd:e(k_count)}}number of count-equation regressors{p_end}
{synopt:{cmd:e(k_inflate)}}number of inflate-equation parameters (incl. constant){p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood, null model{p_end}
{synopt:{cmd:e(lnalpha)}}log of dispersion parameter{p_end}
{synopt:{cmd:e(alpha)}}dispersion parameter{p_end}
{synopt:{cmd:e(chi2)}}Wald chi-squared (count equation){p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(p)}}model p-value{p_end}
{synopt:{cmd:e(converged)}}1 if converged{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:zifenbreg}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(model)}}NB2 or NB1{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(ivar)}}panel variable{p_end}
{synopt:{cmd:e(tvar)}}time variable{p_end}
{synopt:{cmd:e(vce)}}cluster{p_end}
{synopt:{cmd:e(vcetype)}}Sandwich{p_end}
{synopt:{cmd:e(clustvar)}}clustering variable{p_end}
{synopt:{cmd:e(inflate)}}inflate equation variables{p_end}
{synopt:{cmd:e(predict)}}zifenbreg_p{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (count, inflate, lnalpha){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{marker references}{...}
{title:References}

{phang}
Allison, P.D. and R.P. Waterman. 2002. Fixed-effects negative binomial
regression models. {it:Sociological Methodology} 32: 247-265.

{phang}
Bergé, L. 2018. Efficient estimation of maximum likelihood models with
multiple fixed-effects: the R package FENmlm. {it:CREA Discussion Papers}.

{phang}
Guimaraes, P. 2008. The fixed effects negative binomial model revisited.
{it:Economics Letters} 99: 63-66.

{phang}
Lambert, D. 1992. Zero-inflated Poisson regression, with an application
to defects in manufacturing. {it:Technometrics} 34: 1-14.

{phang}
Neyman, J. and E.L. Scott. 1948. Consistent estimates based on partially
consistent observations. {it:Econometrica} 16: 1-32.

{phang}
Vuong, Q.H. 1989. Likelihood ratio tests for model selection and
non-nested hypotheses. {it:Econometrica} 57: 307-333.

{title:Author}

{pstd}
Michael Elwell{p_end}
{pstd}
University of North Texas{p_end}

{title:Also see}

{psee}
{helpb fenbreg}, {helpb zinb}, {helpb xtnbreg}, {helpb xtpoisson}
{p_end}
