{smcl}
{* *! version 1.0.0}{...}
{viewerjumpto "Syntax" "fenbreg##syntax"}{...}
{viewerjumpto "Description" "fenbreg##description"}{...}
{viewerjumpto "Options" "fenbreg##options"}{...}
{viewerjumpto "Remarks" "fenbreg##remarks"}{...}
{viewerjumpto "Examples" "fenbreg##examples"}{...}
{viewerjumpto "Stored results" "fenbreg##results"}{...}
{viewerjumpto "References" "fenbreg##references"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{bf:fenbreg} {hline 2}}Fixed-effects negative binomial regression{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:fenbreg} {depvar} {indepvars} {ifin}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt nb1}}use NB1 parameterization; default is NB2{p_end}
{synopt:{opth exp:osure(varname)}}include ln({it:varname}) as offset{p_end}
{synopt:{opth off:set(varname)}}include {it:varname} as offset{p_end}

{syntab:Bias correction}
{synopt:{opt biascorr(method)}}bias correction method; {it:method} is {opt none} (default) or {opt jackknife}{p_end}

{syntab:SE/Robust}
{synopt:{opt vce(robust)}}sandwich SEs clustered at panel level{p_end}
{synopt:{opt vce(cluster} {it:clustvar}{opt )}}sandwich SEs clustered by {it:clustvar}{p_end}

{syntab:Reporting}
{synopt:{opt irr}}report incidence-rate ratios{p_end}
{synopt:{opt comp:are}}compare with Poisson FE and report Hausman test{p_end}
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
{cmd:fenbreg} fits a fixed-effects negative binomial regression model using
a concentrated (profile) likelihood approach that genuinely absorbs
unit-level heterogeneity from the conditional mean.

{pstd}
Unlike {cmd:xtnbreg, fe} (Hausman, Hall, and Griliches 1984), which only
conditions on the dispersion parameter and does not control for
unit-specific unobserved heterogeneity correlated with regressors,
{cmd:fenbreg} estimates unit-specific effects {it:c_i} that enter the
conditional mean:

{pmore}
E[y_it | x_it, c_i] = c_i * exp(x_it * beta)

{pstd}
The model is estimated by maximizing the concentrated log-likelihood
over (beta, alpha) using BFGS, with unit effects solved via an inner
Newton-Raphson loop at each evaluation. Standard errors are computed
using a panel-clustered sandwich estimator.

{pstd}
Both NB2 (default, quadratic variance) and NB1 (linear variance)
parameterizations are supported.

{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt nb1} specifies the NB1 (linear variance) parameterization where
Var(y|x,c) = mu*(1+alpha). The default is NB2 where
Var(y|x,c) = mu*(1+alpha*mu).

{phang}
{opth exposure(varname)} specifies a variable whose log is included as
an offset. Cannot be combined with {cmd:offset()}.

{phang}
{opth offset(varname)} specifies a variable to be included directly as
an offset. Cannot be combined with {cmd:exposure()}.

{dlgtab:Bias correction}

{phang}
{opt biascorr(method)} specifies the bias-correction method.

{pmore}
{opt none} (the default) applies no bias correction.

{pmore}
{opt jackknife} applies the split-panel jackknife of Dhaene and
Jochmans (2015). Each unit's time series is split into two halves, the
model is re-estimated on each half-panel, and the corrected estimator
is beta_SPJ = 2*beta_full - 0.5*(beta_h1 + beta_h2). Requires T_i >= 4
for all units.

{dlgtab:SE/Robust}

{phang}
{opt vce(robust)} computes sandwich standard errors clustered at the
panel level. This is the default.

{phang}
{opt vce(cluster} {it:clustvar}{opt )} clusters standard errors by
{it:clustvar}.

{dlgtab:Reporting}

{phang}
{opt irr} reports exponentiated coefficients (incidence-rate ratios).

{phang}
{opt compare} estimates a Poisson FE model ({cmd:xtpoisson, fe}) on the
same sample and displays a side-by-side comparison with a Hausman-type
test for the equality of coefficients.

{phang}
{opt level(#)} specifies the confidence level for confidence intervals.

{phang}
{opt nolog} suppresses the iteration log.

{dlgtab:Maximization}

{phang}
{opt tolerance(#)} specifies the convergence tolerance for both the
inner Newton-Raphson loop and the outer BFGS optimization. Default is
1e-8.

{phang}
{opt iterate(#)} specifies the maximum number of outer iterations.
Default is 1000.

{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Why not xtnbreg, fe?}

{pstd}
The Hausman, Hall, and Griliches (1984) "fixed-effects" negative
binomial model implemented in {cmd:xtnbreg, fe} is not a true
fixed-effects estimator. As shown by Allison and Waterman (2002) and
Guimaraes (2008), the HHG model conditions on a nuisance parameter
related to the dispersion rather than the conditional mean, and
estimates are not invariant to the addition of unit dummies.

{pstd}
{cmd:fenbreg} implements a genuine concentrated-likelihood approach
where unit effects c_i enter mu_it = c_i * exp(x_it * beta) and are
profiled out of the likelihood at each evaluation step. This is
analogous to how Poisson FE eliminates unit effects but extended to the
negative binomial.

{pstd}
{bf:Algorithm details}

{pstd}
The estimator uses a nested optimization:

{phang2}1. {it:Inner loop}: For given (beta, alpha), solve for each
unit's log-effect a_i = log(c_i) via Newton-Raphson. This is a scalar
concave problem per unit and typically converges in 2-5 iterations.

{phang2}2. {it:Outer loop}: Maximize the concentrated log-likelihood
over (beta, ln(alpha)) using BFGS with analytically computed gradients.
The gradient of the concentrated likelihood equals the partial
derivative of the full likelihood evaluated at the optimal c_i
(envelope theorem).

{phang2}3. {it:Variance}: Sandwich VCE with panel-clustered scores.

{marker examples}{...}
{title:Examples}

{phang}{cmd:. webuse ships}{p_end}
{phang}{cmd:. xtset ship year}{p_end}

{pstd}Basic NB2 fixed-effects model{p_end}
{phang}{cmd:. fenbreg accident op_75_79 co_65_69 co_70_74 co_75_79}{p_end}

{pstd}With exposure variable{p_end}
{phang}{cmd:. fenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, exposure(service)}{p_end}

{pstd}NB1 parameterization with IRR{p_end}
{phang}{cmd:. fenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, nb1 irr}{p_end}

{pstd}Split-panel jackknife bias correction{p_end}
{phang}{cmd:. fenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, biascorr(jackknife)}{p_end}

{pstd}Compare with Poisson FE{p_end}
{phang}{cmd:. fenbreg accident op_75_79 co_65_69 co_70_74 co_75_79, compare}{p_end}

{pstd}Post-estimation prediction{p_end}
{phang}{cmd:. predict mu_hat, mu}{p_end}
{phang}{cmd:. predict xb_hat, xb}{p_end}
{phang}{cmd:. predict fe_hat, fe}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:fenbreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of groups{p_end}
{synopt:{cmd:e(k)}}number of regressors{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood, constant-only model{p_end}
{synopt:{cmd:e(lnalpha)}}log of dispersion parameter{p_end}
{synopt:{cmd:e(alpha)}}dispersion parameter{p_end}
{synopt:{cmd:e(chi2)}}Wald chi-squared{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(p)}}model p-value{p_end}
{synopt:{cmd:e(converged)}}1 if converged{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:fenbreg}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(model)}}NB2 or NB1{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(ivar)}}panel variable{p_end}
{synopt:{cmd:e(tvar)}}time variable{p_end}
{synopt:{cmd:e(vce)}}cluster{p_end}
{synopt:{cmd:e(vcetype)}}Sandwich{p_end}
{synopt:{cmd:e(clustvar)}}clustering variable{p_end}
{synopt:{cmd:e(biascorr)}}bias correction method{p_end}
{synopt:{cmd:e(predict)}}fenbreg_p{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{marker references}{...}
{title:References}

{phang}
Allison, P.D. and R.P. Waterman. 2002. Fixed-effects negative binomial
regression models. {it:Sociological Methodology} 32: 247-265.

{phang}
Dhaene, G. and K. Jochmans. 2015. Split-panel jackknife estimation of
fixed-effect models. {it:Review of Economic Studies} 82: 991-1030.

{phang}
Guimaraes, P. 2008. The fixed effects negative binomial model revisited.
{it:Economics Letters} 99: 63-66.

{phang}
Hausman, J.A., B.H. Hall, and Z. Griliches. 1984. Econometric models
for count data with an application to the patents-R&D relationship.
{it:Econometrica} 52: 909-938.

{phang}
Wooldridge, J.M. 2010.
{it:Econometric Analysis of Cross Section and Panel Data}. 2nd ed.
Cambridge, MA: MIT Press.

{title:Author}

{pstd}
Michael Elwell{p_end}
{pstd}
University of North Texas{p_end}

{title:Also see}

{psee}
{helpb xtnbreg}, {helpb xtpoisson}, {helpb nbreg}
{p_end}
