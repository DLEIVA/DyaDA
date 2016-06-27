\name{DyaDA-package}
\alias{DyaDA}
\docType{package}
\title{
Dyadic Data Analysis
}
\description{
DyaDA is a package that carries out analyses for dyadic data. It includes several statistical methods for quantifying and testing linearity and steepness of social dominance hierarchies, rank orders in linear hierarchies, matrix correlation methods, social reciprocity based on dyadic discrepancies and dyadic nonindependence in several dyadic designs among others.
}
\details{
\tabular{ll}{
Package: \tab DyaDA\cr
Version: \tab 0.1.1\cr
Date: \tab 2016-23-06\cr
Depends: \tab >= 3.1.0\cr
License: \tab GPL version 2 or newer\cr
}

Index:
\preformatted{
block.SRM.between.groups.tTest				    Between groups t test for SRM block designs
block.SRM.dyadic.reciprocity				      Dyadic reciprocity in SRM block designs
block.SRM.dyadic.reciprocity.tTest		    Dyadic covariance t test in SRM block designs
block.SRM.effects					                SRM effects in block designs
block.SRM.generalized.reciprocity			    Generalized reciprocity in SRM block designs
block.SRM.jackknife					              Jackknife test in SRM block designs
block.SRM.relative.variances				      SRM relative variances in block designs
block.SRM.reliability					            Reliability of SRM estimators in block designs
block.SRM.variances					              SRM absolute variances in block designs
categorical.distinguishable.dyad			    Nonindependence for distinguishable dyads and nominal outcomes
categorical.indistinguishable.dyad			  Nonindependence for indistinguishable dyads and nominal outcomes
dietzRtest						                    Statistical test for Dietz's R statistic
distinguishable.dyad					            Nonindependence for distinguishable dyads and interval outcomes
freq2abil					                        Transformation of a matrix of social interactions
getcirc							                      Number of circular triads
getDs                   				          David's scores based on Dij -DS-
getdc					                            Directional consistency index - DC -
getdelta			                            Delta index - delta -
getepsilon						                    Generalized reciprocity index - epsilon -
getexpeccirc						                  Expected Number of Circular Triads
getexpech						                      Mathematical expectancy for Landau's Linearity Index
getexpecR						                      Mathematical expectancy for Dietz's R statistic
getexpecZ						                      Mathematical expectancy for Mantel's Z statistic
getimplandau						                  Improved Landau's Linearity Index -h'-
getkappa						                      Dyadic reciprocity index - kappa -
getkendall						                    Kendall's Linearity Index
getKr							                        Rowwise Kendall's K statistic
getlandau						                      Landau's Linearity Index -h-
getlambda						                      Matrix of unweighted contributions to symmetry - lambda -
getlambdaj						                    Array of individuals' contribution to symmetry - lambdaj -
getmaxcirc						                    Maximum number of circular triads
getnu						                          Matrix of unweighted contributions to skew-symmetry - nu -
getnuj						                        Array of individuals' contribution to skew-symmetry - nuj -
getomega						                      Matrix of dyadic balanced reciprocity - omega -
getonewayrel						                  Dyads with one way relationships
getpartialZr						                  Partial rowwise Mantel's Z statistic
getpartialRr						                  Partial rowwise Dietz's R statistic
getpartialKr						                  Partial rowwise Kendall's K statistic
getpearson						                    Pearson's Correlation between two matrices
getphi                 				            Skew-symmetry index - Phi -
getphii                 				          Individuals' contributions to asymmetry as actors
getPHIij                				          Dyadic contributions to asymmetry matrix
getphij                 				          Individuals' contributions to asymmetry as partners
getPHIr                 				          Overall asymmetry index
getPHIrexpec            				          Mathematical expectancy of the PHIr statistic
getphirmat              				          Dyadic directional asymmetry matrix
getPHIrSE               				          Standard error of the PHIr statistic
getpsi                 				            Symmetry index - psi -
getpvkendall						                  Statistical Significance for Kendall's Linearity Index
getR							                        Dietz's R statistic
getratiolambda							              Matrix of dyadic contributions to symmetry - ratiolambda -
getrationu							                  Matrix of dyadic contributions to skew-symmetry - rationu -
getRr							                        Rowwise Dietz's R statistic -Rr-
getsignificant							              Dyads with significant relationships
getSER							                      Standard Error for Dietz's R statistic
getSEZ							                      Standard Error for Mantel's Z statistic
getspearman						                    Spearman's Correlation between two matrices
gett						                          Student's t statistic for two matrices
gettaukendall						                  Kendall's statistics proposed by Dietz
gettied							                      Dyads with tied relationships
gettwowayrel						                  Dyads with two way relationships
getunknown						                    Dyads with unknown relationships
getvarh							                      Variance for Landau's Linearity index
getwl							                        Win-loss measures at individual level
getZ							                        Mantel's Z statistic
getZr							                        Rowwise Mantel's Z statistic -Zr-
ISI.method						                    I&SI method in R
indistinguishable.dyad					          Nonindependence for distinguishable dyads and interval outcomes
linear.hierarchy.test					            Testing linearity in dominance hierarchies
mantelZtest						                    Statistical test for Mantel's Z statistic
pairwise.correlation					            Pairwise correlation for measuring nonindependence
reciptest1               				          Statistical significance for social reciprocity at different levels of analysis
reciptest2               				          Statistical significance for social reciprocity at different levels of analysis
rowwiseKrtest						                  Statistical significance for Kr statistic
rowwiseRrtest						                  Statistical significance for Rr statistic
rowwiseZrtest						                  Statistical significance for Zr statistic
SRM.between.groups.tTest				          Between groups t test for round robin designs
SRM.dyadic.covariance.Ftest				        Dyadic covariance F test in SRM round robin designs
SRM.dyadic.reciprocity					          Dyadic reciprocity in SRM round robin designs
SRM.effects						                    SRM effects in round robin designs
SRM.generalized.reciprocity				        Generalized reciprocity in SRM round robin designs
SRM.jackknife						                  Jackknife test in SRM round robin designs
SRM.relative.variances					          SRM relative variances in round robin designs
SRM.reliability						                Reliability of SRM estimators in round robin designs
SRM.variances						                  SRM absolute variances in round robin designs
SRM.within.groups.tTest					          Within groups t test for SRM round robin designs
}
}
\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu>, Han de Vries <J.deVries1@uu.nl>, & David A. Kenny <david.kenny@uconn.edu>.

Maintainer: David Leiva <dleivaur@ub.edu>
}
\references{

de Vries, H. (1993). The rowwise correlation between two proximity matrices and the partial rowwise correlation. \emph{Psychometrika}, \emph{58}, 53-69.
de Vries, H. (1995). An improved test of linearity in dominance hierarchies containing unknown or tied relationships. \emph{Animal Behaviour}, \emph{50}, 1375-1389.
de Vries, H. (1998). Finding a dominance order most consistent with a linear hierarchy: a new procedure and review. \emph{Animal Behaviour}, \emph{55}, 827-843.
de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.
Solanas, A., Leiva, D., Sierra, V., & Salafranca, Ll. (2009). Measuring and making decisions for social reciprocity. \emph{Behavior Research Methods}. \emph{41}, 742-754.
}
\keyword{package}

