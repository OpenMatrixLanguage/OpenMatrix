/**
* @file DistributionFuncs.h
* @date June 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair's dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair's trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/
#ifndef _Stats_DistributionFuncs_h
#define _Stats_DistributionFuncs_h

#include "StatisticsExports.h"

// forward declarations
class hwMathStatus;
class hwMersenneTwisterState;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<int, hwTComplex<int> > hwMatrixI;
typedef hwTMatrix<int64_t, hwTComplex<int64_t> > hwMatrixI64;

//------------------------------------------------------------------------------
//!
//! \brief Distribution functions
//!
//------------------------------------------------------------------------------

//!
//! Generates seed for random number generator using the system clock
//!
unsigned long GetSeedFromClock();

// Uniform distribution methods
//!
//! Uniform distribution probability density function
//! \param x
//! \param a
//! \param b
//! \param density
//!
STATISTICS_DECLS hwMathStatus UnifPDF(double  x, 
                                      double  a, 
                                      double  b, 
                                      double& density);
//!
//! Uniform distribution cumulative density function
//! \param x
//! \param a
//! \param b
//! \param prob
//!
STATISTICS_DECLS hwMathStatus UnifCDF(double  x, 
                                      double  a, 
                                      double  b, 
                                      double& prob);
//!
//! Uniform distribution inverse cumulative density function
//! \param prob
//! \param a
//! \param b
//! \param x
//!
STATISTICS_DECLS hwMathStatus UnifInvCDF(double  prob, 
                                         double  a, 
                                         double  b, 
                                         double& x);
//!
//! Uniform distribution parameter estimation function
//! \param data Data sample
//! \param aHat
//! \param bHat
//! \param aCI  Optional argument
//! \param bCI  Optional argument
//!
STATISTICS_DECLS hwMathStatus UnifFit(const hwMatrix& data, 
                                      double&         aHat, 
                                      double&         bHat,
                                      hwMatrix*       aCI = nullptr, 
                                      hwMatrix*       bCI = nullptr);
//!
//! Uniform distribution random number generation, generating a single value on (0,1)
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus UnifRnd(hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Uniform distribution random number generation, generating a single value on (A,B)
//! \param A       Input
//! \param B       Input
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus UnifRnd(double                  A, 
                                      double                  B, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Uniform distribution random number generation, generating a matrix of values
//! on (A,B)
//! \param A       Input
//! \param B       Input
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus UnifRnd(double                  A, 
                                      double                  B, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);
//!
//! Uniform distribution random number generation, generating a matrix of values
//! on (A,B)
//! \param A       Input
//! \param B       Input
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus UnifRnd(const hwMatrix&         A, 
                                      const hwMatrix&         B,
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);

// Normal distribution methods

//!
//! Normal distribution probability density function
//! \param x       
//! \param mu       
//! \param sigma
//! \param density
//!
STATISTICS_DECLS hwMathStatus NormPDF(double  x, 
                                      double  mu, 
                                      double  sigma, 
                                      double& density);
//!
//! Normal distribution cumulative density function
//! \param x       
//! \param mu       
//! \param sigma
//! \param prob
//!
STATISTICS_DECLS hwMathStatus NormCDF(double  x, 
                                      double  mu, 
                                      double  sigma, 
                                      double& prob);
//!
//! Normal distribution inverse cumulative density function
//! \param prob       
//! \param mu       
//! \param sigma
//! \param x
//!
STATISTICS_DECLS hwMathStatus NormInvCDF(double  prob, 
                                         double  mu, 
                                         double  sigma, 
                                         double& x);
//!
//! Normal distribution parameter estimation function
//! \param data       
//! \param muHat       
//! \param sigmaHat
//! \param muCI
//! \param sigmaCI
//!
STATISTICS_DECLS hwMathStatus NormFit(const hwMatrix& data, 
                                      double&         muHat,
                                      double&         sigmaHat,
                                      hwMatrix*       muCI    = nullptr, 
                                      hwMatrix*       sigmaCI = nullptr);
//!
//! Normal distribution random number generation, generating single value from 
//! N(mu,sigma^2)
//! \param MTstate       
//! \param seed       
//! \param value
//!
STATISTICS_DECLS hwMathStatus NormRnd(hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Normal distribution random number generation, generating single value from 
//! N(mu,sigma^2) 
//! \param mu
//! \param sigma
//! \param MTstate       
//! \param seed       
//! \param value
//!
STATISTICS_DECLS hwMathStatus NormRnd(double                  mu, 
                                      double                  sigma, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Normal distribution random number generation, generating a matrix of values 
//! from N(mu,sigma^2) 
//! \param mu
//! \param sigma
//! \param MTstate       
//! \param seed       
//! \param value
//!
STATISTICS_DECLS hwMathStatus NormRnd(double                  mu, 
                                      double                  sigma, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);
//!
//! Normal distribution random number generation, generating a matrix of values 
//! from N(mu,sigma^2) 
//! \param mu
//! \param sigma
//! \param MTstate       
//! \param seed       
//! \param matrix
//!
STATISTICS_DECLS hwMathStatus NormRnd(const hwMatrix&         mu, 
                                      const hwMatrix&         sigma,
                                      hwMersenneTwisterState* MTstate, 
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);

// Beta distribution methods

//!
//! Beta distribution probability density function
//! \param x
//! \param a
//! \param b
//! \param density
//!
STATISTICS_DECLS hwMathStatus BetaPDF(double  x, 
                                      double  a, 
                                      double  b, 
                                      double& density);
//!
//! Beta distribution cumulative density function
//! \param x
//! \param a
//! \param b
//! \param prob
//!
STATISTICS_DECLS hwMathStatus BetaCDF(double  x, 
                                      double  a, 
                                      double  b, 
                                      double& prob);
//!
//! Beta distribution inverse cumulative density function
//! \param prob
//! \param a
//! \param b
//! \param x
//!
STATISTICS_DECLS hwMathStatus BetaInvCDF(double  prob, 
                                         double  a, 
                                         double  b, 
                                         double& x);
//!
//! Beta distribution parameter estimation function
//! \param data Data sample
//! \param aHat
//! \param bHat
//! \param aCI
//! \param bCI
//!
STATISTICS_DECLS hwMathStatus BetaFit(const hwMatrix& data, 
                                      double&         aHat, 
                                      double&         bHat,
                                      hwMatrix*       aCI = nullptr, 
                                      hwMatrix*       bCI = nullptr);
//!
//! Beta distribution random number generation, generating a single value from
//! Beta(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus BetaRnd(double                  a, 
                                      double                  b, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Beta distribution random number generation, generating a matrix of values 
//! from Beta(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus BetaRnd(double                  a, 
                                      double                  b,
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);
//!
//! Beta distribution random number generation, generating a matrix of values 
//! from Beta(a,b)
//! \param A
//! \param B
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus BetaRnd(const hwMatrix&         A, 
                                      const hwMatrix&         B,
                                      hwMersenneTwisterState* MTstate, 
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);

// Gamma distribution methods

//!
//! Gamma distribution probability density function
//! \param x
//! \param a
//! \param b
//! \param density
//!
STATISTICS_DECLS hwMathStatus GammaPDF(double  x, 
                                       double  a, 
                                       double  b, 
                                       double& density);
//!
//! Gamma distribution cumulative density function
//! \param x
//! \param a
//! \param b
//! \param prob
//!
STATISTICS_DECLS hwMathStatus GammaCDF(double  x, 
                                       double  a, 
                                       double  b, 
                                       double& prob);
//!
//! Gamma distribution inverse cumulative density function
//! \param prob
//! \param a
//! \param b
//! \param x
//!
STATISTICS_DECLS hwMathStatus GammaInvCDF(double  prob, 
                                          double  a, 
                                          double  b, 
                                          double& x);
//!
//! Gamma distribution parameter estimation function
//! \param data
//! \param aHat
//! \param bHat
//! \param aCI
//! \param bCI
//!
STATISTICS_DECLS hwMathStatus GammaFit(const hwMatrix& data, 
                                       double&         aHat, 
                                       double&         bHat,
                                       hwMatrix*       aCI = nullptr, 
                                       hwMatrix*       bCI = nullptr);
//!
//! Gamma distribution random number generation, generating a single value from
//! Gamma(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus GammaRnd(double                  a, 
                                       double                  b, 
                                       hwMersenneTwisterState* MTstate,
                                       unsigned long*          seed,
                                       double&                 value);
//!
//! Gamma distribution random number generation, generating a matrix of values
//! from Gamma(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus GammaRnd(double                  a, 
                                       double                  b, 
                                       hwMersenneTwisterState* MTstate,
                                       unsigned long*          seed, 
                                       hwMatrix&               matrix);
//!
//! Gamma distribution random number generation, generating a matrix of values
//! from Gamma(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value
//!
STATISTICS_DECLS hwMathStatus GammaRnd(const hwMatrix&         A, 
                                       const hwMatrix&         B, 
                                       hwMersenneTwisterState* MTstate,
                                       unsigned long*          seed, 
                                       hwMatrix&               matrix);

// Exponential distribution methods
//!
//! Exponential distribution probability density function
//! \param x
//! \param mu
//! \param density
//!
STATISTICS_DECLS hwMathStatus ExpPDF(double  x, 
                                     double  mu, 
                                     double& density);
//!
//! Exponential distribution cumulative density function
//! \param x
//! \param mu
//! \param prob
//!
STATISTICS_DECLS hwMathStatus ExpCDF(double  x,
                                     double  mu, 
                                     double& prob);
//!
//! Exponential distribution inverse cumulative density function
//! \param prob
//! \param mu
//! \param x
//!
STATISTICS_DECLS hwMathStatus ExpInvCDF(double  prob, 
                                        double  mu, 
                                        double& x);
//!
//! Exponential distribution parameter estimation function
//! \param data
//! \param muHat
//! \param muCI  Optional argument
//!
STATISTICS_DECLS hwMathStatus ExpFit(const hwMatrix& data, 
                                     double&         muHat, 
                                     hwMatrix*       muCI = nullptr);
//!
//! Generates exponential distribution random number with mean mu
//! \param mu
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//!
STATISTICS_DECLS hwMathStatus ExpRnd(double                  mu, 
                                     hwMersenneTwisterState* MTstate,
                                     unsigned long*          seed, 
                                     double&                 value);
//!
//! Generates exponential distribution random numbers with mean mu
//! \param mu
//! \param MTstate
//! \param seed
//! \param value   Generated random numbers
//!
STATISTICS_DECLS hwMathStatus ExpRnd(double                  mu, 
                                     hwMersenneTwisterState* MTstate,
                                     unsigned long*          seed, 
                                     hwMatrix&               matrix);
//!
//! Generates exponential distribution random numbers with mean mu
//! \param Mu
//! \param MTstate
//! \param seed
//! \param value   Generated random numbers
//!
STATISTICS_DECLS hwMathStatus ExpRnd(const hwMatrix&         Mu, 
                                     hwMersenneTwisterState* MTstate,
                                     unsigned long*          seed, 
                                     hwMatrix&               matrix);

// Chi-squared distribution

//!
//! Chi-squared distribution probability density function
//! \param x
//! \param ndof    Degrees of freedom
//! \param density
//!
STATISTICS_DECLS hwMathStatus Chi2PDF(double  x, 
                                      int     ndof, 
                                      double& density);
//!
//! Chi-squared distribution cumulative density function
//! \param x
//! \param ndof Degrees of freedom
//! \param prob
//!
STATISTICS_DECLS hwMathStatus Chi2CDF(double  x, 
                                      int     ndof, 
                                      double& prob);
//!
//! Chi-squared distribution inverse cumulative density function
//! \param prob
//! \param ndof Degrees of freedom
//! \param x
//!
STATISTICS_DECLS hwMathStatus Chi2InvCDF(double  prob, 
                                         int     ndof, 
                                         double& x);
//!
//! Generates Chi-squared distribution random number with ndof degrees of freedom
//! \param ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param value   Generated value
//!
STATISTICS_DECLS hwMathStatus Chi2Rnd(int                     ndof, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      double&                 value);
//!
//! Generates Chi-squared distribution random numbers with ndof degrees of freedom
//! \param ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param matrix  Generated values
//!
STATISTICS_DECLS hwMathStatus Chi2Rnd(int                     ndof, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);
//!
//! Generates Chi-squared distribution random numbers with ndof degrees of freedom
//! \param ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param matrix  Generated values
//!
STATISTICS_DECLS hwMathStatus Chi2Rnd(const hwMatrixI&        Ndof, 
                                      hwMersenneTwisterState* MTstate,
                                      unsigned long*          seed, 
                                      hwMatrix&               matrix);

// T distribution functions

//!
//! Student T distribution probability density function
//! \param x
//! \param ndof    Degrees of freedom
//! \param density
//!
STATISTICS_DECLS hwMathStatus T_PDF(double  x,
                                    int     ndof, 
                                    double& density);
//!
//! Student T distribution cumulative density function
//! \param x
//! \param ndof Degrees of freedom
//! \param prob
//!
STATISTICS_DECLS hwMathStatus T_CDF(double  x, 
                                    int     ndof, 
                                    double& prob);
//!
//! Student T distribution cumulative density function
//! \param prob
//! \param ndof Degrees of freedom
//! \param x
//!
STATISTICS_DECLS hwMathStatus T_InvCDF(double  prob, 
                                       int     ndof, 
                                       double& x);
//!
//! Generates Student T distribution random number with ndof degrees of freedom
//! \param ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//!
STATISTICS_DECLS hwMathStatus T_Rnd(int                     ndof, 
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    double&                 value);
//!
//! Generates Student T distribution random number with ndof degrees of freedom
//! \param ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//!
STATISTICS_DECLS hwMathStatus T_Rnd(int                     ndof, 
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    hwMatrix&               matrix);
//!
//! Generates Student T distribution random number with ndof degrees of freedom
//! \param Ndof    Degrees of freedom
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//!
STATISTICS_DECLS hwMathStatus T_Rnd(const hwMatrixI&        Ndof, 
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    hwMatrix&               matrix);

// F distribution functions

//!
//! F distribution probability density function
//! \param x
//! \param mdof
//! \param ndof
//! \param density
//!
STATISTICS_DECLS hwMathStatus F_PDF(double  x, 
                                    int     mdof, 
                                    int     ndof, 
                                    double& density);
//!
//! F distribution cumulative density function
//! \param x
//! \param mdof
//! \param ndof
//! \param prob
//!
STATISTICS_DECLS hwMathStatus F_CDF(double  x, 
                                    int     mdof, 
                                    int     ndof, 
                                    double& prob);
//!
//! F distribution inverse cumulative density function
//! \param prob
//! \param mdof
//! \param ndof
//! \param x
//!
STATISTICS_DECLS hwMathStatus F_InvCDF(double  prob, 
                                       int     mdof, 
                                       int     ndof, 
                                       double& x);
//!
//! Generates F distribution random number with mdof and ndof degrees of freedom
//! \param mdof
//! \param ndof
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//!
STATISTICS_DECLS hwMathStatus F_Rnd(int                     mdof, 
                                    int                     ndof, 
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    double&                 value);
//!
//! Generates F distribution random number with mdof and ndof degrees of freedom
//! \param mdof
//! \param ndof
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//!
STATISTICS_DECLS hwMathStatus F_Rnd(int                     mdof, 
                                    int                     ndof, 
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    hwMatrix&               matrix);
//!
//! Generates F distribution random number with Mdof and Ndof degrees of freedom
//! \param Mdof
//! \param Ndof
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//!
STATISTICS_DECLS hwMathStatus F_Rnd(const hwMatrixI&        Mdof, 
                                    const hwMatrixI&        Ndof,
                                    hwMersenneTwisterState* MTstate,
                                    unsigned long*          seed, 
                                    hwMatrix&               matrix);

// Pearson distribution functions

//!
//! Compute four moments needed for Pearson distributions
//! \param data   Input
//! \param moment Four moments needed for Pearson distributions
//! 
STATISTICS_DECLS hwMathStatus FourMoments(const hwMatrix& data, 
                                          hwMatrix&       moment);
//!
//! Gets the Pearson distribution family type and parameters
//! \param moment 
//! \param type   Pearson distribution family type
//! \param param  Pearson distribution parameters
//! 
STATISTICS_DECLS hwMathStatus PearsonInfo(const hwMatrix& moment, 
                                          int&            type, 
                                          hwMatrix&       param);
//!
//! Pearson distribution probability density function
//! \param x 
//! \param moment
//! \param density
//! 
STATISTICS_DECLS hwMathStatus PearsonPDF(double          x, 
                                         const hwMatrix& moment, 
                                         double&         density);
//!
//! Pearson distribution cumulative density function
//! \param x 
//! \param moment
//! \param prob
//! 
STATISTICS_DECLS hwMathStatus PearsonCDF(double          x, 
                                         const hwMatrix& moment, 
                                         double&         prob);
//!
//! Pearson distribution inverse cumulative density function
//! \param prob 
//! \param moment
//! \param x
//! 
STATISTICS_DECLS hwMathStatus PearsonInvCDF(double          prob, 
                                            const hwMatrix& moment, 
                                            double&         x);
//!
//! Pearson distribution median
//! \param moment
//! \param median Median
//! 
STATISTICS_DECLS hwMathStatus PearsonMedian(const hwMatrix& moment, 
                                            double&         median);
//!
//! Pearson distribution mode
//! \param moment
//! \param median Median
//! 
STATISTICS_DECLS hwMathStatus PearsonMode(const hwMatrix& moment, 
                                          double&         mode);
//!
//! Generates Pearson distribution random number
//! \param moment
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//! 
STATISTICS_DECLS hwMathStatus PearsonRnd(const hwMatrix&         moment, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         double&                 value);
//!
//! Generates Pearson distribution random numbers
//! \param moment
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//! 
STATISTICS_DECLS hwMathStatus PearsonRnd(const hwMatrix&         moment, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         hwMatrix&               matrix);

// LogNormal distribution functions

//!
//! Lognormal distribution probability density function
//! \param x
//! \param mu
//! \param sigma
//! \param density
//! 
STATISTICS_DECLS hwMathStatus LogNormPDF(double  x, 
                                         double  mu,  
                                         double  sigma, 
                                         double& density);
//!
//! Lognormal distribution cumulative density function
//! \param x
//! \param mu
//! \param sigma
//! \param prob
//! 
STATISTICS_DECLS hwMathStatus LogNormCDF(double  x, 
                                         double  mu, 
                                         double  sigma, 
                                         double& prob);
//!
//! Lognormal distribution inverse cumulative density function
//! \param x
//! \param mu
//! \param sigma
//! \param prob
//! 
STATISTICS_DECLS hwMathStatus LogNormInvCDF(double  prob, 
                                            double  mu, 
                                            double  sigma, 
                                            double& x);
//!
//! Lognormal distribution parameter estimation function
//! \param data
//! \param muHat
//! \param sigmaHat
//! \param muCI     Optional argument
//! \param sigmaCI  Optional argument
//! 
STATISTICS_DECLS hwMathStatus LogNormFit(const hwMatrix& data, 
                                         double&         muHat, 
                                         double&         sigmaHat,
                                         hwMatrix*       muCI    = nullptr, 
                                         hwMatrix*       sigmaCI = nullptr);
//!
//! Generates Lognormal distribution random number using N(mu,sigma^2)
//! \param mu
//! \param sigma
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//! 
STATISTICS_DECLS hwMathStatus LogNormRnd(double                  mu, 
                                         double                  sigma, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed,
                                         double&                 value);
//!
//! Generates Lognormal distribution random numbers using N(mu,sigma^2)
//! \param mu
//! \param sigma
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//! 
STATISTICS_DECLS hwMathStatus LogNormRnd(double                  mu, 
                                         double                  sigma, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed,
                                         hwMatrix&               matrix);
//!
//! Generates Lognormal distribution random numbers using N(mu,sigma^2)
//! \param mu
//! \param sigma
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//! 
STATISTICS_DECLS hwMathStatus LogNormRnd(const hwMatrix&         mu, 
                                         const hwMatrix&         sigma,
                                         hwMersenneTwisterState* MTstate, 
                                         unsigned long*          seed, 
                                         hwMatrix&               matrix);

// Weibull distribution functions

//!
//! Weibull distribution probability density function
//! \param x
//! \param a
//! \param b
//! \param density
//! 
STATISTICS_DECLS hwMathStatus WeibullPDF(double  x, 
                                         double  a,
                                         double  b, 
                                         double& density);
//!
//! Weibull distribution cumulative density function
//! \param x
//! \param a
//! \param b
//! \param prob
//! 
STATISTICS_DECLS hwMathStatus WeibullCDF(double  x, 
                                         double  a, 
                                         double  b, 
                                         double& prob);
//!
//! Weibull distribution inverse cumulative density function
//! \param prob
//! \param a
//! \param b
//! \param x
//! 
STATISTICS_DECLS hwMathStatus WeibullInvCDF(double  prob, 
                                            double  a, 
                                            double  b, 
                                            double& x);
//!
//! Weibull distribution parameter estimation function
//! \param data
//! \param aHat
//! \param bHat
//! \param aCI
//! \param bCI
//! 
STATISTICS_DECLS hwMathStatus WeibullFit(const hwMatrix& data, 
                                         double&         aHat, 
                                         double&         bHat,
                                         hwMatrix*       aCI = nullptr, 
                                         hwMatrix*       bCI = nullptr);
//!
//! Generates Weibull distribution random number from Weibull(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param value   Generated random number
//! 
STATISTICS_DECLS hwMathStatus WeibullRnd(double                  a, 
                                         double                  b, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         double&                 value);
//!
//! Generates Weibull distribution random numbers from Weibull(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//! 
STATISTICS_DECLS hwMathStatus WeibullRnd(double                  a, 
                                         double                  b, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         hwMatrix&               matrix);
//!
//! Generates Weibull distribution random numbers from Weibull(a,b)
//! \param a
//! \param b
//! \param MTstate
//! \param seed
//! \param matrix  Generated random numbers
//! 
STATISTICS_DECLS hwMathStatus WeibullRnd(const hwMatrix&         A, 
                                         const hwMatrix&         B,
                                         hwMersenneTwisterState* MTstate, 
                                         unsigned long*          seed, 
                                         hwMatrix&               matrix);

// Poisson distribution functions

//!
//! Poisson distribution probability density function
//! \param x
//! \param lambda
//! \param density
//! 
STATISTICS_DECLS hwMathStatus PoissonPDF(int     x, 
                                         double  lambda, 
                                         double& density);
//!
//! Poisson distribution cumulative density function
//! \param x
//! \param lambda
//! \param prob
//! 
STATISTICS_DECLS hwMathStatus PoissonCDF(double  x, 
                                         double  lambda, 
                                         double& prob);
//!
//! Poisson distribution inverse cumulative density function
//! \param prob
//! \param lambda
//! \param x
//! 
STATISTICS_DECLS hwMathStatus PoissonInvCDF(double prob, 
                                            double lambda, 
                                            int&   x);
//!
//! Poisson distribution parameter estimation function
//! \param data
//! \param lambdaHat
//! \param lambdaCI
//! 
STATISTICS_DECLS hwMathStatus PoissonFit(const hwMatrix& data, 
                                         double&         lambdaHat,
                                         hwMatrix*       lambdaCI = nullptr);
//!
//! Generates Poisson distribution random number
//! \param lamda
//! \param MTstate
//! \param seed
//! \param value   Generated random value
//!
STATISTICS_DECLS hwMathStatus PoissonRnd(double                  lambda, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         double&                 value);
//!
//! Generates Poisson distribution random numbers
//! \param lamda
//! \param MTstate
//! \param seed
//! \param matrix  Generated random values
//!
STATISTICS_DECLS hwMathStatus PoissonRnd(double                  lambda, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         hwMatrixI&              matrix);
//!
//! Generates Poisson distribution random numbers
//! \param Lamda
//! \param MTstate
//! \param seed
//! \param matrix  Generated random values
//!
STATISTICS_DECLS hwMathStatus PoissonRnd(const hwMatrix&         lambda, 
                                         hwMersenneTwisterState* MTstate,
                                         unsigned long*          seed, 
                                         hwMatrixI&              matrix);

// Misc functions
//!
//! Generate a random permutation vector on [1:max]
//! \param max
//! \param numPts
//! \param permVec
//!
STATISTICS_DECLS hwMathStatus RandPerm(int                     max,
                                       int                     numPts,
                                       hwMersenneTwisterState* pMTState,
                                       hwMatrixI&              permVec);
//!
//! Generate a random permutation vector on [1:max]
//! \param max
//! \param numPts
//! \param permVec
//!
STATISTICS_DECLS hwMathStatus RandPerm(int64_t                 max,
                                       int                     numPts,
                                       hwMersenneTwisterState* pMTState,
                                       hwMatrixI64&            permVec);
//!
//! D'Agostino-Pearson omnibus normality test
//! \param data
//! \param sigLevel
//! \param reject
//!
STATISTICS_DECLS hwMathStatus NormalityTestDP(const hwMatrix& data, 
                                              double          sigLevel, 
                                              bool&           reject);

#endif // _Stats_DistributionFuncs_h
