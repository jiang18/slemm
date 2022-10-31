#ifndef GAMMA_INC_H_
#define GAMMA_INC_H_


/*  == expint: Exponential Integral and Incomplete Gamma Function ==
 *
 *  Declarations for the package and various constant and macro
 *  definitions.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

/* Exported function */
long double kf_gammap(long double s, long double z);
long double kf_gammaq(long double s, long double z);
long double gamma_inc(long double a, long double x);

/* Constants (taken from gsl_machine.h in GSL sources) */
#define EULER_CNST         0.57721566490153286060651209008L
#define LDBL_TINY 1e-1000L

/* Macros */
#define E1_IS_ODD(n)  ((n) & 1)	/* taken from GSL */

#endif
