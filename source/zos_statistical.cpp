#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <iostream>

#include "zos_statistical.h"

#define MAX(x, y) ((x) < (y) ? (y) : (x))
#define MIN(x, y) ((x) > (y) ? (y) : (x))
#define isNonPositive(x) ((x) <= 0.e0 ? 1 : 0)
#define isPositive(x) ((x) > 0.e0 ? 1 : 0)
#define isNegative(x) ((x) < 0.e0 ? 1 : 0)
#define isGreaterThanOne(x) ((x) > 1.e0 ? 1 : 0)
#define isZero(x) ((x) == 0.e0 ? 1 : 0)
#define isOne(x) ((x) == 1.e0 ? 1 : 0)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                         G L O B A L  C O N S T A N T S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

constexpr double ALPHA = 0.01; /* SIGNIFICANCE LEVEL */

static const double rel_error = 1E-12;

constexpr double MACHEP = 1.11022302462515654042E-16;        // 2**-53
constexpr double MAXLOG = 7.09782712893383996732224E2;       // log(MAXNUM)
constexpr double MAXNUM = 1.7976931348623158E308;            // 2**1024*(1-MACHEP)
constexpr double PI = 3.14159265358979323846;                // pi, duh!

static double big = 4.503599627370496e15;
static double biginv = 2.22044604925031308085e-16;

int sgngam = 0;
extern unsigned char *epsilon;

double cephes_igamc (double a, double x)
{
   double ans = 0;
   double ax = 0;
   double c = 0;
   double yc = 0;
   double r = 0;
   double t = 0;
   double y = 0;
   double z = 0;
   double pk = 0;
   double pkm1 = 0;
   double pkm2 = 0;
   double qk = 0;
   double qkm1 = 0;
   double qkm2 = 0;

   if ((x <= 0) || (a <= 0))
   {
      return (1.0);
   }

   if ((x < 1.0) || (x < a))
   {
      return (1.e0 - cephes_igam (a, x));
   }

   ax = a * log (x) - x - cephes_lgam (a);

   if (ax < -MAXLOG)
   {
      std::printf ("igamc: UNDERFLOW\n");
      return 0.0;
   }
   ax = exp (ax);

   /* continued fraction */
   y = 1.0 - a;
   z = x + y + 1.0;
   c = 0.0;
   pkm2 = 1.0;
   qkm2 = x;
   pkm1 = x + 1.0;
   qkm1 = z * x;
   ans = pkm1 / qkm1;

   do
   {
      c += 1.0;
      y += 1.0;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z - pkm2 * yc;
      qk = qkm1 * z - qkm2 * yc;
      if (qk != 0)
      {
         r = pk / qk;
         t = std::fabs ((ans - r) / r);
         ans = r;
      }
      else
      {
         t = 1.0;
      }
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (std::fabs (pk) > big)
      {
         pkm2 *= biginv;
         pkm1 *= biginv;
         qkm2 *= biginv;
         qkm1 *= biginv;
      }
   } while (t > MACHEP);

   return ans * ax;
}

double cephes_igam (double a, double x)
{
   double ans, ax, c, r;

   if ((x <= 0) || (a <= 0))
   {
      return 0.0;
   }

   if ((x > 1.0) && (x > a))
   {
      return 1.e0 - cephes_igamc (a, x);
   }

   /* Compute  x**a * exp(-x) / gamma(a)  */
   ax = a * log (x) - x - cephes_lgam (a);
   if (ax < -MAXLOG)
   {
      std::printf ("igam: UNDERFLOW\n");
      return 0.0;
   }
   ax = exp (ax);

   /* power series */
   r = a;
   c = 1.0;
   ans = 1.0;

   do
   {
      r += 1.0;
      c *= x / r;
      ans += c;
   } while (c / ans > MACHEP);

   return ans * ax / a;
}


/* Logarithm of gamma function */
double cephes_lgam (double x)
{
   return lgamma (x);
}

double cephes_polevl (double x, double *coef, int N)
{
   double ans = 0.0;
   int i = 0;
   double *p = nullptr;

   p = coef;
   ans = *p++;
   i = N;

   do
   {
      ans = ans * x + *p++;
   } while (--i);

   return ans;
}

double cephes_p1evl (double x, double *coef, int N)
{
   double ans;
   double *p;
   int i;

   p = coef;
   ans = x + *p++;
   i = N - 1;

   do
   {
      ans = ans * x + *p++;
   } while (--i);

   return ans;
}

double cephes_erf (double x)
{
   constexpr double two_sqrtpi = 1.128379167095512574;
   double sum = x, term = x, xsqr = x * x;
   int j = 1;

   if (std::fabs (x) > 2.2)
   {
      return 1.0 - cephes_erfc (x);
   }

   do
   {
      term *= xsqr / j;
      sum -= term / (2 * j + 1);
      j++;
      term *= xsqr / j;
      sum += term / (2 * j + 1);
      j++;
   } while (std::fabs (term) / sum > rel_error);

   return two_sqrtpi * sum;
}

double cephes_erfc (double x)
{
   constexpr double one_sqrtpi = 0.564189583547756287;
   double a = 1, b = x, c = x, d = x * x + 0.5;
   double q1, q2 = b / d, n = 1.0, t;

   if (std::fabs (x) < 2.2)
   {
      return 1.0 - cephes_erf (x);
   }
   if (x < 0)
   {
      return 2.0 - cephes_erfc (-x);
   }

   do
   {
      t = a * n + b * x;
      a = b;
      b = t;
      t = c * n + d * x;
      c = d;
      d = t;
      n += 0.5;
      q1 = q2;
      q2 = b / d;
   } while (std::fabs (q1 - q2) / q2 > rel_error);

   return one_sqrtpi * exp (-x * x) * q2;
}


double cephes_normal (double x)
{
   double arg = 0;
   double result = 0;
   constexpr double sqrt2 = 1.414213562373095048801688724209698078569672;

   if (x > 0)
   {
      arg = x / sqrt2;
      result = 0.5 * (1 + erf (arg));
   }
   else
   {
      arg = -x / sqrt2;
      result = 0.5 * (1 - erf (arg));
   }

   return (result);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R U N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Runs (int n)
{
   int S, k;
   double pi, V, erfc_arg, p_value;

   S = 0;
   for (k = 0; k < n; k++)
   {
      if (epsilon[k])
      {
         S++;
      }
   }
   pi = (double) S / (double) n;

   if (std::fabs (pi - 0.5) > (2.0 / std::sqrt (n)))
   {
      std::printf ("    RUNS TEST\n");
      std::printf ("  ------------------------------------------\n");
      std::printf ("  PI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
      p_value = 0.0;
   }
   else
   {

      V = 1;
      for (k = 1; k < n; k++)
      {
         if (epsilon[k] != epsilon[k - 1])
         {
            V++;
         }
      }

      erfc_arg = std::fabs (V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * std::sqrt (2 * n));
      p_value = erfc (erfc_arg);

      std::printf ("    Runs Test\n");
      std::printf ("  ------------------------------------------\n");
      std::printf ("  Computational Information:\n");
      std::printf ("  ------------------------------------------\n");
      std::printf ("  (a) Pi                        = %f\n", pi);
      std::printf ("  (b) V_n_obs (Total # of runs) = %d\n", (int) V);
      std::printf ("  (c) V_n_obs - 2 n pi (1-pi)\n");
      std::printf ("      -----------------------   = %f\n", erfc_arg);
      std::printf ("        2 std::sqrt(2n) pi (1-pi)\n");
      std::printf ("  ------------------------------------------\n");
      if (isNegative (p_value) || isGreaterThanOne (p_value))
      {
         std::printf ("WARNING:  P_VALUE Is Out Of Range.\n");
      }

      std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
   }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     R A N D O M  E X C U R S I O N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandomExcursions (int n)
{
   int b, i, j, k, J, x;
   int cycleStart, cycleStop, *cycle = nullptr, *S_k = nullptr;
   int stateX[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
   int counter[8] = {0, 0, 0, 0, 0, 0, 0, 0};
   double p_value, sum, constraint, nu[6][8];
   double pi[5][6] = {{0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000},
                      {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
                      {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
                      {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
                      {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051}};

   if (((S_k = (int *) std::calloc (n, sizeof (int))) == nullptr) || ((cycle = (int *) std::calloc (MAX (1000, n / 100), sizeof (int))) == nullptr))
   {
      std::printf ("Random Excursions Test:  Insufficient Work Space Allocated.\n");
      if (S_k != nullptr)
      {
         std::free (S_k);
      }
      if (cycle != nullptr)
      {
         std::free (cycle);
      }
      return;
   }

   J = 0; /* DETERMINE CYCLES */
   S_k[0] = 2 * (int) epsilon[0] - 1;
   for (i = 1; i < n; i++)
   {
      S_k[i] = S_k[i - 1] + 2 * epsilon[i] - 1;
      if (S_k[i] == 0)
      {
         J++;
         if (J > MAX (1000, n / 100))
         {
            std::printf ("Error In Function randomExcursions:  Exceeding The Max Number Of Cycles Expected\n.");
            std::free (S_k);
            std::free (cycle);
            return;
         }
         cycle[J] = i;
      }
   }
   if (S_k[n - 1] != 0)
   {
      J++;
   }
   cycle[J] = n;

   std::printf ("     Random Excursions Test\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  (a) Number Of Cycles (J) = %04d\n", J);
   std::printf ("  (b) Sequence Length (n)  = %d\n", n);

   constraint = MAX (0.005 * std::pow (n, 0.5), 500);
   if (J < constraint)
   {
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  Warning:  Test Not Applicable.  There Are An\n");
      std::printf ("     Insufficient Number Of Cycles.\n");
      std::printf ("  ---------------------------------------------\n");
      for (i = 0; i < 8; i++)
      {
         std::printf ("%f\n", 0.0);
      }
   }
   else
   {
      std::printf ("  (c) Rejection Constraint = %f\n", constraint);
      std::printf ("  -------------------------------------------\n");

      cycleStart = 0;
      cycleStop = cycle[1];
      for (k = 0; k < 6; k++)
      {
         for (i = 0; i < 8; i++)
         {
            nu[k][i] = 0.;
         }
      }
      for (j = 1; j <= J; j++)
      { /* FOR EACH CYCLE */
         for (i = 0; i < 8; i++)
         {
            counter[i] = 0;
         }
         for (i = cycleStart; i < cycleStop; i++)
         {
            if ((S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1))
            {
               if (S_k[i] < 0)
               {
                  b = 4;
               }
               else
               {
                  b = 3;
               }
               counter[S_k[i] + b]++;
            }
         }
         cycleStart = cycle[j] + 1;
         if (j < J)
         {
            cycleStop = cycle[j + 1];
         }

         for (i = 0; i < 8; i++)
         {
            if ((counter[i] >= 0) && (counter[i] <= 4))
            {
               nu[counter[i]][i]++;
            }
            else if (counter[i] >= 5)
            {
               nu[5][i]++;
            }
         }
      }

      for (i = 0; i < 8; i++)
      {
         x = stateX[i];
         sum = 0.;
         for (k = 0; k < 6; k++)
         {
            sum += std::pow (nu[k][i] - J * pi[(int) std::fabs (x)][k], 2) / (J * pi[(int) std::fabs (x)][k]);
         }
         p_value = cephes_igamc (2.5, sum / 2.0);

         if (isNegative (p_value) || isGreaterThanOne (p_value))
         {
            std::printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
         }

         std::printf ("%s  x = %2d chi^2 = %9.6f p_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
         std::printf ("%f\n", p_value);
      }
   }
   std::printf ("\n");

   std::free (S_k);
   std::free (cycle);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                      L O N G E S T  R U N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


void LongestRunOfOnes (int n)
{
   double pval, chi2, pi[7];
   int run, v_n_obs, N, i, j, K, M, V[7];
   unsigned int nu[7] = {0, 0, 0, 0, 0, 0, 0};

   if (n < 128)
   {
      std::printf ("     Longest Runs Of Ones Test\n");
      std::printf ("  ---------------------------------------------\n");
      std::printf ("     n=%d is too short\n", n);
      return;
   }
   if (n < 6272)
   {
      K = 3;
      M = 8;
      V[0] = 1;
      V[1] = 2;
      V[2] = 3;
      V[3] = 4;
      pi[0] = 0.21484375;
      pi[1] = 0.3671875;
      pi[2] = 0.23046875;
      pi[3] = 0.1875;
   }
   else if (n < 750000)
   {
      K = 5;
      M = 128;
      V[0] = 4;
      V[1] = 5;
      V[2] = 6;
      V[3] = 7;
      V[4] = 8;
      V[5] = 9;
      pi[0] = 0.1174035788;
      pi[1] = 0.242955959;
      pi[2] = 0.249363483;
      pi[3] = 0.17517706;
      pi[4] = 0.102701071;
      pi[5] = 0.112398847;
   }
   else
   {
      K = 6;
      M = 10000;
      V[0] = 10;
      V[1] = 11;
      V[2] = 12;
      V[3] = 13;
      V[4] = 14;
      V[5] = 15;
      V[6] = 16;
      pi[0] = 0.0882;
      pi[1] = 0.2092;
      pi[2] = 0.2483;
      pi[3] = 0.1933;
      pi[4] = 0.1208;
      pi[5] = 0.0675;
      pi[6] = 0.0727;
   }

   N = n / M;
   for (i = 0; i < N; i++)
   {
      v_n_obs = 0;
      run = 0;
      for (j = 0; j < M; j++)
      {
         if (epsilon[i * M + j] == 1)
         {
            run++;
            if (run > v_n_obs)
            {
               v_n_obs = run;
            }
         }
         else
         {
            run = 0;
         }
      }
      if (v_n_obs < V[0])
      {
         nu[0]++;
      }
      for (j = 0; j <= K; j++)
      {
         if (v_n_obs == V[j])
         {
            nu[j]++;
         }
      }
      if (v_n_obs > V[K])
      {
         nu[K]++;
      }
   }

   chi2 = 0.0;
   for (i = 0; i <= K; i++)
   {
      chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);
   }

   pval = cephes_igamc ((double) (K / 2.0), chi2 / 2.0);

   std::printf ("     Longest Runs Of Ones Test\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  (a) N (# of substrings)  = %d\n", N);
   std::printf ("  (b) M (Substring Length) = %d\n", M);
   std::printf ("  (c) Chi^2                = %f\n", chi2);
   std::printf ("  ---------------------------------------------\n");
   std::printf ("        F R E Q U E N C Y\n");
   std::printf ("  ---------------------------------------------\n");

   if (K == 3)
   {
      std::printf ("    <=1     2     3    >=4   P-value  Assignment");
      std::printf ("\n   %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
   }
   else if (K == 5)
   {
      std::printf ("  <=4  5  6  7  8  >=9 P-value  Assignment");
      std::printf ("\n   %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5]);
   }
   else
   {
      std::printf ("  <=10  11  12  13  14  15 >=16 P-value  Assignment");
      std::printf ("\n   %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], nu[6]);
   }
   if (isNegative (pval) || isGreaterThanOne (pval))
   {
      std::printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
   }

   std::printf ("%s  p_value = %f\n\n", pval < ALPHA ? "FAILURE" : "SUCCESS", pval);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void __ogg_fdrffti (int n, double *wsave, int *ifac);
void __ogg_fdrfftf (int n, double *X, double *wsave, int *ifac);

void drfti1 (int n, double *wa, int *ifac)
{
   static int ntryh[4] = {4, 2, 3, 5};
   constexpr double two_pi = 6.28318530717958647692528676655900577;
   double arg, argh, argld, fi;
   int ntry = 0, i, j = -1;
   int k1, l1, l2, ib;
   int ld, ii, ip, is, nq, nr;
   int ido, ipm, nfm1;
   int nl = n;
   int nf = 0;

L101:
   j++;
   if (j < 4)
   {
      ntry = ntryh[j];
   }
   else
   {
      ntry += 2;
   }

L104:
   nq = nl / ntry;
   nr = nl - ntry * nq;
   if (nr != 0)
   {
      goto L101;
   }

   nf++;
   ifac[nf + 1] = ntry;
   nl = nq;
   if (ntry != 2)
   {
      goto L107;
   }
   if (nf == 1)
   {
      goto L107;
   }

   for (i = 1; i < nf; i++)
   {
      ib = nf - i + 1;
      ifac[ib + 1] = ifac[ib];
   }
   ifac[2] = 2;

L107:
   if (nl != 1)
   {
      goto L104;
   }
   ifac[0] = n;
   ifac[1] = nf;
   argh = two_pi / n;
   is = 0;
   nfm1 = nf - 1;
   l1 = 1;

   if (nfm1 == 0)
   {
      return;
   }

   for (k1 = 0; k1 < nfm1; k1++)
   {
      ip = ifac[k1 + 2];
      ld = 0;
      l2 = l1 * ip;
      ido = n / l2;
      ipm = ip - 1;

      for (j = 0; j < ipm; j++)
      {
         ld += l1;
         i = is;
         argld = (double) ld * argh;
         fi = 0.0;
         for (ii = 2; ii < ido; ii += 2)
         {
            fi += 1.0;
            arg = fi * argld;
            wa[i++] = cos (arg);
            wa[i++] = sin (arg);
         }
         is += ido;
      }
      l1 = l2;
   }
}

void __ogg_fdrffti (int n, double *wsave, int *ifac)
{
   if (n == 1)
   {
      return;
   }
   drfti1 (n, wsave + n, ifac);
}


void dradf2 (int ido, int l1, double *cc, double *ch, double *wa1)
{
   int i, k;
   double ti2, tr2;
   int t0, t1, t2, t3, t4, t5, t6;

   t1 = 0;
   t0 = (t2 = l1 * ido);
   t3 = ido << 1;
   for (k = 0; k < l1; k++)
   {
      ch[t1 << 1] = cc[t1] + cc[t2];
      ch[(t1 << 1) + t3 - 1] = cc[t1] - cc[t2];
      t1 += ido;
      t2 += ido;
   }

   if (ido < 2)
   {
      return;
   }
   if (ido == 2)
   {
      goto L105;
   }

   t1 = 0;
   t2 = t0;
   for (k = 0; k < l1; k++)
   {
      t3 = t2;
      t4 = (t1 << 1) + (ido << 1);
      t5 = t1;
      t6 = t1 + t1;
      for (i = 2; i < ido; i += 2)
      {
         t3 += 2;
         t4 -= 2;
         t5 += 2;
         t6 += 2;
         tr2 = wa1[i - 2] * cc[t3 - 1] + wa1[i - 1] * cc[t3];
         ti2 = wa1[i - 2] * cc[t3] - wa1[i - 1] * cc[t3 - 1];
         ch[t6] = cc[t5] + ti2;
         ch[t4] = ti2 - cc[t5];
         ch[t6 - 1] = cc[t5 - 1] + tr2;
         ch[t4 - 1] = cc[t5 - 1] - tr2;
      }
      t1 += ido;
      t2 += ido;
   }

   if (ido % 2 == 1)
   {
      return;
   }

L105:
   t3 = (t2 = (t1 = ido) - 1);
   t2 += t0;
   for (k = 0; k < l1; k++)
   {
      ch[t1] = -cc[t2];
      ch[t1 - 1] = cc[t3];
      t1 += ido << 1;
      t2 += ido;
      t3 += ido;
   }
}

void dradf4 (int ido, int l1, double *cc, double *ch, double *wa1, double *wa2, double *wa3)
{
   constexpr double hsqt2 = .70710678118654752440084436210485;
   int i, k, t0, t1, t2, t3, t4, t5, t6;
   double ci2, ci3, ci4, cr2, cr3, cr4;
   double ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

   t0 = l1 * ido;
   t1 = t0;
   t4 = t1 << 1;
   t2 = t1 + (t1 << 1);
   t3 = 0;

   for (k = 0; k < l1; k++)
   {
      tr1 = cc[t1] + cc[t2];
      tr2 = cc[t3] + cc[t4];
      ch[t5 = t3 << 2] = tr1 + tr2;
      ch[(ido << 2) + t5 - 1] = tr2 - tr1;
      ch[(t5 += (ido << 1)) - 1] = cc[t3] - cc[t4];
      ch[t5] = cc[t2] - cc[t1];

      t1 += ido;
      t2 += ido;
      t3 += ido;
      t4 += ido;
   }

   if (ido < 2)
   {
      return;
   }
   if (ido == 2)
   {
      goto L105;
   }

   t1 = 0;
   for (k = 0; k < l1; k++)
   {
      t2 = t1;
      t4 = t1 << 2;
      t5 = (t6 = ido << 1) + t4;
      for (i = 2; i < ido; i += 2)
      {
         t3 = (t2 += 2);
         t4 += 2;
         t5 -= 2;

         t3 += t0;
         cr2 = wa1[i - 2] * cc[t3 - 1] + wa1[i - 1] * cc[t3];
         ci2 = wa1[i - 2] * cc[t3] - wa1[i - 1] * cc[t3 - 1];
         t3 += t0;
         cr3 = wa2[i - 2] * cc[t3 - 1] + wa2[i - 1] * cc[t3];
         ci3 = wa2[i - 2] * cc[t3] - wa2[i - 1] * cc[t3 - 1];
         t3 += t0;
         cr4 = wa3[i - 2] * cc[t3 - 1] + wa3[i - 1] * cc[t3];
         ci4 = wa3[i - 2] * cc[t3] - wa3[i - 1] * cc[t3 - 1];

         tr1 = cr2 + cr4;
         tr4 = cr4 - cr2;
         ti1 = ci2 + ci4;
         ti4 = ci2 - ci4;
         ti2 = cc[t2] + ci3;
         ti3 = cc[t2] - ci3;
         tr2 = cc[t2 - 1] + cr3;
         tr3 = cc[t2 - 1] - cr3;


         ch[t4 - 1] = tr1 + tr2;
         ch[t4] = ti1 + ti2;

         ch[t5 - 1] = tr3 - ti4;
         ch[t5] = tr4 - ti3;

         ch[t4 + t6 - 1] = ti4 + tr3;
         ch[t4 + t6] = tr4 + ti3;

         ch[t5 + t6 - 1] = tr2 - tr1;
         ch[t5 + t6] = ti1 - ti2;
      }
      t1 += ido;
   }
   if (ido % 2 == 1)
   {
      return;
   }

L105:

   t2 = (t1 = t0 + ido - 1) + (t0 << 1);
   t3 = ido << 2;
   t4 = ido;
   t5 = ido << 1;
   t6 = ido;

   for (k = 0; k < l1; k++)
   {
      ti1 = -hsqt2 * (cc[t1] + cc[t2]);
      tr1 = hsqt2 * (cc[t1] - cc[t2]);
      ch[t4 - 1] = tr1 + cc[t6 - 1];
      ch[t4 + t5 - 1] = cc[t6 - 1] - tr1;
      ch[t4] = ti1 - cc[t1 + t0];
      ch[t4 + t5] = ti1 + cc[t1 + t0];
      t1 += ido;
      t2 += ido;
      t4 += t3;
      t6 += ido;
   }
}

void dradfg (int ido, int ip, int l1, int idl1, double *cc, double *c1, double *c2, double *ch, double *ch2, double *wa)
{
   constexpr double two_pi = 6.28318530717958647692528676655900577;
   int idij = 0;
   int ipph = 0;
   int i = 0;
   int j = 0;
   int k = 0;
   int l = 0;
   int ic = 0;
   int ik = 0;
   int is = 0;
   int t0 = 0;
   int t1 = 0;
   int t2 = 0;
   int t3 = 0;
   int t4 = 0;
   int t5 = 0;
   int t6 = 0;
   int t7 = 0;
   int t8 = 0;
   int t9 = 0;
   int t10 = 0;
   double dc2 = 0.0;
   double ai1 = 0.0;
   double ai2 = 0;
   double ar1 = 0;
   double ar2 = 0;
   double ds2 = 0;
   int nbd = 0;
   double dcp = 0;
   double arg = 0;
   double dsp = 0;
   double ar1h = 0;
   double ar2h = 0;
   int idp2 = 0;
   int ipp2 = 0;

   arg = two_pi / (double) ip;
   dcp = cos (arg);
   dsp = sin (arg);
   ipph = (ip + 1) >> 1;
   ipp2 = ip;
   idp2 = ido;
   nbd = (ido - 1) >> 1;
   t0 = l1 * ido;
   t10 = ip * ido;

   if (ido == 1)
   {
      goto L119;
   }
   for (ik = 0; ik < idl1; ik++)
   {
      ch2[ik] = c2[ik];
   }

   t1 = 0;
   for (j = 1; j < ip; j++)
   {
      t1 += t0;
      t2 = t1;
      for (k = 0; k < l1; k++)
      {
         ch[t2] = c1[t2];
         t2 += ido;
      }
   }

   is = -ido;
   t1 = 0;
   if (nbd > l1)
   {
      for (j = 1; j < ip; j++)
      {
         t1 += t0;
         is += ido;
         t2 = -ido + t1;
         for (k = 0; k < l1; k++)
         {
            idij = is - 1;
            t2 += ido;
            t3 = t2;
            for (i = 2; i < ido; i += 2)
            {
               idij += 2;
               t3 += 2;
               ch[t3 - 1] = wa[idij - 1] * c1[t3 - 1] + wa[idij] * c1[t3];
               ch[t3] = wa[idij - 1] * c1[t3] - wa[idij] * c1[t3 - 1];
            }
         }
      }
   }
   else
   {
      for (j = 1; j < ip; j++)
      {
         is += ido;
         idij = is - 1;
         t1 += t0;
         t2 = t1;
         for (i = 2; i < ido; i += 2)
         {
            idij += 2;
            t2 += 2;
            t3 = t2;
            for (k = 0; k < l1; k++)
            {
               ch[t3 - 1] = wa[idij - 1] * c1[t3 - 1] + wa[idij] * c1[t3];
               ch[t3] = wa[idij - 1] * c1[t3] - wa[idij] * c1[t3 - 1];
               t3 += ido;
            }
         }
      }
   }

   t1 = 0;
   t2 = ipp2 * t0;
   if (nbd < l1)
   {
      for (j = 1; j < ipph; j++)
      {
         t1 += t0;
         t2 -= t0;
         t3 = t1;
         t4 = t2;
         for (i = 2; i < ido; i += 2)
         {
            t3 += 2;
            t4 += 2;
            t5 = t3 - ido;
            t6 = t4 - ido;
            for (k = 0; k < l1; k++)
            {
               t5 += ido;
               t6 += ido;
               c1[t5 - 1] = ch[t5 - 1] + ch[t6 - 1];
               c1[t6 - 1] = ch[t5] - ch[t6];
               c1[t5] = ch[t5] + ch[t6];
               c1[t6] = ch[t6 - 1] - ch[t5 - 1];
            }
         }
      }
   }
   else
   {
      for (j = 1; j < ipph; j++)
      {
         t1 += t0;
         t2 -= t0;
         t3 = t1;
         t4 = t2;
         for (k = 0; k < l1; k++)
         {
            t5 = t3;
            t6 = t4;
            for (i = 2; i < ido; i += 2)
            {
               t5 += 2;
               t6 += 2;
               c1[t5 - 1] = ch[t5 - 1] + ch[t6 - 1];
               c1[t6 - 1] = ch[t5] - ch[t6];
               c1[t5] = ch[t5] + ch[t6];
               c1[t6] = ch[t6 - 1] - ch[t5 - 1];
            }
            t3 += ido;
            t4 += ido;
         }
      }
   }

L119:
   for (ik = 0; ik < idl1; ik++)
   {
      c2[ik] = ch2[ik];
   }

   t1 = 0;
   t2 = ipp2 * idl1;
   for (j = 1; j < ipph; j++)
   {
      t1 += t0;
      t2 -= t0;
      t3 = t1 - ido;
      t4 = t2 - ido;
      for (k = 0; k < l1; k++)
      {
         t3 += ido;
         t4 += ido;
         c1[t3] = ch[t3] + ch[t4];
         c1[t4] = ch[t4] - ch[t3];
      }
   }

   ar1 = 1.0;
   ai1 = 0.0;
   t1 = 0;
   t2 = ipp2 * idl1;
   t3 = (ip - 1) * idl1;
   for (l = 1; l < ipph; l++)
   {
      t1 += idl1;
      t2 -= idl1;
      ar1h = dcp * ar1 - dsp * ai1;
      ai1 = dcp * ai1 + dsp * ar1;
      ar1 = ar1h;
      t4 = t1;
      t5 = t2;
      t6 = t3;
      t7 = idl1;

      for (ik = 0; ik < idl1; ik++)
      {
         ch2[t4++] = c2[ik] + ar1 * c2[t7++];
         ch2[t5++] = ai1 * c2[t6++];
      }

      dc2 = ar1;
      ds2 = ai1;
      ar2 = ar1;
      ai2 = ai1;

      t4 = idl1;
      t5 = (ipp2 - 1) * idl1;
      for (j = 2; j < ipph; j++)
      {
         t4 += idl1;
         t5 -= idl1;

         ar2h = dc2 * ar2 - ds2 * ai2;
         ai2 = dc2 * ai2 + ds2 * ar2;
         ar2 = ar2h;

         t6 = t1;
         t7 = t2;
         t8 = t4;
         t9 = t5;
         for (ik = 0; ik < idl1; ik++)
         {
            ch2[t6++] += ar2 * c2[t8++];
            ch2[t7++] += ai2 * c2[t9++];
         }
      }
   }

   t1 = 0;
   for (j = 1; j < ipph; j++)
   {
      t1 += idl1;
      t2 = t1;
      for (ik = 0; ik < idl1; ik++)
      {
         ch2[ik] += c2[t2++];
      }
   }

   if (ido < l1)
   {
      goto L132;
   }

   t1 = 0;
   t2 = 0;
   for (k = 0; k < l1; k++)
   {
      t3 = t1;
      t4 = t2;
      for (i = 0; i < ido; i++)
      {
         cc[t4++] = ch[t3++];
      }
      t1 += ido;
      t2 += t10;
   }

   goto L135;

L132:
   for (i = 0; i < ido; i++)
   {
      t1 = i;
      t2 = i;
      for (k = 0; k < l1; k++)
      {
         cc[t2] = ch[t1];
         t1 += ido;
         t2 += t10;
      }
   }

L135:
   t1 = 0;
   t2 = ido << 1;
   t3 = 0;
   t4 = ipp2 * t0;
   for (j = 1; j < ipph; j++)
   {

      t1 += t2;
      t3 += t0;
      t4 -= t0;

      t5 = t1;
      t6 = t3;
      t7 = t4;

      for (k = 0; k < l1; k++)
      {
         cc[t5 - 1] = ch[t6];
         cc[t5] = ch[t7];
         t5 += t10;
         t6 += ido;
         t7 += ido;
      }
   }

   if (ido == 1)
   {
      return;
   }
   if (nbd < l1)
   {
      goto L141;
   }

   t1 = -ido;
   t3 = 0;
   t4 = 0;
   t5 = ipp2 * t0;
   for (j = 1; j < ipph; j++)
   {
      t1 += t2;
      t3 += t2;
      t4 += t0;
      t5 -= t0;
      t6 = t1;
      t7 = t3;
      t8 = t4;
      t9 = t5;
      for (k = 0; k < l1; k++)
      {
         for (i = 2; i < ido; i += 2)
         {
            ic = idp2 - i;
            cc[i + t7 - 1] = ch[i + t8 - 1] + ch[i + t9 - 1];
            cc[ic + t6 - 1] = ch[i + t8 - 1] - ch[i + t9 - 1];
            cc[i + t7] = ch[i + t8] + ch[i + t9];
            cc[ic + t6] = ch[i + t9] - ch[i + t8];
         }
         t6 += t10;
         t7 += t10;
         t8 += ido;
         t9 += ido;
      }
   }
   return;

L141:

   t1 = -ido;
   t3 = 0;
   t4 = 0;
   t5 = ipp2 * t0;
   for (j = 1; j < ipph; j++)
   {
      t1 += t2;
      t3 += t2;
      t4 += t0;
      t5 -= t0;
      for (i = 2; i < ido; i += 2)
      {
         t6 = idp2 + t1 - i;
         t7 = i + t3;
         t8 = i + t4;
         t9 = i + t5;
         for (k = 0; k < l1; k++)
         {
            cc[t7 - 1] = ch[t8 - 1] + ch[t9 - 1];
            cc[t6 - 1] = ch[t8 - 1] - ch[t9 - 1];
            cc[t7] = ch[t8] + ch[t9];
            cc[t6] = ch[t9] - ch[t8];
            t6 += t10;
            t7 += t10;
            t8 += ido;
            t9 += ido;
         }
      }
   }
}

void drftf1 (int n, double *c, double *ch, double *wa, int *ifac)
{
   int i, k1, l1, l2;
   int na, kh, nf;
   int ip, iw, ido, idl1, ix2, ix3;

   nf = ifac[1];
   na = 1;
   l2 = n;
   iw = n;

   for (k1 = 0; k1 < nf; k1++)
   {
      kh = nf - k1;
      ip = ifac[kh + 1];
      l1 = l2 / ip;
      ido = n / l2;
      idl1 = ido * l1;
      iw -= (ip - 1) * ido;
      na = 1 - na;

      if (ip != 4)
      {
         goto L102;
      }

      ix2 = iw + ido;
      ix3 = ix2 + ido;
      if (na != 0)
      {
         dradf4 (ido, l1, ch, c, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);
      }
      else
      {
         dradf4 (ido, l1, c, ch, wa + iw - 1, wa + ix2 - 1, wa + ix3 - 1);
      }
      goto L110;

   L102:
      if (ip != 2)
      {
         goto L104;
      }
      if (na != 0)
      {
         goto L103;
      }

      dradf2 (ido, l1, c, ch, wa + iw - 1);
      goto L110;

   L103:
      dradf2 (ido, l1, ch, c, wa + iw - 1);
      goto L110;

   L104:
      if (ido == 1)
      {
         na = 1 - na;
      }
      if (na != 0)
      {
         goto L109;
      }

      dradfg (ido, ip, l1, idl1, c, c, c, ch, ch, wa + iw - 1);
      na = 1;
      goto L110;

   L109:
      dradfg (ido, ip, l1, idl1, ch, ch, ch, c, c, wa + iw - 1);
      na = 0;

   L110:
      l2 = l1;
   }

   if (na == 1)
   {
      return;
   }

   for (i = 0; i < n; i++)
   {
      c[i] = ch[i];
   }
}

void __ogg_fdrfftf (int n, double *r, double *wsave, int *ifac)
{
   if (n == 1)
   {
      return;
   }
   drftf1 (n, r, wsave, wsave + n, ifac);
}

void DiscreteFourierTransform (int n)
{
   double p_value, upperBound, percentile, N_l, N_o, d, *m = nullptr, *X = nullptr, *wsave = nullptr;

   int i, count, ifac[15];

   if (((X = (double *) std::calloc (n, sizeof (double))) == nullptr) || ((wsave = (double *) std::calloc (2 * n, sizeof (double))) == nullptr) ||
       ((m = (double *) std::calloc (n / 2 + 1, sizeof (double))) == nullptr))
   {
      std::printf ("  Unable to allocate working arrays for the DFT.\n");
      if (X != nullptr)
      {
         std::free (X);
      }
      if (wsave != nullptr)
      {
         std::free (wsave);
      }
      if (m != nullptr)
      {
         std::free (m);
      }
      return;
   }
   for (i = 0; i < n; i++)
   {
      X[i] = 2 * (int) epsilon[i] - 1;
   }

   __ogg_fdrffti (n, wsave, ifac);    /* INITIALIZE WORK ARRAYS */
   __ogg_fdrfftf (n, X, wsave, ifac); /* APPLY FORWARD FFT */

   m[0] = std::sqrt (X[0] * X[0]); /* COMPUTE MAGNITUDE */

   for (i = 0; i < n / 2; i++)
   {
      m[i + 1] = std::sqrt (std::pow (X[2 * i + 1], 2) + std::pow (X[2 * i + 2], 2));
   }
   count = 0; /* Confidence Interval */
   upperBound = std::sqrt (2.995732274 * n);
   for (i = 0; i < n / 2; i++)
   {
      if (m[i] < upperBound)
      {
         count++;
      }
   }
   percentile = (double) count / (n / 2) * 100;
   N_l = (double) count; /* number of peaks less than h = std::sqrt(3*n) */
   N_o = (double) 0.95 * n / 2.0;
   d = (N_l - N_o) / std::sqrt (n / 4.0 * 0.95 * 0.05);
   p_value = erfc (std::fabs (d) / std::sqrt (2.0));

   std::printf ("    FFT Test\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  (a) Percentile = %f\n", percentile);
   std::printf ("  (b) N_l        = %f\n", N_l);
   std::printf ("  (c) N_o        = %f\n", N_o);
   std::printf ("  (d) d          = %f\n", d);
   std::printf ("  -------------------------------------------\n");

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   std::free (X);
   std::free (wsave);
   std::free (m);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                A P P R O X I M A T E  E N T R O P Y   T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void ApproximateEntropy (int m, int n)
{
   int i = 0;
   int j = 0;
   int k = 0;
   int r = 0;
   int blockSize = 0;
   int seqLength = 0;
   int powLen = 0;
   int index = 0;
   double sum = 0;
   double numOfBlocks = 0;
   double ApEn[2] = {0};
   double apen = 0;
   double chi_squared = 0;
   double p_value = 0;
   unsigned int *P = nullptr;

   std::printf ("   Approximate Entropy Test\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  (a) m (block length)    = %d\n", m);

   seqLength = n;
   r = 0;

   for (blockSize = m; blockSize <= m + 1; blockSize++)
   {
      if (blockSize == 0)
      {
         ApEn[0] = 0.00;
         r++;
      }
      else
      {
         numOfBlocks = (double) seqLength;
         powLen = (int) std::pow (2, blockSize + 1) - 1;
         if ((P = (unsigned int *) std::calloc (powLen, sizeof (unsigned int))) == nullptr)
         {
            std::printf ("ApEn:  Insufficient memory available.\n");
            return;
         }
         for (i = 1; i < powLen - 1; i++)
         {
            P[i] = 0;
         }
         for (i = 0; i < numOfBlocks; i++)
         { /* Compute Frequency */
            k = 1;
            for (j = 0; j < blockSize; j++)
            {
               k <<= 1;
               if ((int) epsilon[(i + j) % seqLength] == 1)
               {
                  k++;
               }
            }
            P[k - 1]++;
         }
         /* Display Frequency */
         sum = 0.0;
         index = (int) std::pow (2, blockSize) - 1;
         for (i = 0; i < (int) std::pow (2, blockSize); i++)
         {
            if (P[index] > 0)
            {
               sum += P[index] * std::log (P[index] / numOfBlocks);
            }
            index++;
         }
         sum /= numOfBlocks;
         ApEn[r] = sum;
         r++;
         std::free (P);
      }
   }
   apen = ApEn[0] - ApEn[1];

   chi_squared = 2.0 * seqLength * (log (2) - apen);
   p_value = cephes_igamc (std::pow (2, m - 1), chi_squared / 2.0);

   std::printf ("  (b) n (sequence length) = %d\n", seqLength);
   std::printf ("  (c) Chi^2               = %f\n", chi_squared);
   std::printf ("  (d) Phi(m)	       = %f\n", ApEn[0]);
   std::printf ("  (e) Phi(m+1)	       = %f\n", ApEn[1]);
   std::printf ("  (f) ApEn                = %f\n", apen);
   std::printf ("  (g) Log(2)              = %f\n", std::log (2.0));
   std::printf ("  --------------------------------------------\n");

   if (m > (int) (std::log (seqLength) / std::log (2) - 5))
   {
      std::printf ("  Note: The blockSize = %d exceeds recommended value of %d\n", m, MAX (1, (int) (std::log (seqLength) / std::log (2) - 5)));
      std::printf ("  Results are inaccurate!\n");
      std::printf ("  --------------------------------------------\n");
   }

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
}
/* ======================================================================== *
 * Test harness:  ibm-clang++ -m64 -D ZOS_GETENTROPY_TEST ...
 * ======================================================================== */


void Frequency (int n)
{
   int i = 0;
   double f = 0;
   double s_obs = 0;
   double p_value = 0;
   double sum = 0;
   constexpr double sqrt2 = 1.41421356237309504880;

   sum = 0.0;
   for (i = 0; i < n; i++)
   {
      sum += 2 * (int) epsilon[i] - 1;
   }
   s_obs = std::fabs (sum) / std::sqrt (n);
   f = s_obs / sqrt2;
   p_value = erfc (f);

   std::printf ("         Frequency Test\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  (a) The nth partial sum = %d\n", (int) sum);
   std::printf ("  (b) S_n/n               = %f\n", sum / n);
   std::printf ("  ---------------------------------------------\n");

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
}



void Serial (int m, int n)
{
   double p_value1 = 0;
   double p_value2 = 0;
   double psim0 = 0;
   double psim1 = 0;
   double psim2 = 0;
   double del1 = 0;
   double del2 = 0;

   psim0 = psi2 (m, n);
   psim1 = psi2 (m - 1, n);
   psim2 = psi2 (m - 2, n);
   del1 = psim0 - psim1;
   del2 = psim0 - 2.0 * psim1 + psim2;
   p_value1 = cephes_igamc (std::pow (2, m - 1) / 2, del1 / 2.0);
   p_value2 = cephes_igamc (std::pow (2, m - 2) / 2, del2 / 2.0);

   std::printf ("          Serial Test\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("   Computational Information:		  \n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  (a) Block length    (m) = %d\n", m);
   std::printf ("  (b) Sequence length (n) = %d\n", n);
   std::printf ("  (c) Psi_m               = %f\n", psim0);
   std::printf ("  (d) Psi_m-1             = %f\n", psim1);
   std::printf ("  (e) Psi_m-2             = %f\n", psim2);
   std::printf ("  (f) Del_1               = %f\n", del1);
   std::printf ("  (g) Del_2               = %f\n", del2);
   std::printf ("  ---------------------------------------------\n");

   std::printf ("%s  p_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
   std::printf ("%s  p_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    C U M U L A T I V E  S U M S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void CumulativeSums (int n)
{
   int S, sup, inf, z, zrev, k;
   double sum1, sum2, p_value;

   S = 0;
   sup = 0;
   inf = 0;

   for (k = 0; k < n; k++)
   {
      epsilon[k] ? S++ : S--;
      if (S > sup)
      {
         sup++;
      }
      if (S < inf)
      {
         inf--;
      }
      z = (sup > -inf) ? sup : -inf;
      zrev = (sup - S > S - inf) ? sup - S : S - inf;
   }

   // forward
   sum1 = 0.0;
   for (k = (-n / z + 1) / 4; k <= (n / z - 1) / 4; k++)
   {
      sum1 += cephes_normal (((4 * k + 1) * z) / std::sqrt (n));
      sum1 -= cephes_normal (((4 * k - 1) * z) / std::sqrt (n));
   }
   sum2 = 0.0;
   for (k = (-n / z - 3) / 4; k <= (n / z - 1) / 4; k++)
   {
      sum2 += cephes_normal (((4 * k + 3) * z) / std::sqrt (n));
      sum2 -= cephes_normal (((4 * k + 1) * z) / std::sqrt (n));
   }

   p_value = 1.0 - sum1 + sum2;

   std::printf ("        Cumulative Sums (Forward) Test\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  (a) The maximum partial sum = %d\n", z);
   std::printf ("  -------------------------------------------\n");

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      std::printf ("  Warning:  P_Value Is Out Of Range\n");
   }

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   // backwards
   sum1 = 0.0;
   for (k = (-n / zrev + 1) / 4; k <= (n / zrev - 1) / 4; k++)
   {
      sum1 += cephes_normal (((4 * k + 1) * zrev) / std::sqrt (n));
      sum1 -= cephes_normal (((4 * k - 1) * zrev) / std::sqrt (n));
   }
   sum2 = 0.0;
   for (k = (-n / zrev - 3) / 4; k <= (n / zrev - 1) / 4; k++)
   {
      sum2 += cephes_normal (((4 * k + 3) * zrev) / std::sqrt (n));
      sum2 -= cephes_normal (((4 * k + 1) * zrev) / std::sqrt (n));
   }
   p_value = 1.0 - sum1 + sum2;

   std::printf ("        Cumulative Sums (Reverse) Test\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  -------------------------------------------\n");
   std::printf ("  (a) The maximum partial sum = %d\n", zrev);
   std::printf ("  -------------------------------------------\n");

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      std::printf ("  WARNING:  P_VALUE IS OUT OF RANGE\n");
   }

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
}

double psi2 (int m, int n)
{
   int i, j, k, powLen;
   double sum, numOfBlocks;
   unsigned int *P;

   if ((m == 0) || (m == -1))
   {
      return 0.0;
   }
   numOfBlocks = n;
   powLen = (int) std::pow (2, m + 1) - 1;
   if ((P = (unsigned int *) std::calloc (powLen, sizeof (unsigned int))) == nullptr)
   {
      std::printf ("Serial Test:  Insufficient memory available.\n");
      return 0.0;
   }
   for (i = 1; i < powLen - 1; i++)
   {
      P[i] = 0; /* Initialize Nodes */
   }
   for (i = 0; i < numOfBlocks; i++)
   { /* Compute Frequency */
      k = 1;
      for (j = 0; j < m; j++)
      {
         if (epsilon[(i + j) % n] == 0)
         {
            k *= 2;
         }
         else if (epsilon[(i + j) % n] == 1)
         {
            k = 2 * k + 1;
         }
      }
      P[k - 1]++;
   }
   sum = 0.0;
   for (i = (int) std::pow (2, m) - 1; i < (int) std::pow (2, m + 1) - 1; i++)
   {
      sum += std::pow (P[i], 2);
   }
   sum = (sum * std::pow (2, m) / (double) n) - (double) n;
   std::free (P);

   return sum;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                         U N I V E R S A L  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Universal (int n)
{
   int i, j, p, L, Q, K;
   double arg, sqrt2, sigma, phi, sum, p_value, c;
   long *T, decRep;
   constexpr double expected_value[17] = {
       0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656, 8.1764248, 9.1723243, 10.170032, 11.168765, 12.168070, 13.167693, 14.167488, 15.167379};
   constexpr double variance[17] = {0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384, 3.401, 3.410, 3.416, 3.419, 3.421};

   /* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
    * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
    * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
   L = 5;
   if (n >= 387840)
   {
      L = 6;
   }
   if (n >= 904960)
   {
      L = 7;
   }
   if (n >= 2068480)
   {
      L = 8;
   }
   if (n >= 4654080)
   {
      L = 9;
   }
   if (n >= 10342400)
   {
      L = 10;
   }
   if (n >= 22753280)
   {
      L = 11;
   }
   if (n >= 49643520)
   {
      L = 12;
   }
   if (n >= 107560960)
   {
      L = 13;
   }
   if (n >= 231669760)
   {
      L = 14;
   }
   if (n >= 496435200)
   {
      L = 15;
   }
   if (n >= 1059061760)
   {
      L = 16;
   }

   Q = 10 * (int) std::pow (2, L);
   K = (int) (floor (n / L) - (double) Q); /* BLOCKS TO TEST */

   p = (int) std::pow (2, L);
   if ((L < 6) || (L > 16) || ((double) Q < 10 * std::pow (2, L)) || ((T = (long *) std::calloc (p, sizeof (long))) == nullptr))
   {
      std::printf ("  Universal Statistical Test\n");
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  ERROR:  L is out of range.\n");
      std::printf ("  -OR- :  Q is less than %f.\n", 10 * std::pow (2, L));
      std::printf ("  -OR- :  Unable to allocate T.\n");
      return;
   }

   /* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
   c = 0.7 - 0.8 / (double) L + (4 + 32 / (double) L) * std::pow (K, -3 / (double) L) / 15;
   sigma = c * std::sqrt (variance[L] / (double) K);
   sqrt2 = std::sqrt (2);
   sum = 0.0;
   for (i = 0; i < p; i++)
   {
      T[i] = 0;
   }
   for (i = 1; i <= Q; i++)
   { /* INITIALIZE TABLE */
      decRep = 0;
      for (j = 0; j < L; j++)
      {
         decRep += epsilon[(i - 1) * L + j] * (long) std::pow (2, L - 1 - j);
      }
      T[decRep] = i;
   }
   for (i = Q + 1; i <= Q + K; i++)
   { /* PROCESS BLOCKS */
      decRep = 0;
      for (j = 0; j < L; j++)
      {
         decRep += epsilon[(i - 1) * L + j] * (long) std::pow (2, L - 1 - j);
      }
      sum += log (i - T[decRep]) / log (2);
      T[decRep] = i;
   }
   phi = (double) (sum / (double) K);

   std::printf ("  Universal Statistical Test\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  (a) L         = %d\n", L);
   std::printf ("  (b) Q         = %d\n", Q);
   std::printf ("  (c) K         = %d\n", K);
   std::printf ("  (d) sum       = %f\n", sum);
   std::printf ("  (e) sigma     = %f\n", sigma);
   std::printf ("  (f) variance  = %f\n", variance[L]);
   std::printf ("  (g) exp_value = %f\n", expected_value[L]);
   std::printf ("  (h) phi       = %f\n", phi);
   std::printf ("  (i) WARNING:  %d bits were discarded.\n", n - (Q + K) * L);
   std::printf ("  -----------------------------------------\n");

   arg = std::fabs (phi - expected_value[L]) / (sqrt2 * sigma);
   p_value = std::erfc (arg);

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      std::printf ("  WARNING:  P_VALUE IS OUT OF RANGE\n");
   }

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   std::free (T);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     R A N D O M  E X C U R S I O N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandomExcursions (int n)
{
   int b, i, j, k, J, x;
   int cycleStart, cycleStop, *cycle = nullptr, *S_k = nullptr;
   constexpr int stateX[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
   int counter[8] = {0, 0, 0, 0, 0, 0, 0, 0};
   double p_value, sum, constraint, nu[6][8];
   constexpr double pi[5][6] = {{0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000},
                      {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
                      {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
                      {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
                      {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051}};

   if (((S_k = (int *) std::calloc (n, sizeof (int))) == nullptr) || ((cycle = (int *) std::calloc (MAX (1000, n / 100), sizeof (int))) == nullptr))
   {
      std::printf ("Random Excursions Test:  Insufficient Work Space Allocated.\n");
      if (S_k != nullptr)
      {
         std::free (S_k);
      }
      if (cycle != nullptr)
      {
         std::free (cycle);
      }
      return;
   }

   J = 0; /* DETERMINE CYCLES */
   S_k[0] = 2 * (int) epsilon[0] - 1;
   for (i = 1; i < n; i++)
   {
      S_k[i] = S_k[i - 1] + 2 * epsilon[i] - 1;
      if (S_k[i] == 0)
      {
         J++;
         if (J > MAX (1000, n / 100))
         {
            std::printf ("ERROR IN FUNCTION randomExcursions:  EXCEEDING THE MAX NUMBER OF CYCLES EXPECTED\n.");
            std::free (S_k);
            std::free (cycle);
            return;
         }
         cycle[J] = i;
      }
   }
   if (S_k[n - 1] != 0)
   {
      J++;
   }
   cycle[J] = n;

   std::printf ("     Random Excursions Test\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  (a) Number Of Cycles (J) = %04d\n", J);
   std::printf ("  (b) Sequence Length (n)  = %d\n", n);

   constraint = MAX (0.005 * std::pow (n, 0.5), 500);
   if (J < constraint)
   {
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  WARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
      std::printf ("     INSUFFICIENT NUMBER OF CYCLES.\n");
      std::printf ("  ---------------------------------------------\n");
      for (i = 0; i < 8; i++)
      {
         std::printf (results[TEST_RND_EXCURSION], "%f\n", 0.0);
      }
   }
   else
   {
      std::printf ("  (c) Rejection Constraint = %f\n", constraint);
      std::printf ("  -------------------------------------------\n");

      cycleStart = 0;
      cycleStop = cycle[1];
      for (k = 0; k < 6; k++)
      {
         for (i = 0; i < 8; i++)
         {
            nu[k][i] = 0.;
         }
      }
      for (j = 1; j <= J; j++)
      { /* FOR EACH CYCLE */
         for (i = 0; i < 8; i++)
         {
            counter[i] = 0;
         }
         for (i = cycleStart; i < cycleStop; i++)
         {
            if ((S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1))
            {
               if (S_k[i] < 0)
               {
                  b = 4;
               }
               else
               {
                  b = 3;
               }
               counter[S_k[i] + b]++;
            }
         }
         cycleStart = cycle[j] + 1;
         if (j < J)
         {
            cycleStop = cycle[j + 1];
         }

         for (i = 0; i < 8; i++)
         {
            if ((counter[i] >= 0) && (counter[i] <= 4))
            {
               nu[counter[i]][i]++;
            }
            else if (counter[i] >= 5)
            {
               nu[5][i]++;
            }
         }
      }

      for (i = 0; i < 8; i++)
      {
         x = stateX[i];
         sum = 0.;
         for (k = 0; k < 6; k++)
         {
            sum += std::pow (nu[k][i] - J * pi[(int) std::fabs (x)][k], 2) / (J * pi[(int) std::fabs (x)][k]);
         }
         p_value = cephes_igamc (2.5, sum / 2.0);

         if (isNegative (p_value) || isGreaterThanOne (p_value))
         {
            std::printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
         }

         std::printf ("%s  x = %2d chi^2 = %9.6f p_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
      }
   }
   std::printf ("\n");

   std::free (S_k);
   std::free (cycle);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            R A N D O M  E X C U R S I O N S  V A R I A N T  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void RandomExcursionsVariant (int n)
{
   int i, p, J, x, constraint, count, *S_k;
   constexpr int stateX[18] = {-9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   double p_value;

   if ((S_k = (int *) std::calloc (n, sizeof (int))) == nullptr)
   {
      std::printf ("  RANDOM EXCURSIONS VARIANT: Insufficient memory allocated.\n");
      return;
   }
   J = 0;
   S_k[0] = 2 * (int) epsilon[0] - 1;
   for (i = 1; i < n; i++)
   {
      S_k[i] = S_k[i - 1] + 2 * epsilon[i] - 1;
      if (S_k[i] == 0)
      {
         J++;
      }
   }
   if (S_k[n - 1] != 0)
   {
      J++;
   }

   std::printf ("   Random Excursions Variant Test\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  Computational Information:\n");
   std::printf ("  --------------------------------------------\n");
   std::printf ("  (a) Number Of Cycles (J) = %d\n", J);
   std::printf ("  (b) Sequence Length (n)  = %d\n", n);
   std::printf ("  --------------------------------------------\n");

   constraint = (int) MAX (0.005 * std::pow (n, 0.5), 500);
   if (J < constraint)
   {
      std::printf ("\n  WARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
      std::printf ("     INSUFFICIENT NUMBER OF CYCLES.\n");
      std::printf ("  ---------------------------------------------\n");
      for (i = 0; i < 18; i++)
      {
         std::printf ("%f\n", 0.0);
      }
   }
   else
   {
      for (p = 0; p <= 17; p++)
      {
         x = stateX[p];
         count = 0;
         for (i = 0; i < n; i++)
         {
            if (S_k[i] == x)
            {
               count++;
            }
         }
         p_value = std::erfc (std::fabs (count - J) / (std::sqrt (2.0 * J * (4.0 * std::fabs (x) - 2))));

         if (isNegative (p_value) || isGreaterThanOne (p_value))
         {
            std::printf ("  (b) WARNING: P_VALUE IS OUT OF RANGE.\n");
         }
         std::printf ("%s  ", p_value < ALPHA ? "FAILURE" : "SUCCESS");
         std::printf ("(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
      }
   }
   std::printf ("\n");

   std::free (S_k);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


void OverlappingTemplateMatchings (int m, int n)
{
   int i, k, match;
   double W_obs, eta, sum, chi2, p_value, lambda;
   int M, N, j, K = 5;
   unsigned int nu[6] = {0, 0, 0, 0, 0, 0};
   // double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };
   double pi[6] = {0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865};
   unsigned char *sequence;

   M = 1032;
   N = n / M;

   if ((sequence = (unsigned char *) std::calloc (m, sizeof (unsigned char))) == nullptr)
   {
      std::printf ("      OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  TEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
   }
   else
   {
      for (i = 0; i < m; i++)
      {
         sequence[i] = 1;
      }
   }

   lambda = (double) (M - m + 1) / std::pow (2, m);
   eta = lambda / 2.0;
   sum = 0.0;
   for (i = 0; i < K; i++)
   { /* Compute Probabilities */
      pi[i] = Pr (i, eta);
      sum += pi[i];
   }
   pi[K] = 1 - sum;

   for (i = 0; i < N; i++)
   {
      W_obs = 0;
      for (j = 0; j < M - m + 1; j++)
      {
         match = 1;
         for (k = 0; k < m; k++)
         {
            if (sequence[k] != epsilon[i * M + j + k])
            {
               match = 0;
            }
         }
         if (match == 1)
         {
            W_obs++;
         }
      }
      if (W_obs <= 4)
      {
         nu[(int) W_obs]++;
      }
      else
      {
         nu[K]++;
      }
   }
   sum = 0;
   chi2 = 0.0; /* Compute Chi Square */
   for (i = 0; i < K + 1; i++)
   {
      chi2 += std::pow ((double) nu[i] - (double) N * pi[i], 2) / ((double) N * pi[i]);
      sum += nu[i];
   }
   p_value = cephes_igamc (K / 2.0, chi2 / 2.0);

   std::printf ("      OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
   std::printf ("  -----------------------------------------------\n");
   std::printf ("  COMPUTATIONAL INFORMATION:\n");
   std::printf ("  -----------------------------------------------\n");
   std::printf ("  (a) n (sequence_length)      = %d\n", n);
   std::printf ("  (b) m (block length of 1s)   = %d\n", m);
   std::printf ("  (c) M (length of substring)  = %d\n", M);
   std::printf ("  (d) N (number of substrings) = %d\n", N);
   std::printf ("  (e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
   std::printf ("  (f) eta                      = %f\n", eta);
   std::printf ("  -----------------------------------------------\n");
   std::printf ("     F R E Q U E N C Y\n");
   std::printf ("    0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
   std::printf ("  -----------------------------------------------\n");
   std::printf ("  %3d %3d %3d %3d %3d %3d  %f ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      std::printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
   }

   std::free (sequence);
   std::printf ("%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS");
}

double Pr (int u, double eta)
{
   int l;
   double sum, p;

   if (u == 0)
   {
      p = exp (-eta);
   }
   else
   {
      sum = 0.0;
      for (l = 1; l <= u; l++)
      {
         sum += std::exp (-eta - u * std::log (2) + l * std::log (eta) - cephes_lgam (l + 1) + cephes_lgam (u) - cephes_lgam (l) - cephes_lgam (u - l + 1));
      }
      p = sum;
   }
   return p;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void NonOverlappingTemplateMatchings (int m, int n)
{
   int numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116, 2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
   /*----------------------------------------------------------------------------
   NOTE:  Should additional templates lengths beyond 21 be desired, they must
   first be constructed, saved into files and then the corresponding
   number of non periodic templates for that file be stored in the m-th
   position in the numOfTemplates variable.
   ----------------------------------------------------------------------------*/
   unsigned int bit, W_obs, nu[6], *Wj = nullptr;
   FILE *fp = nullptr;
   double sum, chi2, p_value, lambda, pi[6], varWj;
   int i, j, jj, k, match, SKIP, M, N, K = 5;
   char directory[100];
   unsigned char *sequence = nullptr;

   N = 8;
   M = n / N;

   if ((Wj = (unsigned int *) std::calloc (N, sizeof (unsigned int))) == nullptr)
   {
      std::printf (" NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
      std::printf (" Insufficient memory for required work space.\n");
      return;
   }
   lambda = (M - m + 1) / std::pow (2, m);
   varWj = M * (1.0 / std::pow (2.0, m) - (2.0 * m - 1.0) / std::pow (2.0, 2.0 * m));
   std::sprintf (directory, "templates/template%d", m);

   if (((isNegative (lambda)) || (isZero (lambda))) || ((fp = fopen (directory, "r")) == nullptr) ||
       ((sequence = (unsigned char *) std::calloc (m, sizeof (unsigned char))) == nullptr))
   {
      std::printf (" NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
      std::printf (" Lambda (%f) not being positive!\n", lambda);
      std::printf (" Template file <%s> not existing\n", directory);
      std::printf (" Insufficient memory for required work space.\n");
      if (sequence != nullptr)
      {
         std::free (sequence);
      }
   }
   else
   {
      std::printf ("    NONPERIODIC TEMPLATES TEST\n");
      std::printf ("-------------------------------------------------------------------------------------\n");
      std::printf ("    COMPUTATIONAL INFORMATION\n");
      std::printf ("-------------------------------------------------------------------------------------\n");
      std::printf (" LAMBDA = %f M = %d N = %d m = %d n = %d\n", lambda, M, N, m, n);
      std::printf ("-------------------------------------------------------------------------------------\n");
      std::printf ("  F R E Q U E N C Y\n");
      std::printf ("Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
      std::printf ("-------------------------------------------------------------------------------------\n");

      if (numOfTemplates[m] < MAXNUMOFTEMPLATES)
      {
         SKIP = 1;
      }
      else
      {
         SKIP = (int) (numOfTemplates[m] / MAXNUMOFTEMPLATES);
      }
      numOfTemplates[m] = (int) numOfTemplates[m] / SKIP;

      sum = 0.0;
      for (i = 0; i < 2; i++)
      { /* Compute Probabilities */
         pi[i] = std::exp (-lambda + i * std::log (lambda) - cephes_lgam (i + 1));
         sum += pi[i];
      }
      pi[0] = sum;
      for (i = 2; i <= K; i++)
      { /* Compute Probabilities */
         pi[i - 1] = std::exp (-lambda + i * std::log (lambda) - cephes_lgam (i + 1));
         sum += pi[i - 1];
      }
      pi[K] = 1 - sum;

      for (jj = 0; jj < MIN (MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++)
      {
         sum = 0;

         for (k = 0; k < m; k++)
         {
            std::fscanf (fp, "%d", &bit);
            sequence[k] = bit;
            std::printf ("%d", sequence[k]);
         }
         std::printf (" ");
         for (k = 0; k <= K; k++)
         {
            nu[k] = 0;
         }
         for (i = 0; i < N; i++)
         {
            W_obs = 0;
            for (j = 0; j < M - m + 1; j++)
            {
               match = 1;
               for (k = 0; k < m; k++)
               {
                  if ((int) sequence[k] != (int) epsilon[i * M + j + k])
                  {
                     match = 0;
                     break;
                  }
               }
               if (match == 1)
               {
                  W_obs++;
                  j += m - 1;
               }
            }
            Wj[i] = W_obs;
         }
         sum = 0;
         chi2 = 0.0; /* Compute Chi Square */
         for (i = 0; i < N; i++)
         {
            if (m == 10)
            {
               std::printf ("%3d  ", Wj[i]);
            }
            else
            {
               std::printf ("%4d ", Wj[i]);
            }
            chi2 += std::pow (((double) Wj[i] - lambda) / std::pow (varWj, 0.5), 2);
         }
         p_value = cephes_igamc (N / 2.0, chi2 / 2.0);

         if (isNegative (p_value) || isGreaterThanOne (p_value))
         {
            std::printf ("  WARNING:  P_VALUE IS OUT OF RANGE.\n");
         }

         std::printf ("%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
         if (SKIP > 1)
         {
            std::fseek (fp, (long) (SKIP - 1) * 2 * m, SEEK_CUR);
         }
      }
   }

   std::printf ("\n");
   if (sequence != nullptr)
   {
      std::free (sequence);
   }

   std::free (Wj);
   if (fp != nullptr)
   {
      std::fclose (fp);
   }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    B L O C K  F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void BlockFrequency (int M, int n)
{
   int i, j, N, blockSum;
   double p_value, sum, pi, v, chi_squared;

   N = n / M; /* # OF SUBSTRING BLOCKS      */
   sum = 0.0;

   for (i = 0; i < N; i++)
   {
      blockSum = 0;
      for (j = 0; j < M; j++)
      {
         blockSum += epsilon[j + i * M];
      }
      pi = (double) blockSum / (double) M;
      v = pi - 0.5;
      sum += v * v;
   }
   chi_squared = 4.0 * M * sum;
   p_value = cephes_igamc (N / 2.0, chi_squared / 2.0);

   std::printf ("   BLOCK FREQUENCY TEST\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  COMPUTATIONAL INFORMATION:\n");
   std::printf ("  ---------------------------------------------\n");
   std::printf ("  (a) Chi^2           = %f\n", chi_squared);
   std::printf ("  (b) # of substrings = %d\n", N);
   std::printf ("  (c) block length    = %d\n", M);
   std::printf ("  (d) Note: %d bits were discarded.\n", n % M);
   std::printf ("  ---------------------------------------------\n");

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R A N K  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void Rank (int n)
{
  
    int N, i, k, r;
   double p_value, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
   unsigned char **matrix = create_matrix (32, 32);

   N = n / (32 * 32);
   if (isZero (N))
   {
      std::printf ("    RANK TEST\n");
      std::printf ("  Error: Insufficient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
      p_value = 0.00;
   }
   else
   {
      r = 32; /* COMPUTE PROBABILITIES */
      product = 1;
      for (i = 0; i <= r - 1; i++)
      {
         product *= ((1.e0 - std::pow (2, i - 32)) * (1.e0 - std::pow (2, i - 32))) / (1.e0 - std::pow (2, i - r));
      }
      p_32 = std::pow (2, r * (32 + 32 - r) - 32 * 32) * product;

      r = 31;
      product = 1;
      for (i = 0; i <= r - 1; i++)
      {
         product *= ((1.e0 - std::pow (2, i - 32)) * (1.e0 - std::pow (2, i - 32))) / (1.e0 - std::pow (2, i - r));
      }
      p_31 = std::pow (2, r * (32 + 32 - r) - 32 * 32) * product;

      p_30 = 1 - (p_32 + p_31);

      F_32 = 0;
      F_31 = 0;
      for (k = 0; k < N; k++)
      { /* FOR EACH 32x32 MATRIX   */
         def_matrix (32, 32, matrix, k);
#if (DISPLAY_MATRICES == 1)
         display_matrix (32, 32, matrix);
#endif
         R = computeRank (32, 32, matrix);
         if (R == 32)
         {
            F_32++; /* DETERMINE FREQUENCIES */
         }
         if (R == 31)
         {
            F_31++;
         }
      }
      F_30 = (double) N - (F_32 + F_31);

      chi_squared = (std::pow (F_32 - N * p_32, 2) / (double) (N * p_32) + std::pow (F_31 - N * p_31, 2) / (double) (N * p_31) +
                     std::pow (F_30 - N * p_30, 2) / (double) (N * p_30));

      arg1 = -chi_squared / 2.e0;

      std::printf ("    RANK TEST\n");
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  COMPUTATIONAL INFORMATION:\n");
      std::printf ("  ---------------------------------------------\n");
      std::printf ("  (a) Probability P_%d = %f\n", 32, p_32);
      std::printf ("  (b)             P_%d = %f\n", 31, p_31);
      std::printf ("  (c)             P_%d = %f\n", 30, p_30);
      std::printf ("  (d) Frequency   F_%d = %d\n", 32, (int) F_32);
      std::printf ("  (e)             F_%d = %d\n", 31, (int) F_31);
      std::printf ("  (f)             F_%d = %d\n", 30, (int) F_30);
      std::printf ("  (g) # of matrices    = %d\n", N);
      std::printf ("  (h) Chi^2            = %f\n", chi_squared);
      std::printf ("  (i) NOTE: %d BITS WERE DISCARDED.\n", n % (32 * 32));
      std::printf ("  ---------------------------------------------\n");

      p_value = std::exp (arg1);

      if (isNegative (p_value) || isGreaterThanOne (p_value))
      {
         std::printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
      }

      for (i = 0; i < 32; i++) /* DEALLOCATE MATRIX  */
      {
         std::free (matrix[i]);
      }
      std::free (matrix);
   }

   std::printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

}
