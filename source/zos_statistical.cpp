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

constexpr double ALPHA=0.01;            /* SIGNIFICANCE LEVEL */

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
      printf ("igamc: UNDERFLOW\n");
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
         t = fabs ((ans - r) / r);
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
      if (fabs (pk) > big)
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
      printf ("igam: UNDERFLOW\n");
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

   if (fabs (x) > 2.2)
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
   } while (fabs (term) / sum > rel_error);

   return two_sqrtpi * sum;
}

double cephes_erfc (double x)
{
   constexpr double one_sqrtpi = 0.564189583547756287;
   double a = 1, b = x, c = x, d = x * x + 0.5;
   double q1, q2 = b / d, n = 1.0, t;

   if (fabs (x) < 2.2)
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
   } while (fabs (q1 - q2) / q2 > rel_error);

   return one_sqrtpi * exp (-x * x) * q2;
}


double cephes_normal (double x)
{
   double arg=0;double result=0;
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

   if (fabs (pi - 0.5) > (2.0 / sqrt (n)))
   {
      printf ("    RUNS TEST\n");
      printf ("  ------------------------------------------\n");
      printf ("  PI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
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

      erfc_arg = fabs (V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt (2 * n));
      p_value = erfc (erfc_arg);

      printf ("    Runs Test\n");
      printf ("  ------------------------------------------\n");
      printf ("  Computational Information:\n");
      printf ("  ------------------------------------------\n");
      printf ("  (a) Pi                        = %f\n", pi);
      printf ("  (b) V_n_obs (Total # of runs) = %d\n", (int) V);
      printf ("  (c) V_n_obs - 2 n pi (1-pi)\n");
      printf ("      -----------------------   = %f\n", erfc_arg);
      printf ("        2 sqrt(2n) pi (1-pi)\n");
      printf ("  ------------------------------------------\n");
      if (isNegative (p_value) || isGreaterThanOne (p_value))
      {
         printf ("WARNING:  P_VALUE Is Out Of Range.\n");
      }

      printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
   }

   printf ("%f\n", p_value);
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

   if (((S_k = (int *) calloc (n, sizeof (int))) == nullptr) || ((cycle = (int *) calloc (MAX (1000, n / 100), sizeof (int))) == nullptr))
   {
      printf ("Random Excursions Test:  Insufficient Work Space Allocated.\n");
      if (S_k != nullptr)
      {
         free (S_k);
      }
      if (cycle != nullptr)
      {
         free (cycle);
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
            printf ("Error In Function randomExcursions:  Exceeding The Max Number Of Cycles Expected\n.");
            free (S_k);
            free (cycle);
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

   printf ("     Random Excursions Test\n");
   printf ("  --------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  --------------------------------------------\n");
   printf ("  (a) Number Of Cycles (J) = %04d\n", J);
   printf ("  (b) Sequence Length (n)  = %d\n", n);

   constraint = MAX (0.005 * pow (n, 0.5), 500);
   if (J < constraint)
   {
      printf ("  ---------------------------------------------\n");
      printf ("  Warning:  Test Not Applicable.  There Are An\n");
      printf ("     Insufficient Number Of Cycles.\n");
      printf ("  ---------------------------------------------\n");
      for (i = 0; i < 8; i++)
      {
         printf ("%f\n", 0.0);
      }
   }
   else
   {
      printf ("  (c) Rejection Constraint = %f\n", constraint);
      printf ("  -------------------------------------------\n");

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
            sum += pow (nu[k][i] - J * pi[(int) fabs (x)][k], 2) / (J * pi[(int) fabs (x)][k]);
         }
         p_value = cephes_igamc (2.5, sum / 2.0);

         if (isNegative (p_value) || isGreaterThanOne (p_value))
         {
            printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
         }

         printf ("%s  x = %2d chi^2 = %9.6f p_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
         printf ("%f\n", p_value);
      }
   }
   printf ("\n");

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
      printf ("  ---------------------------------------------\n");
      printf ("     n=%d is too short\n", n);
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

   printf ("     Longest Runs Of Ones Test\n");
   printf ("  ---------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  ---------------------------------------------\n");
   printf ("  (a) N (# of substrings)  = %d\n", N);
   printf ("  (b) M (Substring Length) = %d\n", M);
   printf ("  (c) Chi^2                = %f\n", chi2);
   printf ("  ---------------------------------------------\n");
   printf ("        F R E Q U E N C Y\n");
   printf ("  ---------------------------------------------\n");

   if (K == 3)
   {
      printf ("    <=1     2     3    >=4   P-value  Assignment");
      printf ("\n   %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
   }
   else if (K == 5)
   {
      printf ("  <=4  5  6  7  8  >=9 P-value  Assignment");
      printf ("\n   %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5]);
   }
   else
   {
      printf ("  <=10  11  12  13  14  15 >=16 P-value  Assignment");
      printf ("\n   %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], nu[6]);
   }
   if (isNegative (pval) || isGreaterThanOne (pval))
   {
      printf ("WARNING:  P_VALUE IS OUT OF RANGE.\n");
   }

   printf ("%s  p_value = %f\n\n", pval < ALPHA ? "FAILURE" : "SUCCESS", pval);
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
      printf ("  Unable to allocate working arrays for the DFT.\n");
      if (X != nullptr)
      {
         free (X);
      }
      if (wsave != nullptr)
      {
         free (wsave);
      }
      if (m != nullptr)
      {
         free (m);
      }
      return;
   }
   for (i = 0; i < n; i++)
   {
      X[i] = 2 * (int) epsilon[i] - 1;
   }

   __ogg_fdrffti (n, wsave, ifac);    /* INITIALIZE WORK ARRAYS */
   __ogg_fdrfftf (n, X, wsave, ifac); /* APPLY FORWARD FFT */

   m[0] = sqrt (X[0] * X[0]); /* COMPUTE MAGNITUDE */

   for (i = 0; i < n / 2; i++)
   {
      m[i + 1] = sqrt (pow (X[2 * i + 1], 2) + pow (X[2 * i + 2], 2));
   }
   count = 0; /* Confidence Interval */
   upperBound = sqrt (2.995732274 * n);
   for (i = 0; i < n / 2; i++)
   {
      if (m[i] < upperBound)
      {
         count++;
      }
   }
   percentile = (double) count / (n / 2) * 100;
   N_l = (double) count; /* number of peaks less than h = sqrt(3*n) */
   N_o = (double) 0.95 * n / 2.0;
   d = (N_l - N_o) / sqrt (n / 4.0 * 0.95 * 0.05);
   p_value = erfc (fabs (d) / sqrt (2.0));

   printf ("    FFT Test\n");
   printf ("  -------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  -------------------------------------------\n");
   printf ("  (a) Percentile = %f\n", percentile);
   printf ("  (b) N_l        = %f\n", N_l);
   printf ("  (c) N_o        = %f\n", N_o);
   printf ("  (d) d          = %f\n", d);
   printf ("  -------------------------------------------\n");

   printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   free (X);
   free (wsave);
   free (m);
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

   printf ("   Approximate Entropy Test\n");
   printf ("  --------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  --------------------------------------------\n");
   printf ("  (a) m (block length)    = %d\n", m);

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
         powLen = (int) pow (2, blockSize + 1) - 1;
         if ((P = (unsigned int *) calloc (powLen, sizeof (unsigned int))) == nullptr)
         {
            printf ("ApEn:  Insufficient memory available.\n");
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
         index = (int) pow (2, blockSize) - 1;
         for (i = 0; i < (int) pow (2, blockSize); i++)
         {
            if (P[index] > 0)
            {
               sum += P[index] * log (P[index] / numOfBlocks);
            }
            index++;
         }
         sum /= numOfBlocks;
         ApEn[r] = sum;
         r++;
         free (P);
      }
   }
   apen = ApEn[0] - ApEn[1];

   chi_squared = 2.0 * seqLength * (log (2) - apen);
   p_value = cephes_igamc (pow (2, m - 1), chi_squared / 2.0);

   printf ("  (b) n (sequence length) = %d\n", seqLength);
   printf ("  (c) Chi^2               = %f\n", chi_squared);
   printf ("  (d) Phi(m)	       = %f\n", ApEn[0]);
   printf ("  (e) Phi(m+1)	       = %f\n", ApEn[1]);
   printf ("  (f) ApEn                = %f\n", apen);
   printf ("  (g) Log(2)              = %f\n", log (2.0));
   printf ("  --------------------------------------------\n");

   if (m > (int) (log (seqLength) / log (2) - 5))
   {
      printf ("  Note: The blockSize = %d exceeds recommended value of %d\n", m, MAX (1, (int) (log (seqLength) / log (2) - 5)));
      printf ("  Results are inaccurate!\n");
      printf ("  --------------------------------------------\n");
   }

   printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
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
   s_obs = fabs (sum) / sqrt (n);
   f = s_obs / sqrt2;
   p_value = erfc (f);

   printf ("         Frequency Test\n");
   printf ("  ---------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  ---------------------------------------------\n");
   printf ("  (a) The nth partial sum = %d\n", (int) sum);
   printf ("  (b) S_n/n               = %f\n", sum / n);
   printf ("  ---------------------------------------------\n");

   printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
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
   p_value1 = cephes_igamc (pow (2, m - 1) / 2, del1 / 2.0);
   p_value2 = cephes_igamc (pow (2, m - 2) / 2, del2 / 2.0);

   printf ("          Serial Test\n");
   printf ("  ---------------------------------------------\n");
   printf ("   Computational Information:		  \n");
   printf ("  ---------------------------------------------\n");
   printf ("  (a) Block length    (m) = %d\n", m);
   printf ("  (b) Sequence length (n) = %d\n", n);
   printf ("  (c) Psi_m               = %f\n", psim0);
   printf ("  (d) Psi_m-1             = %f\n", psim1);
   printf ("  (e) Psi_m-2             = %f\n", psim2);
   printf ("  (f) Del_1               = %f\n", del1);
   printf ("  (g) Del_2               = %f\n", del2);
   printf ("  ---------------------------------------------\n");

   printf ("%s  p_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
   printf ("%s  p_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
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
      sum1 += cephes_normal (((4 * k + 1) * z) / sqrt (n));
      sum1 -= cephes_normal (((4 * k - 1) * z) / sqrt (n));
   }
   sum2 = 0.0;
   for (k = (-n / z - 3) / 4; k <= (n / z - 1) / 4; k++)
   {
      sum2 += cephes_normal (((4 * k + 3) * z) / sqrt (n));
      sum2 -= cephes_normal (((4 * k + 1) * z) / sqrt (n));
   }

   p_value = 1.0 - sum1 + sum2;

   printf ("        Cumulative Sums (Forward) Test\n");
   printf ("  -------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  -------------------------------------------\n");
   printf ("  (a) The maximum partial sum = %d\n", z);
   printf ("  -------------------------------------------\n");

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      printf ("  Warning:  P_Value Is Out Of Range\n");
   }

   printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   // backwards
   sum1 = 0.0;
   for (k = (-n / zrev + 1) / 4; k <= (n / zrev - 1) / 4; k++)
   {
      sum1 += cephes_normal (((4 * k + 1) * zrev) / sqrt (n));
      sum1 -= cephes_normal (((4 * k - 1) * zrev) / sqrt (n));
   }
   sum2 = 0.0;
   for (k = (-n / zrev - 3) / 4; k <= (n / zrev - 1) / 4; k++)
   {
      sum2 += cephes_normal (((4 * k + 3) * zrev) / sqrt (n));
      sum2 -= cephes_normal (((4 * k + 1) * zrev) / sqrt (n));
   }
   p_value = 1.0 - sum1 + sum2;

   printf ("        Cumulative Sums (Reverse) Test\n");
   printf ("  -------------------------------------------\n");
   printf ("  Computational Information:\n");
   printf ("  -------------------------------------------\n");
   printf ("  (a) The maximum partial sum = %d\n", zrev);
   printf ("  -------------------------------------------\n");

   if (isNegative (p_value) || isGreaterThanOne (p_value))
   {
      printf ("  WARNING:  P_VALUE IS OUT OF RANGE\n");
   }

   printf ("%s  p_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);

   printf ("%f\n", p_value);
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
   powLen = (int) pow (2, m + 1) - 1;
   if ((P = (unsigned int *) calloc (powLen, sizeof (unsigned int))) == nullptr)
   {
      printf ("Serial Test:  Insufficient memory available.\n");
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
   for (i = (int) pow (2, m) - 1; i < (int) pow (2, m + 1) - 1; i++)
   {
      sum += pow (P[i], 2);
   }
   sum = (sum * pow (2, m) / (double) n) - (double) n;
   free (P);

   return sum;
}
