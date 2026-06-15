#pragma once

double cephes_erf (double x);
double cephes_erfc (double x);
double cephes_igam (double a, double x);
double cephes_igamc (double a, double x);
double cephes_lgam (double x);
double cephes_normal (double x);
double cephes_p1evl (double x, double *coef, int N);
double cephes_polevl (double x, double *coef, int N);
double Pr (int u, double eta);
double psi2 (int m, int n);
void __ogg_fdrfftf (int n, double *X, double *wsave, int *ifac);
void __ogg_fdrffti (int n, double *wsave, int *ifac);
void ApproximateEntropy (int m, int n);
void BlockFrequency (int M, int n);
void CumulativeSums (int n);
void DiscreteFourierTransform (int n);
void dradf2 (int ido, int l1, double *cc, double *ch, double *wa1);
void dradf4 (int ido, int l1, double *cc, double *ch, double *wa1, double *wa2, double *wa3);
void dradfg (int ido, int ip, int l1, int idl1, double *cc, double *c1, double *c2, double *ch, double *ch2, double *wa);
void drftf1 (int n, double *c, double *ch, double *wa, int *ifac);
void drfti1 (int n, double *wa, int *ifac);
void Frequency (int n);
void LongestRunOfOnes (int n);
void NonOverlappingTemplateMatchings (int m, int n);
void OverlappingTemplateMatchings (int m, int n);
void RandomExcursions (int n);
void Rank (int n);
void Runs (int n);
void Serial (int m, int n);
void Universal (int n);
unsigned char *epsilon;

