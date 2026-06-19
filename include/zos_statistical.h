#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


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
void dradf2 (int ido, int l1, double *cc, double *ch, double *wa1);
void dradf4 (int ido, int l1, double *cc, double *ch, double *wa1, double *wa2, double *wa3);
void dradfg (int ido, int ip, int l1, int idl1, double *cc, double *c1, double *c2, double *ch, double *ch2, double *wa);
void drftf1 (int n, double *c, double *ch, double *wa, int *ifac);
void drfti1 (int n, double *wa, int *ifac);

double ApproximateEntropy (int m, int n);
double BlockFrequency (int M, int n);
std::pair<double, double> CumulativeSums (int n);
double DiscreteFourierTransform (int n);
double Frequency (int n);
double LinearComplexity (int M, int n);
double LongestRunOfOnes (int n);
std::vector<double> NonOverlappingTemplateMatchings (int m, int n);
double OverlappingTemplateMatchings (int m, int n);
std::vector<double> RandomExcursions (int n);
std::vector<double> RandomExcursionsVariant (int n);
double Rank (int n);
double Runs (int n);
std::pair<double, double> Serial (int m, int n);
double Universal (int n);

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN unsigned char *epsilon;

constexpr double ALPHA = 0.01; /* SIGNIFICANCE LEVEL */
