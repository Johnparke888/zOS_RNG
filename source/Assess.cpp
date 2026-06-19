#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Assess.h"
#include "zos_statistical.h"

template <typename T> bool compareInsensitive (const std::string &str1, const T &str2)
{
   std::string s2 (str2);

   if (str1.length () != s2.length ())
   {
      return false;
   }

   return std::equal (
       str1.begin (), str1.end (), s2.begin (), [] (unsigned char a, unsigned char b) { return std::tolower (a) == std::tolower (b); });
}

struct PassCountRange
{
   int minPass = 0;
   int maxPass = 0;
   double minRate = 0.0;
   double maxRate = 0.0;

   std::string toString () const
   {
      std::ostringstream os;
      os << minPass << ".." << maxPass;
      return os.str ();
   }
};

static PassCountRange calculate_acceptable_pass_count_range (int sampleSize, double alpha)
{
   PassCountRange range{};

   if (sampleSize <= 0)
   {
      return range;
   }

   const double p_hat = 1.0 - alpha;
   const double margin = 3.0 * std::sqrt ((p_hat * alpha) / static_cast<double> (sampleSize));

   range.minRate = std::max (0.0, p_hat - margin);
   range.maxRate = std::min (1.0, p_hat + margin);

   /*
    * Match the NIST STS style: the floating-point threshold is stored
    * in an int, effectively truncating toward zero.  For alpha=0.01
    * and sampleSize=100, the lower threshold is about 96.015, so
    * minPass becomes 96.
    */
   range.minPass = static_cast<int> (range.minRate * static_cast<double> (sampleSize));
   range.maxPass = static_cast<int> (range.maxRate * static_cast<double> (sampleSize));

   range.minPass = std::clamp (range.minPass, 0, sampleSize);
   range.maxPass = std::clamp (range.maxPass, 0, sampleSize);

   return range;
}

struct SeriesAssessment
{
   std::string testName;

   int sampleSize = 0;
   int passed = 0;
   int failed = 0;

   int minPass = 0;
   int maxPass = 0;

   double failurePercent = 0.0;
   double uniformityPValue = 0.0;

   bool enoughSamples = false;
   bool proportionOK = false;
   bool uniformityOK = false;
   bool seriesOK = false;
};

static SeriesAssessment assess_one_series (const std::string &testName, const std::vector<double> &values, double ALPHA)
{
   SeriesAssessment a;
   a.testName = testName;
   a.sampleSize = static_cast<int> (values.size ());

   if (a.sampleSize == 0)
   {
      return a;
   }

   a.failed = static_cast<int> (std::count_if (values.begin (), values.end (), [ALPHA] (double p) { return p < ALPHA; }));

   a.passed = a.sampleSize - a.failed;
   a.failurePercent = 100.0 * static_cast<double> (a.failed) / static_cast<double> (a.sampleSize);

   const PassCountRange passRange = calculate_acceptable_pass_count_range (a.sampleSize, ALPHA);
   a.minPass = passRange.minPass;
   a.maxPass = passRange.maxPass;
   a.proportionOK = (a.passed >= a.minPass && a.passed <= a.maxPass);

   /*
    * Uniformity needs enough p-values to be meaningful. Ten is the bare
    * minimum because there are ten histogram bins; 100+ is much better.
    */
   a.enoughSamples = (a.sampleSize >= 10);

   if (a.enoughSamples)
   {
      UniformityResult u = calculate_uniformity_p_value (values);
      a.uniformityPValue = u.uniformity_p_value;

      /* NIST STS marks p-value uniformity below 0.0001 as suspicious. */
      a.uniformityOK = (a.uniformityPValue >= 0.0001);
   }
   else
   {
      a.uniformityPValue = 0.0;
      a.uniformityOK = false;
   }

   a.seriesOK = a.proportionOK && a.uniformityOK;

   return a;
}

void Evaluate_Results (const std::string &generatorName, const std::map<std::string, std::vector<double>> &pvalues, double ALPHA)
{
   std::vector<SeriesAssessment> assessments;

   int totalSeries = 0;
   int assessedSeries = 0;
   int passedSeries = 0;
   int proportionFailures = 0;
   int uniformityFailures = 0;
   int bothFailures = 0;
   int insufficientSeries = 0;

   int totalPValues = 0;
   int totalPassedPValues = 0;
   int totalFailedPValues = 0;

   for (const auto &[testName, values] : pvalues)
   {
      ++totalSeries;

      SeriesAssessment a = assess_one_series (testName, values, ALPHA);
      assessments.push_back (a);

      if (a.sampleSize == 0)
      {
         ++insufficientSeries;
         continue;
      }

      ++assessedSeries;

      totalPValues += a.sampleSize;
      totalPassedPValues += a.passed;
      totalFailedPValues += a.failed;

      if (!a.enoughSamples)
      {
         ++insufficientSeries;
         continue;
      }

      if (a.proportionOK && a.uniformityOK)
      {
         ++passedSeries;
      }

      if (!a.proportionOK)
      {
         ++proportionFailures;
      }

      if (!a.uniformityOK)
      {
         ++uniformityFailures;
      }

      if (!a.proportionOK && !a.uniformityOK)
      {
         ++bothFailures;
      }
   }

   const int badSeries = proportionFailures + uniformityFailures - bothFailures;

   double seriesScore = 0.0;
   if (assessedSeries > 0)
   {
      seriesScore = 100.0 * static_cast<double> (passedSeries) / static_cast<double> (assessedSeries);
   }

   double totalFailureRate = 0.0;
   if (totalPValues > 0)
   {
      totalFailureRate = 100.0 * static_cast<double> (totalFailedPValues) / static_cast<double> (totalPValues);
   }

   const double expectedFailedPValues = static_cast<double> (totalPValues) * ALPHA;

   std::string verdict;

   /*
    * Uniformity failures are more concerning than a small number of isolated
    * proportion failures among many subtests.  With 148 non-overlapping
    * template subtests, a few borderline proportion warnings can occur by
    * chance even when the generator is healthy.
    */
   const int minorWarningLimit = std::max (2, assessedSeries / 50);       // about 2%
   const int suspiciousLimit = std::max (4, assessedSeries / 20);         // about 5%

   if (assessedSeries == 0)
   {
      verdict = "INCONCLUSIVE";
   }
   else if (badSeries == 0 && insufficientSeries == 0)
   {
      verdict = "GOOD statistical profile";
   }
   else if (uniformityFailures == 0 && badSeries <= minorWarningLimit)
   {
      verdict = "ACCEPTABLE with minor warnings";
   }
   else if (badSeries <= suspiciousLimit)
   {
      verdict = "SUSPICIOUS - investigate failed test series";
   }
   else
   {
      verdict = "POOR statistical profile";
   }

   std::cout << "\n\n";
   std::cout << "============================================================\n";
   std::cout << "RNG QUALITY SUMMARY\n";
   std::cout << "============================================================\n";
   std::cout << "Generator                 : " << generatorName << "\n";
   std::cout << "Alpha                     : " << ALPHA << "\n";
   std::cout << "Test/subtest series        : " << totalSeries << "\n";
   std::cout << "Assessed series            : " << assessedSeries << "\n";
   std::cout << "Series passed              : " << passedSeries << "\n";
   std::cout << "Series score               : " << std::fixed << std::setprecision (2) << seriesScore << "%\n";
   std::cout << "Total p-values             : " << totalPValues << "\n";
   std::cout << "Total passed p-values      : " << totalPassedPValues << "\n";
   std::cout << "Total failed p-values      : " << totalFailedPValues << "\n";
   std::cout << "Observed p-value fail rate : " << std::fixed << std::setprecision (3) << totalFailureRate << "%\n";
   std::cout << "Expected failed p-values   : " << std::fixed << std::setprecision (2) << expectedFailedPValues << "\n";
   std::cout << "Proportion-check failures  : " << proportionFailures << "\n";
   std::cout << "Uniformity-check failures  : " << uniformityFailures << "\n";
   std::cout << "Both-check failures        : " << bothFailures << "\n";
   std::cout << "Insufficient/empty series   : " << insufficientSeries << "\n";
   std::cout << "Verdict                   : " << verdict << "\n";
   std::cout << "============================================================\n\n";

   std::cout << std::left << std::setw (38) << "Test" << std::right << std::setw (8) << "N" << std::setw (8) << "Pass" << std::setw (8) << "Fail"
             << std::setw (14) << "PassRange" << std::setw (14) << "Uniformity" << std::setw (14) << "Result" << "\n";

   std::cout << std::string (104, '-') << "\n";

   for (const SeriesAssessment &a : assessments)
   {
      std::string result;

      if (a.sampleSize == 0)
      {
         result = "NO DATA";
      }
      else if (!a.enoughSamples)
      {
         result = "LOW N";
      }
      else if (a.proportionOK && a.uniformityOK)
      {
         result = "PASS";
      }
      else if (!a.proportionOK && !a.uniformityOK)
      {
         result = "FAIL BOTH";
      }
      else if (!a.proportionOK)
      {
         result = "FAIL PROP";
      }
      else
      {
         result = "FAIL UNIF";
      }

      std::ostringstream range;
      range << a.minPass << ".." << a.maxPass;

      std::cout << std::left << std::setw (38) << a.testName << std::right << std::setw (8) << a.sampleSize << std::setw (8) << a.passed
                << std::setw (8) << a.failed << std::setw (14) << range.str () << std::setw (14) << std::fixed << std::setprecision (6)
                << a.uniformityPValue << std::setw (14) << result << "\n";
   }

   std::cout << "\nNotes:\n";
   std::cout << "  PASS means both the pass-count proportion and p-value uniformity checks passed.\n";
   std::cout << "  FAIL PROP means too many or too few sequences passed for that test series.\n";
   std::cout << "  FAIL UNIF means the p-values were not uniformly distributed across 0..1.\n";
   std::cout << "  Random Excursions p-values equal to 0.0 should be omitted before this function is called.\n";
   std::cout << "  This is a statistical quality summary, not a proof of cryptographic security.\n";
}

void Summarize_pvalue_Series (const std::string &testName, const std::vector<double> &values, double ALPHA)
{
   const int n = static_cast<int> (values.size ());

   std::cout << "\n\n";
   std::cout << "Test: " << testName << "\n";
   std::cout << "Alpha: " << ALPHA << "\n";
   std::cout << "Total p-values: " << n << "\n";

   if (n == 0)
   {
      std::cout << "No valid p-values for this test series.\n";
      return;
   }

   const int total_failed = static_cast<int> (std::count_if (values.begin (), values.end (), [ALPHA] (double p) { return p < ALPHA; }));
   const int total_passed = n - total_failed;
   const double proportion_below_alpha = static_cast<double> (total_failed) / static_cast<double> (n);
   const PassCountRange passRange = calculate_acceptable_pass_count_range (n, ALPHA);
   const bool proportionOK = (total_passed >= passRange.minPass && total_passed <= passRange.maxPass);

   std::cout << "Total passed: " << total_passed << "\n";
   std::cout << "Total failed: " << total_failed << "\n";
   std::cout << "Proportion failed: " << proportion_below_alpha * 100.0 << "%\n";
   std::cout << "Acceptable pass count: " << passRange.toString () << "\n";
   std::cout << "Proportion result: " << (proportionOK ? "PASS" : "FAIL") << "\n";

   if (n >= 10)
   {
      UniformityResult result = calculate_uniformity_p_value (values);
      std::cout << "Uniformity Results\n";
      std::cout << result.toString () << std::endl;
   }
   else
   {
      std::cout << "Uniformity Results\n";
      std::cout << "Not calculated: fewer than 10 p-values.\n";
   }
}

/*
 * uniformity_p_value >= 0.0001  => acceptable distribution
 * uniformity_p_value <  0.0001  => suspicious distribution
 */
UniformityResult calculate_uniformity_p_value (const std::vector<double> &pvalues)
{
   if (pvalues.empty ())
   {
      throw std::invalid_argument ("No p-values supplied");
   }

   UniformityResult result{};

   const int sampleSize = static_cast<int> (pvalues.size ());
   const double expected = static_cast<double> (sampleSize) / 10.0;

   for (double p : pvalues)
   {
      if (p < 0.0 || p > 1.0)
      {
        if (p > 1.00)
        {
          std::cout << "warning pvalue of " << p << std::endl;
          p = 1.00;
        }
        if (p < 0.00)
        {
          std::cout << "warning pvalue of " << p << std::endl;
          p = 0.00;
        }        // throw std::invalid_argument ("p-value out of range");
      }

      int bin = static_cast<int> (std::floor (p * 10.0));
      if (bin == 10)
      {
         bin = 9;       // p == 1.0 belongs in the last bin
      }

      result.bins[bin]++;
   }

   for (int i = 0; i < 10; ++i)
   {
      const double diff = static_cast<double> (result.bins[i]) - expected;
      result.chi_square += (diff * diff) / expected;
   }

   result.uniformity_p_value = cephes_igamc (9.0 / 2.0, result.chi_square / 2.0);

   return result;
}
