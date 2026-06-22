#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <vector>

#include "Assess.h"
#include "matrix.h"
#include "UniformityResult.h"
#include "zos_entropy.h"
#include "zos_generators.h"
#define EXTERN
#include "zos_statistical.h"

static int floor_log2_int (int n)
{
   int r = 0;
   while (n > 1)
   {
      n >>= 1;
      ++r;
   }
   return r;
}

static int choose_serial_m (int n)
{
   /* Serial uses all overlapping m-bit patterns. */
   int m = floor_log2_int (n) - 5;
   return std::clamp (m, 2, 16);
}

static int choose_approximate_entropy_m (int n)
{
   /* ApproximateEntropy warns when m > log2(n) - 5. */
   int m = floor_log2_int (n) - 5;
   return std::clamp (m, 2, 10);
}

static int choose_block_frequency_M (int n)
{
   int M = n / 100;
   return std::max (20, M);
}

static bool universal_is_applicable (int n)
{
   return n >= 387840;
}

static bool random_excursions_likely_applicable (int n)
{
   return n >= 387840;
}

static bool rank_is_reasonable (int n)
{
   return n >= 38 * 32 * 32;
}

static int choose_linear_complexity_M (int n)
{
   if (n >= 500 * 100)
   {
      return 500;
   }

   int M = n / 100;
   return std::clamp (M, 100, 500);
}

static bool linear_complexity_is_reasonable (int n, int M)
{
   return M > 0 && (n / M) >= 100;
}

enum class GeneratorType
{
   TRNO,
   JITTER,
   DEVURANDOM,
   NAIVE,
   BADRAWCLOCK,
   BADHASHCOUNTER
};

/* ------------------------------------------------------------------------ *
 * Public entry point.
 * ------------------------------------------------------------------------ */
int zos_getentropy (void *output_buffer_ptr, size_t size, const GeneratorType &generatorType)
{
   if (output_buffer_ptr == nullptr)
   {
      errno = EINVAL;
      return -1;
   }

   unsigned char *out = static_cast<unsigned char *> (output_buffer_ptr);

   switch (generatorType)
   {
      case GeneratorType::TRNO:
      {
         prno_trng_generate (out, size);
         break;
      }
      case GeneratorType::JITTER:
      {
         jitter_fill (out, size);
         break;
      }
      case GeneratorType::DEVURANDOM:
      {
         dev_urandom_generate (out, size);
         break;
      }
      case GeneratorType::NAIVE:
      {
         naive_prng_generate (out, size);
         break;
      }
      case GeneratorType::BADRAWCLOCK:
      {
         bad_raw_clock_generate (out, size);
         break;
      }
      case GeneratorType::BADHASHCOUNTER:
      {
         bad_hash_counter_generate (out, size);
         break;
      }
   }

   return 0;
}

static std::string format_elapsed_time (std::chrono::steady_clock::duration elapsed)
{
   using namespace std::chrono;

   const auto total_ms = duration_cast<milliseconds> (elapsed).count ();

   const auto hours = total_ms / (1000 * 60 * 60);
   const auto minutes = (total_ms / (1000 * 60)) % 60;
   const auto seconds = (total_ms / 1000) % 60;
   const auto millis = total_ms % 1000;

   std::ostringstream os;

   if (hours > 0)
   {
      os << hours << "h ";
   }

   if (hours > 0 || minutes > 0)
   {
      os << minutes << "m ";
   }

   os << seconds << "." << std::setfill ('0') << std::setw (3) << millis << "s";

   return os.str ();
}

[[maybe_unused]] static void hexdump (const unsigned char *p, size_t n)
{
   for (size_t i = 0; i < n; ++i)
   {
      printf ("%02x", p[i]);
      printf (((i & 31) == 31) ? "\n" : " ");
   }
   if (n & 31)
   {
      printf ("\n");
   }
}

static constexpr int randomExcursionStates[8] = {-4, -3, -2, -1, 1, 2, 3, 4};

static constexpr int randomExcursionVariantStates[18] = {-9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9};

int main ()
{
   int rc = 0;
   int numberOfRuns = 10;

   constexpr std::array<GeneratorType, 6> allGenerators = {GeneratorType::TRNO,
                                                           GeneratorType::JITTER,
                                                           GeneratorType::DEVURANDOM,
                                                           GeneratorType::NAIVE,
                                                           GeneratorType::BADRAWCLOCK,
                                                           GeneratorType::BADHASHCOUNTER};

   constexpr std::array<const char *, 6> generatorNames = {
       "PRNO-TRNG", "CPU jitter", "/dev/urandom", "Naive PRNG", "Bad Raw Clock", "Bad Hash Counter"};

   constexpr std::size_t sample_size = 16 * 1024;
   std::vector<unsigned char> random_data (sample_size);

   epsilon = static_cast<unsigned char *> (std::malloc (sample_size * 8));

   if (epsilon == nullptr)
   {
      perror ("malloc");
      return 1;
   }

   std::cout << "Start" << std::endl;
   std::map<std::string, std::vector<double>> pvalues;

   for (std::size_t generator = 0; generator < allGenerators.size (); ++generator)
   {
      pvalues.clear ();
      bool haveReported = false;
      const auto generatorStart = std::chrono::steady_clock::now ();

      for (int run = 0; run < numberOfRuns; ++run)
      {
         const auto runStart = std::chrono::steady_clock::now ();
         rc = zos_getentropy (random_data.data (), sample_size, allGenerators[generator]);
         if (rc != 0)
         {
            std::cerr << "zos_getentropy failed for generator " << generatorNames[generator] << " on run " << (run + 1) << "\n";
            continue;
         }

         for (std::size_t i = 0; i < sample_size; ++i)
         {
            for (int b = 0; b < 8; ++b)
            {
               epsilon[i * 8 + b] = (random_data[i] >> (7 - b)) & 1;
            }
         }

         const int n = static_cast<int> (sample_size * 8);

         const int blockFrequencyM = choose_block_frequency_M (n);
         const int serialM = choose_serial_m (n);
         const int approximateM = choose_approximate_entropy_m (n);
         const int templateM = 9;
         const int linearComplexityM = choose_linear_complexity_M (n);

         if (!haveReported)
         {
            std::cout << "\n============================================================\n";
            std::cout << "Generator: " << generatorNames[generator] << std::endl;
            haveReported = true;
            std::cout << "Run number           : " << run + 1 << " of " << numberOfRuns << '\n';
            std::cout << "Sample size in bytes : " << sample_size << '\n';
            std::cout << "Sample size in bits  : " << n << '\n';
            std::cout << "BlockFrequency M     : " << blockFrequencyM << '\n';
            std::cout << "Serial m             : " << serialM << '\n';
            std::cout << "ApproxEntropy m      : " << approximateM << '\n';
            std::cout << "Template m           : " << templateM << '\n';
            std::cout << "LinearComplexity M   : " << linearComplexityM << '\n';
            std::cout << "============================================================\n\n";
         }
         auto approximateEntropyResult = ApproximateEntropy (approximateM, n);
         if (EvaluateApproximateEntropyResult (approximateEntropyResult))
         {
             pvalues["ApproximateEntropy"].push_back (approximateEntropyResult);
         }
         pvalues["Frequency"].push_back (Frequency (n));

         pvalues["Runs"].push_back (Runs (n));

         if (rank_is_reasonable (n))
         {
            pvalues["Rank"].push_back (Rank (n));
         }

         auto cusum = CumulativeSums (n);
         pvalues["CumulativeSums.Forward"].push_back (cusum.first);
         pvalues["CumulativeSums.Reverse"].push_back (cusum.second);
         pvalues["DiscreteFourierTransform"].push_back (DiscreteFourierTransform (n));
         pvalues["LongestRunOfOnes"].push_back (LongestRunOfOnes (n));

         if (universal_is_applicable (n))
         {
            pvalues["Universal"].push_back (Universal (n));
         }
         if (CanRunSerial (serialM, n))
         {
            auto serial = Serial (serialM, n);
            auto serialResults = EvaluateSerialResults (serial.first, serial.second);
            if (serialResults.first)
            {
               pvalues["Serial.m=" + std::to_string (serialM)].push_back (serial.first);
            }
            if (serialResults.second)
            {
               pvalues["Serial.2"].push_back (serial.second);
            }
         }

         auto templates = NonOverlappingTemplateMatchings (templateM, n);
         for (std::size_t i = 0; i < templates.size (); ++i)
         {
            pvalues["NonOverlappingTemplate." + std::to_string (i)].push_back (templates[i]);
         }

         if (random_excursions_likely_applicable (n))
         {
            auto excursions = RandomExcursions (n);
            for (std::size_t i = 0; i < excursions.size () && i < std::size (randomExcursionStates); ++i)
            {
               if (excursions[i] > 0.0)
               {
                  pvalues["RandomExcursions.x=" + std::to_string (randomExcursionStates[i])].push_back (excursions[i]);
               }
            }

            auto excursionsVariant = RandomExcursionsVariant (n);
            for (std::size_t i = 0; i < excursionsVariant.size () && i < std::size (randomExcursionVariantStates); ++i)
            {
               if (excursionsVariant[i] > 0.0)
               {
                  pvalues["RandomExcursionsVariant.x=" + std::to_string (randomExcursionVariantStates[i])].push_back (excursionsVariant[i]);
               }
            }
         }
         if (CanRunBlockFrequencyTest (n, blockFrequencyM))
         {
            pvalues["BlockFrequency"].push_back (BlockFrequency (blockFrequencyM, n));
         }
         pvalues["OverlappingTemplate"].push_back (OverlappingTemplateMatchings (templateM, n));

         if (linear_complexity_is_reasonable (n, linearComplexityM))
         {
            pvalues["LinearComplexity"].push_back (LinearComplexity (linearComplexityM, n));
         }

         const auto runEnd = std::chrono::steady_clock::now ();
         const auto runElapsed = runEnd - runStart;

         std::cout << "Completed run " << run + 1 << " of " << numberOfRuns << " for generator: " << generatorNames[generator] << " in "
                   << format_elapsed_time (runElapsed) << std::endl;
      }

      for (const auto &[testName, values] : pvalues)
      {
         Summarize_pvalue_Series (testName, values, ALPHA);
      }

      Evaluate_Results (generatorNames[generator], pvalues, ALPHA);

      const auto generatorEnd = std::chrono::steady_clock::now ();
      std::cout << "Completed generator: " << generatorNames[generator] << " in " << format_elapsed_time (generatorEnd - generatorStart) << "\n";
   }

   std::free (epsilon);
   epsilon = nullptr;

   return 0;
}
