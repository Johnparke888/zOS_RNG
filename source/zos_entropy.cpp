#include <cstdio>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "zos_entropy.h"
#include "zos_statistical.h"
#include "zos_generators.h"
#include "matrix.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

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
   /*
    * Serial uses all overlapping m-bit patterns.
    * Keep average observations per pattern reasonable.
    */
   int m = floor_log2_int (n) - 5;
   return std::clamp (m, 2, 16);
}

static int choose_approximate_entropy_m (int n)
{
   /*
    * Your ApproximateEntropy() function itself warns when:
    *
    *     m > log2(n) - 5
    *
    * Cap at 10 to keep the table size modest.
    */
   int m = floor_log2_int (n) - 5;
   return std::clamp (m, 2, 10);
}

static int choose_block_frequency_M (int n)
{
   /*
    * M is the block size, not a template length.
    * Prefer 128 if it gives at least 100 blocks.
    */
   if (n / 128 >= 100)
   {
      return 128;
   }

   int M = n / 100;
   return std::clamp (M, 20, 128);
}

static bool universal_is_applicable (int n)
{
   return n >= 387840;
}

static bool random_excursions_likely_applicable (int n)
{
   /*
    * The actual applicability depends on the generated walk's cycle count J.
    * This avoids running it on obviously too-small samples.
    */
   return n >= 387840;
}

static bool rank_is_reasonable (int n)
{
   /*
    * Rank uses 32x32 matrices, 1024 bits each.
    * This gives at least 38 matrices.
    */
   return n >= 38 * 32 * 32;
}
enum class GeneratorType
{
   TRNO,
   JITTER,
   DEVURANDOM,
   NAIVE
};


/* ------------------------------------------------------------------------ *
 * Public entry point.
 *
 *
 * ------------------------------------------------------------------------ */
int zos_getentropy (void *output_buffer_ptr, size_t size, const GeneratorType &generatorType)
{
   if (output_buffer_ptr == nullptr)
   {
      errno = EINVAL;
      return -1;
   }

   unsigned char *out = (unsigned char *) output_buffer_ptr;

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
   }
   return 0;
}


// unsigned char *epsilon;


static void hexdump (const unsigned char *p, size_t n)
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


int main ()
{

   int rc = 0;

   constexpr std::array<GeneratorType, 4> allGenerators = {
       GeneratorType::TRNO, GeneratorType::JITTER, GeneratorType::DEVURANDOM, GeneratorType::NAIVE};

   constexpr std::array<const char *, 4> generatorNames = {"PRNO-TRNG", "CPU jitter", "/dev/urandom", "Naive PRNG"};

   constexpr std::size_t sample_size = 128 * 1024;
   std::vector<unsigned char> random_data (sample_size);

   epsilon = static_cast<unsigned char *> (std::malloc (sample_size * 8));
   
   if (epsilon == nullptr)
   {
      perror ("malloc");
      return 1;
   }

   for (std::size_t g = 0; g < allGenerators.size (); ++g)
   {
      rc = zos_getentropy (random_data.data (), random_data.size (), allGenerators[g]);
      if (rc != 0)
      {
         perror ("zos_getentropy");
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

      std::cout << "\n============================================================\n";
      std::cout << generatorNames[g] << '\n';
      std::cout << "Sample size in bytes : " << sample_size << '\n';
      std::cout << "Sample size in bits  : " << n << '\n';
      std::cout << "BlockFrequency M     : " << blockFrequencyM << '\n';
      std::cout << "Serial m             : " << serialM << '\n';
      std::cout << "ApproxEntropy m      : " << approximateM << '\n';
      std::cout << "Template m           : " << templateM << '\n';
      std::cout << "============================================================\n\n";

      Frequency (n);
      Runs (n);
      CumulativeSums (n);
      LongestRunOfOnes (n);
      DiscreteFourierTransform (n);

      if (rank_is_reasonable (n))
      {
         Rank (n);
      }
      else
      {
         std::cout << "Skipping Rank: sample is too small for a useful 32x32 matrix count.\n\n";
      }

      if (universal_is_applicable (n))
      {
         Universal (n);
      }
      else
      {
         std::cout << "Skipping Universal: requires at least 387840 bits.\n\n";
      }

      if (random_excursions_likely_applicable (n))
      {
         RandomExcursions (n);
         RandomExcursionsVariant (n);
      }
      else
      {
         std::cout << "Skipping Random Excursions: sample is probably too small.\n\n";
      }

      ApproximateEntropy (approximateM, n);
      Serial (serialM, n);

      OverlappingTemplateMatchings (templateM, n);
      NonOverlappingTemplateMatchings (templateM, n);

      BlockFrequency (blockFrequencyM, n);
      std::cout << "============================================================\n\n";
   }

   std::free (epsilon);
 
  
   return 0;
}
