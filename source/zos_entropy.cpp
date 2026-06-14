#include <cstdio>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "zos_entropy.h"
#include "zos_statistical.h"
#include "zos_generators.h"

enum class GeneratorType
{
   TRNO,
   JITTER,
   DEVURANDOM
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
         jitter_fill (out, size);
         break;
      }
      case GeneratorType::JITTER:
      {
         prno_trng_generate (out, size);
         break;
      }
      case GeneratorType::DEVURANDOM:
      {
         dev_urandom_generate (out, size);
         break;
      }
   }
   return 0;
}


//unsigned char *epsilon;


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
   const int sample_size = 4096;
   unsigned char random_data[4096];

constexpr std::array<GeneratorType, 3> allGenerators = {
    GeneratorType::TRNO, GeneratorType::JITTER, GeneratorType::DEVURANDOM
};
constexpr std::array<const char*, 3> generatorNames = {
    "PRNO-TRNG", "CPU jitter", "/dev/urandom"
};
 epsilon = (unsigned char *) malloc (sample_size * 8);

for (std::size_t i = 0; i < allGenerators.size(); ++i)
{

   rc = zos_getentropy (random_data, sample_size, allGenerators[i]);

   if (epsilon == nullptr)
   {
      perror ("malloc");
      return 1;
   }

   for (int i = 0; i < sample_size; i++)
   {
      for (int b = 0; b < 8; b++)
      {
         epsilon[i * 8 + b] = (random_data[i] >> (7 - b)) & 1;
      }
   }
   std::cout << generatorNames[i] << std::endl << "Sample Size in bytes: " << sample_size << std::endl << std::endl;

   LongestRunOfOnes (sample_size * 8);
   Runs (sample_size * 8);
   CumulativeSums (sample_size * 8);
   ApproximateEntropy (6, sample_size * 8);
   Frequency (sample_size * 8);
   Serial (4, 2048);                  /* very conservative */
   Serial (5, 2048);                  /* good */
   Serial (6, 2048);                  /* good choice */
   Serial (7, 2048);                  /* still okay, but thinner counts */
   Serial (6, sample_size * 8);       // m = 6 means the test counts all overlapping 6-bit patterns:
   DiscreteFourierTransform (sample_size * 8);
}
   std::free (epsilon);

   return 0;
}
