#include "zos_entropy.h"
#include "zos_statistical.h"
#include "zos_generators.h"
/*******************************************************************************
 * Sources of Information and Code
 *
 * 1. z/Architecture Principles of Operation
 *    Fourteenth Edition (May, 2022)
 *    (C) Copyright International Business Machines Corporation 2000, 2022.
 *    All rights reserved.
 *
 * 2. National Institute of Standards and Technology
 *    Special Publication 800-22 revision 1a
 *    Natl. Inst. Stand. Technol. Spec. Publ. 800-22rev1a, 131 pages
 *    (April 2010)
 *    "A Statistical Test Suite for Random and Pseudorandom Number Generators
 *    for Cryptographic Applications"
 *
 *    This is free and unencumbered software released into the public domain.
 *
 *    The following applies to the NIST statistical code:
 *
 *    This code was accessed on 2019-01-05 during the US government shutdown,
 *    while the NIST government servers were unavailable. Code was located on
 *    the Internet Archive at:
 *
 *      https://web.archive.org/web/20180720195312/https://csrc.nist.gov/projects/random-bit-generation/documentation-and-software
 *
 *    which includes the following notice:
 *
 *      Software disclaimer: "This software was developed at the National
 *      Institute of Standards and Technology by employees of the Federal
 *      Government in the course of their official duties. Pursuant to title 17
 *      Section 105 of the United States Code this software is not subject to
 *      copyright protection and is in the public domain. The NIST Statistical
 *      Test Suite is an experimental system. NIST assumes no responsibility
 *      whatsoever for its use by other parties, and makes no guarantees,
 *      expressed or implied, about its quality, reliability, or any other
 *      characteristic. We would appreciate acknowledgment if the software is
 *      used."
 *
 *    End of NIST statement
 *
 * 3. ZOSLIB
 *    License: Apache-2.0
 *    ZOSLIB is a z/OS C/C++ library, available at:
 *      https://github.com/ibmruntimes/zoslib
 *    It is an extended implementation of the z/OS LE C Runtime Library.
 ******************************************************************************/

enum class GeneratorType
{
   TRNO,
   JITTER,
   DEVURANDOM
};


/* ------------------------------------------------------------------------ *
 * Public entry point.
 *
 * Contract mirrors POSIX/BSD getentropy(): fill "size" random_data, return 0, or -1
 * with errno set. The documented maximum for getentropy() is 256 random_data (the
 * original used 257 here; 256 is kept below, with the discrepancy noted).
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
        // prno_trng_generate (out, size);
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
   const int sample_size = 256;
   unsigned char random_data[256];
   GeneratorType generatorType = GeneratorType::JITTER;


       rc = zos_getentropy (random_data, sample_size, generatorType);

   generatorType = GeneratorType::TRNO;

   epsilon = (unsigned char *) malloc (sample_size * 8);

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

   std::free (epsilon);

   return 0;
}
