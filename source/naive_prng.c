/*
 * naive_prng.c
 *
 * Naive z/OS USS pseudo/random-ish generator:
 *   - sample STCKF
 *   - mix in a counter and a little process-local variation
 *   - hash that seed material with KLMD-SHA-512
 *
 * Return:
 *   0  success
 *  -1  error, errno set
 *
 * This is NOT a cryptographic RNG. Prefer PRNO function code 114 when
 * true random data is available.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#if !defined(__LP64__)
#error "This sample expects 64-bit mode. Compile with ibm-clang64 or -m64."
#endif

#define KLMD_FC_SHA512 3u
#define SHA512_DIGEST_LEN 64u

/*
 * KLMD-SHA-512 parameter block:
 *   64 bytes: SHA-512 chaining value
 *   16 bytes: 128-bit bit count for data already processed by prior KIMD calls
 *
 * For a one-shot KLMD call, initialize the chaining value to the SHA-512 IV
 * and set the prior bit count to zero.
 */
typedef struct
{
   uint64_t h[8];
   unsigned long long block_length_high;
   unsigned long long block_length_low;
} sha512_klmd_parm_t __attribute__ ((aligned (16)));

static const sha512_klmd_parm_t sha512_initial_parm = {{0x6a09e667f3bcc908ULL,
                                                        0xbb67ae8584caa73bULL,
                                                        0x3c6ef372fe94f82bULL,
                                                        0xa54ff53a5f1d36f1ULL,
                                                        0x510e527fade682d1ULL,
                                                        0x9b05688c2b3e6c1fULL,
                                                        0x1f83d9abfb41bd6bULL,
                                                        0x5be0cd19137e2179ULL},
                                                       {0}};

static int naive_counter = 0;

/*
 * Store Clock Fast: returns the 8-byte TOD clock value.
 */
static inline uint64_t z_stckf64 (void)
{
   uint64_t v = 0;

   asm volatile (" STCKF %0\n" : "=m"(v) : : "cc", "memory");

   return v;
}

/*
 * One-shot SHA-512 using KLMD.
 *
 * R0 = 3 (SHA-512 function code)
 * R1 = pointer to parm block
 *
The R2 field designates an even-odd pair of general registers and must designate an even-numbered register
other than general register 0.
The location of the leftmost byte of the second operand is specified by the contents of the R2 general
register. The number of bytes in the second-operand location is specified in general register R2 + 1.
As part of the operation, the address in general register R2 is incremented by the number of bytes processed
from the second operand, and the length in general register R2 + 1 is decremented by the same
number. The formation and updating of the address and length is dependent on the addressing mode.

 *
 * KLMD may be interruptible, so reissue until R3 reaches zero.
 */
static int z_sha512_klmd (const void *input_buffer_ptr, size_t input_length, unsigned char digest[SHA512_DIGEST_LEN])
{
   sha512_klmd_parm_t parameter_block;
   unsigned char dummy = 0;

   memcpy (&parameter_block, &sha512_initial_parm, sizeof (parameter_block));
   parameter_block.block_length_high = 0;
   parameter_block.block_length_low = input_length * 8ull;

   if (input_buffer_ptr == NULL && input_length != 0)
   {
      errno = EINVAL;
      return -1;
   }
   // R0 = Bit positions 57-63 of general register 0 contain the function code.
   // Register pairs: R4/R5
   // r2 - Operand register 1 is ignored.
   // R1 - parameter block address
   // R4 = input_buffer_ptr, R5 = input_length
   unsigned long r0 = KLMD_FC_SHA512;
   unsigned long r1 = (unsigned long) (uintptr_t) &parameter_block;
   unsigned long r2 = 0;
   unsigned long r4 = (unsigned long) (uintptr_t) (input_length ? input_buffer_ptr : &dummy);
   unsigned long r5 = (unsigned long) input_length;


   parameter_block.block_length_low = r4 * 8ull;

   while (r5 != 0 || input_length == 0)
   {
      /*
       * Format of any asm statement is:
       * asm volatile ("instruction" : output_operands : input_operands : clobbers);
       */
      asm volatile (" KLMD 2,4\n"
                    " jo *-4\n" /* CC==3 (partial completion) -> retry */
                    : "+{r4}"(r4), "+{r5}"(r5)
                    : "{r0}"(r0), "{r1}"(r1), "{r4}"(r4), "{r5}"(r5)
                    : "cc", "memory");
       
   }

   memcpy (digest, parameter_block.h, SHA512_DIGEST_LEN);
   return 0;
}

typedef struct
{
   char domain[16];
   uint64_t stckf_before;
   uint64_t stckf_after;
   uint64_t counter;
   uint64_t pid;
   uint64_t requested_len;
   uint64_t offset;
   uintptr_t output_addr;
   uintptr_t stack_addr;
} naive_seed_material_t;

/*
 * Public function.
 */
int naive_prng_generate (unsigned char *output, size_t length)
{
   if (length == 0)
   {
      return 0;
   }

   if (output == NULL)
   {
      errno = EINVAL;
      return -1;
   }

   size_t produced = 0;

   while (produced < length)
   {
      naive_seed_material_t seed;
      unsigned char digest[SHA512_DIGEST_LEN];

      memset (&seed, 0, sizeof (seed));
      memcpy (seed.domain, "naive-prng-v1", 13);

      seed.stckf_before = z_stckf64 ();
      seed.counter ++;
      seed.pid = (uint64_t) getpid ();
      seed.requested_len = (uint64_t) length;
      seed.offset = (uint64_t) produced;
      seed.output_addr = (uintptr_t) (output + produced);
      seed.stack_addr = (uintptr_t) &seed;
      seed.stckf_after = z_stckf64 ();

      if (z_sha512_klmd (&seed, sizeof (seed), digest) != 0)
      {
         return -1;
      }

      size_t remaining = length - produced;
      size_t n = remaining < SHA512_DIGEST_LEN ? remaining : SHA512_DIGEST_LEN;

      memcpy (output + produced, digest, n);
      produced += n;

      /*
       * Not secret in the strong sense, but clear the temporary digest anyway.
       */
      memset (digest, 0, sizeof (digest));
   }

   return 0;
}

#ifdef TEST_NAIVE_PRNG
#include <stdio.h>

int main (void)
{
   unsigned char buf[128];

   if (naive_prng_generate (buf, sizeof (buf)) != 0)
   {
      perror ("naive_prng_generate");
      return 1;
   }

   for (size_t i = 0; i < sizeof (buf); ++i)
   {
      printf ("%02x", buf[i]);
      if ((i + 1) % 32 == 0)
      {
         putchar ('\n');
      }
   }

   return 0;
}
#endif
