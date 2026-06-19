#define _XOPEN_SOURCE 600
#define _OPEN_SYS_FILE_EXT 1
#define _OPEN_MSGQ_EXT 1

#include <errno.h>
#include <fcntl.h>        // Required for O_RDONLY
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>       // Required for close()


#ifdef __MVS__
#include <builtins.h> /* __stckf */
#include <ctest.h>
#include <psa.h>
#endif

#ifndef __MVS__
// #error "This file targets z/OS USS only."
#define __ptr32
#endif

#include "zos_generators.h"
       // clock_gettime, struct timespec

static int naive_counter = 0;


// XL-specific NR parameter constraint:
// https://www.ibm.com/docs/en/zos/2.4.0?topic=statements-inline-assembly-extension

/* PRNO function codes (bits 57-63 of GR0 select the function). */
enum
{
   PRNO_QUERY = 0, /* report installed functions into a 16-byte parm parm_block */
   PRNO_TRNG = 114 /* 0x72: true random generate; no parameter parm_block       */
};
static int cached = -1; /* Cache the result of the PRNO-TRNG check. -1 = not yet checked. */

/*
 * read() may return fewer bytes than requested, or be interrupted by a
 * signal. Loop until the whole buffer is filled. EINTR is retried; a return
 * of 0 (EOF) from /dev/urandom is abnormal and is treated as a failure.
 */
int read_full (int fd, unsigned char *output_buffer, size_t number_of_bytes_requested)
{
   size_t total_retrieved = 0;
   while (total_retrieved < number_of_bytes_requested)
   {
      ssize_t bytes_read = read (fd, output_buffer + total_retrieved, number_of_bytes_requested - total_retrieved);
      if (bytes_read <= 0)
      {
         if (bytes_read < 0 && errno == EINTR)
         {
            continue;
         }
         return -1;
      }
      total_retrieved += (size_t) bytes_read;
   }
   return 0;
}

/*
 */
int dev_urandom_generate (unsigned char *output, size_t length)
{

   int fd = open ("/dev/urandom", O_RDONLY);

   if (fd < 0)
   {
      return -1;
   }

   if (read_full (fd, output, length) != 0)
   {
      close (fd);
      return -2;
   }

   close (fd);

   return 0;
}


int test_function_code (struct Parm_Block *parm_block, int function)
{
   if (parm_block == NULL || function < 0 || function > 127)
   {
      return false;
   }
   if (function < 64)
   {
      return (parm_block->word1 & (1ULL << (63 - function))) != 0;
   }
   else
   {
      return (parm_block->word2 & (1ULL << (127 - function))) != 0;
   }
}

/* ------------------------------------------------------------------------ *
 * Feature detection: is PRNO-TRNG (function 114) installed?
 *
 * Two gates:
 *   (1) A cheap low-storage facility-indicator check (carried over from the
 *       original; verify on your box -- see file header).
 *   (2) PRNO Query (function 0), which writes a bitmap of installed function
 *       codes into a 16-byte parameter parm_block addressed by GR1. Installed-
 *       function bit 114 sits in the second doubleword of that bitmap, at
 *       1ULL << (127 - 114) == 0x2000.
 *
 * ------------------------------------------------------------------------ */
int prno_trng_installed ()
{

#ifdef __linux__
   return 1;
#endif
#ifdef __MVS__
   struct psa *__ptr32 psa_ptr = 0; /* PSA is always at virtual address 0. */
   /* Gate 1: facility indicator in the PSA. 0x40 at byte 207. */
   if (!(0x40 & psa_ptr->flcfacl7))
   {
      return false;
   }
#endif
   /* Gate 2: PRNO Query. Parameter parm_block is 16 random_data for function 0. */

   /* For Query, GR0 = function code 0 and GR1 = parm-parm_block address; the
    * R1/R2 register fields (here GR8/GR9 and GR10/GR11) are not used to move
    * data, so nothing meaningful is bound to them. The "memory" clobber tells
    * the compiler the parm parm_block is written.
    *
    * PRNO-Query (PRNO Function Code 0)
    *
    * The contents of general registers R1, R1 + 1, R2, and R2 + 1 are ignored by the query function.
    * A 128-bit status word is stored in the parameter parm_block. Bits 0-127 of this field correspond
    * to function codes 0-127, respectively, of the PERFORM RANDOM NUMBER OPERATION instruction.
    * When a bit is one, the corresponding function is installed; otherwise, the function is not installed.
    * Condition code 0 is set when execution of the PRNOQuery function completes; condition code 3 is not applicable to this function

    */
   if (cached != -1)
   {
      return cached;
   }
   struct Parm_Block parm_block = {0, 0};
#ifdef __MVS__
   asm volatile (" prno 8,10\n"
                 " jo *-4\n" /* CC3 => operation incomplete; reissue - instruction length is 4   */
                 :
                 : "{r0}"((unsigned long) PRNO_QUERY), "{r1}"(&parm_block)
                 : "memory");
#endif
   cached = test_function_code (&parm_block, PRNO_TRNG) ? 1 : 0;

   return cached;
}


/*
 * PRNO-TRNG (PRNO Function Code 114)
 *
 * A series of hardware-generated random numbers is stored at either or both the first- and second-operand locations.
 * A parameter parm_block is not used by the PRNO-TRNG function, and general register 1 is ignored.
 *
 *The R1 field designates an even-odd pair of general registers. The even-numbered register contains the address of the leftmost
 * byte of the first operand, and the odd-numbered register contains the length of the first operand. The first operand comprises random numbers in
 * the form of raw entropy, produced directly by a hardware source.
 *
 * The R2 field designates an even-odd pair of general registers. The even-numbered register contains the address of the leftmost byte of the second
 * operand, and the odd-numbered register contains the length of the second operand. The second operand comprises random numbers extracted from the
 * raw entropy source and then conditioned by an approved algorithm.
 *
 * When the length of the first operand is nonzero but the length of the second operand is zero, the process continues with storing the raw entropy in
 * the first operand only. Similarly, when the length of the second operand is nonzero but the length of the first operand is zero, the process
 * continues with storing the conditioned entropy in the second operand only.
 *
 * The process continues either until both operand lengths are zero, called normal completion, or until a CPU-determined number of random_data have
 * been stored. called partial completion, whichever occurs first. When the operation ends due to normal completion, condition code 0 is set. When the
 * operation ends due to partial completion, condition code 3 is set.
 *
 */

#ifdef __MVS__
void prno_trng_generate (unsigned char *output_buffer_ptr, size_t size)
{
   /* First operand (raw random_data) is unused; pin it to defined values. */
   unsigned char *raw_addr = 0;
   unsigned long raw_length = 0;

#ifdef __XPLINK__

   __asm__ volatile (" prno 10,2\n"
                     " jo *-4\n"
                     : "+{r2}"(output_buffer_ptr), "+{r3}"(size)
                     : "{r0}"((unsigned long) PRNO_TRNG), "{r10}"(raw_addr), "{r11}"(raw_length)
                     : "memory");
#else
   /* R1 = GR8/GR9 (raw, len 0), R2 = GR10/GR11 (output_buffer, size).
    */

   __asm__ volatile (" prno 8,10\n"
                     " jo *-4\n"
                     : "+{r10}"(output_buffer_ptr), "+{r11}"(size)
                     : "{r0}"((unsigned long) PRNO_TRNG), "{r8}"(raw_addr), "{r9}"(raw_length)
                     : "memory");

#endif
}
#endif
/* ------------------------------------------------------------------------ *
 * Fallback path: CPU timing-jitter entropy (only used pre-z14).
 *
 * jitter_sample_byte() returns one byte folded from many clock readings taken
 * across supervisor transitions. SVC 137 is issued purely to force a
 * problem-state -> supervisor-state round trip whose latency varies
 * unpredictably; the STCKF read that follows captures that variance. "shift"
 * discards the high (predictable, monotonic) clock bits so only the noisy
 * low-order bits are accumulated. The loop stops once enough wall-clock time
 * has elapsed (delta > 0xfffff) to bound how long sampling can run.
 *
 *
 * CALLDISP â Pass control to another ready task
 *
 * Description
 * The CALLDISP macro saves the caller's status in the current TCB/RB, and passes
 * control to another ready task. The task with the highest priority is
 * the one that receives control. When the original task is re-dispatched,
 * control is returned to the next sequential instruction.
 *
 * These are the requirements for the caller:
 *  Environment
 *  When BRANCH=NO
 *  Environmental factor
 *  Requirement
 *  Minimum authorization:      None.
 *  Dispatchable unit mode:     Task
 *  Cross memory mode:  PASN=HASN=SASN
 *  AMODE:      24- or 31- or 64-bit
 *  ASC mode:   Primary
 *  Interrupt status: Enabled for I/O and external interrupts
 *  Locks:      No locks held
 *  Control parameters: None.
 *
 * Input register information
 *
 * Before issuing the CALLDISP macro, the caller does not have to place any
 * information into any register unless using it in register notation for a
 * particular parameter, or using it as a base register.
 *
 * Output register information
 *
 * When control returns to the caller, the GPRs contain:
 * Register   Contents
 * 0-13      Unchanged
 * 14        Unchanged when BRANCH=NO, used as a work register by the system when BRANCH=YES
 * 15        Used as a work register by the system
 *
 * ------------------------------------------------------------------------ */

unsigned char jitter_sample_byte (int shift)
{
   unsigned long long int start_time;
   unsigned long long int t0;
   unsigned long long int accumulator = 0;
   int i = 0;
#ifdef __MVS__
   //__asm__ [volatile] ( template : outputs : inputs : clobbers );

   // CALLDISP branch=no, sets r15 to 0 then issues SVC 137
   __asm__ volatile (" la 15,0\n svc 137\n" ::: "r15", "r6");
#endif
   __stckf (&start_time);

   start_time >>= shift;
   // Up to 400 iterations of SVC+STCKF
   // Stops early if too much shifted-clock time has passed since start

   for (i = 0; i < 400; ++i)
   {
#ifdef __MVS__
      __asm__ volatile (" la 15,0\n svc 137\n" ::: "r15", "r6");       // r15 and r6 are clobber register
#endif
      __stckf (&t0);
      t0 >>= shift;
      if ((t0 - start_time) > 0xfffff)
      {
         break;
      }
      accumulator ^= t0;
   }
   return (unsigned char) accumulator;
}

/*
 * jitter_fill() first calibrates the timer's effective resolution: it takes
 * the minimum of several XOR'd STCKF deltas, isolates that value's lowest set
 * bit (-m & m), and looks its position up in a de Bruijn table. That position
 * ("bits") estimates how many low-order clock bits are pure noise, and is used
 * as the shift for jitter_sample_byte(). Bytes are then chained through an
 * XOR-accumulator so each output depends on the running state.
 *
 * clock_reading_1 and clock_reading_2 hold consecutive STCKF clock reads.
 * r holds their XOR difference.
 * m tracks the smallest observed XOR value.
 * byte_accumulator is a running byte accumulator initialized to 0xaa.
 * the table is a compact way to say: âfind the index of the least-significant 1 bit.â
 *
 * Theory of operation
 *
 * These routines provide a fallback entropy source for z/OS USS systems when a true hardware
 * or operating-system random number generator is not available.
 *
 * The design is based on timing jitter. It repeatedly samples the z/Architecture high-resolution clock using STCKF and mixes small variations in
 * execution timing into output random_data. The timing variation is influenced by dispatcher activity, interrupt timing, CPU scheduling, cache
 * effects, pipeline state, and other system noise.
 *
 * Before producing random_data, jitter_fill() performs a short calibration step. It reads the clock back-to-back several times and examines the
 * low-order changing bits. From this it estimates how many low-order clock bits are unstable or noisy. That value is used as a shift count when
 * collecting samples, so the sampler avoids relying directly on clock bits that are too fine or too deterministic.
 *
 * For each output byte, jitter_sample_byte() performs repeated dispatcher calls using SVC 137, reads the clock after each call, shifts the timestamp
 * by the calibrated amount, and XORs the sampled values into an accumulator. The loop is bounded by both a maximum iteration count and an
 * elapsed-time cutoff so it cannot run indefinitely.
 *
 * jitter_fill() then chains the generated random_data through a running XOR accumulator. Each output byte depends on the previous accumulator state
 * and the latest jitter sample, rather than being a raw timestamp byte.
 *
 * Conceptually, the flow is:
 *
 * Calibrate STCKF timing behavior
 *         â
 * Estimate useful timing-jitter bit position
 *         â
 * For each requested byte:
 *     force dispatcher interaction
 *    sample STCKF repeatedly
 *     mix shifted timestamps with XOR
 *     chain result into output accumulator
 *         â
 * Return fallback random-looking random_data
 *
 * This mechanism should be considered a best-effort fallback entropy source, not a replacement for a real TRNG, ICSF, CPACF/PRNO, /dev/random, or a
 * cryptographic DRBG seeded from trusted entropy. Its output should preferably be mixed into a cryptographic hash or used only to seed a standard
 * DRBG,
 * rather than used directly as security-critical random data.
 *  */

void jitter_fill (unsigned char *output_buffer_ptr, size_t size)
{

#ifdef _LP64
   static int bits = 0;
#else
   int bits = 0;
#endif
   unsigned long long clock_reading_1 = 0;
   unsigned long long clock_reading_2 = 0;
   unsigned long long r = 0;
   unsigned long long m = 0xffff;
   unsigned int byte_accumulator = 0xaa;
   size_t i = 0;

   /* de Bruijn-style lowest-set-bit position lookup (index = value % 37). */
   const unsigned int zbitcnt[] = {0xffffffff, 0,  1,  26, 2,  23, 27, 0,  3, 16, 24, 30, 28, 11, 0,  13, 4,  7, 17,
                                   0,          25, 22, 31, 15, 29, 19, 12, 6, 0,  21, 14, 9,  5,  20, 8,  19, 18};

   /* Calibrate: find a noise-floor bit position in (1, 11]. */
   while (bits == 0 || bits > 11)
   {
      for (int k = 0; k < 10; ++k)
      {
         __stckf (&clock_reading_1);
         __stckf (&clock_reading_2);
         r = clock_reading_1 ^ clock_reading_2;
         if (r < m)
         {
            m = r;
         }
      }
      bits = zbitcnt[(-m & m) % 37];
   }

   for (i = 0; i < size; ++i)
   {
      byte_accumulator ^= jitter_sample_byte (bits);
      output_buffer_ptr[i] = (unsigned char) byte_accumulator;
   }
}



/*
 * KLMD-SHA-512 parameter block:
 *   64 bytes: SHA-512 chaining value
 *   16 bytes: 128-bit bit count for data already processed by prior KIMD calls
 *
 * For a one-shot KLMD call, initialize the chaining value to the SHA-512 IV
 * and set the prior bit count to zero.
 */
struct sha512_klmd_parm
{
   unsigned long long H[8];
   unsigned long long block_length_high;       // High 64 bits of total bit count
   unsigned long long block_length_low;        // bit length
};

static const struct sha512_klmd_parm sha512_initial_parm = {{0x6a09e667f3bcc908ULL,
                                                      0xbb67ae8584caa73bULL,
                                                      0x3c6ef372fe94f82bULL,
                                                      0xa54ff53a5f1d36f1ULL,
                                                      0x510e527fade682d1ULL,
                                                      0x9b05688c2b3e6c1fULL,
                                                      0x1f83d9abfb41bd6bULL,
                                                      0x5be0cd19137e2179ULL}};



/*
 * Store Clock Fast: returns the 8-byte TOD clock value.
 */
unsigned long long z_stckf64 ()
{
   unsigned long long start_time = 0;

   (void) __stckf (&start_time);

   return start_time;
}

/*
*
*
* COMPUTE LAST MESSAGE DIGEST
* KLMD R1,R2 [RRE]
* R1 and R2 do not represent specific general-purpose registers; they are positional notations for the first
* and second register operands.
* Here R1 is a positional designation, while GR1 is a specific register. The same applies to R2 and GR2.

*
* Operand register 1 (R1)is ignored.
* The R2 field designates an even-odd pair of general registers and must designate an even-numbered register
* other than general register 0.
*
* The location of the leftmost byte of the second operand is specified by the contents of the R2 general
* register. The number of bytes in the second-operand location is specified in general register R2 + 1.
* As part of the operation, the address in general register R2 is incremented by the number of bytes processed
* from the second operand, and the length in general register R2 + 1 is decremented by the same
* number. The formation and updating of the address and length is dependent on the addressing mode.
*
* The second operand is processed as specified by the function code using an initial
* chaining value in the parameter block, and the result replaces the chaining value.
* For the SHA-1,SHA-256, and SHA-512 functions, the operation also uses a message bit length
* in the parameter block.
*
* The message digest for the message (M) in operand 2 is generated using the SHA-512 algorithm with the
* chaining value (called H fields) and message-bit length information in the parameter block.
* If the message in operand 2 is equal to or greater than 128 bytes, an intermediate message digest is
* generated for each 128-byte message block using the SHA-512 block digest algorithm with the 64-byte
* chaining value in the parameter block, and the generated intermediate message digest, also called the
* output chaining value (OCV), is stored into the chaining-value field of the parameter block.
* This operation repeats until the remaining message is less than 128 bytes or until a CPU-determined
* number of blocks have been stored.
*
* For the KLMD-SHA-1, KLMD-SHA-256, and KLMD-SHA-512 functions, when processing the last message part,
* the program must compute the length of the original message in bits and place this length value in the
* message-bit-length field of the parameter block, and use the COMPUTE LAST MESSAGE
* DIGEST instruction.
*
* The COMPUTE LAST MESSAGE DIGEST instruction does not require the second operand to be a multiple of the data block size. It
* first processes complete blocks, and may set condition code 3 before processing all blocks. After processing all complete blocks,
* it then performs the padding operation including the remaining portion of the second operand. This may require one or two
* iterations of the designated block digest algorithm.
* KIMD

 * Code   Function       Parm.          Data
                         Block Size Block Size
                         (bytes)     (bytes)
 0        KLMD-Query       16            â
 1        KLMD-SHA-1       28           64
 2        KLMD-SHA-256     40           64
 3        KLMD-SHA-512     80          128
32        KLMD-SHA3-224   200          144
33        KLMD-SHA3-256   200          136
34        KLMD-SHA3-384   200          104
35        KLMD-SHA3-512   200           72
36        KLMD-SHAKE-128  200          168
37        KLMD-SHAKE-256  200          136

* Register Use:
* GR0 = 3 (SHA-512 function code)
* GR2 = ignored
* GR1 = pointer to parm block
* GR4 = pointer to data block
* GR5 = length of data block
*
* condition code 0 normal completion
* condition code 3 partial completion
*
*/

#ifdef __MVS__
int z_sha512_klmd (DataBlock *data_block, unsigned char digest[SHA512_DIGEST_LENGTH])
{

   sha512_klmd_parm parameter_block;
   unsigned char dummy = 0;

   if (data_block == NULL)
   {
      errno = EINVAL;
      printf ("KLMD-SHA-512: invalid input buffer pointer\n");
      return -1;
   }

   memcpy (&parameter_block.H, &sha512_initial_parm, sizeof (sha512_initial_parm));
   parameter_block.block_length_high = 0;
   parameter_block.block_length_low = sizeof (DataBlock) * 8ull;

   // printf ("KLMD-SHA-512: input length = %zu bytes\n", sizeof (DataBlock));


   // R0 = Bit positions 57-63 of general register 0 contain the function code.
   // Register pairs: R4/R5
   // r2 - Operand register 1 is ignored.
   // R1 - parameter block address
   // R4 = data_block, R5 = input_length
   unsigned long long int r0 = KLMD_FC_SHA512;
   unsigned long long int r1 = (unsigned long long int) (uintptr_t) &parameter_block;
   unsigned long long int r2 = 0;
   unsigned long long int r4 = (unsigned long long int) (uintptr_t) data_block;
   unsigned long long int r5 = sizeof (DataBlock);


   //  print register values for debugging
   // printf ("KLMD-SHA-512: R0 = %lu, R1 = 0x%lx, R2 = %lu, R4 = 0x%lx, R5 = %lu\n", r0, r1, r2, r4, r5);

   /*
    * Format of any asm statement is:
    * asm volatile ("instruction" : output_operands : input_operands : clobbers);
    */

   /* GR0=FC, GR1=parm (architecture-mandated); operand-2 pair held in
      GR2/GR3 â volatile under XPLINK and usable under standard linkage.
      R1 field is ignored, so 0 is fine. */
   asm volatile (" KLMD 0,2\n"
                 " jo   *-4\n" /* CC3 = partial completion; reissue */
                 : "+{r2}"(r4), "+{r3}"(r5)
                 : "{r0}"(r0), "{r1}"(r1)
                 : "memory", "cc");


   memcpy (digest, parameter_block.H, SHA512_DIGEST_LENGTH);
   return 0;
}
#endif

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

   printf ("Generating %zu bytes of random data using naive_prng_generate...\n", length);

   struct DataBlock data_block;
   unsigned char digest[SHA512_DIGEST_LENGTH];

   while (produced < length)
   {
      /* Re-seed every iteration: z_stckf64() returns a fresh store-clock
         value on each call, so each block (and thus each digest) differs. */
      memset (&data_block, 0, sizeof (data_block));

      for (size_t i = 0; i < 16; ++i)
      {
         data_block.data[i] = z_stckf64 ();
      }

      z_sha512_klmd (&data_block, digest);

      size_t remaining = length - produced;
      size_t chunk = (remaining < SHA512_DIGEST_LENGTH) ? remaining : SHA512_DIGEST_LENGTH;

      memcpy (output + produced, digest, chunk);
      produced += chunk;
   }

   return 0;
}
int bad_raw_clock_generate (unsigned char *output, size_t length)
{
   size_t produced = 0;

   while (produced < length)
   {
      unsigned long long time = z_stckf64 ();
      size_t chunk = length - produced < sizeof (time) ? length - produced : sizeof (time);

      memcpy (output + produced, &time, chunk);
      produced += chunk;
   }

   return 0;
}

int bad_hash_counter_generate (unsigned char *output, size_t length)
{
   static unsigned long long counter = 0;
   size_t produced = 0;

   while (produced < length)
   {
      struct DataBlock data_block;
      unsigned char digest[SHA512_DIGEST_LENGTH];

      data_block.data[0] = counter++;

      z_sha512_klmd (&data_block, digest);

      size_t remaining = length - produced;
      size_t chunk = remaining < SHA512_DIGEST_LENGTH ? remaining : SHA512_DIGEST_LENGTH;

      memcpy (output + produced, digest, chunk);
      produced += chunk;
   }

   return 0;
}
