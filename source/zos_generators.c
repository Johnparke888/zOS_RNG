
/*
 * zos_getentropy.cpp
 *
 * WHAT IT DOES
 * ------------
 * getentropy()-style fill with two backends, chosen at run time:
 *
 *   1. CPACF True Random Number Generator (the fast/strong path).
 *      Uses the PERFORM RANDOM NUMBER OPERATION (PRNO) instruction with
 *      function code 114 (0x72) = "PRNO-TRNG", a hardware noise source added
 *      with Message-Security-Assist Extension 7 (z14 and later). This is the
 *      same facility that backs /dev/[u]random on z14+; it produces conditioned
 *      true-random random_data and uses no parameter parm_block.
 *
 *   2. Timing-jitter fallback (the slow path, pre-z14 only).
 *      Harvests entropy from the nondeterministic latency of repeated
 *      supervisor transitions, sampled with STORE CLOCK FAST (__stckf). Note
 *      this measures timing *variance*, not a raw clock value -- a single
 *      STCKF reading would be predictable and worthless; the jitter between
 *      many readings is the actual entropy source.
 *
 * BUILD (test harness)
 * --------------------
 *   ibm-clang++ -m64 -D ZOS_GETENTROPY_TEST -o zostest zos_getentropy.cpp
 *   ./zostest
 *
 */

#ifndef __MVS__
// #error "This file targets z/OS USS only."
#define __ptr32
#endif
#include <errno.h>
#include <stddef.h>
#include <builtins.h> /* __stckf */
#include <psa.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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



/* PRNO function codes (bits 57-63 of GR0 select the function). */
enum
{
   PRNO_QUERY = 0, /* report installed functions into a 16-byte parm parm_block */
   PRNO_TRNG = 114 /* 0x72: true random generate; no parameter parm_block       */
};
static int cached = -1; /* Cache the result of the PRNO-TRNG check. -1 = not yet checked. */
struct Parm_Block
{
   unsigned long long word1;
   unsigned long long word2;
};

int test_function_code (const struct Parm_Block *parm_block, int function)
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


   psa *__ptr32 psa_ptr = 0; /* PSA is always at virtual address 0. */
   /* Gate 1: facility indicator in the PSA. 0x40 at byte 207. */
   if (!(0x40 & psa_ptr->flcfacl7))
   {
      return false;
   }

   /* Gate 2: PRNO Query. Parameter parm_block is 16 random_data for function 0. */

   /* For Query, GR0 = function code 0 and GR1 = parm-parm_block address; the
    * R1/R2 register fields (here GR8/GR9 and GR10/GR11) are not used to move
    * data, so nothing meaningful is bound to them. The "memory" clobber tells
    * the compiler the parm parm_block is written.
    *
    * PRNO-Query (PRNO Function Code 0)
    *
    * The contents of general registers R1, R1 + 1, R2, and R2 + 1 are ignored by the query function.
    * A 128-bit status word is stored in the parameter parm_block. Bits 0-127 of this field correspond to function codes 0-127, respectively, of the
    PERFORM
    * RANDOM NUMBER OPERATION instruction. When a bit is one, the corresponding function is installed; otherwise, the function is not installed.
    Condition
    * code 0 is set when execution of the PRNOQuery function completes; condition code 3 is not applicable to this function

    */
   if (cached != -1)
   {
      return cached;
   }
   struct Parm_Block parm_block = {0, 0};

   asm volatile (" prno 8,10\n"
                 " jo *-4\n" /* CC3 => operation incomplete; reissue - instruction length is 4   */
                 :
                 : "{r0}"((unsigned long) PRNO_QUERY), "{r1}"(&parm_block)
                 : "memory");

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
   unsigned long long start_time;
   unsigned long long t0;
   unsigned long long accumulator = 0;
   int i = 0;
   //__asm__ [volatile] ( template : outputs : inputs : clobbers );
   // CALLDISP branch=no, sets r15 to 0 then issues SVC 137
   __asm__ volatile (" la 15,0\n svc 137\n" ::: "r15", "r6");
   (void) __stckf (&start_time);

   start_time >>= shift;
   // Up to 400 iterations of SVC+STCKF
   // Stops early if too much shifted-clock time has passed since start

   for (i = 0; i < 400; ++i)
   {
      __asm__ volatile (" la 15,0\n svc 137\n" ::: "r15", "r6");       // r15 and r6 are clobber register
      (void) __stckf (&t0);
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
* low-order changing bits. From this it estimates how many low-order clock bits are unstable or noisy. That value is used as a shift count when collecting
* samples, so the sampler avoids relying directly on clock bits that are too fine or too deterministic.
*
* For each output byte, jitter_sample_byte() performs repeated dispatcher calls using SVC 137, reads the clock after each call, shifts the timestamp
* by the calibrated amount, and XORs the sampled values into an accumulator. The loop is bounded by both a maximum iteration count and an elapsed-time
* cutoff so it cannot run indefinitely.
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
   constexpr unsigned int zbitcnt[] = {0xffffffff, 0,  1,  26, 2,  23, 27, 0,  3, 16, 24, 30, 28, 11, 0,  13, 4,  7, 17,
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
 