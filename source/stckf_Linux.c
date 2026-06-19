#include <errno.h>
#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>       // clock_gettime, struct timespec

#include "zos_generators.h"


//============================================================================
// Emulation of the z/Architecture STORE CLOCK FAST (STCKF) instruction.
//
// STCKF stores the 64-bit Time-Of-Day (TOD) clock. In the short (64-bit)
// architectural format the clock is an unsigned counter with:
//
//   * epoch  = 1900-01-01 00:00:00 UTC
//   * bit 51 = incremented once per microsecond
//
// Because bit 51 is worth 1 us, the LSB (bit 63) is worth 2^-12 us, so the
// whole value counts in units of 1/4096 of a microsecond. Converting a Unix
// CLOCK_REALTIME value (epoch 1970-01-01 UTC) therefore means:
//
//   1. add the fixed 70-year offset (2,208,988,800 s, includes 17 leap days)
//   2. rescale into TOD units (* 4096, i.e. << 12)
//
// Unlike STCK, STCKF does NOT guarantee a unique/monotonic value across rapid
// successive calls -- it just reads the clock. This emulation matches that:
// two calls within the same resolution tick can return the same value.
//
//============================================================================

unsigned long long int __stckf (unsigned long long int *tod)
{
   long long int stckf_value = 0;
   // Seconds between 1900-01-01 and 1970-01-01 (70 years, 17 leap days).
   const long long int kEpochDelta = 2208988800ULL;

   struct timespec ts;
   clock_gettime (CLOCK_REALTIME, &ts);

   const long long int sec_since_1900 = (long long int) (ts.tv_sec) + kEpochDelta;
   const long long int nsec = (long long int) (ts.tv_nsec);

   // Whole microseconds since 1900, scaled to TOD units (<<12 == * 4096).
   *tod = (sec_since_1900 * 1000000ULL + nsec / 1000ULL) << 12;

   // Fold the sub-microsecond remainder (0..999 ns) into the low TOD bits.
   *tod += (nsec % 1000ULL) * 4096ULL / 1000ULL;
   stckf_value = *tod;
   return stckf_value;
}


//----------------------------------------------------------------------------
// Inverse helper: turn a TOD value back into a Unix timespec. Handy in tests
// for asserting that a captured TOD decodes to the expected wall-clock time.
//----------------------------------------------------------------------------
struct timespec tod_to_timespec (long long int tod)
{
   const long long int kEpochDelta = 2208988800ULL;

   const long long int usec_since_1900 = tod >> 12;         // whole us
   const long long int sub_us_units = tod & 0xFFFULL;       // 1/4096 us units
   const long long int sec_since_1900 = usec_since_1900 / 1000000ULL;
   const long long int usec = usec_since_1900 % 1000000ULL;

   struct timespec ts;
   ts.tv_sec = (time_t) (sec_since_1900 - kEpochDelta);
   ts.tv_nsec = (long) (usec * 1000ULL + (sub_us_units * 1000ULL) / 4096ULL);
   return ts;
}
