#include <cstddef>           // size_t
#include <cstdint>
#include <cstring>           // std::memcpy
#include <ctime>
#include <errno.h>
#include <immintrin.h>
#include <iostream>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>       // clock_gettime, struct timespec

#include "zos_generators.h"
#include <botan/hash.h>
#include <botan/hex.h>


constexpr int RDRAND_BUFFER_SIZE = 16;
constexpr int MAX_RDRAND_ATTEMPTS = 10;

int z_sha512_klmd (DataBlock *data_block, unsigned char digest[SHA512_DIGEST_LENGTH])
{

   // 2. Create the SHA-512 hash object
   std::unique_ptr<Botan::HashFunction> sha512 = Botan::HashFunction::create ("SHA-512");
   if (!sha512)
   {
      std::cerr << "Error: SHA-512 algorithm not available in this build." << std::endl;
      return 1;
   }

   // 3. Feed data into the hash function
   sha512->update (reinterpret_cast<const uint8_t *> (data_block), sizeof (DataBlock));

   // 4. Retrieve the output buffer size and allocate memory
   std::vector<uint8_t> digest_result (sha512->output_length ());

   // 5. Finalize the computation and write to the buffer
   sha512->final (digest_result.data ());
   // copy the digest to the output parameter
   std::memcpy (digest, digest_result.data (), SHA512_DIGEST_LENGTH);
   return 0;
}



int prno_trng_generate (unsigned char *output_buffer_ptr, size_t size)
{
   if (output_buffer_ptr == nullptr)
   {
      errno = EINVAL;
      return -1;
   }

   unsigned char *out = output_buffer_ptr;
   size_t bytesRemaining = size;

   while (bytesRemaining > 0)
   {
      unsigned long long value = 0;
      bool retrievedData = false;

      for (int attempt = 0; attempt < MAX_RDRAND_ATTEMPTS; ++attempt)
      {
         if (_rdrand64_step (&value) != 0)
         {
            retrievedData = true;
            break;
         }
      }

      if (!retrievedData)
      {
         errno = EIO;
         std::cerr << "Error: RDRAND failed after " << MAX_RDRAND_ATTEMPTS << " attempts." << std::endl;
         return -1;
      }

      size_t retrievedSize = std::min (bytesRemaining, sizeof (value));
      std::memcpy (out, &value, retrievedSize);

      out += retrievedSize;
      bytesRemaining -= retrievedSize;
   }

   return 0;
}
