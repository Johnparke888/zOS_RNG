#pragma once

struct Parm_Block
{
   unsigned long long word1;
   unsigned long long word2;
};

#ifdef __cplusplus
extern "C"
{
#endif
   int test_function_code (struct Parm_Block *parm_block, int function);
   void prno_trng_generate (unsigned char *output_buffer_ptr, size_t size);
   unsigned char jitter_sample_byte (int shift);
   int prno_trng_installed ();
   void jitter_fill (unsigned char *output_buffer_ptr, size_t size);
   int dev_urandom_generate (unsigned char *output, size_t length);
   int naive_prng_generate (unsigned char *output, size_t length);
#ifdef __cplusplus
}
#endif
