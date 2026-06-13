#pragma once
#ifdef __cplusplus
extern "C"
{
#endif
   int test_function_code (const struct Parm_Block *parm_block, int function);
   static void prno_trng_generate (unsigned char *output_buffer_ptr, size_t size);
   static unsigned char jitter_sample_byte (int shift);
   static bool prno_trng_installed ();
#ifdef __cplusplus
}
#endif
