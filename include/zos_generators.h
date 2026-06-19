
#define KLMD_FC_SHA512 3u
#define SHA512_DIGEST_LENGTH 64u

struct Parm_Block
{
   unsigned long long word1;
   unsigned long long word2;
};
struct DataBlock
{
   unsigned long long data[16];
};

#ifdef __cplusplus
extern "C"
{
#endif
   int test_function_code (struct Parm_Block *parm_block, int function);
   int prno_trng_generate (unsigned char *output_buffer_ptr, size_t size);
   unsigned char jitter_sample_byte (int shift);
   int prno_trng_installed ();
   void jitter_fill (unsigned char *output_buffer_ptr, size_t size);
   int dev_urandom_generate (unsigned char *output, size_t length);
   int naive_prng_generate (unsigned char *output, size_t length);
   int bad_raw_clock_generate (unsigned char *output, size_t length);
   int bad_hash_counter_generate (unsigned char *output, size_t length);
   unsigned char jitter_sample_byte (int shift);
   int read_full (int fd, unsigned char *output_buffer, size_t number_of_bytes_requested);
   int dev_urandom_generate (unsigned char *output, size_t length);
   int prno_trng_installed ();
   int test_function_code (struct Parm_Block *parm_block, int function);
   unsigned long long z_stckf64 ();
   int z_sha512_klmd (struct DataBlock *data_block, unsigned char digest[SHA512_DIGEST_LENGTH]);
#ifdef __linux__
   unsigned long long __stckf (unsigned long long *result);
#endif
#ifdef __cplusplus
}

#endif
