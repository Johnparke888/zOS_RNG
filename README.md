# z/OS Random Number Generation Test Suite

This project was created solely as a way to experiment with and test various
methods of random number generation on z/OS Unix System Services. It should
also work on z/OS.

It was compiled with the IBM Open XL C/C++ 2.1 compiler for z/OS
(Â© Copyright International Business Machines Corporation 2024).

## Random Number Generators

Currently it tests three methods of random number generation:

1. **PRNO-TRNG (function 114)** â uses `asm` statements to access the 
    Perform Random Number Operation PRNO-TRNG (PRNO Function Code 114) instruction
2. **Jitter** â CPU timing-jitter entropy
3. **`/dev/urandom`** â provides cryptographically secure random output
4. **Naive PRNG** â a very simple and insecurepseudo-random number generator that uses the system clock and SHA-512 to generate random output.
5. **Bad Raw Clock** â a very simple and insecure generator that uses the
	 system clock as a source of entropy
6. **Bad Hash Counter** â a very simple and insecure generator that uses a
	 hash function on a counter as a source of entropy


The first two generators were adapted from code found in ZOSLIB. The
statistical code was taken from the NIST Statistical Test Suite.
PRNG generators 4-6 were added for testing purposes, to verify that the statistical tests can detect poor generators.

## Sources of Information and Code

### 1. z/Architecture Principles of Operation

Fourteenth Edition (May, 2022)\
Â© Copyright International Business Machines Corporation 2000, 2022.
All rights reserved.

### 2. NIST Special Publication 800-22 revision 1a

National Institute of Standards and Technology\
Natl. Inst. Stand. Technol. Spec. Publ. 800-22rev1a, 131 pages (April 2010)\
*A Statistical Test Suite for Random and Pseudorandom Number Generators for
Cryptographic Applications*

This is free and unencumbered software released into the public domain.

The following applies to the NIST statistical code:

This code was accessed on 2019-01-05 during the US government shutdown, while
the NIST government servers were unavailable. Code was located on the Internet
Archive at:

<https://web.archive.org/web/20180720195312/https://csrc.nist.gov/projects/random-bit-generation/documentation-and-software>

which includes the following notice:

> **Software disclaimer:** "This software was developed at the National
> Institute of Standards and Technology by employees of the Federal Government
> in the course of their official duties. Pursuant to title 17 Section 105 of
> the United States Code this software is not subject to copyright protection
> and is in the public domain. The NIST Statistical Test Suite is an
> experimental system. NIST assumes no responsibility whatsoever for its use by
> other parties, and makes no guarantees, expressed or implied, about its
> quality, reliability, or any other characteristic. We would appreciate
> acknowledgment if the software is used."

*End of NIST statement*

### 3. ZOSLIB

**License:** Apache-2.0

ZOSLIB is a z/OS C/C++ library, available at
<https://github.com/ibmruntimes/zoslib>. It is an extended implementation of
the z/OS LE C Runtime Library.
