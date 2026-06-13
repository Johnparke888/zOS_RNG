#pragma once
#ifndef __MVS__
// #error "This file targets z/OS USS only."
#define __ptr32
#endif
#include <errno.h>
#include <stddef.h>
#ifdef __MVS__
#include <builtins.h> /* __stckf */
#include <psa.h>
#endif

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
