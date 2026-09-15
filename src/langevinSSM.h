#ifndef LANGEVIN_SSM_H
#define LANGEVIN_SSM_H

// Only load x86 SIMD intrinsics if we are compiling for an x86/x64 architecture
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <immintrin.h>
#endif

#endif
