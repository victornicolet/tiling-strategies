#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdlib.h>

// Dimanesions
#define SMALL (1 << 5)
#define DEFAULT ( 1 << 10)
#define BIG (1 << 13)
#define HUGE (1 << 16)

int main(int argc, char ** argv){
  return 0;
}

void jacobi1d_kernel(){
  int i,t;
  static int t_max, i_max;

  t_max = DEFAULT;
  i_max = BIG;

  double * m1[i_max];
  double * m0[i_max];

/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5;
 register int lbv, ubv;
/* Start of CLooG code */
if ((i_max >= 1) && (t_max >= 1)) {
  for (t1=0;t1<=floord(t_max-1,32);t1++) {
    for (t2=2*t1;t2<=min(floord(64*t1+i_max+62,32),floord(2*t_max+i_max-2,32));t2++) {
      if (t1 <= floord(32*t2-i_max,64)) {
        if (i_max%2 == 0) {
          m0[(i_max-1)] = m1[(i_max-1)];;
        }
      }
      if (i_max == 1) {
        for (t3=16*t2;t3<=min(16*t2+15,t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          m0[0] = m1[0];;
        }
      }
      for (t3=max(ceild(32*t2-i_max+1,2),32*t1);t3<=min(min(min(floord(32*t2-i_max+31,2),32*t1+31),16*t2-1),t_max-1);t3++) {
        for (t4=32*t2;t4<=2*t3+i_max-1;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
        m0[(i_max-1)] = m1[(i_max-1)];;
      }
      for (t3=max(ceild(32*t2-i_max+32,2),32*t1);t3<=min(min(32*t1+31,16*t2-1),t_max-1);t3++) {
        for (t4=32*t2;t4<=32*t2+31;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
      }
      if (i_max >= 2) {
        for (t3=16*t2;t3<=min(min(floord(32*t2-i_max+31,2),32*t1+31),t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          for (t4=2*t3+1;t4<=2*t3+i_max-1;t4++) {
            m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
          m0[(i_max-1)] = m1[(i_max-1)];;
        }
      }
      for (t3=max(ceild(32*t2-i_max+32,2),16*t2);t3<=min(min(32*t1+31,16*t2+15),t_max-1);t3++) {
        m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
        for (t4=2*t3+1;t4<=32*t2+31;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
      }
    }
  }
}
/* End of CLooG code */
}
