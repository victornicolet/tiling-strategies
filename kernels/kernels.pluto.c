#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdlib.h>

// Dimensions
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
  int t1, t2, t3, t4;
 register int lbv, ubv;
/* Start of CLooG code */
if ((i_max >= 1) && (t_max >= 1)) {
  for (t1=ceild(-i_max-30,32);t1<=floord(t_max-1,16);t1++) {
    for (t2=max(t1,-t1-1);t2<=min(min(floord(-8*t1+t_max-1,8),floord(16*t1+i_max+14,16)),floord(2*t_max+i_max-2,32));t2++) {
      if ((t1 <= floord(16*t2-i_max,16)) && (t2 >= ceild(i_max,32))) {
        if (i_max%2 == 0) {
          m0[(i_max-1)] = m1[(i_max-1)];;
        }
      }
      if ((i_max >= 2) && (t1 == t2)) {
        m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
        for (t4=32*t1+1;t4<=32*t1+2;t4++) {
          m0[(-32*t1+t4-1)] = m1[(-32*t1+t4-1)];;
        }
      }
      if ((i_max == 1) && (t1 == t2)) {
        m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
        m0[0] = m1[0];;
      }
      for (t3=max(max(0,8*t1+8*t2),16*t1+8);t3<=min(min(floord(32*t1+i_max-1,2),t_max-1),8*t1+8*t2+7);t3++) {
        for (t4=32*t2;t4<=-32*t1+4*t3;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
        for (t4=-32*t1+4*t3+1;t4<=min(2*t3+i_max,-32*t1+4*t3+2);t4++) {
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
      }
      if (t1 == t2) {
        for (t3=16*t1+1;t3<=min(min(floord(32*t1+i_max-1,2),16*t1+7),t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          for (t4=2*t3+1;t4<=-32*t1+4*t3;t4++) {
            m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
          for (t4=-32*t1+4*t3+1;t4<=min(2*t3+i_max,-32*t1+4*t3+2);t4++) {
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
        }
      }
      if (t1 == t2) {
        for (t3=16*t1+8;t3<=min(min(floord(32*t1+i_max-1,2),16*t1+15),t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          for (t4=2*t3+1;t4<=32*t1+31;t4++) {
            m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
        }
      }
      if ((i_max == 1) && (t1 == t2)) {
        for (t3=16*t1+1;t3<=min(16*t1+15,t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          m0[0] = m1[0];;
        }
      }
      for (t3=max(max(max(0,ceild(32*t1+i_max,2)),ceild(32*t1-i_max+33,2)),ceild(32*t2-i_max+1,2));t3<=min(t_max-1,8*t1+8*t2+7);t3++) {
        for (t4=32*t2;t4<=2*t3+i_max-1;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
        m0[(i_max-1)] = m1[(i_max-1)];;
      }
      for (t3=max(16*t1+16,8*t1+8*t2+8);t3<=min(min(floord(32*t1+i_max+28,2),floord(32*t2-i_max+31,2)),t_max-1);t3++) {
        for (t4=-32*t1+4*t3-31;t4<=-32*t1+4*t3-30;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
        }
        for (t4=-32*t1+4*t3-29;t4<=2*t3+i_max-1;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
        m0[(i_max-1)] = m1[(i_max-1)];;
      }
      for (t3=max(max(ceild(32*t2-i_max+32,2),16*t1+16),8*t1+8*t2+8);t3<=min(t_max-1,8*t1+8*t2+15);t3++) {
        for (t4=-32*t1+4*t3-31;t4<=-32*t1+4*t3-30;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
        }
        for (t4=-32*t1+4*t3-29;t4<=32*t2+31;t4++) {
          m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
          m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
        }
      }
      if ((i_max >= 3) && (t1 <= min(floord(16*t2-i_max+1,16),floord(2*t_max-i_max-31,32)))) {
        if ((i_max+1)%2 == 0) {
          for (t4=32*t1+2*i_max+27;t4<=32*t1+2*i_max+28;t4++) {
            m1[(-32*t1+t4-i_max-29)] = (m0[(-32*t1+t4-i_max-29) - 1] + m0[(-32*t1+t4-i_max-29)] + m0[(-32*t1+t4-i_max-29) + 1]) / 3.0;;
          }
          m0[(i_max-1)] = m1[(i_max-1)];;
        }
      }
      if ((i_max >= 2) && (t1 == t2)) {
        for (t3=ceild(32*t1+i_max,2);t3<=min(floord(32*t1-i_max+31,2),t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          for (t4=2*t3+1;t4<=2*t3+i_max-1;t4++) {
            m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
          m0[(i_max-1)] = m1[(i_max-1)];;
        }
      }
      if (t1 == t2) {
        for (t3=max(ceild(32*t1+i_max,2),ceild(32*t1-i_max+32,2));t3<=min(16*t1+15,t_max-1);t3++) {
          m1[0] = (m0[0 - 1] + m0[0] + m0[0 + 1]) / 3.0;;
          for (t4=2*t3+1;t4<=32*t1+31;t4++) {
            m1[(-2*t3+t4)] = (m0[(-2*t3+t4) - 1] + m0[(-2*t3+t4)] + m0[(-2*t3+t4) + 1]) / 3.0;;
            m0[(-2*t3+t4-1)] = m1[(-2*t3+t4-1)];;
          }
        }
      }
      if (t1 <= min(floord(16*t2-i_max,16),floord(2*t_max-i_max-32,32))) {
        if (i_max%2 == 0) {
          m1[(i_max-1)] = (m0[(i_max-1) - 1] + m0[(i_max-1)] + m0[(i_max-1) + 1]) / 3.0;;
        }
      }
    }
  }
}
/* End of CLooG code */
}

void jacobi2d_kernel(){
  int i, j, k, l, t;
  int t_max, i_max, j_max;

  t_max = DEFAULT;
  i_max = BIG;
  j_max = BIG;

  double ** m1;
  double ** m0;

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
  int t1, t2, t3, t4, t5, t6;
 register int lbv, ubv;
/* Start of CLooG code */
if ((i_max >= 3) && (t_max >= 1)) {
  for (t1=0;t1<=floord(t_max-1,32);t1++) {
    for (t2=2*t1;;t2++) {
      if ((j <= j_max-2) && (j_max >= 3)) {
        for (t3=32*t1;t3<=min(min(floord(32*t2-i_max,2),32*t1+31),t_max-1);t3++) {
          for (t4=32*t2;t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max <= 2)) {
        for (t3=32*t1;t3<=min(min(32*t1+31,16*t2+15),t_max-1);t3++) {
          for (t4=max(32*t2,2*t3+1);t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
          }
        }
      }
      if ((j >= j_max-1) && (j_max >= 3)) {
        for (t3=max(ceild(32*t2-i_max+1,2),32*t1);t3<=min(min(32*t1+31,16*t2+14),t_max-1);t3++) {
          for (t4=max(32*t2,2*t3+2);t4<=min(32*t2+31,2*t3+i_max-1);t4++) {
            lbv=1;
            ubv=j_max-2;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              m0[(-2*t3+t4-1)][t6] = m1[(-2*t3+t4-1)][t6];;
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max >= 3)) {
        for (t3=max(ceild(32*t2-i_max+1,2),32*t1);t3<=min(min(min(floord(32*t2-i_max+31,2),32*t1+31),16*t2-1),t_max-1);t3++) {
          for (t4=32*t2;t4<=2*t3+i_max-1;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
            lbv=1;
            ubv=j_max-2;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              m0[(-2*t3+t4-1)][t6] = m1[(-2*t3+t4-1)][t6];;
            }
          }
          for (t4=2*t3+i_max;t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max >= 3)) {
        for (t3=max(ceild(32*t2-i_max+32,2),32*t1);t3<=min(min(32*t1+31,16*t2-1),t_max-1);t3++) {
          for (t4=32*t2;t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
            lbv=1;
            ubv=j_max-2;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              m0[(-2*t3+t4-1)][t6] = m1[(-2*t3+t4-1)][t6];;
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max >= 3)) {
        for (t3=16*t2;t3<=min(min(floord(32*t2-i_max+31,2),32*t1+31),t_max-1);t3++) {
          for (t6=1;t6<=i_max-2;t6++) {
            S1(t2,t1,t3,t6,1);
          }
          for (t4=2*t3+2;t4<=2*t3+i_max-1;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
            lbv=1;
            ubv=j_max-2;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              m0[(-2*t3+t4-1)][t6] = m1[(-2*t3+t4-1)][t6];;
            }
          }
          for (t4=2*t3+i_max;t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max >= 3)) {
        for (t3=max(ceild(32*t2-i_max+32,2),16*t2);t3<=min(min(32*t1+31,16*t2+14),t_max-1);t3++) {
          for (t6=1;t6<=i_max-2;t6++) {
            S1(t2,t1,t3,t6,1);
          }
          for (t4=2*t3+2;t4<=32*t2+31;t4++) {
            for (t6=1;t6<=i_max-2;t6++) {
              S1(t2,t1,t3,t6,(-2*t3+t4));
            }
            lbv=1;
            ubv=j_max-2;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              m0[(-2*t3+t4-1)][t6] = m1[(-2*t3+t4-1)][t6];;
            }
          }
        }
      }
      if ((j <= j_max-2) && (j_max >= 3) && (t1 >= ceild(t2-1,2)) && (t2 <= floord(t_max-16,16))) {
        for (t6=1;t6<=i_max-2;t6++) {
          S1(t2,t1,(16*t2+15),t6,1);
        }
      }
    }
  }
}
/* End of CLooG code */
}
