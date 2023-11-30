// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef A_TRANSPOSE_INCLUDED
#define A_TRANSPOSE_INCLUDED

#include "parallel.h"

#define _TRANS_THRESHHOLD 64

template <class E, class intV> struct transpose {
  E *A, *B;
  transpose(E *AA, E *BB) : A(AA), B(BB) {}

  void transR(intV rStart, intV rCount, intV rLength, intV cStart, intV cCount,
              intV cLength) {
    // cout << "cc,rc: " << cCount << "," << rCount << endl;
    if (cCount < _TRANS_THRESHHOLD && rCount < _TRANS_THRESHHOLD) {
      for (intV i = rStart; i < rStart + rCount; i++)
        for (intV j = cStart; j < cStart + cCount; j++)
          B[j * cLength + i] = A[i * rLength + j];
    } else if (cCount > rCount) {
      intV l1 = cCount / 2;
      intV l2 = cCount - cCount / 2;
      cilk_spawn this->transR(rStart, rCount, rLength, cStart, l1, cLength);
      transR(rStart, rCount, rLength, cStart + l1, l2, cLength);
      cilk_sync;
    } else {
      intV l1 = rCount / 2;
      intV l2 = rCount - rCount / 2;
      cilk_spawn this->transR(rStart, l1, rLength, cStart, cCount, cLength);
      transR(rStart + l1, l2, rLength, cStart, cCount, cLength);
      cilk_sync;
    }
  }

  void trans(intV rCount, intV cCount) {
    transR(0, rCount, cCount, 0, cCount, rCount);
  }
};

template <class E, class intV> struct blockTrans {
  E *A, *B;
  intV *OA, *OB, *L;

  blockTrans(E *AA, E *BB, intV *OOA, intV *OOB, intV *LL)
      : A(AA), B(BB), OA(OOA), OB(OOB), L(LL) {}

  void transR(intV rStart, intV rCount, intV rLength, intV cStart, intV cCount,
              intV cLength) {
    // cout << "cc,rc: " << cCount << "," << rCount << endl;
    if (cCount < _TRANS_THRESHHOLD && rCount < _TRANS_THRESHHOLD) {
      for (intV i = rStart; i < rStart + rCount; i++)
        for (intV j = cStart; j < cStart + cCount; j++) {
          E *pa = A + OA[i * rLength + j];
          E *pb = B + OB[j * cLength + i];
          intV l = L[i * rLength + j];
          // cout << "pa,pb,l: " << pa << "," << pb << "," << l << endl;
          for (intV k = 0; k < l; k++)
            *(pb++) = *(pa++);
        }
    } else if (cCount > rCount) {
      intV l1 = cCount / 2;
      intV l2 = cCount - cCount / 2;
      cilk_spawn this->transR(rStart, rCount, rLength, cStart, l1, cLength);
      transR(rStart, rCount, rLength, cStart + l1, l2, cLength);
      cilk_sync;
    } else {
      intV l1 = rCount / 2;
      intV l2 = rCount - rCount / 2;
      cilk_spawn this->transR(rStart, l1, rLength, cStart, cCount, cLength);
      transR(rStart + l1, l2, rLength, cStart, cCount, cLength);
      cilk_sync;
    }
  }

  void trans(intV rCount, intV cCount) {
    transR(0, rCount, cCount, 0, cCount, rCount);
  }
};

#endif // A_TRANSPOSE_INCLUDED
