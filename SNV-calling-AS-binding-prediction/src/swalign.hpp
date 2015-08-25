/*
 *  Copyright (c) 2010 Nicolaus Lance Hepler
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cfloat>
#include <iostream>

using namespace std;

#define GAP -1.0
#define MATCH 2.0
#define MISMATCH -0.5

typedef struct {
  string a;
  string b;
  double score;
} seq_pair;

typedef struct {
  double score;
  unsigned int prev[2];
} entry;

typedef struct {
  int m;
  int n;
  vector<vector<entry> > mat;
} matrix;

//static void reverse(string &str);

//static void traceback(seq_pair &problem, matrix &S, bool local, seq_pair &result);

//static void create_matrix(int m, int n, matrix &S);

void destroy_matrix(matrix &S);

void destroy_seq_pair(seq_pair &pair);

void smith_waterman(seq_pair &problem, bool local, seq_pair &result);
