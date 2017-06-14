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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum { false, true } bool;

typedef struct {
  char *a;
  unsigned int alen;
  char *b;
  unsigned int blen;
} seq_pair_t;

// An entry is a cell in the matrix.
// prev holds the coordinates of the previous cell in the matrix.
typedef struct {
  double score;
  unsigned int prev[2];
} entry_t;

typedef struct {
  unsigned int m;
  unsigned int n;
  entry_t **mat;
} matrix_t;

typedef struct {
  seq_pair_t *seqs;
  int start_a;
  int start_b;
  int end_a;
  int end_b;
  int matches;
  double score;
} align_t;

static char* reverse(char *str);

static char get_char_comp(char c);

char* revcomp(char *str);

static align_t *traceback(seq_pair_t *problem, matrix_t *S, bool local);

static matrix_t *create_matrix(unsigned int m, unsigned int n);

void destroy_matrix(matrix_t *S);

void print_matrix(matrix_t *matrix, seq_pair_t *seq_pair);

void destroy_seq_pair(seq_pair_t *pair);

align_t *smith_waterman(seq_pair_t *problem, bool local);

void print_alignment(align_t *result, int target_len, int query_len);
