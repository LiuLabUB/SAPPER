/*
 *  Copyright (c) 2010 Nicolaus Lance Hepler
 *  
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
// Note: That's an MIT license.
// All double-commented comments below are from Nicolaus Lance Hepler.
// Original repository: https://code.google.com/archive/p/swalign/

#include "swalign.h"

#define GAP -1.0
#define MATCH 2.0
#define MISMATCH -0.5
//             ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz
#define TRANS "TVGHEFCDIJMLKNOPQYWAABSXRZ[\\]^_`tvghefcdijmlknopqywaabsxrz"
#define TRANS_OFFSET 65


// /* reverse a string in place, return str */
static char* reverse(char *str) {
  char *left  = str;
  char *right = left + strlen(str) - 1;
  char tmp;

  while (left < right) {
    tmp        = *left;
    *(left++)  = *right;
    *(right--) = tmp;
  }

  return str;
}

// Return the reverse complement of a sequence.
char* revcomp(char *str) {
  char *left = str;
  char *right = left + strlen(str) - 1;
  char tmp;

  while (left < right) {
    tmp        = get_char_comp(*left);
    *(left++)  = get_char_comp(*right);
    *(right--) = tmp;
  }

  return str;
}

// Return the complement of a base.
// Uses a simple lookup table: a string with the complements of all possible sequence characters.
static char get_char_comp(char c) {
  int i = c - TRANS_OFFSET;
  if (i < 0 || i > 57) {
    return c;
  } else {
    return TRANS[i];
  }
}

// // works globally
// Note: Currently the "local" flag isn't functional. It seems to always do a local alignment.
static align_t *traceback(seq_pair_t *problem, matrix_t *S, bool local) {
  align_t *result = malloc(sizeof(align_t));
  seq_pair_t *seqs = malloc(sizeof(seq_pair_t));
  unsigned int i    = S->m - 1;
  unsigned int j    = S->n - 1;
  unsigned int k    = 0;
  // Create output strings. Allocate maximum potential length.
  char c[S->m + S->n + 1];
  char d[S->m + S->n + 1];

  memset(c, '\0', sizeof(c));
  memset(d, '\0', sizeof(d));

  // This wasn't finished by NLH. Not functioning correctly yet.
  // It seems the purpose is to start the traceback from the place where the score reaches its
  // maximum instead of the very end (set i and j to those coordinates).
  if (local == true) {
    unsigned int l, m;
    double max = FLT_MIN;

    for (l = 0; l < S->m; l++) {
      for (m = 0; m < S->n; m++) {
        if (S->mat[l][m].score > max) {
          i = l;
          j = m;
          max = S->mat[l][m].score;
        } 
      } 
    }
  }

  double score = DBL_MIN;
  int matches = 0;
  int start_a = 0;
  int start_b = 0;
  int end_a = 0;
  int end_b = 0;
  bool move_i = false;
  bool move_j = false;
  // Walk back through the matrix from the end, taking the path determined by the "prev" values of
  // each cell. Assemble the sequence along the way.
  if (S->mat[i][j].prev[0] != 0 && S->mat[i][j].prev[1] != 0) {
    while (i > 0 || j > 0) {
      unsigned int new_i = S->mat[i][j].prev[0];
      unsigned int new_j = S->mat[i][j].prev[1];
  
      // If we've moved in the i axis, add the new base to the sequence. Otherwise, it's a gap.
      if (new_i < i) {
        *(c+k) = *(problem->a+i-1);
        move_i = true;
      } else {
        *(c+k) = '-';
        move_i = false;
      }
  
      // If we've moved in the j axis, add the new base to the sequence. Otherwise, it's a gap.
      if (new_j < j) {
        *(d+k) = *(problem->b+j-1);
        move_j = true;
      } else {
        *(d+k) = '-';
        move_j = false;
      }

      if (S->mat[i][j].score > score) {
        score = S->mat[i][j].score;
      }

      if (move_i && move_j) {
        if (*(c+k) == *(d+k)) {
          matches++;
        }
        // Start and end of each sequence are in the coordinates of the other sequence.
        start_a = new_j + 1;
        start_b = new_i + 1;
        if (! end_a) {
          end_a = new_j + 1;
        }
        if (! end_b) {
          end_b = new_i + 1;
        }
      }
  
      k++;
  
      i = new_i;
      j = new_j;
    }
  }

  seqs->a = malloc(sizeof(char) * k + 1);
  seqs->b = malloc(sizeof(char) * k + 1);

  memset(seqs->a, '\0', sizeof(seqs->a));
  memset(seqs->b, '\0', sizeof(seqs->b));

  reverse(c);
  reverse(d);

  strcpy(seqs->a, c);
  strcpy(seqs->b, d);

  seqs->alen = k;
  seqs->blen = k;

  result->seqs = seqs;
  result->score = score;
  result->matches = matches;
  result->start_a = start_a;
  result->start_b = start_b;
  result->end_a = end_a;
  result->end_b = end_b;

  return result; 
}

static matrix_t *create_matrix(unsigned int m, unsigned int n) {
  matrix_t *S = malloc(sizeof(matrix_t));
  unsigned int i;

  S->m = m;
  S->n = n;

  S->mat = malloc(sizeof(entry_t) * m * n);

  for (i = 0; i < m; i++) {
    S->mat[i] = malloc(sizeof(entry_t) * n);
  }

  return S;
}

void destroy_matrix(matrix_t *S) {
  unsigned int i;
  for (i = 0; i < S->m; i++) {
    free(S->mat[i]);
  }
  free(S->mat);
  free(S);
  return;
}

// Print a visual representation of the path through the matrix.
void print_matrix(matrix_t *matrix, seq_pair_t *seq_pair) {
  int i, j;
  for (i = 0; i < matrix->m; i++) {
    if (i == 0) {
      printf("\t\t");
      for (j = 0; j < seq_pair->blen; j++) {
        printf("%c\t", seq_pair->b[j]);
      }
      printf("\n");
      printf("        ");
      for (j = 0; j < matrix->n; j++) {
        printf("%d\t", j);
      }
      printf("\n");
    }
    if (i == 0) {
      printf("     0  ");
    } else {
      printf("%c %4d  ", seq_pair->a[i-1], i);
    }
    for (j = 0; j < matrix->n; j++) {
      printf("%d,%d|%0.0f\t", matrix->mat[i][j].prev[0], matrix->mat[i][j].prev[1], matrix->mat[i][j].score);
    }
    printf("\n");
  }
}

void destroy_seq_pair(seq_pair_t *pair) {
  free(pair->a);
  free(pair->b);
  free(pair);
  return;
}

align_t *smith_waterman(seq_pair_t *problem, bool local) {
  unsigned int m = problem->alen + 1;
  unsigned int n = problem->blen + 1;
  matrix_t *S = create_matrix(m, n);
  align_t *result;
  unsigned int i, j, k, l;

  S->mat[0][0].score   = 0;
  S->mat[0][0].prev[0] = 0;
  S->mat[0][0].prev[1] = 0;

  for (i = 1; i <= problem->alen; i++) {
    S->mat[i][0].score   = 0.0;
    S->mat[i][0].prev[0] = i-1;
    S->mat[i][0].prev[1] = 0;
  }

  for (j = 1; j <= problem->blen; j++) {
    S->mat[0][j].score   = 0.0;
    S->mat[0][j].prev[0] = 0;
    S->mat[0][j].prev[1] = j-1;
  }

  for (i = 1; i <= problem->alen; i++) {
    for (j = 1; j <= problem->blen; j++) {
      int nw_score = (strncmp(problem->a+(i-1), problem->b+(j-1), 1) == 0) ? MATCH : MISMATCH;

      S->mat[i][j].score   = DBL_MIN;
      S->mat[i][j].prev[0] = 0;
      S->mat[i][j].prev[1] = 0;

      for (k = 0; k <= 1; k++) {
        for (l = 0; l <= 1; l++) {
          int val;

          if (k == 0 && l == 0) {
            continue;
          } else if (k > 0 && l > 0) {
            val = nw_score; 
          } else if (k > 0 || l > 0) {
            if ((i == problem->alen && k == 0) ||
                (j == problem->blen && l == 0))
              val = 0.0;
            else
              val = GAP;
          } else {
            // do nothing..
          }

          val += S->mat[i-k][j-l].score;

          if (val > S->mat[i][j].score) {
            S->mat[i][j].score   = val;
            S->mat[i][j].prev[0] = i-k;
            S->mat[i][j].prev[1] = j-l;
          }
        }
      }
    }
  }

  result = traceback(problem, S, local);

  // print_matrix(S, problem);

  destroy_matrix(S);

  return result;
}

void print_alignment(align_t *result, int target_len, int query_len) {
  printf("Score: %0.0f  Matches: %d\n", result->score, result->matches);
  printf("Target: %3d %s %-3d\n", result->start_a, result->seqs->a, result->end_a);
  printf("Query:  %3d %s %-3d\n", result->start_b, result->seqs->b, result->end_b);
}

int main(int argc, const char **argv) {

  if (argc != 3) {
    printf("usage: swalign TARGET_SEQ QUERY_SEQ\n");
    exit(1);
  }

  {
    seq_pair_t problem;
    align_t *result;
    char c[strlen(argv[1])], d[strlen(argv[2])];
  
    strcpy(c, argv[1]);
    strcpy(d, argv[2]);
  
    problem.a = c;
    problem.alen = strlen(problem.a);
    problem.b = d;
    problem.blen = strlen(problem.b);
  
    result = smith_waterman(&problem, false);
  
    print_alignment(result, problem.alen, problem.blen);
  }

  exit(0);
} 



