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

#include "swalign.hpp"

/* reverse a string in place, return str */
static void reverse(string &str) {
  string stemp(str);
  str.clear();
  for(int i=stemp.size()-1;i>=0;i--) str.push_back(stemp[i]);
}

// works globally
static void traceback(seq_pair &problem, matrix &S, bool local, seq_pair &result) {
  int i    = S.m - 1;
  int j    = S.n - 1;
  int k    = 0;
  string c,d;
  unsigned long t;
  
  for(int m=0;i<S.m+S.n+1;m++) {c.push_back(' ');d.push_back(' ');}

  if (local == true) {
    int l, m;
    double max = FLT_MIN;

    for (l = 0; l < S.m; l++) {
      for (m = 0; m < S.m; m++) {
        if (S.mat[l][m].score > max) {
          i = l;
          j = m;
        } 
      } 
    }
  }

  if (S.mat[i][j].prev[0] != 0 && S.mat[i][j].prev[1] != 0) {
    while (i > 0 || j > 0) {
      //printf("%d %d %f\n", i, j, S.mat[i][j].score);

      int new_i = S.mat[i][j].prev[0];
      int new_j = S.mat[i][j].prev[1];
  
      if (new_i < i)
        c[k] = problem.a[i-1];
      else
    	c[k] = '-';
  
      if (new_j < j)
    	d[k] = problem.b[j-1];
      else
    	d[k] = '-';
  
      k++;
  
      i = new_i;
      j = new_j; 
    }
  }


  reverse(c);
  reverse(d);

  for(t=0;t<c.size();t++)
  {
	  if(c[t]!=' ') result.a.push_back(c[t]);
  }
  for(t=0;t<d.size();t++)
  {
  	  if(d[t]!=' ') result.b.push_back(d[t]);
   }
}

static void create_matrix(int m, int n, matrix &S) {
   int i, j;
   entry myentry;

   S.m=m;
   S.n=n;

  S.mat.resize(m);
  for(i=0;i<m;i++)
  {
	  S.mat[i].resize(n);
	  for(j=0;j<n;j++) S.mat[i][j]=myentry;
  }
}

void destroy_matrix(matrix &S) {
  vector<vector<entry> > v;
  v.swap(S.mat);

  return; 
}

void destroy_seq_pair(seq_pair &pair) {
  pair.a="";
  pair.b="";

  return;
}

void smith_waterman(seq_pair &problem, bool local, seq_pair &result) {
  int m = problem.a.length() + 1;
  int n = problem.b.length() + 1;

  unsigned long i, j;
  int k, l;
  int val;
  int itemp=0.0;
  int nw_score;
  matrix S;

  create_matrix(m, n, S);


  S.mat[0][0].score   = 0;
  S.mat[0][0].prev[0] = 0;
  S.mat[0][0].prev[1] = 0;

  for (i = 1; i <= problem.a.length(); i++) {
    S.mat[i][0].score   = 0.0;
    S.mat[i][0].prev[0] = i-1;
    S.mat[i][0].prev[1] = 0;
  }

  for (j = 1; j <= problem.b.length(); j++) {
    S.mat[0][j].score   = 0.0;
    S.mat[0][j].prev[0] = 0;
    S.mat[0][j].prev[1] = j-1;
  }

  for (i = 1; i <= problem.a.length(); i++) {
    for (j = 1; j <= problem.b.length(); j++) {
      if((problem.a[i-1] == problem.b[j-1]) && (problem.a[i-1] != 'N')) nw_score=(int)MATCH;
      else nw_score=(int)MISMATCH;

      S.mat[i][j].score   = DBL_MIN;
      S.mat[i][j].prev[0] = 0;
      S.mat[i][j].prev[1] = 0;

      for (k = 0; k <= 1; k++) {
        for (l = 0; l <= 1; l++) {

          if (k == 0 && l == 0) {
            continue;
          } else if (k > 0 && l > 0) {
            val = nw_score; 
          } else if (k > 0 || l > 0) {
            if ((i == problem.a.length() && k == 0) ||
                (j == problem.b.length() && l == 0))
              val = 0.0;
            else
              val = GAP;
          } else {
            // do nothing..
          }

          val += S.mat[i-k][j-l].score;

          if (val > S.mat[i][j].score) {
            S.mat[i][j].score   = val;
            S.mat[i][j].prev[0] = i-k;
            S.mat[i][j].prev[1] = j-l;
          }
        }
      }

      itemp=S.mat[i][j].score;
    }
  }

  traceback(problem, S, local, result);
  result.score=itemp;

  destroy_matrix(S);
}


