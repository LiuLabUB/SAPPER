#include "math.hpp"

typedef struct {
    int n11, n1_, n_1, n;
    double p;
} hgacc_t;

// log\binom{n}{k}
static double lbinom(int n, int k)
{
    if (k == 0 || n == k) return 0;
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n
// hypergeometric distribution
static double hypergeo(int n11, int n1_, int n_1, int n)
{
    return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux)
{
    if (n1_ || n_1 || n) {
        aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
    } else { // then only n11 changed; the rest fixed
        if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
            if (n11 == aux->n11 + 1) { // incremental
                aux->p *= (double)(aux->n1_ - aux->n11) / n11
                    * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
                aux->n11 = n11;
                return aux->p;
            }
            if (n11 == aux->n11 - 1) { // incremental
                aux->p *= (double)aux->n11 / (aux->n1_ - n11)
                    * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
                aux->n11 = n11;
                return aux->p;
            }
        }
        aux->n11 = n11;
    }
    aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
    return aux->p;
}

double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two)
{
    int i, j, max, min;
    double p, q, left, right;
    hgacc_t aux;
    int n1_, n_1, n;

    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
    max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
    min = n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
    if (min < 0) min = 0; // min n11, for left tail
    *two = *_left = *_right = 1.;
    if (min == max) return 1.; // no need to do test
    q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
    // left tail
    p = hypergeo_acc(min, 0, 0, 0, &aux);
    for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
        left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
    --i;
    if (p < 1.00000001 * q) left += p;
    else --i;
    // right tail
    p = hypergeo_acc(max, 0, 0, 0, &aux);
    for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
        right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
    ++j;
    if (p < 1.00000001 * q) right += p;
    else ++j;
    // two-tail
    *two = left + right;
    if (*two > 1.) *two = 1.;
    // adjust left and right
    if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
    else left = 1.0 - right + q;
    *_left = left; *_right = right;
    return q;
}

