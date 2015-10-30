#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

//fisher exact p_value and OR
double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);


#endif
