#ifndef __TESTS_H
#define __TESTS_H

void test_jacobi_sequential(int N,double tol,int maxiter);
void test_jacobi_openmp(int N,double tol,int maxiter);
void test_gs_sequential(int N,double tol,int maxiter);
void test_jacobi_openmpft(int N,double tol,int maxiter);
void test_jacobi_openmpv2(int N,double tol,int maxiter);

#endif // __TESTS_H
