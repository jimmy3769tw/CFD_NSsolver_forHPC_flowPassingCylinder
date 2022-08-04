
#pragma once 

#include <vector>
#include <tuple>
#include <omp.h>

#include<iostream>
namespace math{

    using namespace std;
    // --------------------------
    inline double 
    dot(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};

        #pragma omp parallel reduction(+:r) default(none) shared(a, b)
        for(size_t i=0; i<a.size() ; ++i)
        {  r += a[i] * b[i]; }

        return r;
    }
    // --------------------------

    // --------------------------
    inline double 
    L2Norm(const std::vector<double> & a)
    {
        double r{0.0};

        for(size_t i=0; i<a.size(); ++i)
        {  r += a[i] * a[i]; }

        return sqrt( r );
    }
    // --------------------------


    // --------------------------
    inline void 
    zero(vector<double> & a){
        #pragma omp parallel  default(none) shared(a)
        {
            #pragma omp for 
            for (size_t i = 0; i < a.size() ; ++i)
            { a[i] = 0.0; }
        }
    }
    // --------------------------



    // --------------------------
    inline vector<double> 
    copy(vector<double> & a, vector <double> & r)
    {
        #pragma omp parallel default(none) shared(a,r)
        {
            #pragma omp for 
            for (size_t i = 0; i < a.size()  ; ++i)
            {  r[i] = a[i];}
        }
        return r;
    }
    // --------------------------


    // --------------------------
    inline double 
    getError(vector<double> &a, vector<double> &b){
        double r = 0.0;
        #pragma omp parallel default(none) shared(a,b) reduction(+:r)
        {
            #pragma omp for 
            for (size_t i = 0; i < a.size() ; ++i)
                r += std::abs(a[i]-b[i]);
        }
        return r;
    }
    // --------------------------
}