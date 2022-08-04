#pragma once
#include <vector>
#include <cmath>


namespace solver{
    using namespace std;
    using namespace math;

    template<typename matrixT>
    class bicgstabRe2
    {
        public:

        bicgstabRe2(){ }

        bicgstabRe2(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~bicgstabRe2(){ }

        // ! init
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();
            r0_.resize(length_, 0.0);
            r1_.resize(length_, 0.0);
            p1_.resize(length_, 0.0);
            s1_.resize(length_, 0.0);
            ap_.resize(length_, 0.0);
            as_.resize(length_, 0.0);
            // PA_.resize(length_, 0.0);
            px1_.resize(length_, 0.0);
            px2_.resize(length_, 0.0);
            return true;
        }

        void setTolerance(double tolerance) { pResNormTol_ = tolerance;}

        std::pair<int, double> solve(const vector<double> & rhs, vector<double> & x); 

        private:
            // ----------------------
            int length_;
            // ----------------------

            int iters_max_{5000};

            // ----------------------
            double pResNormTol_      = 1e-3;
            double // residual criteria for P
                   pRestartFactor_   = 1e-7; 
            // ----------------------

            vector<double> r0_;
            vector<double> r1_;
            vector<double> p1_;
            vector<double> s1_;
            vector<double> ap_;
            vector<double> as_;
            // vector<double> PA_;

            // ----------------------------
            vector<double> px1_;
            vector<double> px2_;
            // ----------------------------

                
            // ----------------------------
            int timestep_ = 0;
            int pIterTotal_ = 0, 
                pResetMembers_ = 0,
                pResets_ = 0,
                pRestarts_  = 0;
            // ----------------------------

            matrixT lhs_mat_;
    };


    // ! main
    template<typename matrixT>
    inline std::pair<int, double> 
    bicgstabRe2<matrixT>::solve(
        const vector<double> & rhs, vector<double> & x)
    {

        timestep_ ++;
        double
            r1r0 = 0, r1r0l = 0, a = 0, w = 0, b = 0, norm0 = L2Norm(rhs);

        lhs_mat_.multiply_omp(x, r0_);
        // ---------------------------
        copy(px1_, px2_);
        copy(x, px1_);

        #pragma omp parallel for default(none) shared(rhs, r0_)
        for (int i=0; i<length_; ++i)
        { r0_[i] = rhs[i] - r0_[i]; }

        copy(r0_, r1_);
        copy(r0_, p1_);
        // ---------------------------

        r1r0 = dot(r1_, r0_);

        int pIter = 0;

        while (1) {

            lhs_mat_.multiply_omp(p1_, ap_);

            a  = r1r0 / dot(ap_, r0_);

            for (int i=0; i<length_; ++i) 
            { s1_[i] = r1_[i] - a*ap_[i]; }

            lhs_mat_.multiply_omp(s1_, as_);
        
            w = dot(as_, s1_) / dot(as_, as_); 

            for (int i=0; i<length_; ++i)
            {
                x[i] += a*p1_[i] + w*s1_[i];
                r1_[i] = s1_[i] - w*as_[i];
            }

            // Check for terminating conditions
            if ((L2Norm(r1_) / norm0) < pResNormTol_) { break; }


            // ! ## Check for reset of bCGSTAB
            if (pIter >= iters_max_) {
                int checkReset = 0;

                for (int i=0; i<length_; ++i){
                    if ( std::isnan(x[i]) || std::isinf(x[i]) ) 
                    { x[i] = px1_[i], ++checkReset;}
                }

                if (checkReset > length_*0.1) {

                    for (int i=0; i<length_; ++i) 
                    { x[i] = px2_[i]; }

                    pResetMembers_ += length_, ++pResets_;

                } 
                else if (checkReset != 0)
                { pResetMembers_ += checkReset, ++pResets_; }

                break;
            }
            // *------------------------------------

            r1r0l = r1r0;
            r1r0 = dot(r1_, r0_);
            b  = (a / w) * (r1r0 / r1r0l);

            // Check rho for restart
            if (timestep_ < 10 || abs(r1r0l) > pResNormTol_*pRestartFactor_)
            {
                for (int i=0; i<length_; ++i) 
                { p1_[i] = r1_[i] + b*(p1_[i] - w*ap_[i]); }
            }
            else
            {
                for (int i = 0; i < length_; ++i)
                {
                    r0_[i] = r1_[i];
                    p1_[i] = r1_[i];
                }
                ++pRestarts_;
            }
            ++pIter;
            ++pIterTotal_;
        }

        return std::make_pair( pIter, L2Norm(r1_) / norm0);
    }
}


