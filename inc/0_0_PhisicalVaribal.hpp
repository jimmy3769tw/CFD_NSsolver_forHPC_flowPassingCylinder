#pragma once
#include <vector>
#include <string>

class velocity{

	public:
		// -------------------------------
		std::vector <double> u, v, w,
							 su, sv, sw,
							 Viseff, c ;
		// -------------------------------
		

		bool resize_s(size_t n){

			su.resize(n);

			sv.resize(n);

			sw.resize(n);			

			return true;
		}


		bool resize(size_t n){

			u.resize(n);
			
			v.resize(n);
			
			w.resize(n);
			
			su.resize(1);		
			
			sv.resize(1);
			
			sw.resize(1);		

			return true;
		}

		bool iniU(double ui, double vi, double wi){
			std::fill(u.begin(), u.end(), ui);
			std::fill(v.begin(), v.end(), vi);
			std::fill(w.begin(), w.end(), wi);
			return true;
		}


		inline void fill_(std::vector<double> &v, double val){

			#pragma omp for schedule(static)
			for(size_t i = 0; i < v.size() ; ++i){
				v[i] = val;
			}
		}

		bool iniU_omp(double ui, double vi, double wi){

				fill_(u,ui);
				fill_(v,vi);
				fill_(w,wi);
			return true;
		}



	private:


};

class pressure{
	public:
		std::vector<double> p;
		std::vector<double> T;

		bool init_p(double pi){

			std::fill(p.begin(), p.end(), pi);
			return true;
		}


	private:


};

struct DfibArray{
    std::vector<double> f, eta, 
		Val_sumz, Val_sumz_sumy, ValSum;

    double cylinderDimension;
	std::vector<double> cylinderCenter;
};


struct shareMenory{

	double pChangeMax;
	
	double Out;

	// * Check Max val 
	double uDif_Max;
	
	double vDif_Max;
	
	double wDif_Max;
};
