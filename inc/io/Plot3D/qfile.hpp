#pragma once

#include <sys/types.h>
#include <dirent.h>

auto OutputPLOT3D_Qfile(
    DfibArray& Dfib,
    simpulationVariable& simu,
    pressure& t1,
    velocity& T3,
    divideLocal& Lo,
    grid& gA
)
{
	int icel;
	auto [nx, ny, nz , gC] = gA.nxyzgC;
	auto [nxCal, nyCal, nzCal] = gA.nxyzCal;

	// double temp = 1.0;
    float mach = 0;    
    float alpha = 0;    
    float reyn = simu.Re; 
    float time = simu.get_time();    



	int io = simu.currentFile;


	ofstream file;

	int Nblock = 1,
		tempNX = gA.nxCal,
		tempNY = gA.nyCal,
		tempNZ = gA.nzCal;



	if(opendir("../mx_out") == NULL)
	{ if (system("mkdir ../mx_out") != 0){ return 1; } }

	string filename = "mx_out/P3D";

	if (io < 10){	filename += "0";}

	filename += std::to_string(io);
	filename += ".q";	
	file.open(filename,ofstream::binary);
	cout << "(Output Plaot3D):" << io << endl;

	file.write((char *)(&Nblock), sizeof(int));  
	file.write((char *)(&tempNX), sizeof(int));
	file.write((char *)(&tempNY), sizeof(int));
	file.write((char *)(&tempNZ), sizeof(int));

	file.write((char *)(&mach), sizeof(float));
	file.write((char *)(&alpha), sizeof(float));
	file.write((char *)(&reyn), sizeof(float));
	file.write((char *)(&time), sizeof(float));

	auto ii = [&](auto i, auto j, auto k ){ return gA.icel(i,j,k);};

	for (size_t k = gC; k < nzCal+gC; ++k) 
	for (size_t j = gC; j < nyCal+gC; ++j) 
	for (size_t i = gC; i < nxCal+gC; ++i)
	{
		float out = t1.p[ii(i,j,k)];
		file.write((char *)(&out), sizeof(float));
	}

	for (size_t k = gC; k < nzCal+gC; ++k) 
	for (size_t j = gC; j < nyCal+gC; ++j) 
	for (size_t i = gC; i < nxCal+gC; ++i)
	{
		float out = 0.5 * (T3.u[ii(i-1,j,k)] + T3.u[ii(i,j,k)]);
		file.write((char *)(&out), sizeof(float));
	}


	for (size_t k = gC; k < nzCal+gC; ++k) 
	for (size_t j = gC; j < nyCal+gC; ++j) 
	for (size_t i = gC; i < nxCal+gC; ++i)
	{
		float out = 0.5 * (T3.v[ii(i,j-1,k)] + T3.v[ii(i,j,k)]);
		file.write((char *)(&out), sizeof(float));
	}

	for (size_t k = gC; k < nzCal+gC; ++k) 
	for (size_t j = gC; j < nyCal+gC; ++j) 
	for (size_t i = gC; i < nxCal+gC; ++i)
	{
		float out = 0.5 * (T3.w[ii(i,j,k-1)] + T3.w[ii(i,j,k)]);
		file.write((char *)(&out), sizeof(float));
	}


	for (size_t k = gC; k < nzCal+gC; ++k) 
	for (size_t j = gC; j < nyCal+gC; ++j) 
	for (size_t i = gC; i < nxCal+gC; ++i)
	{
		float out =Dfib.eta[ii(i,j,k)];;
		file.write((char *)(&out), sizeof(float));
	}


	return 0;

}



