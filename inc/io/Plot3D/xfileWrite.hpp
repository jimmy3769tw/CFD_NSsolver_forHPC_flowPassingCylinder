#pragma once

#include"0_General.hpp"


auto OutputPlot3D_Xfile(
    simpulationVariable& simu,
    grid& gridA
){

	auto [nxCal, nyCal, nzCal] = gridA.nxyzCal;

    vector<vector<vector<float> > > Xout (nxCal,vector<vector<float> >(nyCal,vector<float> (nzCal)) );
    vector<vector<vector<float> > > Yout (nxCal,vector<vector<float> >(nyCal,vector<float> (nzCal)) );
    vector<vector<vector<float> > > Zout (nxCal,vector<vector<float> >(nyCal,vector<float> (nzCal)) );

	int Nblock = 1,
		tempNX = nxCal,
		tempNY = nyCal,
		tempNZ = nzCal;


    for(size_t i = 0; i < nxCal; ++i )
    for(size_t j = 0; j < nyCal; ++j )
    for(size_t k = 0; k < nzCal; ++k )
    {
        Xout[i][j][k] = gridA.Xc.at(i);
        Yout[i][j][k] = gridA.Yc.at(j);
        Zout[i][j][k] = gridA.Zc.at(k);
    }

	FILE *fptr;

	if(opendir("../mx_out") == NULL)
	{ if (system("mkdir ../mx_out") != 0){ return 1; } }

	fptr = fopen("mx_out/P3D.x","wb");

	fwrite(&Nblock, sizeof(int), 1,fptr);
	fwrite(&tempNX, sizeof(int), 1,fptr);
	fwrite(&tempNY, sizeof(int), 1,fptr);
	fwrite(&tempNZ, sizeof(int), 1,fptr);


    for(size_t k = 0; k < nzCal; ++k )
    for(size_t j = 0; j < nyCal; ++j )
    for(size_t i = 0; i < nxCal; ++i )
	    fwrite(&Xout[i][j][k],sizeof(float),1,fptr);

    for(size_t k = 0; k < nzCal; ++k )
    for(size_t j = 0; j < nyCal; ++j )
    for(size_t i = 0; i < nxCal; ++i )
	    fwrite(&Yout[i][j][k],sizeof(float),1,fptr);

    for(size_t k = 0; k < nzCal; ++k )
    for(size_t j = 0; j < nyCal; ++j )
    for(size_t i = 0; i < nxCal; ++i )
	    fwrite(&Zout[i][j][k],sizeof(float),1,fptr);



    fclose(fptr);
    return 0; 
}
