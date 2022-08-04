#pragma once

void OutputDataFile(
    DfibArray& Dfib,
    simpulationVariable& simu,
    pressure& t,
    velocity& T,
    divideLocal& Lo,
    grid& gridA	
){
	cout << "mx_out/DATA" << endl;
	auto [nx, ny, nz ] = gridA.nxyz;
    const int gC = gridA.gC;

	file_(T.u, "mx_out/u", nx, ny, nz,0);
	file_(T.v, "mx_out/v", nx, ny, nz,0);
	file_(T.w, "mx_out/w", nx, ny, nz,0);
	file_(t.p, "mx_out/p", nx, ny, nz,0);
}
