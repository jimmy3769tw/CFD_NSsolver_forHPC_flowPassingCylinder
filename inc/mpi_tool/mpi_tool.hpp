#pragma once


#include "0_General.hpp"


// true
//   auto [mpi_word_size, mpi_word_rank, mpi_coord, mpi_neighborhood] = mpi_init(reorder, XYZ, argc, argv);

std::tuple< MPI_Comm, int, int, std::vector<int>, std::vector<int> > mpi_init(
    bool &reorder,
    std::vector<int> &XYZ,
    int argc, 
    char **argv
)
{

    MPI_Init(&argc, &argv);

    int mpi_word_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_word_rank);

    int mpi_word_size;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_word_size);

    // int Check = 1;

    // for (auto dims:XYZ){
    //   Check *= dims;
    // }
   
    // if (Check!=mpi_word_size) {
    //     MPI_Barrier(comm_world);
    //     MPI_Finalize();
    //    throw std::invalid_argument("MPI size doesn't correct.");
    // }

    // * virtual process topology functions 

    // * number of dimensions of cartesian grid
    const int ndims = XYZ.size();

    MPI_Comm comm_world;

    int periods[ndims] = {0}; 

    MPI_Cart_create(MPI_COMM_WORLD, ndims, &XYZ[0], periods, reorder, &comm_world);

    std::vector<int> mpi_coord(ndims);
    
    MPI_Cart_coords(comm_world, mpi_word_rank, mpi_coord.size(), &mpi_coord[0]);

    std::vector<int> mpi_neighborhood(2*ndims);

    auto disp = 1;
    for (size_t i = 0; i < ndims ; ++i){
      MPI_Cart_shift(comm_world, i, disp, 
        &mpi_neighborhood[i*2], &mpi_neighborhood[i*2+1]);
    }

    cout << "[ank, coord ] " << mpi_word_rank <<", {";
    for(size_t i = 0; i < mpi_coord.size(); i++){
      cout << mpi_coord.at(i) ;
      if (i != mpi_coord.size()-1)
        cout << ", ";
    }
    for (size_t i = 0 ; i < ndims ;++i){

      std::cout << "}>> Neberhood | " 
        << mpi_neighborhood[i*2]<< ", "
        << mpi_neighborhood[i*2+1] << endl;
    }

    // ! -------------------------------------------------
    char name[1024];
    int length=1024, minor;
    MPI_Get_processor_name(name, &length);

    int major;
    MPI_Get_version(&major, &minor);
    
    // MPI_Request req[10];

    if (mpi_word_rank == 0){
        cout << "\nMPI Version " << major << "." << minor << endl;
        cout << "This Project is from " << name << endl;
    }

    // ! -------------------------------------------------

    return std::tie(comm_world, mpi_word_size, mpi_word_rank, mpi_coord, mpi_neighborhood);
}


// #define MPI_NONBLOCKING_SR
#define MPI_BLOCKING_SR


// mx_comm_world, t1.p, Lo, gridA)

void mpi_iSR_double_x(
  int nolayers,
  MPI_Comm &comm_world,
  std::vector<int> &mpi_neighborhood,
  std::vector<double> &v,
  divideLocal& Lo,
  grid& gridA
){

  if(nolayers < 0){
    throw std::invalid_argument("nolayers < 0 is not suporting!!");
  }

  MPI_Barrier(comm_world);

  MPI_Request requests[4];

  int nynz = gridA.nz * gridA.ny;

  int itag[2];
  
  int nolayers_nynz = nolayers * nynz;
  
  int shift_recv[2];

  int shift_send[2];

  itag[0] = 110;
  itag[1] = 220;
  // make a tuple 
  
  shift_send[0] = (Lo.i_begin) * nynz;
  shift_recv[0] = (Lo.i_endof) * nynz;// !`+`

    // const int icel_endof = (Lo.i_endof_table.at(i_rank)-1)*nz*ny + (ny-1)*nz + (nz-1);
    // const int count = 1- icel_begin + icel_endof;

  shift_send[1] = (Lo.i_endof - nolayers) * nynz;
  shift_recv[1] = (Lo.i_begin - nolayers) * nynz;// !`-`

  #ifdef MPI_NONBLOCKING_SR
    // !`+`
    MPI_Isend(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, requests+1);
    
    // !`-`
    MPI_Isend(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world, requests+2);
    MPI_Irecv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, requests+3);
   
    MPI_Status status[4];

    MPI_Waitall(4, requests, status);
    
  #endif

  #ifdef MPI_BLOCKING_SR
    MPI_Status status[2];
    // !`+`
    MPI_Send(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world);
    MPI_Recv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, status);
    
    // !`-`
    MPI_Send(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world);
    MPI_Recv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, status+1);
  #endif

  MPI_Barrier(comm_world);
}



void mpi_iSR_double_x_half(
  int nolayers,
  MPI_Comm &comm_world,
  std::vector<int> &mpi_neighborhood,
  std::vector<double> &v,
  divideLocal& Lo,
  grid& gridA
){

  MPI_Barrier(comm_world);

  int itag[2];

  const auto abs_nolayers = std::abs(nolayers);
  
  const int nolayers_nynz = abs_nolayers* gridA.nz *gridA.ny;
  
  const int nynz = gridA.nz * gridA.ny;

  int shift_recv[2];

  int shift_send[2];

  itag[0] = 100;
  itag[1] = 200;

  shift_send[0] = (Lo.i_begin) * nynz;
  shift_recv[0] = (Lo.i_endof) * nynz;// !`+`
  shift_send[1] = (Lo.i_endof-abs_nolayers) * nynz;
  shift_recv[1] = (Lo.i_begin-abs_nolayers) * nynz;// !`-`

  if (nolayers > 0){

  #ifdef MPI_NONBLOCKING_SR
    MPI_Request requests[2];
    // * --------------------- `+`*
    MPI_Isend(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, requests+1);
    }
    else{
    // * --------------------- `-`*
    MPI_Request requests[2];
    MPI_Isend(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world, requests+0);
    MPI_Irecv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, requests+1);
  #endif



  #ifdef MPI_BLOCKING_SR
    MPI_Status status[1];
    // * --------------------- `+`*
    MPI_Send(&v[shift_send[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[0], comm_world);
    MPI_Recv(&v[shift_recv[0]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[0], comm_world, status);
    }
    else{
    // * --------------------- `-`*
    MPI_Status status[1];
    MPI_Send(&v[shift_send[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(1), itag[1], comm_world);
    MPI_Recv(&v[shift_recv[1]], nolayers_nynz, MPI_DOUBLE, mpi_neighborhood.at(0), itag[1], comm_world, status);
  #endif


  }
  #ifdef MPI_NONBLOCKING_SR
  MPI_Status status[2];

  MPI_Waitall(2, requests, status);
  #endif
  
  MPI_Barrier(comm_world);
}



void mpi_iSR_double_x_Collect_to_Master(
  MPI_Comm &comm_world,
  int rank,
  const int rank_size,
  std::vector<double> &v,
  divideLocal& Lo,
  grid& gridA
){

  MPI_Barrier(comm_world);
  
  MPI_Request requests[rank_size];

  MPI_Status status[rank_size];
  
  int itag = 100;
  
  size_t nynz = gridA.nz *gridA.ny;

  const int master = 0;

  for (size_t i_rank = 0; i_rank < rank_size ; ++i_rank){

    itag += 10;
    
    int shift = (Lo.i_begin_table.at(i_rank)) * nynz ;
    
    int count = nynz*(Lo.i_length_table.at(i_rank));

    #ifdef MPI_NONBLOCKING_SR
    // TODO: FIXME: 
      if (rank != master){
        MPI_Isend(&v[shift], count, MPI_DOUBLE, master, itag, comm_world, requests+0*i_rank*2);
      }
      else{
        MPI_Irecv(&v[shift], count, MPI_DOUBLE, i_rank, itag, comm_world, requests+1*i_rank*2);
      }
    #endif

    #ifdef MPI_BLOCKING_SR
      if (rank != master){
        if (rank == i_rank){

          MPI_Send((void *)&v[shift], count, MPI_DOUBLE, master, itag, comm_world);
        }
      }
      else{
        MPI_Recv((void *)&v[shift], count, MPI_DOUBLE, i_rank, itag, comm_world, status+i_rank);
      }
    #endif

  }

  MPI_Barrier(comm_world);

  #ifdef MPI_NONBLOCKING_SR
    MPI_Waitall(rank_size, requests, status);
  #endif

}



void mpi_iSR_double_x_debugger(
  int size,
  int rank,
  MPI_Comm &comm_world,
  std::vector<int> &mpi_neighborhood,
  std::vector<double> &v,
  divideLocal& Lo,
  grid& gridA
){

  MPI_Barrier(comm_world);
  std::vector<std::vector<double> > Buffer(size,v);
  // TODO: check this is correct:??
  // TODO: Unit test

  int nynz = gridA.nz *gridA.ny;

  int nz = gridA.nz;

  int nx = gridA.nx;

  int ny = gridA.ny;

  for(size_t i_rank = 0; i_rank < size ; ++i_rank){
    MPI_Bcast((void *)&Buffer[i_rank][0], v.size(), MPI_DOUBLE, i_rank, comm_world);
    MPI_Barrier(comm_world);
  }

  MPI_Barrier(comm_world);

  for(size_t i_rank = 0; i_rank < size ; ++i_rank){
    const int icel_begin =  Lo.i_begin_table.at(i_rank)*nz*ny + 0*nz + 0;
    const int icel_endof = (Lo.i_endof_table.at(i_rank)-1)*nz*ny + (ny-1)*nz + (nz-1);
    const int count = 1- icel_begin + icel_endof;
    // cout << "count " << count << ", " << Lo.i_length_table.at(i_rank)*nynz << endl;
    // Lo.i_endof_table.at(i_rank);


    // if (i_rank != rank){
    //   for (size_t i = Lo.i_begin_table.at(i_rank) ; i < Lo.i_endof_table.at(i_rank) ; ++i )
    //   for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    //   for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    //   {
    //     const int icel = i*nz*ny + j*nz + k;
    //     v.at(icel) = Buffer[i_rank].at(icel);
    //   }
    // }

    if (i_rank != rank){
      for (size_t j = 0; j < count ;++j){
        v.at(icel_begin+j) = Buffer[i_rank].at(icel_begin+j);
      }
    }

  }

  MPI_Barrier(comm_world);
}




void mpi_debug_file(
  divideLocal& Lo,
  simpulationVariable& simu,
  int vec_size,
  int mpi_word_size,
  MPI_Comm &mx_comm_world,
  std::vector<int> &mpi_neighborhood,
  grid& gridA
)
{
  std::vector<double> deVec1(vec_size);

  std::fill(deVec1.begin(), deVec1.end(),simu.PID);

  for (size_t i = 0; i < gridA.nx ; ++i)
  for (size_t j = 0; j < gridA.ny ; ++j)
  for (size_t k = 0; k < gridA.nz ; ++k)
  {
    const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
    deVec1.at(icel) *= 10000;
    deVec1.at(icel) += icel;
  } 

  auto deVec2 = deVec1;
  mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, deVec1, Lo, gridA);

  // mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, deVec2, Lo, gridA );

  std::cout <<  std::flush ;
  MPI_Barrier(mx_comm_world);


	ofstream file;
  std::string filename = "Information/MPI_file";
	filename += std::to_string(simu.PID);
	filename += ".dat";	

	file.open(filename);

	for (size_t i = 0; i < gridA.nx; ++i)
	for (size_t j = 0; j < gridA.ny; ++j) 
	for (size_t k = 0; k < gridA.nz; ++k) 
	{
    const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
		file <<"(" <<  icel << ")\ti\t" << i << "\tj\t" << j << "\tk\t" << k 
				<< "\t" <<  deVec1.at(icel) << "\t" << deVec2.at(icel)<< "\t" << deVec1.at(icel) -  deVec2.at(icel) <<endl;

		// file.write((char *)(&T3.u.at(icel)), sizeof(double));
	}
	file.close();

      // for (size_t i = 0; i < gridA.nx ; ++i)
      // for (size_t j = 0; j < gridA.ny ; ++j)
      // for (size_t k = 0; k < gridA.nz ; ++k)
      // {
      //   const int icel = i*gridA.nz*gridA.ny + j*gridA.nz + k;
      //   // if (abs(deVec1.at(icel) - deVec2.at(icel)) > 1.0e-4){
      //   cout << std::flush;
      //   MPI_Barrier(mx_comm_world);
      //   if (simu.PID == 0 ){
      //   cout << endl;
      //   cout << std::flush  << i << ", " 
      //                       << j << ", " 
      //                       << k << "::";
      //   }

      //   for(size_t ii = 0; ii < mpi_word_size ; ii++ ){
      //     cout << std::flush ;
      //     MPI_Barrier(mx_comm_world);
      //     if (simu.PID == ii ){
      //       cout <<  " | " << deVec1.at(icel) << " | " << deVec2.at(icel);
      //     } 
      //   }
      // }

}





