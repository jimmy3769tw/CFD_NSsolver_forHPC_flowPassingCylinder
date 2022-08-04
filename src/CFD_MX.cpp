#include "0_General.hpp"
#include "1_Run.hpp"

using std::string, 
      std::vector, 
      std::abs, 
      std::pow, 
      std::sqrt, 
      std::cout, 
      std::endl; 

int main(int argc, char **argv){

    // * CFD_MX_struct
    clockstruct timer;


    simpulationVariable simu;

    grid gA;

    // *   parameter  ------------
    simu.set_time_max(40);

    simu.set_dt(0.002); 

    simu.set_Re(40.0);

    simu.init_fileStep(1);
    // *   parameter  ------------


    // *   Accuracy  ------------
    std::cout << "\n----------------------------------------------------" 
              << "\n| Gridder :                       | "             << gA.Gridder
              << "\n| [Lx:Ly:Lz]" << gA.lx << " : " << gA.ly << " : " << gA.lz
              << "\n----------------------------------------------------" ;

    simu.DfibMethod = "DFIB_Cylinder-Z";//DFIB_Cylinder-X, DFIB_Cylinder-Z, "DFIB_Cylinder-Y", "OFF"

    // *  >>>>>>>>>>>>>----- SETTING (RUNTIME)

    simu.Locality = 1;

    // *  >>>>>>>>>>>>>----- SETTING (RUNTIME)


    // ! RUN ============
    int mpi_word_rank = 0;

    if (mpi_word_rank == 0){

        std::cout << "\n# MainLoop :" << simu.loop_max << std::endl;

        std::cout << "----------------------------------------------------" << std::endl;
    }
    


    // ! TEXT ============
    #if defined (PC_SEQ) || defined (PC_OMP)

        cout << simu.Sequential_TEXT << endl;

    #elif defined (PC_HYBRID_MPI_OMP)

    #elif defined (PC_OCL)

    #endif

    runSeqOmp(gA, timer, simu, argc, argv);

    // ! RUN RUN RUN ============

    std::cout << "Finish program !\n";

    return 0;
}
