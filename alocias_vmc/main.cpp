#ifdef TESTS
    #include "../tests/test_main.cpp"
#endif

#include <mpi.h>

#include "examples/run.h"

int main(int argc, char *argv[]) {
    /* main function */
    #ifdef TESTS
        test_main();
        return 0;
    #endif

    // print usage and exit if not enough arguments
    if (argc < 1) {
        /* print error message and exit if number of arguments are too few */
        std::cout << "USAGE: mpirun np 'P' ./main 'inputfile'" << std::endl;
        exit(1);
    } // end if

    // initialize MPI
    int myRank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // let eigen use threads
    Eigen::initParallel();

    // parse input, initialize wavefunction compiled for and vmc and run
    // sampling. Minimize(find optimal parameters) if parameters are not given.
    run(argv[1]);

    // clean up MPI
    MPI_Finalize();

    return 0;
} // end if
