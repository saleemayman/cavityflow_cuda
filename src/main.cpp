/*
 * Copyright
 * 2010 Martin Schreiber
 * 2013 Arash Bakhtiari
 * 2016 Christoph Riesinger, Ayman Saleem
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cassert>
#include <fstream>
#include <limits>
#include <sstream>
#include <typeinfo>

#include <mpi.h>

#include "libmath/CMath.hpp"
#include "CConfiguration.hpp"
#include "CManager.hpp"

CVector<3,int> lbm_units[] = {
        E0, E1, E2, E3,
        E4, E5, E6, E7,
        E8, E9, E10, E11,
        E12, E13, E14, E15,
        E16, E17, E18
};

int main(int argc, char** argv)
{
    /*
     * TODO
     */
    if (argc != 2)
    {
        std::cerr << "----- main() -----" << std::endl;
        std::cerr << "Exactly one parameter has to be passed to this executable which specifies the location of the XML configuration file." << std::endl;
        std::cerr << "Instead, " << argc << " parameters have been passed." << std::endl;
        for (int i = 0; i < argc; i++) {
            std::cerr << "argv[" << i << "] = " << argv[i] << std::endl;
        }
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "------------------" << std::endl;

        exit (EXIT_FAILURE);
    }

    CConfiguration<TYPE>* configuration = new CConfiguration<TYPE>(argv[1]);

    if (configuration->doLogging)
    {
        configuration->print();
    }

    /*
     * Setup MPI environment
     */
    int rank, numOfRanks, nodeNameLength;
    char nodeName[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfRanks);
    MPI_Get_processor_name(nodeName, &nodeNameLength);

    if (configuration->doLogging)
    {
        std::cout << "----- main() -----" << std::endl;
        std::cout << "MPI has been successfully initialized." << std::endl;
        std::cout << "------------------" << std::endl;
        std::cout << "local rank number:     " << rank << std::endl;
        std::cout << "total number of ranks: " << numOfRanks << std::endl;
        std::cout << "local node name:       " << nodeName << std::endl;
        std::cout << "------------------" << std::endl;
    }

    /*
     * Make sure there's a dedicated MPI rank/process for every subdomain.
     */
    if(numOfRanks != configuration->numOfSubdomains.elements())
    {
        std::cerr << "----- main() -----" << std::endl;
        std::cerr << "Number of launched MPI ranks/processes not equal to the number of specified subdomains." << std::endl;
        std::cerr << "number of launched MPI ranks/processes: " << numOfRanks << std::endl;
        std::cerr << "number of specified subdomains:         " << configuration->numOfSubdomains.elements() << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "------------------" << std::endl;

        exit (EXIT_FAILURE);
    }

    CManager<TYPE>* manager = new CManager<TYPE>(rank, configuration);

    manager->run();

    /*
     * Validation performs a comparison between the velocities results computed
     * on multiple ranks and the results computed on one single rank. Results
     * are compared after all iteration steps have finished.
     */
    if (configuration->doValidation)
    {
        TYPE* velocities = NULL;
        TYPE* velocitiesLocal = NULL;
        TYPE* velocitiesValidation = NULL;
        CConfiguration<TYPE>* configurationValidation;
        CManager<TYPE>* managerValidation;

        int domainSize[3] = {
                manager->getDomain()->getSize()[0],
                manager->getDomain()->getSize()[1],
                manager->getDomain()->getSize()[2]};
        int subdomainSize[3] = {
                manager->getController()->getDomain()->getSize()[0],
                manager->getController()->getDomain()->getSize()[1],
                manager->getController()->getDomain()->getSize()[2]};
        int numOfDomainCells = manager->getDomain()->getSize().elements();
        int numOfSubomainCells = manager->getController()->getDomain()->getSize().elements();

        /*
         * All ranks compute their local result.
         */
        velocitiesLocal = new TYPE[3 * numOfSubomainCells];

        manager->getController()->getSolver()->getVelocities(velocitiesLocal);

        /*
         * Rank 0 allocates additional arrays:
         * - velocities stores the gathered data from all ranks
         * - velocitiesValidation stores data computed on one single rank
         * A global solution using only rank 0 is computed.
         */
        if (rank == 0)
        {
            configurationValidation = new CConfiguration<TYPE>(argv[1]);
            configurationValidation->serialize();
            managerValidation = new CManager<TYPE>(rank, configurationValidation);

            managerValidation->run();

            velocities = new TYPE[3 * numOfDomainCells];
            velocitiesValidation = new TYPE[3 * numOfDomainCells];

            managerValidation->getController()->getSolver()->getVelocities(velocitiesValidation);

            delete managerValidation;
            delete configurationValidation;
        }

        /*
         * Gather results.
         */
        MPI_Datatype sendArray;
        MPI_Datatype recvArray;
        MPI_Datatype recvArrayResized;
        int start[3] = {0, 0, 0};
        int sendcounts[configuration->numOfSubdomains.elements()];
        int displacements[configuration->numOfSubdomains.elements()];

        MPI_Type_create_subarray(3, subdomainSize, subdomainSize, start, MPI_ORDER_FORTRAN, ((typeid(TYPE) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE), &sendArray);
        MPI_Type_commit(&sendArray);

        if (rank == 0)
        {
            MPI_Type_create_subarray(3, domainSize, subdomainSize, start, MPI_ORDER_FORTRAN, ((typeid(TYPE) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE), &recvArray);
            MPI_Type_create_resized(recvArray, 0, sizeof(TYPE), &recvArrayResized);
            MPI_Type_commit(&recvArrayResized);

            int displacementIdx;
            int subdomainIdx;

            for (int i = 0; i < configuration->numOfSubdomains[0]; i++)
            {
                for (int j = 0; j < configuration->numOfSubdomains[1]; j++)
                {
                    for (int k = 0; k < configuration->numOfSubdomains[2]; k++)
                    {
                        subdomainIdx = k * configuration->numOfSubdomains[0] * configuration->numOfSubdomains[1] + j * configuration->numOfSubdomains[0] + i;
                        displacementIdx = k * subdomainSize[2] * domainSize[0] * domainSize[1] + j * subdomainSize[1] * domainSize[0] + i * subdomainSize[0];
                        sendcounts[subdomainIdx] = 1;
                        displacements[subdomainIdx] = displacementIdx;
                    }
                }
            }
        }

        MPI_Gatherv(
                &velocitiesLocal[0], 1, sendArray,
                &velocities[0], sendcounts, displacements, recvArrayResized,
                0, MPI_COMM_WORLD);
        MPI_Gatherv(
                &velocitiesLocal[numOfSubomainCells], 1, sendArray,
                &velocities[numOfDomainCells], sendcounts, displacements, recvArrayResized,
                0, MPI_COMM_WORLD);
        MPI_Gatherv(
                &velocitiesLocal[2 * numOfSubomainCells], 1, sendArray,
                &velocities[2 * numOfDomainCells], sendcounts, displacements, recvArrayResized,
                0, MPI_COMM_WORLD);

        MPI_Type_free(&sendArray);
        if (rank == 0)
            MPI_Type_free(&recvArrayResized);
        delete[] velocitiesLocal;

        /*
         * Finally, rank 0 compares the parallel result in velocities with the
         * serial result in velocitiesValidation and outputs information if an
         * inequality is detected.
         */
        if (rank == 0)
        {
            int numOfInequalities = 0;
            int globalIdx, velocitiesX, velocitiesY, velocitiesZ;

            std::stringstream validationFileName;
            validationFileName << configuration->validationOutputDir << "/validation.txt";
            std::ofstream validationFile(validationFileName.str().c_str(), std::ios::out);

            for (int i = 0; i < domainSize[0]; i++)
                {
                    for (int j = 0; j < domainSize[1]; j++)
                    {
                        for (int k = 0; k < domainSize[2]; k++)
                        {
                            globalIdx = k * domainSize[1] * domainSize[2] + j * domainSize[1] + i;
                            velocitiesX = globalIdx;
                            velocitiesY = numOfDomainCells + globalIdx;
                            velocitiesZ = 2 * numOfDomainCells + globalIdx;

                            if (CMath<TYPE>::abs(velocities[velocitiesX] - velocitiesValidation[velocitiesX]) / velocitiesValidation[velocitiesX] > std::numeric_limits<TYPE>::epsilon() ||
                                    CMath<TYPE>::abs(velocities[velocitiesY] - velocitiesValidation[velocitiesY]) / velocitiesValidation[velocitiesY] > std::numeric_limits<TYPE>::epsilon() ||
                                    CMath<TYPE>::abs(velocities[velocitiesZ] - velocitiesValidation[velocitiesZ]) / velocitiesValidation[velocitiesZ] > std::numeric_limits<TYPE>::epsilon())
                            {
                                if (validationFile.is_open())
                                {
                                    validationFile << "[" << i << "," << j << "," << k << "]: ";
                                    validationFile << "[" << velocities[velocitiesX] << ", " << velocities[velocitiesY] << ", " << velocities[velocitiesZ] << "] vs. ";
                                    validationFile << "[" << velocitiesValidation[velocitiesX] << ", " << velocitiesValidation[velocitiesY] << ", " << velocitiesValidation[velocitiesZ] << "]" << std::endl;
                                } else {
                                    std::cerr << "----- main() -----" << std::endl;
                                    std::cerr << "There is no open file to write validation results." << std::endl;
                                    std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                                    std::cerr << "------------------" << std::endl;

                                    exit (EXIT_FAILURE);
                                }

                                numOfInequalities++;
                        }
                    }
                }
            }

            if (validationFile.is_open())
            {
                validationFile << "number of inequalities: " << numOfInequalities << std::endl;
                validationFile.close();
            } else {
                std::cerr << "----- main() -----" << std::endl;
                std::cerr << "There is no open validation file to close." << std::endl;
                std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                std::cerr << "------------------" << std::endl;

                exit (EXIT_FAILURE);
            }

            delete[] velocitiesValidation;
            delete[] velocities;
        }
    }

    delete manager;

    /*
     * Tear down MPI environment
     */
    MPI_Finalize();

    delete configuration;

    return 0;
}
