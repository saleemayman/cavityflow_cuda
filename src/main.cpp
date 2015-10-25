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

#include <mpi.h>

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

    CManager<TYPE> *manager = new CManager<TYPE>(rank, configuration);

    if (configuration->doValidation)
    {
    } else {
        manager->run();
    }

    delete manager;

    /*
	 * Tear down MPI environment
	 */
	MPI_Finalize();

	delete configuration;

	return 0;
}
