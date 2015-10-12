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

#include <getopt.h>

// MPI
#include <mpi.h>
// TinyXML2
#include <tinyxml2.h>
// UnitTest++
/*
#include <TestReporterStdout.h>
#include <UnitTest++.h>
*/

// internals
#include "libmath/CMath.hpp"
#include "CConfiguration.hpp"
#include "CController.hpp"
#include "CDomain.hpp"
#include "CManager.hpp"

CVector<3,int> lbm_units[] = {
		E0, E1, E2, E3,
        E4, E5, E6, E7,
        E8, E9, E10, E11,
        E12, E13, E14, E15,
        E16, E17, E18
};

#define VALIDATION_RANK 0
void extract_comma_separated_integers(std::vector<int> &int_list, std::string &int_string)
{
    size_t start_pos = 0;

    size_t comma_pos;

    comma_pos = int_string.find_first_of(',', start_pos);

    while (comma_pos != std::string::npos)
    {
        std::string substring = int_string.substr(start_pos, comma_pos-start_pos);

        int number;

        if (substring.empty())  number = 0;
        else                    number = atoi(substring.c_str());

        int_list.push_back(number);
        start_pos = comma_pos+1;
        comma_pos = int_string.find_first_of(',', start_pos);
    }
    int_list.push_back(atoi(int_string.substr(start_pos).c_str()));
}

int main(int argc, char** argv)
{
	/*
     * TODO
     * Check if the following variables should be moved to constants.h
     */
    bool debug = false;
    bool unit_test = false;
    std::string test_suite;
    std::string conf_file;

    int device_nr = 0;
    char optchar;
    while ((optchar = getopt(argc, argv, "u:c:")) > 0)
    {
        switch(optchar)
        {
        case 'u':
            unit_test = true;
            test_suite = optarg;
            break;
        case 'c':
            conf_file = optarg;
            break;
        default:
            std::cout << "Not enough input arguments. Exiting!" << std::endl;
			return -1;
        }
    }

	// load the config file and read the configuration parameters
    CConfiguration<TYPE>* configuration = new CConfiguration<TYPE>(conf_file);
    std::vector<dim3> lbm_opencl_number_of_threads_list;

#if DEBUG
    configuration->debug_mode = true;
    configuration->printMe();
#endif

    /* if(unit_test)
    {
        if(strcmp( "all", test_suite.c_str() ) == 0 )
        {
            if( debug)
                std::cout << "running all test." << std::endl;
            return UnitTest::RunAllTests();
        } else {
            const UnitTest::TestList& allTests(UnitTest::Test::GetTestList());
            UnitTest::TestList selectedTests;
            UnitTest::Test* p = allTests.GetHead();
            while( p )
            {
                //            for( int i = 1 ; i < argc ; ++i )
                if(strcmp( p->m_details.suiteName , test_suite.c_str()) == 0) {
                    selectedTests.Add( p );
                    if(debug)
                        std::cout << "Added test " << p->m_details.testName << "from suite " <<  p->m_details.suiteName << " to tes list." << std::endl;
                }
                p = p->next;
            }
            //run selected test(s) only
            UnitTest::TestReporterStdout reporter;
            UnitTest::TestRunner runner(reporter);
            return runner.RunTestsIf(selectedTests, 0, UnitTest::True(), 0 );
        }
    } else */ if(configuration->do_validate) {
        CVector<3,int> origin(0,0,0);
    	/*
        CDomain<TYPE> domain(-1, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size, origin, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length);
        */

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);    /// Start MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// Get number of processes

        CManager<TYPE> *manager = new CManager<TYPE>(my_rank, configuration);

        /*
         * TODO
         * Should be obsolete since parameter checks are done in CConfig and
         * the debug commands in CManager
         */
        /*
        if(num_procs != CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements()+1)
        {
            std::cout << "Number of allocated processors should be one more than of the number of subdomains! The last processor is used to generate validation data." << std::endl;
            std::cout << "Current number of subdomains is : " << CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }
        */

        if(my_rank != num_procs - 1)
        {
            // manager->initSimulation(my_rank);
            manager->run();

            // getting the local data
            /*
            if(my_rank == VALIDATION_RANK)
            {
                CController<TYPE>* controller = manager->getController();
                CVector<3,int> local_size_with_halo = controller->getDomain().getSize();
                CVector<3,int> local_size_without_halo(local_size_with_halo[0]-2, local_size_with_halo[1]-2, local_size_with_halo[2]-2 );
                int local_size[] = {local_size_without_halo[0], local_size_without_halo[1], local_size_without_halo[2]};
                TYPE* local_data = new TYPE[local_size_without_halo.elements()*3];
                CVector<3,int> local_origin(1,1,1);
                controller->getSolver()->getVelocities(local_origin, local_size_without_halo, local_data);
                MPI_Send(local_size, 3, MPI_INT, num_procs - 1, 0, MPI_COMM_WORLD);
                MPI_Send(local_data, local_size_without_halo.elements()*3, MPI_FLOAT, num_procs - 1, 1, MPI_COMM_WORLD );

                delete[] local_data;
            }
            */

            delete manager;
        }
        if (my_rank == num_procs - 1)
        {
            MPI_Status stat1;
            MPI_Status stat2;
            int local_size[3];
            MPI_Recv(local_size, 3, MPI_INT, VALIDATION_RANK, 0, MPI_COMM_WORLD, &stat1);
            CVector<3,int> local_size_without_halo(local_size[0], local_size[1], local_size[2]);
            TYPE* local_data = new TYPE[local_size_without_halo.elements()*3];
            MPI_Recv(local_data,local_size_without_halo.elements()*3,MPI_FLOAT, VALIDATION_RANK, 1, MPI_COMM_WORLD, &stat2);

            //  until here is correct
            //  creating the correct domain size for the case that there is only one domain.
            //  the halo regions is subtracted form each direction.
            //  The size of halo regions is dependent to the number of subdomains in each direction
            CVector<3,int> validation_domain_size(
                    configuration->domain_size[0] - 2*(configuration->subdomain_num[0] - 1),
                    configuration->domain_size[1] - 2*(configuration->subdomain_num[1] - 1),
                    configuration->domain_size[2] - 2*(configuration->subdomain_num[2] - 1));
            TYPE cell_length_x = configuration->domain_length[0] / (TYPE) configuration->domain_size[0];
            TYPE cell_length_y = configuration->domain_length[1] / (TYPE) configuration->domain_size[1];
            TYPE cell_length_z = configuration->domain_length[2] / (TYPE) configuration->domain_size[2];
            CVector<3,TYPE> validation_domain_length(validation_domain_size[0]*cell_length_x, validation_domain_size[1]*cell_length_y, validation_domain_size[2]*cell_length_z);
            CDomain<TYPE> validatiaon_domain(-2, validation_domain_size, origin, validation_domain_length);
            CManager<TYPE> validataion_manager(my_rank, configuration);
            std::cout <<  my_rank << "--> Compute the results for one domain." << std::endl;

            // validataion_manager.initSimulation(-1);
            validataion_manager.run();
            TYPE* sub_global_data = new TYPE[local_size_without_halo.elements()*3];
            // getting the global data
            int id = VALIDATION_RANK;
            int tmpid = id;
            int nx, ny, nz;
            nx = tmpid % configuration->subdomain_num[0];
            tmpid /= configuration->subdomain_num[0];
            ny = tmpid % configuration->subdomain_num[1];
            tmpid /= configuration->subdomain_num[1];
            nz = tmpid;
            CVector<3,int> sub_origin(
                    1+nx*(local_size_without_halo[0]),
                    1+ny*(local_size_without_halo[1]),
                    1+nz*(local_size_without_halo[2])
            );
            validataion_manager.getController()->getSolver()->getVelocities(sub_origin, local_size_without_halo, sub_global_data);

            std::cout << "PROC. RANK: " << my_rank << " VALIDATION SIZE: " << validation_domain_size << std::endl;
            // comparing local and global data
            double tolerance = 1.0e-15;
            int error_counter = 0;
            for ( int i = 0; i < local_size_without_halo.elements()*3; i++)
            {
                if (CMath<TYPE>::abs(sub_global_data[i] - local_data[i]) > tolerance )
                {
#if 0
                    if ((i % 30) == 0)
                    {
                        std::cout << "PROC. RANK: " << my_rank << " VALIDATION FAILED at Index: " << i << std::endl;
                        std::cout << "GLOBAL_DATA[" << i <<  "] = " << sub_global_data[i] << std::endl;
                        std::cout << "LOCAL_DATA[" << i << "] = " << local_data[i] << std::endl;
                        std::cout << "PROC. RANK: " << my_rank << " DIFF: " << fabs((double)sub_global_data[i] - (double)local_data[i]) << std::endl;
                    }
#endif
                    error_counter++;
                }
            }
            std::cout << "--> PROC. RANK: " << my_rank << " TOLERANCE: "<< tolerance << " NUMBER OF FAILED CELLS/TOTAL NUMBER OF CELLS: " << error_counter << "/" << local_size_without_halo.elements() << std::endl;
        }
        MPI_Finalize();    /// Cleanup MPI
    } else {
    	/*
        CVector<3,int> origin(0,0,0);
        CDomain<TYPE> domain(-1, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size, origin, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length);
        */
        // printf("main() CDomain initialized. \n");

        // printf("main() CManager initialized. \n");

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        CManager<TYPE> *manager = new CManager<TYPE>(my_rank, configuration);

        /*
         * TODO
         * Should be obsolete since parameter checks are done in CConfig and
         * the debug commands in CManager
         */
        /*
        if(num_procs != CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements()+1)
        {
            std::cout << "Number of allocated processors should be one more than of the number of subdomains! The last processor is used to generate validation data." << std::endl;
            std::cout << "Current number of subdomains is : " << CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }
        */

        manager->run();

        printf("Rank: %i done.\n", my_rank);
        MPI_Finalize();    /// Cleanup MPI
    }

    delete configuration;

    return 0;
}
