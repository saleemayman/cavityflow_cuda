/*
 * Copyright 2010 Martin Schreiber
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

#include "mpi.h"
// standard
#include <stdlib.h>
#include <iostream>
#include <list>
#include <string>
#include <cmath>

// externals
#include <UnitTest++.h>
#include <../externals/unittest-cpp/UnitTest++/src/UnitTest++.h>
#include "TestReporterStdout.h"

// internals
#include "CDomain.hpp"
#include "CController.hpp"
#include "common.h"
#include "CManager.hpp"
#include "CConfiguration.hpp"
#include "Singleton.hpp"
#include "tinyxml2.h"

CVector<3,int> E0(1,0,0)    ;
CVector<3,int> E1(-1,0,0)   ;
CVector<3,int> E2(0,1,0)    ;
CVector<3,int> E3(0,-1,0)   ;

CVector<3,int> E4(1,1,0)    ;
CVector<3,int> E5(-1,-1,0)  ;
CVector<3,int> E6(1,-1,0)   ;
CVector<3,int> E7(-1,1,0)   ;

CVector<3,int> E8(1,0,1)    ;
CVector<3,int> E9(-1,0,-1)  ;
CVector<3,int> E10(1,0,-1)  ;
CVector<3,int> E11(-1,0,1)  ;

CVector<3,int> E12(0,1,1)   ;
CVector<3,int> E13(0,-1,-1) ;
CVector<3,int> E14(0,1,-1)  ;
CVector<3,int> E15(0,-1,1)  ;

CVector<3,int> E16(0,0,1)   ;
CVector<3,int> E17(0,0,-1)  ;
CVector<3,int> E18(0,0,0)   ;
CVector<3,int> lbm_units[] = {  E0,E1,E2,E3,
        E4,E5,E6,E7,
        E8,E9,E10,E11,
        E12,E13,E14,E15,
        E16,E17,E18
};

// simulation type
// typedef float T;
//typedef double T;
#define VALIDATION_RANK 0
void extract_comma_separated_integers(std::list<int> &int_list, std::string &int_string)
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
	ConfigSingleton::Instance()->loadFile(conf_file);

#if DEBUG
    ConfigSingleton::Instance()->debug_mode = true;
    ConfigSingleton::Instance()->printMe();
#endif

    if (unit_test) {

        if ( strcmp(  "all", test_suite.c_str() ) == 0 ) {
            if( debug)
                std::cout << "running all test." << std::endl;
            return UnitTest::RunAllTests();
        }
        else {
            const UnitTest::TestList& allTests( UnitTest::Test::GetTestList() );
            UnitTest::TestList selectedTests;
            UnitTest::Test* p = allTests.GetHead();
            while( p )
            {
                //            for( int i = 1 ; i < argc ; ++i )
                if( strcmp(  p->m_details.suiteName , test_suite.c_str() ) == 0 ) {
                    selectedTests.Add( p );
                    if( debug)
                        std::cout << "Added test " << p->m_details.testName << "from suite " <<  p->m_details.suiteName << " to tes list." << std::endl;
                }
                p = p->next;
            }
            //run selected test(s) only
            UnitTest::TestReporterStdout reporter;
            UnitTest::TestRunner runner( reporter );
            return runner.RunTestsIf( selectedTests, 0, UnitTest::True(), 0 );
        }
    } else if ( ConfigSingleton::Instance()->do_validate ) {
        CVector<3,int> origin(0,0,0);
        CDomain<T> domain(-1, ConfigSingleton::Instance()->domain_size, origin, ConfigSingleton::Instance()->domain_length);
        CManager<T> *manager = new CManager<T>(domain, ConfigSingleton::Instance()->subdomain_num);

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);    /// Start MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// Get number of processes

        if ( num_procs != ConfigSingleton::Instance()->subdomain_num.elements()+1) {
            std::cout << "Number of allocated processors should be one more than of the number of subdomains! The last processor is used to generate validation data." << std::endl;
            std::cout << "Current number of subdomains is : " << ConfigSingleton::Instance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }

        if (my_rank != num_procs - 1) {
            manager->initSimulation(my_rank);
            manager->startSimulation();

            // getting the local data
            if ( my_rank == VALIDATION_RANK) {
                CController<T>* controller = manager->getController();
                CVector<3,int> local_size_with_halo = controller->getDomain().getSize();
                CVector<3,int> local_size_without_halo(local_size_with_halo[0]-2, local_size_with_halo[1]-2, local_size_with_halo[2]-2 );
                int local_size[] = {local_size_without_halo[0], local_size_without_halo[1], local_size_without_halo[2]};
                T* local_data = new T[local_size_without_halo.elements()*3];
                CVector<3,int> local_origin(1,1,1);
                controller->getSolver()->storeVelocity(local_data, local_origin, local_size_without_halo);
                MPI_Send(local_size, 3, MPI_INT, num_procs - 1, 0, MPI_COMM_WORLD);
                MPI_Send(local_data, local_size_without_halo.elements()*3, MPI_FLOAT, num_procs - 1, 1, MPI_COMM_WORLD );

                delete[] local_data;
            }

            delete manager;
        }
        if (my_rank == num_procs - 1) {
            MPI_Status stat1;
            MPI_Status stat2;
            int local_size[3];
            MPI_Recv(local_size, 3, MPI_INT, VALIDATION_RANK, 0, MPI_COMM_WORLD, &stat1);
            CVector<3,int> local_size_without_halo(local_size[0], local_size[1], local_size[2]);
            T* local_data = new T[local_size_without_halo.elements()*3];
            MPI_Recv(local_data,local_size_without_halo.elements()*3,MPI_FLOAT, VALIDATION_RANK, 1, MPI_COMM_WORLD, &stat2 );

            //  until here is correct
            //  creating the correct domain size for the case that there is only one domain.
            //  the halo regions is subtracted form each direction.
            //  The size of halo regions is dependent to the number of subdomains in each direction
            CVector<3,int> validation_domain_size(
                    ConfigSingleton::Instance()->domain_size[0] - 2*(ConfigSingleton::Instance()->subdomain_num[0] - 1),
                    ConfigSingleton::Instance()->domain_size[1] - 2*(ConfigSingleton::Instance()->subdomain_num[1] - 1),
                    ConfigSingleton::Instance()->domain_size[2] - 2*(ConfigSingleton::Instance()->subdomain_num[2] - 1));
            T cell_length_x = ConfigSingleton::Instance()->domain_length[0] / (T)ConfigSingleton::Instance()->domain_size[0];
            T cell_length_y = ConfigSingleton::Instance()->domain_length[1] / (T)ConfigSingleton::Instance()->domain_size[1];
            T cell_length_z = ConfigSingleton::Instance()->domain_length[2] / (T)ConfigSingleton::Instance()->domain_size[2];
            CVector<3,T> validation_domain_length(validation_domain_size[0]*cell_length_x, validation_domain_size[1]*cell_length_y, validation_domain_size[2]*cell_length_z);
            CDomain<T> validatiaon_domain(-2, validation_domain_size, origin, validation_domain_length);
            CManager<T> validataion_manager(validatiaon_domain, CVector<3,int>(1,1,1));
            std::cout <<  my_rank << "--> Compute the results for one domain." << std::endl;

            validataion_manager.initSimulation(-1);
            validataion_manager.startSimulation();
            T* sub_global_data = new T[local_size_without_halo.elements()*3];
            // getting the global data
            int id = VALIDATION_RANK;
            int tmpid = id;
            int nx, ny, nz;
            nx = tmpid % ConfigSingleton::Instance()->subdomain_num[0];
            tmpid /= ConfigSingleton::Instance()->subdomain_num[0];
            ny = tmpid % ConfigSingleton::Instance()->subdomain_num[1];
            tmpid /= ConfigSingleton::Instance()->subdomain_num[1];
            nz = tmpid;
            CVector<3,int> sub_origin(
                    1+nx*(local_size_without_halo[0]),
                    1+ny*(local_size_without_halo[1]),
                    1+nz*(local_size_without_halo[2])
            );
            validataion_manager.getController()->getSolver()->storeVelocity(sub_global_data, sub_origin, local_size_without_halo);

            std::cout << "PROC. RANK: " << my_rank << " VALIDATION SIZE: " << validation_domain_size << std::endl;
            // comparing local and global data
            double tolerance = 1.0e-15;
            int error_counter = 0;
            for ( int i = 0; i < local_size_without_halo.elements()*3; i++)
            {
                if (fabs(sub_global_data[i] - local_data[i]) > tolerance ){
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
    }
    else {
        CVector<3,int> origin(0,0,0);
        CDomain<T> domain(-1, ConfigSingleton::Instance()->domain_size, origin, ConfigSingleton::Instance()->domain_length);
        // printf("main() CDomain initialized. \n");

        CManager<T> *manager = new CManager<T>(domain, ConfigSingleton::Instance()->subdomain_num);
        // printf("main() CManager initialized. \n");

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);    /// Start MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// Get number of processes

        if ( num_procs != ConfigSingleton::Instance()->subdomain_num.elements()) {
            std::cout << "Number of allocated processors should be equal to the number of subdomains!" << std::endl;
            std::cout << "Current number of subdomains is : " << ConfigSingleton::Instance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }

        //printf("main() before initSimulation. \n");
        manager->initSimulation(my_rank);

        //printf("main() start simulation. \n");
        manager->startSimulation();
//        delete manager;

        printf("Rank: %i done.\n", my_rank);
        MPI_Finalize();    /// Cleanup MPI
    }

    return 0;

}
