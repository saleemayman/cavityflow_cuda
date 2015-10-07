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
#include <TestReporterStdout.h>
#include <UnitTest++.h>

// internals
// #include "helper.h"
#include "CConfiguration.hpp"
#include "CController.hpp"
#include "CDomain.hpp"
#include "CManager.hpp"
#include "CSingleton.hpp"

CVector<3,int> lbm_units[] = {
		E0, E1, E2, E3,
        E4, E5, E6, E7,
        E8, E9, E10, E11,
        E12, E13, E14, E15,
        E16, E17, E18
};

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
    /*
     * TODO
     * Check if the following variables should be moved to constants.h
     */
    bool debug = false;

    CVector<3,int> domain_size(32,32,32);
    CVector<3,int> subdomain_nums(1,1,1);
    CVector<3,TYPE> domain_length(0.1,0.1,0.1);

    CVector<3,TYPE> gravitation(0,-9.81,0);
    CVector<4,TYPE> drivenCavityVelocity(100,0,0,1);
    TYPE viscosity = 0.001308;
    TYPE timestep = -1.0;
    int loops = -1;

    bool do_validate = false;
    bool do_visualisation = false;
    bool take_frame_screenshots = false;
    bool unit_test = false;
    size_t computation_kernel_count = 128;
    size_t threads_per_dimension = 32;

    std::string number_of_registers_string; ///< string storing the number of registers for opencl threads separated with comma
    std::string number_of_threads_string;   ///< string storing the number of threads for opencl separated with comma
    std::string test_suite;
    bool use_config_file = false;
    std::string conf_file;

    int device_nr = 0;

    char optchar;
    while ((optchar = getopt(argc, argv, "x:y:z:d:vr:k:f:gG:t:sl:R:T:X:Y:Z:S:u:c:n:m:p:")) > 0)
    {
        switch(optchar)
        {
        case 'd':
            device_nr = atoi(optarg);
            break;

        case 'v':
            debug = true;
            break;

        case 'R':
            number_of_registers_string = optarg;
            break;

        case 'T':
            number_of_threads_string = optarg;
            break;

        case 'x':
            domain_size[0] = atoi(optarg);
            break;

        case 'X':
            subdomain_nums[0] = atoi(optarg);
            break;
        case 'Y':
            subdomain_nums[1] = atoi(optarg);
            break;
        case 'Z':
            subdomain_nums[2] = atoi(optarg);
            break;
        case 'S':
            domain_size[0] = atoi(optarg);
            domain_size[1] = atoi(optarg);
            domain_size[2] = atoi(optarg);
            break;
        case 'l':
            loops = atoi(optarg);
            break;

        case 'y':
            domain_size[1] = atoi(optarg);
            break;

        case 'k':
            computation_kernel_count = atoi(optarg);
            break;
        
        case 'f':
            threads_per_dimension = atoi(optarg);
            break;

        case 'z':
            domain_size[2] = atoi(optarg);
            break;

        case 'r':
            viscosity = atof(optarg);
            break;
        case 'n':
          domain_length[0] = atof(optarg);
          break;
        case 'm':
          domain_length[1] = atof(optarg);
          break;
        case 'p':
          domain_length[2] = atof(optarg);
          break;
        case 'g':
            do_visualisation = true;
            break;

        case 'G':
            gravitation[1] = atof(optarg);
            break;

        case 't':
            timestep = atof(optarg);
            break;

        case 's':
            take_frame_screenshots = true;
            break;

        case 'u':
            unit_test = true;
            test_suite = optarg;
            break;
        case 'c':
            use_config_file = true;
            conf_file = optarg;
            break;
        default:
            goto parameter_error;
        }
    }

    goto parameter_error_ok;

    parameter_error:
    std::cout << "usage: " << argv[0] << std::endl;
    std::cout << "      [-x resolution_x, default: 32]" << std::endl;
    std::cout << "      [-y resolution_y, default: 32]" << std::endl;
    std::cout << "      [-z resolution_z, default: 32]" << std::endl;
    std::cout << "      [-X subdomain resolution_x, default: 1]" << std::endl;
    std::cout << "      [-Y subdomain resolution_y, default: 1]" << std::endl;
    std::cout << "      [-Z subdomain resolution_z, default: 1]" << std::endl;
    std::cout << std::endl;
    std::cout << "      [-G gravitation in down direction, default: -9.81]" << std::endl;
    std::cout << std::endl;
    std::cout << "      [-r viscosity, default: 0.001308]" << std::endl;
    std::cout << "      [-v]    (debug mode on, be verbose)" << std::endl;
    std::cout << "      [-d device_num] (-1: list available devices, or select device)" << std::endl;
    std::cout << "      [-k number of kernels to run ]  (default: 0 - autodetect maximum)" << std::endl;
    std::cout << "      [-g]    (activate gui, default:disabled)" << std::endl;
    std::cout << "      [-p]    (pause simulation at start, default:disabled)" << std::endl;
    std::cout << std::endl;
    std::cout << "      [-R registers for OpenCL]   (comma separated list of number of registers per threads in OpenCL)" << std::endl;
    std::cout << "      [-T work group threads]     (comma separated list of number of threads in OpenCL)" << std::endl;
    std::cout << std::endl;
    std::cout << "      [-t timestep]   (default: -1 for automatic detection)" << std::endl;
    std::cout << "      [-s]    (take a screenshot every frame - default: disabled)" << std::endl;
    std::cout << "      [-u] run unit tests" << std::endl;
    std::cout << "      [-c] read the configuration file" << std::endl;
    return -1;

    parameter_error_ok:
    // std::list<int> lbm_opencl_number_of_registers_list;
    std::array<dim3, 3> lbm_opencl_number_of_threads_list;

    /*
     * TODO
     * Implement new way to manage grid configurations
     */
    /*
    if (!number_of_threads_string.empty())
        extract_comma_separated_integers(lbm_opencl_number_of_threads_list, number_of_threads_string);

    if (!number_of_registers_string.empty())
        extract_comma_separated_integers(lbm_opencl_number_of_registers_list, number_of_registers_string);
    */

    if(use_config_file)
    {
        CSingleton<CConfiguration<TYPE> >::getInstance()->loadFile(conf_file);
    } else {
        CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size = domain_size;
        CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num = subdomain_nums;
        CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length = domain_length;
        CSingleton<CConfiguration<TYPE> >::getInstance()->gravitation = gravitation;
        CSingleton<CConfiguration<TYPE> >::getInstance()->viscosity = viscosity;
        CSingleton<CConfiguration<TYPE> >::getInstance()->computation_kernel_count = computation_kernel_count;
        CSingleton<CConfiguration<TYPE> >::getInstance()->threads_per_dimension = threads_per_dimension;
        CSingleton<CConfiguration<TYPE> >::getInstance()->device_nr = device_nr;
        CSingleton<CConfiguration<TYPE> >::getInstance()->do_visualization = do_visualisation;
        CSingleton<CConfiguration<TYPE> >::getInstance()->timestep = timestep;
        CSingleton<CConfiguration<TYPE> >::getInstance()->loops = loops;
        // CSingleton<CConfiguration<TYPE> >::getInstance()->lbm_opencl_number_of_registers_list = lbm_opencl_number_of_registers_list;
        CSingleton<CConfiguration<TYPE> >::getInstance()->lbm_opencl_number_of_threads_list = lbm_opencl_number_of_threads_list;
        CSingleton<CConfiguration<TYPE> >::getInstance()->do_validate = do_validate;
        CSingleton<CConfiguration<TYPE> >::getInstance()->drivenCavityVelocity = drivenCavityVelocity;
    }
#if DEBUG
    CSingleton<CConfiguration<TYPE> >::getInstance()->debug_mode = true;
    CSingleton<CConfiguration<TYPE> >::getInstance()->printMe();
#endif

    if(unit_test)
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
    } else if(CSingleton<CConfiguration<TYPE> >::getInstance()->do_validate) {
        CVector<3,int> origin(0,0,0);
        CDomain<TYPE> domain(-1, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size, origin, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length);
        CManager<TYPE> *manager = new CManager<TYPE>(domain, CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num);

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);    /// Start MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// Get number of processes

        if(num_procs != CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements()+1)
        {
            std::cout << "Number of allocated processors should be one more than of the number of subdomains! The last processor is used to generate validation data." << std::endl;
            std::cout << "Current number of subdomains is : " << CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }

        if(my_rank != num_procs - 1)
        {
            manager->initSimulation(my_rank);
            manager->startSimulation();

            // getting the local data
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
                    CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[0] - 2*(CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[0] - 1),
                    CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[1] - 2*(CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[1] - 1),
                    CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[2] - 2*(CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[2] - 1));
            TYPE cell_length_x = CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length[0] / (TYPE) CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[0];
            TYPE cell_length_y = CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length[1] / (TYPE) CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[1];
            TYPE cell_length_z = CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length[2] / (TYPE) CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size[2];
            CVector<3,TYPE> validation_domain_length(validation_domain_size[0]*cell_length_x, validation_domain_size[1]*cell_length_y, validation_domain_size[2]*cell_length_z);
            CDomain<TYPE> validatiaon_domain(-2, validation_domain_size, origin, validation_domain_length);
            CManager<TYPE> validataion_manager(validatiaon_domain, CVector<3,int>(1,1,1));
            std::cout <<  my_rank << "--> Compute the results for one domain." << std::endl;

            validataion_manager.initSimulation(-1);
            validataion_manager.startSimulation();
            TYPE* sub_global_data = new TYPE[local_size_without_halo.elements()*3];
            // getting the global data
            int id = VALIDATION_RANK;
            int tmpid = id;
            int nx, ny, nz;
            nx = tmpid % CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[0];
            tmpid /= CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[0];
            ny = tmpid % CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[1];
            tmpid /= CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num[1];
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
                if (fabs(sub_global_data[i] - local_data[i]) > tolerance )
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
        CVector<3,int> origin(0,0,0);
        CDomain<TYPE> domain(-1, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_size, origin, CSingleton<CConfiguration<TYPE> >::getInstance()->domain_length);
        // printf("main() CDomain initialized. \n");

        CManager<TYPE> *manager = new CManager<TYPE>(domain, CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num);
        // printf("main() CManager initialized. \n");

        int my_rank, num_procs;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        if(num_procs != CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements())
        {
            std::cout << "Number of allocated processors should be equal to the number of subdomains!" << std::endl;
            std::cout << "Current number of subdomains is : " << CSingleton<CConfiguration<TYPE> >::getInstance()->subdomain_num.elements() << std::endl;
            exit(EXIT_FAILURE);
        }

        // printf("main() before initSimulation. \n");
        manager->initSimulation(my_rank);

        // printf("main() start simulation. \n");
        manager->startSimulation();
        // delete manager;

        printf("Rank: %i done.\n", my_rank);
        MPI_Finalize();    /// Cleanup MPI
    }

    return 0;
}
