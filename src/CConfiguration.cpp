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

#include "CConfiguration.hpp"

template <class T>
CConfiguration<T>::CConfiguration() : block_threads_per_dim(3)
{
}

template <class T>
CConfiguration<T>::CConfiguration(std::string file_name) : block_threads_per_dim(3)
{
	int loading_file_res = load_xml(file_name);
	if (loading_file_res != txml::XML_SUCCESS)
		throw "Loading XML file failed";
	interpret_xml_doc();
	check_parameters();
}

template <class T>
CConfiguration<T>::~CConfiguration()
{
}

template <class T>
int CConfiguration<T>::load_xml(std::string file_name)
{
	doc.LoadFile( file_name.c_str() );
	return doc.ErrorID();
}

template <class T>
void CConfiguration<T>::interpret_physiscs_data(const txml::XMLNode* root)
{
	// viscosity
	const txml::XMLNode* child_one = root->FirstChildElement(TAG_NAME_CHILD_ONE);
	viscosity = atof(child_one->FirstChildElement( "viscosity" )->GetText());

	// gravitation
	gravitation[0] = atof(child_one->FirstChildElement( "gravitation" )->FirstChildElement( "x" )->GetText());
	gravitation[1] = atof(child_one->FirstChildElement( "gravitation" )->FirstChildElement( "y" )->GetText());
	gravitation[2] = atof(child_one->FirstChildElement( "gravitation" )->FirstChildElement( "z" )->GetText());

	// cavity velocity
	drivenCavityVelocity[0] = atof(child_one->FirstChildElement( "cavity-velocity" )->FirstChildElement( "x" )->GetText());
	drivenCavityVelocity[1] = atof(child_one->FirstChildElement( "cavity-velocity" )->FirstChildElement( "y" )->GetText());
	drivenCavityVelocity[2] = atof(child_one->FirstChildElement( "cavity-velocity" )->FirstChildElement( "z" )->GetText());
	drivenCavityVelocity[3] = atof(child_one->FirstChildElement( "cavity-velocity" )->FirstChildElement( "w" )->GetText());
}

template <class T>
void CConfiguration<T>::interpret_grid_data(const txml::XMLNode* root)
{
	const txml::XMLNode* child_two = root->FirstChildElement(TAG_NAME_CHILD_TWO);
	domain_size[0] = atoi(child_two->FirstChildElement( "domain-size" )->FirstChildElement( "x" )->GetText());
	domain_size[1] = atoi(child_two->FirstChildElement( "domain-size" )->FirstChildElement( "y" )->GetText());
	domain_size[2] = atoi(child_two->FirstChildElement( "domain-size" )->FirstChildElement( "z" )->GetText());

	subdomain_num[0] = atoi(child_two->FirstChildElement( "subdomain-num" )->FirstChildElement( "x" )->GetText());
	subdomain_num[1] = atoi(child_two->FirstChildElement( "subdomain-num" )->FirstChildElement( "y" )->GetText());
	subdomain_num[2] = atoi(child_two->FirstChildElement( "subdomain-num" )->FirstChildElement( "z" )->GetText());

	domain_length[0] = atof(child_two->FirstChildElement( "domian-length" )->FirstChildElement( "x" )->GetText());
	domain_length[1] = atof(child_two->FirstChildElement( "domian-length" )->FirstChildElement( "y" )->GetText());
	domain_length[2] = atof(child_two->FirstChildElement( "domian-length" )->FirstChildElement( "z" )->GetText());
}

template <class T>
void CConfiguration<T>::interpret_simulation_data(const txml::XMLNode* root)
{
	const txml::XMLNode* child_three = root->FirstChildElement(TAG_NAME_CHILD_THREE);
	loops = atoi(child_three->FirstChildElement( "loops" )->GetText());
	timestep = atof(child_three->FirstChildElement( "timestep" )->GetText());
	do_visualization = atoi(child_three->FirstChildElement( "visualization" )->FirstChildElement("VTK")->GetText());
	do_validate = atoi(child_three->FirstChildElement( "validate" )->GetText());
}

template <class T>
void CConfiguration<T>::interpret_device_data(const txml::XMLNode* root)
{
    const txml::XMLNode* child_four = root->FirstChildElement(TAG_NAME_CHILD_FOUR);

	int threads_per_block[3];
    threads_per_block[0] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
    threads_per_block[1] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
    threads_per_block[2] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "z" )->GetText());
	block_threads_per_dim[0] = dim3(threads_per_block[0], threads_per_block[1], threads_per_block[2]);

    threads_per_block[0] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
    threads_per_block[1] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
    threads_per_block[2] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "z" )->GetText());
	block_threads_per_dim[1] = dim3(threads_per_block[0], threads_per_block[1], threads_per_block[2]);

    threads_per_block[0] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
    threads_per_block[1] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
    threads_per_block[2] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "z" )->GetText());  
	block_threads_per_dim[2] = dim3(threads_per_block[0], threads_per_block[1], threads_per_block[2]);
    device_nr = atoi(child_four->FirstChildElement( "device-number" )->GetText());                                                   
}

template <class T>
void CConfiguration<T>::check_parameters()
{
    // check parameter correctness
    if (viscosity < 0)
        throw "Loading XML file failed. Invalid value: viscosity must be positive.";

    if (domain_size[0] <= 0 || domain_size[1] <= 0 || domain_size[2] <= 0)                                                              
        throw "Loading XML file failed. Invalid value: domain size must be greater than zero.";
        
    if (subdomain_num[0] <= 0 || subdomain_num[1] <= 0 || subdomain_num[2] <= 0)                                                     
        throw "Loading XML file failed. Invalid value: number of sub-domains must be greater than zero.";

    if (domain_length[0] <= 0 || domain_length[1] <= 0 || domain_length[2] <= 0)                                                        
        throw "Loading XML file failed. Invalid value: domain-length must be greater than positive.";

    // TODO: parameter check for the kernels' block size to be done in the solver class.
}


template <class T>
void CConfiguration<T>::interpret_xml_doc()
{
	const txml::XMLNode* root = doc.FirstChildElement(TAG_NAME_ROOT);
	interpret_device_data(root);
	interpret_grid_data(root);
	interpret_physiscs_data(root);
	interpret_simulation_data(root);

}

template <class T>
void CConfiguration<T>::loadFile(std::string file_name)
{
#if DEBUG
	std::cout << "Loading XML File: " << file_name << std::endl;
#endif
	int loading_file_res = load_xml(file_name);
	if (loading_file_res != txml::XML_SUCCESS)
		throw "Loading XML file failed";
	interpret_xml_doc();
	check_parameters();
}

template <class T>
void CConfiguration<T>::printMe()
{
	std::cout << "################" << std::endl;
	std::cout << "# CONFIGURATION " << std::endl;
	std::cout << "################" << std::endl;
	std::cout <<  "PHYSICS:" << std::endl;
	std::cout <<  "       VISCOSITY: " << viscosity << std::endl;
	std::cout <<  "     GRAVITATION: " << gravitation<< std::endl;
	std::cout <<  "      CAVITY VEL: " << drivenCavityVelocity << std::endl;
	std::cout <<  "GRID:" << std::endl;
	std::cout <<  "     DOMAIN_SIZE: " << domain_size << std::endl;
	std::cout <<  "   SUBDOMIAN_NUM: " << subdomain_num<< std::endl;
	std::cout <<  "SIMULATION: " << std::endl;
	std::cout <<  "           LOOPS: " << loops << std::endl;
	std::cout <<  "        TIMESTEP: " << timestep << std::endl;
	std::cout <<  "             VTK: " << do_visualization << std::endl;
	std::cout <<  "        VALIDATE: " << do_validate << std::endl;
	std::cout <<  "DEVICE:" << std::endl;
//	std::cout <<  " INIT BLOCK DIM: " <<  block_threads_per_dim[0] << std::endl;
//	std::cout <<  "ALPHA BLOCK DIM: " <<  block_threads_per_dim[1] << std::endl;
//	std::cout <<  " BETA BLOCK DIM: " <<  block_threads_per_dim[2] << std::endl;
	std::cout <<  "        DEVICE_NR: " << device_nr << std::endl;
}

template class CConfiguration<double>;
template class CConfiguration<float>;
