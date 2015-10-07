#ifndef __CCONFIGURATION_HPP__
#define __CCONFIGURATION_HPP__

#include "tinyxml2.h"
#include <iostream>
#include <iomanip>
#include <list>

/*
 * Class CConfiguration stores the necessary information for the simulation process.
 *
 */
namespace txml = tinyxml2;

#define TAG_NAME_ROOT           "lbm-configuration"
#define TAG_NAME_CHILD_ONE      "physics"
#define TAG_NAME_CHILD_TWO      "grid"
#define TAG_NAME_CHILD_THREE    "simulation"
#define TAG_NAME_CHILD_FOUR     "device"

template <typename T>
class CConfiguration
{
public:
	// physics parameters
    T viscosity;
    CVector<3, T> gravitation;
    CVector<4, T> drivenCavityVelocity;

    // grid data
    CVector<3, int> domain_size;
    CVector<3, int> subdomain_num;
    CVector<3, T> domain_length;

    // simulation configuration data
    int loops;
    T timestep;
    bool do_visualization;
    bool do_validate;

    // device configuration data
    CVector<3, int> init_kernel_block_dim;
    CVector<3, int> alpha_kernel_block_dim;
    CVector<3, int> beta_kernel_block_dim;
    int device_nr;
	
    bool debug_mode;

    CConfiguration()
    {
    }

    CConfiguration(std::string file_name)
    {
        int loading_file_res = load_xml(file_name);
        if (loading_file_res != txml::XML_SUCCESS)
            throw "Loading XML file failed";
        interpret_xml_doc();
		check_parameters();
    }

    ~CConfiguration() {

    }

    void loadFile(std::string file_name)
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

    void printMe() {
        std::cout << "################" << std::endl;
        std::cout << "# CONFIGURATION " << std::endl;
        std::cout << "################" << std::endl;
        std::cout <<  "PHYSICS: " << std::endl;
        std::cout <<  "      VISCOSITY: " <<  viscosity << std::endl;
        std::cout <<  "    GRAVITATION: " <<  gravitation<< std::endl;
        std::cout <<  "     CAVITY VEL: " <<  drivenCavityVelocity << std::endl;
        std::cout <<  "GRID: " << std::endl;
        std::cout <<  "    DOMAIN_SIZE: " <<  domain_size << std::endl;
        std::cout <<  "  SUBDOMIAN_NUM: " <<  subdomain_num<< std::endl;
        std::cout <<  "SIMULATION: " << std::endl;
        std::cout <<  "          LOOPS: " <<  loops << std::endl;
        std::cout <<  "       TIMESTEP: " <<  timestep << std::endl;
        std::cout <<  "            VTK: " <<  do_visualization << std::endl;
        std::cout <<  "       VALIDATE: " <<  do_validate << std::endl;
        std::cout <<  "DEVICE: " << std::endl;
        std::cout <<  " INIT BLOCK DIM: " <<  init_kernel_block_dim << std::endl;
        std::cout <<  "ALPHA BLOCK DIM: " <<  alpha_kernel_block_dim << std::endl;
        std::cout <<  " BETA BLOCK DIM: " <<  beta_kernel_block_dim << std::endl;
        std::cout <<  "      DEVICE_NR: " <<  device_nr << std::endl;
    }
private:
    txml::XMLDocument doc;

    int load_xml(std::string file_name)
    {
        doc.LoadFile( file_name.c_str() );
        return doc.ErrorID();
    }

    void interpret_physiscs_data(const txml::XMLNode* root) {
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
    void interpret_grid_data(const txml::XMLNode* root) {
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
    void interpret_simulation_data(const txml::XMLNode* root) {
        const txml::XMLNode* child_three = root->FirstChildElement(TAG_NAME_CHILD_THREE);
        loops = atoi(child_three->FirstChildElement( "loops" )->GetText());
        timestep = atof(child_three->FirstChildElement( "timestep" )->GetText());
        do_visualization = atoi(child_three->FirstChildElement( "visualization" )->FirstChildElement("VTK")->GetText());
        do_validate = atoi(child_three->FirstChildElement( "validate" )->GetText());
    }
    void interpret_device_data(const txml::XMLNode* root) {
        const txml::XMLNode* child_four = root->FirstChildElement(TAG_NAME_CHILD_FOUR);
        init_kernel_block_dim[0] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
        init_kernel_block_dim[1] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
        init_kernel_block_dim[2] = atoi(child_four->FirstChildElement( "init-kernel-block-dim" )->FirstChildElement( "z" )->GetText());
        alpha_kernel_block_dim[0] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
        alpha_kernel_block_dim[1] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
        alpha_kernel_block_dim[2] = atoi(child_four->FirstChildElement( "alpha-kernel-block-dim" )->FirstChildElement( "z" )->GetText());
        beta_kernel_block_dim[0] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "x" )->GetText());
        beta_kernel_block_dim[1] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "y" )->GetText());
        beta_kernel_block_dim[2] = atoi(child_four->FirstChildElement( "beta-kernel-block-dim" )->FirstChildElement( "z" )->GetText());
        device_nr = atoi(child_four->FirstChildElement( "device-number" )->GetText());
    }

	void check_parameters()
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

    void interpret_xml_doc() {
        const txml::XMLNode* root = doc.FirstChildElement(TAG_NAME_ROOT);
        interpret_device_data(root);
        interpret_grid_data(root);
        interpret_physiscs_data(root);
        interpret_simulation_data(root);
    }

};

#endif
