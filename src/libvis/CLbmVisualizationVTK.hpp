#ifndef CLBMVISUALIZATIONVTK_HPP
#define CLBMVISUALIZATIONVTK_HPP

#include "ILbmVisualization.hpp"
#include "../libmath/CVector.hpp"
#include "VTK_Common.hpp"

template <typename T>
class CLbmVisualizationVTK : virtual public ILbmVisualization<T>
{
	int _UID;
	std::string  _file_name;
	int _timeStepNumber;

public:
	CLbmVisualizationVTK(int UID,std::string file_name): ILbmVisualization<T>()
	{
		_UID = UID;
		_file_name = file_name;
		_timeStepNumber = -1;
	}

	void setup(CLbmSolver<T> *p_cLbmOpencl) {
		ILbmVisualization<T>::setup(p_cLbmOpencl);
	}

	void render(int increment = -1)
	{
		// printf(" ->CLbmVisualization::render()");
		// using the increment in VTK visualization as the time step value
		_timeStepNumber = increment;

		// number of grid cells in each dimension (don't confuse it with the number of grid points!)
		int imax = this->cLbmOpencl->getDomain()->getSize()[0];
		int jmax = this->cLbmOpencl->getDomain()->getSize()[1];
		int kmax = this->cLbmOpencl->getDomain()->getSize()[2];

		int total_el = imax*jmax*kmax;

		// load array from device
		CVector<3, int> origin(0, 0, 0);
		CVector<3, int> size(this->cLbmOpencl->getDomain()->getSize());
		this->cLbmOpencl->getVelocities(origin, size, this->velocity);
		this->cLbmOpencl->getDensities(origin, size, this->density);
		this->cLbmOpencl->getFlags(origin, size, this->flags);

		T dx = this->cLbmOpencl->getDomain()->getLength()[0] / (T)this->cLbmOpencl->getDomain()->getSize()[0];
		T dy = dx;
		T dz = dx;
		CDomain<T> domain = *(this->cLbmOpencl->getDomain());
		char szFileName[80];
		FILE *fp = NULL;
		sprintf( szFileName, "%s.%i.%i.vtk",  _file_name.c_str(), _UID,_timeStepNumber );
		fp = fopen( szFileName, "w");
		if( fp == NULL )
		{
			fprintf(stderr, "Failed to open %s", szFileName );
			return;
		}

		write_vtkHeader( fp );
	    fprintf(fp,"DATASET STRUCTURED_GRID\n");
	    fprintf(fp,"DIMENSIONS  %i %i %i \n", imax+1, jmax+1, kmax+1);
	    fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1)*(kmax+1) );
	    fprintf(fp,"\n");
		write_vtkPointCoordinates<T>(fp, imax, jmax, kmax, dx, dy, dz, domain.getOrigin() );

		T *velx = this->velocity;
		T *vely = this->velocity+total_el;
		T *velz = this->velocity+total_el*2;
//		fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1)*(kmax+1) );

	    fprintf(fp,"\n");
	    fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)*(kmax)) );
	    fprintf(fp, "SCALARS density float 1 \n");
	    fprintf(fp, "LOOKUP_TABLE default \n");
		for (int a = 0; a < total_el ; a++)
		{
			 fprintf(fp, "%f\n", this->density[a] );
		}
		fprintf(fp,"\n");
		fprintf(fp, "SCALARS flag INT 1 \n");
		fprintf(fp, "LOOKUP_TABLE default \n");
		for (int a = 0; a < total_el ; a++)
		{
			fprintf(fp, "%i\n", this->flags[a] );
		}

		fprintf(fp,"\n");
		fprintf(fp, "VECTORS velocity float\n");
		for (int a = 0; a < total_el ; a++)
		{
			fprintf(fp, "%f %f %f\n", velx[a], vely[a], velz[a]);
		}


		if( fclose(fp) )
		{
			fprintf(stderr, "Failed to close %s", szFileName);
		}

	}
};
#endif
