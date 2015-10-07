#ifndef ILBMVISUALIZATION_HPP
#define ILBMVISUALIZATION_HPP

#include "../CLbmSolverGPU.cuh"
#include "../libmath/CVector.hpp"

/**
 * \brief Visualization Abstraction Class
 */
template <typename T>
class ILbmVisualization
{
protected:
	T *velocity;
	T *density;
	Flag *flags;
	CLbmSolverGPU<T> *cLbmOpencl;

public:
	ILbmVisualization()
		: velocity(NULL),
		  density(NULL),
		  flags(NULL)
	{
	}

    virtual ~ILbmVisualization() {
		if (velocity)
			delete[] velocity;

		if (density)
			delete[] density;

		if (flags)
			delete[] flags;
    };

	virtual void setup(CLbmSolverGPU<T> *p_cLbmOpencl) {
		cLbmOpencl = p_cLbmOpencl;

		CVector<3,int> domain_cells = cLbmOpencl->getDomain()->getNumOfCells();

		delete [] velocity;
		// printf(" ->ILbmVisualization::setup()");
		velocity = new T[domain_cells.elements()*3];

		delete [] density;
		density = new T[domain_cells.elements()];

		delete [] flags;
		flags = new Flag[domain_cells.elements()];

	}

	virtual void render(int increment = -1) = 0;
};

#endif
