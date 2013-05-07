/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef JDFTX_FLUID_FEX_H2O_SCALAREOS_H
#define JDFTX_FLUID_FEX_H2O_SCALAREOS_H
#include <fluid/Fex.h>

class Fex_H2O_ScalarEOS : public Fex
{
public:
	//! Create water with the ScalarEOS functional (can choose soft or hard sphere version)
	Fex_H2O_ScalarEOS(const FluidMixture*, const FluidComponent*);
    virtual ~Fex_H2O_ScalarEOS();
	
	double compute(const DataGptr* Ntilde, DataGptr* Phi_Ntilde) const;
	double computeUniform(const double* N, double* Phi_N) const;

private:
	RadialFunctionG fex_LJatt;
	struct ScalarEOS_eval* eval;
};


#endif // JDFTX_FLUID_FEX_H2O_SCALAREOS_H
