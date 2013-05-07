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

#include <fluid/Fex_H2O_ScalarEOS_internal.h>
#include <fluid/Fex_H2O_ScalarEOS.h>
#include <fluid/Fex_LJ.h>
#include <core/Units.h>
#include <core/Operators.h>
#include <electronic/operators.h>

string rigidMoleculeCDFT_ScalarEOSpaper = "R. Sundararaman and T.A. Arias, arXiv:1302.0026";

Fex_H2O_ScalarEOS::Fex_H2O_ScalarEOS(const FluidMixture* fluidMixture, const FluidComponent* comp)
: Fex(fluidMixture, comp), eval(new ScalarEOS_eval(T))
{
	//Initialize the kernels:
	setLJatt(fex_LJatt, gInfo, -9.0/(32*sqrt(2)*M_PI*pow(2*eval->sphereRadius,3)), 2*eval->sphereRadius);
	Citations::add("Scalar-EOS water functional", rigidMoleculeCDFT_ScalarEOSpaper);
}
Fex_H2O_ScalarEOS::~Fex_H2O_ScalarEOS()
{	fex_LJatt.free();
}

#ifdef GPU_ENABLED
void Fex_H20_ScalarEOS_gpu(int nr, const double* Nbar, double* Fex, double* Phi_Nbar, ScalarEOS_eval eval);
#endif
double Fex_H2O_ScalarEOS::compute(const DataGptr* Ntilde, DataGptr* Phi_Ntilde) const
{	//Compute LJatt weighted density:
	DataRptr Nbar = I(fex_LJatt*Ntilde[0]), Phi_Nbar; nullToZero(Phi_Nbar, gInfo);
	//Evaluated weighted density functional:
	DataRptr Aex(DataR::alloc(gInfo,isGpuEnabled()));
	#ifdef GPU_ENABLED
	Fex_H20_ScalarEOS_gpu(gInfo.nr, Nbar->dataGpu(), Aex->dataGpu(), Phi_Nbar->dataGpu(), *eval);
	#else
	threadedLoop(eval, gInfo.nr, Nbar->data(), Aex->data(), Phi_Nbar->data());
	#endif
	//Convert gradients:
	DataRptr NO = I(Ntilde[0]);
	Phi_Ntilde[0] += fex_LJatt*Idag(NO*Phi_Nbar) + Idag(Aex);
	return gInfo.dV*dot(NO,Aex);
}

double Fex_H2O_ScalarEOS::computeUniform(const double* N, double* Phi_N) const
{	double AexPrime, Aex;
	(*eval)(0, &N[0], &Aex, &AexPrime);
	Phi_N[0] += Aex + N[0]*AexPrime;	
	return N[0]*Aex;
}
