/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <fluid/IdealGasPomega.h>
#include <fluid/Euler.h>

IdealGasPomega::IdealGasPomega(const FluidMixture* fluidMixture, const FluidComponent* comp, const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(quad.nOrientations(),fluidMixture,comp), quad(quad), trans(trans), pMol(molecule.getDipole())
{
}

void IdealGasPomega::initState(const DataRptr* Vex, DataRptr* logPomega, double scale, double Elo, double Ehi) const
{	DataRptrCollection Veff(molecule.sites.size()); nullToZero(Veff, gInfo);
	for(unsigned i=0; i<molecule.sites.size(); i++)
	{	Veff[i] += V[i];
		Veff[i] += Vex[i];
	}
	double Emin=+DBL_MAX, Emax=-DBL_MAX, Emean=0.0;
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr Emolecule;
		//Sum the potentials collected over sites for each orientation:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-(rot*pos), 1.0, Veff[i], Emolecule);
		//Accumulate stats and cap:
		Emean += quad.weight(o) * sum(Emolecule)/gInfo.nr;
		double Emin_o, Emax_o;
		callPref(eblas_capMinMax)(gInfo.nr, Emolecule->dataPref(), Emin_o, Emax_o, Elo, Ehi);
		if(Emin_o<Emin) Emin=Emin_o;
		if(Emax_o>Emax) Emax=Emax_o;
		//Set contributions to the state (with appropriate scale factor):
		logPomega[o] = (-scale/T) * Emolecule;
	}
	//Print stats:
	logPrintf("\tIdealGasPomega[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule.name.c_str(), Emin, Emax, Emean);
}

void IdealGasPomega::getDensities(const DataRptr* logPomega, DataRptr* N, vector3<>& P) const
{	for(unsigned i=0; i<molecule.sites.size(); i++) N[i]=0;
	P = vector3<>(0,0,0);
	double& S = ((IdealGasPomega*)this)->S;
	S=0.0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr N_o = (quad.weight(o) * Nbulk) * exp(logPomega[o]); //contribution form this orientation
		//Accumulate N_o to each site dneisty with appropriate translations:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(rot*pos, 1.0, N_o, N[i]);
		//Accumulate contributions to the entropy:
		S += gInfo.dV*dot(N_o, logPomega[o]);
		//Accumulate the cell dipole moments:
		P += (rot * pMol) * integral(N_o);
	}
}

double IdealGasPomega::compute(const DataRptr* logPomega, const DataRptr* N, DataRptr* grad_N,
	const vector3<>& P, vector3<>& grad_P, const double Nscale, double& grad_Nscale) const
{	double PhiNI = 0.0;
	//Add contributions due to external potentials:
	for(unsigned i=0; i<molecule.sites.size(); i++)
		if(V[i])
		{	grad_N[i] += V[i];
			PhiNI += gInfo.dV*dot(N[i], V[i]);
		}
	//Contributions due to uniform electric field:
	grad_P -= Eexternal;
	PhiNI -= dot(Eexternal, P);
	//KE and mu:
	double invSite0mult = 1./molecule.sites[0]->positions.size();
	grad_N[0] -= mu * invSite0mult;
	PhiNI -= (T+mu)*integral(N[0])* invSite0mult;
	//Entropy (this part deals with Nscale explicitly, so need to increment grad_Nscale):
	grad_Nscale += T*S;
	PhiNI += Nscale*T*S;
	return PhiNI;
}

void IdealGasPomega::convertGradients(const DataRptr* logPomega, const DataRptr* N,
	const DataRptr* grad_N, vector3<> grad_P, DataRptr* grad_logPomega, const double Nscale) const
{	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr grad_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each grad_N in grad_N_o
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1.0, grad_N[i], grad_N_o);
		//Collect the contributions the entropy:
		grad_N_o += T*logPomega[o];
		//Collect the contribution from grad_P:
		grad_N_o += dot(grad_P, rot * pMol);
		//Propagate grad_N_o to grad_logPomega[o]:
		grad_logPomega[o] = (quad.weight(o) * Nbulk * Nscale) * exp(logPomega[o]) * grad_N_o;
	}
}

