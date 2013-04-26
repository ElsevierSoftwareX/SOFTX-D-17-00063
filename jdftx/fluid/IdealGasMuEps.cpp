/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <fluid/IdealGasMuEps.h>
#include <fluid/Euler.h>

IdealGasMuEps::IdealGasMuEps(const FluidMixture* fluidMixture, const FluidComponent* comp,  const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(4,fluidMixture,comp), quad(quad), trans(trans), pMol(molecule.getDipole())
{
}

void IdealGasMuEps::initState(const DataRptr* Vex, DataRptr* mueps, double scale, double Elo, double Ehi) const
{	DataRptrCollection Veff(molecule.sites.size()); nullToZero(Veff, gInfo);
	for(int k=0; k<4; k++) mueps[k]=0;
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
		//Add contributions to the state:
		vector3<> pVec = rot * pMol;
		mueps[0] += quad.weight(o) * Emolecule;
		for(int k=0; k<3; k++)
			mueps[k+1] += (pVec[k]*quad.weight(o)) * Emolecule;
	}
	//Set the scale factor:
	for(int k=0; k<4; k++) mueps[k] *= (-scale/T);
	//Print stats:
	logPrintf("\tIdealGasMuEps[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule.name.c_str(), Emin, Emax, Emean);
}

void IdealGasMuEps::getDensities(const DataRptr* mueps, DataRptr* N, vector3<>& P) const
{	for(unsigned i=0; i<molecule.sites.size(); i++) N[i]=0;
	P = vector3<>(0,0,0);
	double& S = ((IdealGasMuEps*)this)->S;
	S=0.0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		vector3<> pVec = rot * pMol;
		DataRptr sum_mueps; //the exponent in the boltzmann factor
		sum_mueps += mueps[0];
		for(int k=0; k<3; k++)
			sum_mueps += pVec[k]*mueps[k+1];
		DataRptr N_o = quad.weight(o) * Nbulk * exp(sum_mueps); //contribution from this orientation
		//Accumulate N_o to each site density with appropriate translations:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(rot*pos, 1.0, N_o, N[i]);
		//Accumulate contributions to the entropy:
		S += gInfo.dV*dot(N_o, mueps[0]);
		for(int k=0; k<3; k++)
			S += pVec[k] * gInfo.dV*dot(N_o, mueps[k+1]);
		//Accumulate the cell dipole moments:
		P += pVec * integral(N_o);
	}
}

double IdealGasMuEps::compute(const DataRptr* mueps, const DataRptr* N, DataRptr* Phi_N,
	const vector3<>& P, vector3<>& Phi_P, const double Nscale, double& Phi_Nscale) const
{	double PhiNI = 0.0;
	//Add contributions due to external potentials:
	for(unsigned i=0; i<molecule.sites.size(); i++)
		if(V[i])
		{	Phi_N[i] += V[i];
			PhiNI += gInfo.dV*dot(N[i], V[i]);
		}
	//Contributions due to uniform electric field:
	Phi_P -= Eexternal;
	PhiNI -= dot(Eexternal, P);
	//KE and mu:
	double invSite0mult = 1./molecule.sites[0]->positions.size();
	Phi_N[0] -= mu * invSite0mult;
	PhiNI -= (T+mu)*integral(N[0]) * invSite0mult;
	//Entropy (this part deals with Nscale explicitly, so need to increment Phi_Nscale):
	Phi_Nscale += T*S;
	PhiNI += Nscale*T*S;
	return PhiNI;
}

void IdealGasMuEps::convertGradients(const DataRptr* mueps, const DataRptr* N,
	const DataRptr* Phi_N, vector3<> Phi_P, DataRptr* Phi_mueps, const double Nscale) const
{	for(int k=0; k<4; k++) Phi_mueps[k]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		vector3<> pVec = rot*vector3<>(0,0,1);
		DataRptr Phi_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each Phi_N in Phi_N_o
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1.0, Phi_N[i], Phi_N_o);
		//Collect the contributions the entropy:
		Phi_N_o += T*mueps[0];
		for(int k=0; k<3; k++)
			Phi_N_o += (T*pVec[k])*mueps[k+1];
		//Collect the contribution from Phi_P:
		Phi_N_o += dot(Phi_P, pVec);
		//Calculate N_o again:
		DataRptr sum_mueps;
		sum_mueps += mueps[0];
		for(int k=0; k<3; k++)
			sum_mueps += pVec[k]*mueps[k+1];
		DataRptr N_o = (quad.weight(o) * Nbulk * Nscale) * exp(sum_mueps);
		//Accumulate N_o * Phi_N_o into each component of Phi_mueps with appropriate weights:
		DataRptr Phi_mueps_term = N_o * Phi_N_o;
		Phi_mueps[0] += Phi_mueps_term;
		for(int k=0; k<3; k++)
			Phi_mueps[k+1] += pVec[k] * Phi_mueps_term;
	}
}
