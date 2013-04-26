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

#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/FluidComponent.h>
#include <fluid/Euler.h>


IdealGasPsiAlpha::IdealGasPsiAlpha(const FluidMixture* fluidMixture, const FluidComponent* comp, const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(comp->molecule.sites.size(),fluidMixture,comp), quad(quad), trans(trans)
{
}

void IdealGasPsiAlpha::initState(const DataRptr* Vex, DataRptr* psi, double scale, double Elo, double Ehi) const
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
	}
	//Print stats:
	logPrintf("\tIdealGasPsiAlpha[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule.name.c_str(), Emin, Emax, Emean);
	//Initialize the state (simply a constant factor times the potential):
	for(unsigned i=0; i<molecule.sites.size(); i++) psi[i] = (-scale/T)*Veff[i];
}

void IdealGasPsiAlpha::getDensities(const DataRptr* psi, DataRptr* N, vector3<>& P) const
{	for(unsigned i=0; i<molecule.sites.size(); i++) N[i]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr sum_psi; //the exponent in the boltzmann factor
		//Collect psi's from each site in this orientation:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1., psi[i], sum_psi);
		DataRptr N_o = quad.weight(o) * Nbulk * exp(sum_psi); //contribution from this orientation
		//Accumulate N_o to each site density with appropriate translations:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(rot*pos, 1., N_o, N[i]);
	}
}

double IdealGasPsiAlpha::compute(const DataRptr* psi, const DataRptr* N, DataRptr* Phi_N,
	const vector3<>& P, vector3<>& Phi_P, const double Nscale, double& Phi_Nscale) const
{	double PhiNI = 0.0;
	for(unsigned i=0; i<molecule.sites.size(); i++)
	{	DataRptr PhiNI_Ni = T*psi[i] + V[i];
		if(i==0) //Add terms associated with whole molecule (mu and KE) to first site
		{	double invSite0mult = 1./molecule.sites[0]->positions.size();
			PhiNI_Ni -= mu * invSite0mult;
			PhiNI -= T*integral(N[0]) * invSite0mult;
		}
		if(PhiNI_Ni)
		{	Phi_N[i] += PhiNI_Ni;
			PhiNI += gInfo.dV*dot(N[i], PhiNI_Ni);
		}
	}
	return PhiNI;
}

void IdealGasPsiAlpha::convertGradients(const DataRptr* psi, const DataRptr* N,
	const DataRptr* Phi_N, vector3<> Phi_P, DataRptr* Phi_psi, const double Nscale) const
{
	for(unsigned i=0; i<molecule.sites.size(); i++) Phi_psi[i]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr Phi_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each Phi_N in Phi_N_o
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1., Phi_N[i], Phi_N_o);
		//Calculate N_o again (with Nscale this time):
		DataRptr sum_psi;
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1., psi[i], sum_psi);
		DataRptr N_o = (quad.weight(o) * Nbulk * Nscale) * exp(sum_psi); //contribution from this orientation
		//Accumulate N_o * Phi_N_o into each Phi_psi with appropriate translations:
		DataRptr Phi_psi_term = N_o * Phi_N_o;
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(rot*pos, 1., Phi_psi_term, Phi_psi[i]);
	}
}
