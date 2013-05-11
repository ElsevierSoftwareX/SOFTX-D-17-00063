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

#include <fluid/FluidMixture.h>
#include <fluid/MixedFMT.h>
#include <fluid/IdealGas.h>
#include <fluid/Fex.h>
#include <electronic/operators.h>
#include <core/EnergyComponents.h>
#include <core/DataMultiplet.h>

extern string rigidMoleculeCDFT_ScalarEOSpaper;

FluidMixture::FluidMixture(const GridInfo& gInfo, const double T)
: gInfo(gInfo), T(T), verboseLog(false), Qtol(1e-12), nIndepIdgas(0), nDensities(0), polarizable(false)
{
	logPrintf("Initializing fluid mixture at T=%lf K ...\n", T/Kelvin);
	Citations::add("Rigid-molecule density functional theory framework", rigidMoleculeCDFT_ScalarEOSpaper);
}

FluidMixture::~FluidMixture()
{	
}

void FluidMixture::initialize(double p, double epsBulkOverride)
{	logPrintf("Adjusting fluid pressure to p=%lf bar\n", p/Bar);
	//Compute the maximum possible density (core packed limit)
	double Nguess=0., n3=0.;
	assert(component.size());
	for(const FluidComponent* c: component)
	{	Nguess += c->Nbulk;
		n3 += c->Nbulk * c->molecule.getVhs();
	}
	const double mulStep = 0.99;
	double Nstart = Nguess;
	if(n3 > mulStep) Nstart *= mulStep/n3; //ensure that mixtyure doesn;t exceed packing limit
	double pTest = compute_p(Nstart);
	//Find an interval of N that brackets P:
	double Nlo, Nhi;
	if(pTest > p)
	{	Nlo = Nstart;
		do
		{	Nhi = Nlo;
			Nlo = mulStep*Nhi;
			pTest = compute_p(Nlo);
		}
		while(pTest>p);
	}
	else
	{	Nhi = Nstart;
		do
		{	Nlo = Nhi;
			Nhi = Nlo/mulStep;
			pTest = compute_p(Nhi);
		}
		while(pTest<p);
	}
	//Bisect on N to get the pressure:
	double Ntot;
	do
	{	Ntot = 0.5*(Nlo + Nhi);
		pTest = compute_p(Ntot);
		if(pTest<p) Nlo = Ntot;
		else Nhi = Ntot;
	}
	while(Nhi-Nlo>1e-12*Nhi);
	//Set the chemical potentials and bulk densities for each component:
	std::vector<double> Nmol(component.size()), Phi_Nmol(component.size());
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = (Ntot/Nguess)*component[ic]->Nbulk;
	computeUniformEx(Nmol, Phi_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		c.idealGas->Nbulk = Nmol[ic];
		c.idealGas->mu = Phi_Nmol[ic];
		logPrintf("   Component '%s' at bulk density %le bohr^-3\n", c.molecule.name.c_str(), Nmol[ic]);
	}
	this->p = p;
	
	//Determine dipole correlation factor:
	double epsBulk = epsBulkOverride;
	if(!epsBulk)
	{	epsBulk = 1.;
		for(const auto& c: component)
			epsBulk += (c->idealGas->get_Nbulk()/c->pureNbulk(T)) * (c->epsBulk - 1.);
	}
	double chiMF = 0.;
	for(const FluidComponent* c: component)
		chiMF += c->idealGas->get_Nbulk() * (c->molecule.getDipole().length_squared()/(3.*T) + c->molecule.getAlphaTot());
	Ceps = (epsBulk>1.) ? (1./(4.*M_PI*chiMF) - 1./(epsBulk-1.)) : 0.;
	logPrintf("   Local polarization-density correlation factor, Ceps = %lg\n", Ceps);
	
	//Initialize preconditioners:
	Kindep.resize(component.size());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	//Determine second derivative of total excess functional:
		double Nbulk = Nmol[ic];
		const double dNfac = 1e-4;
		Nmol[ic]=Nbulk*(1.+dNfac); computeUniformEx(Nmol, Phi_Nmol); double Phi_Np = Phi_Nmol[ic];
		Nmol[ic]=Nbulk*(1.-dNfac); computeUniformEx(Nmol, Phi_Nmol); double Phi_Nm = Phi_Nmol[ic];
		Nmol[ic]=Nbulk;
		double NNPhi_NN = Nbulk * (Phi_Np - Phi_Nm) / (2.*dNfac);
		if(NNPhi_NN < -0.5*Nbulk*T) NNPhi_NN = -0.5*Nbulk*T;
		//Set preconditioner:
		Kindep[ic] = 1./(gInfo.dV * (Nbulk*T + NNPhi_NN));
	}
	if(polarizable)
	{	double chiPol = 0.;
		for(const FluidComponent* c: component)
			chiPol += c->idealGas->get_Nbulk() * c->molecule.getAlphaTot();
		Keps = 1./(gInfo.dV * chiPol);
	}
}

const std::vector<const FluidComponent*>& FluidMixture::getComponents() const
{	return component;
}


void FluidMixture::addComponent(FluidComponent* comp)
{	component.push_back(comp);
	//Set the offsets for this component:
	comp->offsetIndep = nIndepIdgas;
	comp->offsetDensity = nDensities;
	//Update the totals, which become the offset for the next component
	nIndepIdgas += comp->idealGas->nIndep;
	nDensities += comp->molecule.sites.size();
	//Update the polarizable flag:
	polarizable |= bool(comp->molecule.getAlphaTot());
}

void FluidMixture::addFmix(const Fmix* fmix)
{	fmixArr.push_back(fmix);
}


void FluidMixture::initState(double scale, double Elo, double Ehi)
{	//Compute the effective nonlinear coupling potential for the uniform fluid:
	DataRptrCollection Vex(nDensities);
	{	//Get the uniform fluid site densities:
		std::vector<double> Nbulk(nDensities);
		for(const FluidComponent* c: component)
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
				Nbulk[c->offsetDensity + i] = c->idealGas->get_Nbulk() * c->molecule.sites[i]->positions.size();
		DataGptrCollection Ntilde(nDensities);
		nullToZero(Ntilde, gInfo);
		for(unsigned i=0; i<nDensities; i++) Ntilde[i]->setGzero(Nbulk[i]);
		//Set Vex to the difference between the potentials returned by Fmix::compute() and Fmix::computeUniform():
		DataGptrCollection Phi_Ntilde(nDensities);
		std::vector<double> Phi_Nbulk(nDensities);
		for(const Fmix* fmix: fmixArr)
		{	fmix->compute(Ntilde, Phi_Ntilde);
			fmix->computeUniform(Nbulk, Phi_Nbulk);
		}
		for(unsigned i=0; i<nDensities; i++)
			if(Phi_Ntilde[i]) Vex[i] = Jdag(Phi_Ntilde[i]) - Phi_Nbulk[i];
		//Add electrostatic coupling:
		if(rhoExternal)
		{	DataGptr dExternal = (-4*M_PI)*Linv(O(rhoExternal));
			for(const FluidComponent* c: component)
				for(unsigned i=0; i<c->molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c->molecule.sites[i]);
					if(s.chargeKernel) Vex[c->offsetDensity+i] += I(s.chargeKernel * dExternal);
				}
		}
	}
	state.assign(get_nIndep(), 0);
	//Call initState for each component
	logPrintf("\n----- FluidMixture::initState() -----\n");
	for(const FluidComponent* c: component)
		c->idealGas->initState(&Vex[c->offsetDensity], &state[c->offsetIndep], scale, Elo, Ehi);
	//Initialize polarizability state if necessary:
	if(polarizable)
	{	DataGptrVec PMFtilde; DataGptr rhoMF;
		for(const FluidComponent* c: component)
		{	DataRptrCollection N(c->molecule.sites.size()); DataRptrVec P;
			c->idealGas->getDensities(&state[c->offsetIndep], &N[0], P);
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c->molecule.sites[i]);
				if(s.chargeKernel) rhoMF += s.chargeKernel(0) * (c->molecule.mfKernel * J(N[i]));
			}
			if(P) PMFtilde += c->molecule.mfKernel * J(P);
		}
		DataGptrVec epsMFtilde;
		if(rhoMF) epsMFtilde += gradient(Linv((4*M_PI)*O(rhoMF)));
		if(PMFtilde) epsMFtilde += (4*M_PI*Ceps) * PMFtilde;
		nullToZero(epsMFtilde, gInfo);
		for(unsigned k=nIndepIdgas; k<get_nIndep(); k++) state[k] = scale * I(epsMFtilde[k-nIndepIdgas]);
	}
	logPrintf("\n");
}

void FluidMixture::loadState(const char* filename)
{	nullToZero(state, gInfo, get_nIndep());
	loadFromFile(state, filename);
}

void FluidMixture::saveState(const char* filename) const
{	saveToFile(state, filename);
}

FluidMixture::Outputs::Outputs(DataRptrCollection* N, vector3<>* electricP,
	DataGptr* Phi_rhoExternal, DataRptrCollection* psiEff, EnergyComponents* Phi)
:N(N),electricP(electricP),Phi_rhoExternal(Phi_rhoExternal),psiEff(psiEff),Phi(Phi)
{
}


//! Compute the total charge of a set of components: original number of molecules N0 and charge per molecule Q
//! given as the vector of pairs N0Q, where the actual number of molecules of each component is N = N0 exp(-Q betaV)
//! Also return the gradient w.r.t betaV
double Qtot(double betaV, double& Qtot_betaV, const std::vector<std::pair<double,double> > N0Q,
			const std::vector<string>* names=0, const GridInfo* gInfo=0, const bool verboseLog=0)
{	double Qsum=0.0, Qsum_betaV=0.0;
	for(unsigned i=0; i<N0Q.size(); i++)
	{	
		double N0 = N0Q[i].first;
		double Q = N0Q[i].second;
		double N = N0*exp(-Q*betaV);
		if(verboseLog) logPrintf("%s N0: %le Q: %le N: %le\n", (*names)[i].c_str(), N0, Q, N);
		Qsum += Q*N;
		Qsum_betaV += -Q*Q*N;
	}
	Qtot_betaV = Qsum_betaV;
	return Qsum;
}


double FluidMixture::operator()(const DataRptrCollection& indep, DataRptrCollection& Phi_indep, Outputs outputs) const
{	
	//logPrintf("indep.size: %d nIndep: %d\n",indep.size(),nIndep);
	assert(indep.size()==get_nIndep());

	//---------- Compute site densities from the independent variables ---------
	DataGptrCollection Ntilde(nDensities); //site densities (in reciprocal space)
	std::vector<DataGptrVec> Ptilde(component.size()); //polarization densities
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		DataRptrCollection N(c.molecule.sites.size()); DataRptrVec P;
		c.idealGas->getDensities(&indep[c.offsetIndep], &N[0], P);
		if(P) { Ptilde[ic] = J(P); P = 0; } //Collect mean field polarization density
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
		{	//Replace negative densities with 0:
			double Nmin, Nmax;
			callPref(eblas_capMinMax)(gInfo.nr, N[i]->dataPref(), Nmin, Nmax, 0.);
			//store site densities in fourier space
			Ntilde[c.offsetDensity+i] = J(N[i]);
			N[i] = 0; //Really skimp on memory!
		}
	}

	//----------- Handle density constraints ------------
	std::vector<double> Nscale(component.size(), 1.0); //density scale factor that satisfies the constraint
	std::vector<double> Nscale_Qfixed(component.size(), 0.); //derivative of Nscale w.r.t the fixed charge
	std::vector<std::vector<double> > Nscale_N0(component.size(), std::vector<double>(component.size(),0.0)); //jacobian of Nscale w.r.t the uncorrected molecule counts
	std::vector<string> names; //list of molecule names
	
	//Find fixed N and charged species:
	double Qfixed = 0.0;
	if(rhoExternal) Qfixed += integral(rhoExternal);
	std::vector<std::pair<double,double> > N0Q;
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		double Qmolecule = c.molecule.getCharge();
		if(c.Nnorm>0 || Qmolecule)
		{	double N0 = integral(Ntilde[c.offsetDensity])/c.molecule.sites[0]->positions.size();
			if(c.Nnorm>0)
			{	Nscale[ic] = c.Nnorm/N0;
				Nscale_N0[ic][ic] = -c.Nnorm/pow(N0,2);
				Qfixed += Qmolecule*c.Nnorm;
			}
			else
			{	N0Q.push_back(std::make_pair(N0, Qmolecule));
				names.push_back(c.molecule.name);
			}
		}
	}
	//Find the betaV (see Qtot()) that makes the unit cell neutral
	if(N0Q.size()==0)
	{	if(fabs(Qfixed)>fabs(Qtol))
			die("Unit cell has a fixed net charge %le,"
				"and there are no free charged species to neutralize it.\n", Qfixed);
	}
	else
	{	double Qprime, Qcell, betaV=0.0;
		if(Qfixed+Qtot(-HUGE_VAL,Qprime,N0Q)<0)
			die("Unit cell will always have a net negative charge (no free positive charges).\n")
		if(Qfixed+Qtot(+HUGE_VAL,Qprime,N0Q)>0)
			die("Unit cell will always have a net positive charge (no free negative charges).\n")

		for(int iter=0; iter<10; iter++) //while(1)
		{	Qcell = Qfixed+Qtot(betaV,Qprime,N0Q,&names,&gInfo,verboseLog);
			if(verboseLog) logPrintf("betaV = %le, Qcell = %le, Qprime = %le\n", betaV, Qcell, Qprime);
			if(fabs(Qcell)<fabs(Qtol)) break;
			if(std::isnan(Qcell)) die("NaN encountered in Q convergence.\n")
			betaV -= Qcell/Qprime; //Newton-Raphson update
		}
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& ci = *component[ic];
			double Qi = ci.molecule.getCharge();
			if(Qi && ci.Nnorm<=0)
			{	Nscale[ic] = exp(-Qi*betaV);
				Nscale_Qfixed[ic] = exp(-Qi*betaV) * Qi / Qprime;
				for(unsigned jc=0; jc<component.size(); jc++)
				{	const FluidComponent& cj = *component[jc];
					double Qj = cj.molecule.getCharge();
					if(Qj && cj.Nnorm<=0)
						Nscale_N0[ic][jc] += Qi*Qj*exp(-(Qi+Qj)*betaV)/Qprime;
				}
			}
		}
	}
	std::vector<double> Phi_Nscale(component.size(), 0.0); //accumulate explicit derivatives w.r.t Nscale here

	//Apply the scale factors to the site densities
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
			Ntilde[c.offsetDensity+i] *= Nscale[ic];
		if(Ptilde[ic]) Ptilde[ic] *= Nscale[ic];
	}

	EnergyComponents Phi; //the grand free energy (with component information)
	DataGptrCollection Phi_Ntilde(nDensities); //gradients (functional derivative) w.r.t reciprocal space site densities
	std::vector<DataGptrVec> Phi_Ptilde(component.size()); //functional derivative w.r.t polarization density
	DataGptrVec Phi_epsMF; //functional derivative w.r.t mean field electric field
	
	//--------- Compute the (scaled) mean field coulomb interaction --------
	{	DataGptr rho; //total charge density
		DataGptr rhoMF; //effective charge density for mean-field term
		DataGptrVec PMFtilde; //effective polarization density for mean field correction
		bool needRho = rhoExternal || outputs.Phi_rhoExternal;
		
		DataGptrVec epsMF = polarizable ? J(DataRptrVec(&indep[nIndepIdgas])) : 0; //mean field electric field
		
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *component[ic];
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.chargeKernel)
				{	if(needRho) rho += s.chargeKernel * Ntilde[c.offsetDensity+i];
					rhoMF += s.chargeKernel(0) * (c.molecule.mfKernel * Ntilde[c.offsetDensity+i]);
				}
				//Polarization contributions:
				if(s.polKernel)
				{	
					#define Polarization_Compute_Pi_Ni \
						DataRptrVec Pi = I(s.polKernel(0)*(c.molecule.mfKernel*epsMF) + (rhoExternal ? gradient(s.polKernel*Linv((gInfo.detR*4*M_PI)*rhoExternal)) : 0)); \
						DataRptr Ni = I(Ntilde[c.offsetDensity+i]);
					Polarization_Compute_Pi_Ni
					
					DataRptr Phi_Ni = 0.5*lengthSquared(Pi);
					Phi["Apol"] += gInfo.dV * dot(Ni, Phi_Ni);
					//Derivative contribution to site densities:
					Phi_Ntilde[c.offsetDensity+i] += Idag(Phi_Ni); Phi_Ni=0;
					//Update contributions to bound charge:
					DataGptrVec NPtilde = J(Ni * Pi); Pi=0; Ni=0;
					DataGptr divNPbar;
					if(needRho)
					{	divNPbar = s.polKernel*divergence(NPtilde);
						rho -= divNPbar;
					}
					DataGptrVec NPbarMF = s.polKernel(0)*(c.molecule.mfKernel*NPtilde); NPtilde=0;
					rhoMF -= divergence(NPbarMF);
					PMFtilde += NPbarMF;
					Phi_epsMF += gInfo.nr * NPbarMF;
				}
			}
			if(Ptilde[ic]) PMFtilde += c.molecule.mfKernel * Ptilde[ic];
		}
		
		if(rhoMF)
		{	//External charge interaction:
			DataGptr Phi_rho;
			if(needRho)
			{	if(rhoExternal)
				{	DataGptr OdExternal = O(-4*M_PI*Linv(O(rhoExternal)));
					Phi["ExtCoulomb"] += dot(rho, OdExternal);
					Phi_rho += OdExternal;
				}
				if(outputs.Phi_rhoExternal) *outputs.Phi_rhoExternal = -4*M_PI*Linv(O(rho));
			}
		
			//Mean field contributions:
			DataGptr Phi_rhoMF;
			{	DataGptr OdMF = O(-4*M_PI*Linv(O(rhoMF))); //mean-field electrostatic potential
				Phi["Coulomb"] += 0.5*dot(rhoMF, OdMF);
				Phi_rhoMF += OdMF;
			}
			
			//Polarization density interactions:
			DataGptrVec Phi_PMFtilde;
			if(PMFtilde)
			{	Phi_PMFtilde = (-4*M_PI*Ceps) * O(PMFtilde);
				Phi["LPDA"] += 0.5 * dot(Phi_PMFtilde, PMFtilde);
				
				//Corrections for net dipole in cell:
				vector3<> P0; if(PMFtilde) for(int k=0; k<3; k++) P0[k] = PMFtilde[k]->getGzero();
				if(outputs.electricP) *outputs.electricP = P0 * gInfo.detR;
				vector3<> Phi_P0 = (4*M_PI*gInfo.detR) * P0;
				Phi["PsqCell"] += 0.5 * dot(Phi_P0, P0);
				
				//External electric field interactions:
				Phi["ExtCoulomb"] -= gInfo.detR * dot(Eexternal, P0); //external uniform electric field
				Phi_P0 -= gInfo.detR * Eexternal;

				for(int k=0; k<3; k++) Phi_PMFtilde[k]->setGzero(Phi_PMFtilde[k]->getGzero() + Phi_P0[k]);
			}
			
			//Propagate gradients:
			for(unsigned ic=0; ic<component.size(); ic++)
			{	const FluidComponent& c = *component[ic];
				for(unsigned i=0; i<c.molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c.molecule.sites[i]);
					if(s.chargeKernel)
					{	if(Phi_rho) Phi_Ntilde[c.offsetDensity+i] += (1./gInfo.dV) * (s.chargeKernel * Phi_rho);
						Phi_Ntilde[c.offsetDensity+i] += (s.chargeKernel(0)/gInfo.dV) * (c.molecule.mfKernel * Phi_rhoMF);
					}
					//Polarization contributions:
					if(s.polKernel)
					{	DataRptrVec Phi_NP = Jdag( s.polKernel(0)*(c.molecule.mfKernel*(Phi_PMFtilde + gradient(Phi_rhoMF)))
							+ (needRho ? gradient(s.polKernel * Phi_rho) : 0) );
						//propagate gradients from NP to N, epsMF and rhoExternal:
						Polarization_Compute_Pi_Ni
						#undef Polarization_Compute_Pi_Ni
						// --> via Ni
						DataRptr Phi_Ni; for(int k=0; k<3; k++) Phi_Ni += Phi_NP[k]*Pi[k];
						Phi_Ntilde[c.offsetDensity+i] += (1./gInfo.dV) * Idag(Phi_Ni); Phi_Ni=0;
						// --> via Pi
						DataGptrVec Phi_PiTilde = Idag(Phi_NP * Ni); Phi_NP=0;
						Phi_epsMF += (s.polKernel(0)/gInfo.dV)*(c.molecule.mfKernel*Phi_PiTilde);
					}
				}
				if(Ptilde[ic]) Phi_Ptilde[ic] += (1./gInfo.dV) * (c.molecule.mfKernel * Phi_PMFtilde);
			}
		}
	}
	
	//--------- Hard sphere mixture and bonding -------------
	{	//Compute the FMT weighted densities:
		DataGptr n0tilde, n1tilde, n2tilde, n3tilde, n1vTilde, n2mTilde;
		std::vector<DataRptr> n0mol(component.size(), 0); //partial n0 for molecules that need bonding corrections
		std::vector<int> n0mult(component.size(), 0); //number of sites which contribute to n0 for each molecule
		std::vector<std::map<double,int> > bond(component.size()); //sets of bonds for each molecule
		bool bondsPresent = false; //whether bonds are present for any molecule
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *component[ic];
			bond[ic] = c.molecule.getBonds();
			DataGptr n0molTilde;
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.Rhs)
				{	const DataGptr& Nsite = Ntilde[c.offsetDensity+i];
					n0mult[ic] += s.positions.size();
					n0molTilde += s.w0  * Nsite;
					n1tilde    += s.w1  * Nsite;
					n2tilde    += s.w2  * Nsite;
					n3tilde    += s.w3  * Nsite;
					n1vTilde   += s.w1v * Nsite;
					n2mTilde   += s.w2m * Nsite;
				}
			}
			if(n0molTilde) n0tilde += n0molTilde;
			if(bond[ic].size())
			{	n0mol[ic] = I(n0molTilde);
				bondsPresent = true;
			}
		}
		if(n0tilde) //at least one sphere in the mixture
		{	DataRptr n0 = I(n0tilde); n0tilde=0;
			DataRptr n1 = I(n1tilde); n1tilde=0;
			DataRptr n2 = I(n2tilde); n2tilde=0;
			DataRptr Phi_n0, Phi_n1, Phi_n2; DataGptr Phi_n3tilde, Phi_n1vTilde, Phi_n2mTilde;
			//Compute the sphere mixture free energy:
			Phi["MixedFMT"] += T * PhiFMT(n0, n1, n2, n3tilde, n1vTilde, n2mTilde,
				Phi_n0, Phi_n1, Phi_n2, Phi_n3tilde, Phi_n1vTilde, Phi_n2mTilde);
			//Bonding corrections if required
			if(bondsPresent)
			{	for(unsigned ic=0; ic<component.size(); ic++)
				{	const FluidComponent& c = *component[ic];
					DataRptr Phi_n0mol;
					for(const auto& b: bond[ic])
						Phi["Bonding"] += T * PhiBond(b.first, b.second*1./n0mult[ic],
							n0mol[ic], n2, n3tilde, Phi_n0mol, Phi_n2, Phi_n3tilde);
					if(Phi_n0mol)
					{	//Propagate gradient w.r.t n0mol[ic] to the site densities:
						DataGptr Phi_n0molTilde = Idag(Phi_n0mol);
						for(unsigned i=0; i<c.molecule.sites.size(); i++)
						{	const Molecule::Site& s = *(c.molecule.sites[i]);
							if(s.Rhs)
								Phi_Ntilde[c.offsetDensity+i] += T * (s.w0  * Phi_n0molTilde);
						}
					}
				}
			}
			//Accumulate gradients w.r.t weighted densities to site densities:
			DataGptr Phi_n0tilde = Idag(Phi_n0); Phi_n0=0;
			DataGptr Phi_n1tilde = Idag(Phi_n1); Phi_n1=0;
			DataGptr Phi_n2tilde = Idag(Phi_n2); Phi_n2=0;
			for(const FluidComponent* c: component)
			{	for(unsigned i=0; i<c->molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c->molecule.sites[i]);
					if(s.Rhs)
					{	DataGptr& Phi_Nsite = Phi_Ntilde[c->offsetDensity+i];
						Phi_Nsite += T * (s.w0  * Phi_n0tilde);
						Phi_Nsite += T * (s.w1  * Phi_n1tilde);
						Phi_Nsite += T * (s.w2  * Phi_n2tilde);
						Phi_Nsite += T * (s.w3  * Phi_n3tilde);
						Phi_Nsite += T * (s.w1v * Phi_n1vTilde);
						Phi_Nsite += T * (s.w2m * Phi_n2mTilde);
					}
				}
			}
		}
	}

	//---------- Excess functionals --------------
	for(const FluidComponent* c: component) if(c->fex)
		Phi["Fex("+c->molecule.name+")"] += c->fex->compute(&Ntilde[c->offsetDensity], &Phi_Ntilde[c->offsetDensity]);

	//--------- Mixing functionals --------------
	for(const Fmix* fmix: fmixArr)
		Phi["Fmix("+fmix->getName()+")"] += fmix->compute(Ntilde, Phi_Ntilde);

	//--------- PhiNI ---------
	nullToZero(Phi_Ntilde, gInfo);
	if(outputs.N) outputs.N->resize(nDensities);
	//Put the site densities and gradients back in real space
	DataRptrCollection N(nDensities);
	DataRptrCollection Phi_N(nDensities);
	for(unsigned i=0; i<nDensities; i++)
	{	N[i] = I(Ntilde[i]); Ntilde[i]=0;
		Phi_N[i] = Jdag(Phi_Ntilde[i]); Phi_Ntilde[i] = 0;
		if(outputs.N) (*outputs.N)[i] = N[i]; //Link site-density to return pointer if necessary
	}
	//Put the polarization densities and gradients back in real space:
	std::vector<DataRptrVec> P(component.size());
	std::vector<DataRptrVec> Phi_P(component.size());
	for(unsigned ic=0; ic<component.size(); ic++)
		if(Ptilde[ic])
		{	P[ic] = I(Ptilde[ic]); Ptilde[ic] = 0;
			Phi_P[ic] = Jdag(Phi_Ptilde[ic]); Phi_Ptilde[ic] = 0;
		}
	//Estimate psiEff based on gradients, if requested
	if(outputs.psiEff)
	{	outputs.psiEff->resize(nDensities);
		for(const FluidComponent* c: component)
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
			{	DataRptr& psiCur = outputs.psiEff->at(c->offsetDensity+i);
				psiCur = Phi_N[c->offsetDensity+i] + c->idealGas->V[i];
				if(i==0) psiCur -= c->idealGas->mu / c->molecule.sites[0]->positions.size();
				psiCur *= (-1./T);
			}
	}
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		Phi["PhiNI("+c.molecule.name+")"] +=
			c.idealGas->compute(&indep[c.offsetIndep], &N[c.offsetDensity], &Phi_N[c.offsetDensity], Nscale[ic], Phi_Nscale[ic]);

		//Fixed N correction to entropy:
		if(Nscale[ic]!=1.0)
		{	double deltaTs = T*log(Nscale[ic]) / c.molecule.sites[0]->positions.size();
			Phi_N[c.offsetDensity] += deltaTs;
			Phi["PhiNI("+c.molecule.name+")"] += integral(N[c.offsetDensity])*deltaTs;
		}
	}
	//Add in the implicit contributions to Phi_Nscale
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& ci = *component[ic];
		bool anyNonzero=false;
		for(unsigned jc=0; jc<component.size(); jc++)
			if(Nscale_N0[ic][jc])
				anyNonzero=true;
		if(anyNonzero)
		{	if(P[ic] && Phi_P[ic]) Phi_Nscale[ic] += gInfo.dV * dot(P[ic], Phi_P[ic])/ Nscale[ic];
			for(unsigned i=0; i<ci.molecule.sites.size(); i++)
				Phi_Nscale[ic] += gInfo.dV*dot(N[ci.offsetDensity+i], Phi_N[ci.offsetDensity+i])/ Nscale[ic];
		}
	}
	//Propagate gradients from Nscale to N:
	for(unsigned jc=0; jc<component.size(); jc++)
	{	const FluidComponent& cj = *component[jc];
		double Phi_Ncontrib = 0.0;
		for(unsigned ic=0; ic<component.size(); ic++)
			if(Nscale_N0[ic][jc])
				Phi_Ncontrib += Phi_Nscale[ic] * Nscale_N0[ic][jc];
		if(Phi_Ncontrib)
			Phi_N[cj.offsetDensity] += Phi_Ncontrib / (Nscale[jc] * cj.molecule.sites[0]->positions.size());
	}

	//Propagate gradients from Phi_N and Phi_P to Phi_indep
	Phi_indep.resize(get_nIndep());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		c.idealGas->convertGradients(&indep[c.offsetIndep], &N[c.offsetDensity],
			&Phi_N[c.offsetDensity], Phi_P[ic], &Phi_indep[c.offsetIndep], Nscale[ic]);
	}
	for(unsigned k=nIndepIdgas; k<get_nIndep(); k++) Phi_indep[k] = Jdag(Phi_epsMF[k-nIndepIdgas]);
	
	//Propagate gradients from Nscale to Qfixed / rhoExternal (Natural G=0 solution)
	if(outputs.Phi_rhoExternal)
	{	double Phi_Qfixed = 0.;
		for(unsigned ic=0; ic<component.size(); ic++)
			Phi_Qfixed += Phi_Nscale[ic] * Nscale_Qfixed[ic];
		nullToZero(*outputs.Phi_rhoExternal, gInfo);
		(*outputs.Phi_rhoExternal)->setGzero(Phi_Qfixed);
	}
	
	Phi["+pV"] += p * gInfo.detR; //background correction

	if(verboseLog) Phi.print(globalLog, true, "\t\t\t\t%15s = %25.16lf\n");
	if(outputs.Phi) *(outputs.Phi) = Phi;
	
	Phi_indep *= gInfo.dV; //convert functional derivative to partial derivative
	return Phi;
}

double FluidMixture::getFreeEnergy(Outputs outputs) const
{	DataRptrCollection Phi_state;
	return (*this)(state, Phi_state, outputs);
	//calls operator() and provides whichever outputs are requested
	//before calling must create an Outputs structure with the required quantities allocated.
}

void FluidMixture::step(const DataRptrCollection& dir, double alpha)
{	axpy(alpha, dir, state);
}

double FluidMixture::compute(DataRptrCollection* grad)
{	DataRptrCollection tempGrad;
	return (*this)(state, grad ? *grad : tempGrad, Outputs());
}

DataRptrCollection FluidMixture::precondition(const DataRptrCollection& grad)
{	DataRptrCollection Kgrad(get_nIndep());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned k=c.offsetIndep; k<c.offsetIndep+c.idealGas->nIndep; k++)
			Kgrad[k] = Kindep[ic]*grad[k];
	}
	for(unsigned k=nIndepIdgas; k<get_nIndep(); k++)
		Kgrad[k] = Keps * grad[k];
	return Kgrad;
}


double FluidMixture::computeUniformEx(const std::vector<double>& Nmol, std::vector<double>& Phi_Nmol) const
{	//--------- Compute the site densities ----------
	std::vector<double> N(nDensities);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
			N[c.offsetDensity + i] += Nmol[ic] * c.molecule.sites[i]->positions.size();
	}

	EnergyComponents phi; std::vector<double> Phi_N(nDensities);

	//-------- Mean field coulomb -------- (No contribution in uniform fluid)

	//-------- Hard sphere/bonding -----------
	double n0=0.0, n1=0.0, n2=0.0, n3=0.0;
	std::vector<double> n0mol(component.size(), 0.0); //partial n0 for molecules that need bonding corrections
	std::vector<int> n0mult(component.size(), 0); //number of sites which contribute to n0 for each molecule
	std::vector<std::map<double,int> > bond(component.size()); //sets of bonds for each molecule
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		bond[ic] = c.molecule.getBonds();
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
		{	const Molecule::Site& s = *(c.molecule.sites[i]);
			if(s.Rhs)
			{	double Nsite = N[c.offsetDensity+i];
				n0mult[ic] += s.positions.size();
				n0mol[ic]  += s.w0(0)  * Nsite;
				n1         += s.w1(0)  * Nsite;
				n2         += s.w2(0)  * Nsite;
				n3         += s.w3(0)  * Nsite;
			}
		}
		n0 += n0mol[ic];
	}
	if(n0)
	{	double Phi_n0=0.0, Phi_n1=0.0, Phi_n2=0.0, Phi_n3=0.0;
		phi["MixedFMT"] = T * phiFMTuniform(n0, n1, n2, n3, Phi_n0, Phi_n1, Phi_n2, Phi_n3);
		//Bonding corrections:
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *component[ic];
			double Phi_n0mol=0.0;
			for(const auto& b: bond[ic])
				phi["Bonding"] += T * phiBondUniform(b.first, b.second*1.0/n0mult[ic],
					n0mol[ic], n2, n3, Phi_n0mol, Phi_n2, Phi_n3);
			if(Phi_n0mol)
			{	//Propagate gradient w.r.t n0mol[ic] to the site densities:
				for(unsigned i=0; i<c.molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c.molecule.sites[i]);
					if(s.Rhs) Phi_N[c.offsetDensity+i] += T * (s.w0(0)  * Phi_n0mol);
				}
			}
		}
		//Convert FMT weighted gradients to site density gradients:
		for(const FluidComponent* c: component)
		{	for(unsigned i=0; i<c->molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c->molecule.sites[i]);
				if(s.Rhs)
				{	double& Phi_Nsite = Phi_N[c->offsetDensity+i];
					Phi_Nsite += T * (s.w0(0)  * Phi_n0);
					Phi_Nsite += T * (s.w1(0)  * Phi_n1);
					Phi_Nsite += T * (s.w2(0)  * Phi_n2);
					Phi_Nsite += T * (s.w3(0)  * Phi_n3);
				}
			}
		}
	}

	//---------- Excess functionals --------------
	for(const FluidComponent* c: component) if(c->fex)
		phi["Fex("+c->molecule.name+")"] += c->fex->computeUniform(&N[c->offsetDensity], &Phi_N[c->offsetDensity]);

	//--------- Mixing functionals --------------
	for(const Fmix* fmix: fmixArr)
		phi["Fmix("+fmix->getName()+")"] += fmix->computeUniform(N, Phi_N);

	//--------- Convert site density gradients to molecular density gradients ---------
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		Phi_Nmol[ic] = 0.0;
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
			Phi_Nmol[ic] += Phi_N[c.offsetDensity+i] * c.molecule.sites[i]->positions.size();;
	}

	//phi.print(gInfo.fpLog, true, "\t\t\t\t%15s = %12.4le\n"); //Uncomment when debugging to get contributions
	return phi;
}

double FluidMixture::compute_p(double Ntot) const
{	std::vector<double> Nmol(component.size()), Phi_Nmol(component.size());
	double Nguess=0.;
	for(const FluidComponent* c: component) Nguess += c->Nbulk;
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = (Ntot/Nguess)*component[ic]->Nbulk;
	//p = N T + Sum_ic Nic d(aEx)/dNic - aEx where aEx is Helmholtz energy density
	double p = Ntot*T - computeUniformEx(Nmol, Phi_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
		p += Nmol[ic]*Phi_Nmol[ic];
	return p;
}
