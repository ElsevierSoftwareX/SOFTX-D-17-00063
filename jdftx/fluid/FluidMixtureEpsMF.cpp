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

#include <fluid/FluidMixture.h>
#include <fluid/IdealGas.h>
#include <electronic/operators.h>
#include <core/DataMultiplet.h>

struct EpsMFsolver : public LinearSolvable<DataGptrVec>
{
	const GridInfo& gInfo;
	const std::vector<const FluidComponent*>& component; //!< fluid components
	const double Ceps;
	RadialFunctionG Kkernel; //!< preconditioner
	DataRptrCollection alphaN; //!< sum alpha_i N_i for each fluid component
	MinimizeParams mp;
	
	EpsMFsolver(const GridInfo& gInfo, const std::vector<const FluidComponent*>& component, const double Ceps)
	: gInfo(gInfo), component(component), Ceps(Ceps)
	{
		//Initialize preconditioner:
		unsigned nGradial = unsigned(ceil(gInfo.GmaxGrid/gInfo.dGradial))+5;
		std::vector<double> KkernelSamples(nGradial);
		for(unsigned i=0; i<nGradial; i++)
		{	double G = i*gInfo.dGradial;
			//Compute diagonal part of the hessian (K^-1 - (1-C)chi):
			double diagChi = 0.;
			for(const FluidComponent* c: component)
				diagChi += c->idealGas->get_Nbulk() * c->molecule.getAlphaTot() * pow(c->molecule.mfKernel(G), 2);
			//Compute diagonal part of hessian:
			double diagH = diagChi * (1. + 4*M_PI*(1.-Ceps)*diagChi);
			//Set its inverse as the preconditioner:
			KkernelSamples[i] = (diagH>1e-8) ? 1./(gInfo.dV*diagH) : 0.;
		}
		Kkernel.init(0, KkernelSamples, gInfo.dGradial);
		
		//Initialize mnimize parameters:
		mp.nDim = gInfo.nr;
		mp.knormThreshold = 1e-12;
		mp.fpLog = globalLog;
		mp.linePrefix = "\tPolMinimize: ";
	}
	
	~EpsMFsolver()
	{	Kkernel.free();
	}
	
	static DataGptr coulomb(const DataGptr& rho) { return (-4*M_PI) * Linv(O(rho)); }
	
	DataGptrVec chi(const DataGptrVec& eps) const
	{	DataGptrVec chiEps;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const RadialFunctionG& mfKernel = component[ic]->molecule.mfKernel;
			chiEps += mfKernel*J(alphaN[ic] * I(mfKernel*eps));
		}
		return chiEps;
	}
	
	DataGptrVec hessian(const DataGptrVec& eps) const
	{	DataGptrVec chiEps = chi(eps);
		return chi(eps - (4*M_PI*Ceps)*chiEps - gradient(coulomb(divergence(chiEps))));
	}
	
//     DataGptrVec precondition(const DataGptrVec& v) const
//     {	return Kkernel * v;
// 	}
	
	void solve(DataRptrCollection& indep, const DataGptr& rhoExternal)
	{	
		//Compute site densities and RHS of equation:
		alphaN.assign(component.size(), 0);
		DataGptrVec PMF; DataGptr rhoMF;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *(component[ic]);
			DataRptrCollection N(c.molecule.sites.size()); DataRptrVec P;
			c.idealGas->getDensities(&indep[c.offsetIndep], &N[0], P);
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.chargeKernel)
					rhoMF += s.chargeKernel(0) * (c.molecule.mfKernel * J(N[i]));
				if(s.polKernel)
				{	alphaN[ic] += N[i] * pow(s.polKernel(0),2);
					//Polarizability contributions to rhoMF and PMF:
					if(rhoExternal)
					{	DataGptrVec NPbarMF = s.chargeKernel(0) * (c.molecule.mfKernel * J(N[i]
							* I(gradient(s.polKernel*Linv((gInfo.detR*4*M_PI)*rhoExternal))) ) );
						rhoMF -= divergence(NPbarMF);
						PMF += NPbarMF;
					}
				}
			}
			if(P) PMF += c.molecule.mfKernel * J(P);
		}
		DataGptrVec RHS = (4*M_PI*Ceps) * PMF - gradient(coulomb(rhoMF));
		
		//Linear solve:
		unsigned nIndepIdgas = component.back()->offsetIndep + component.back()->idealGas->nIndep;
		nullToZero(indep, gInfo);
		DataGptrVec& epsMF = state;
		epsMF = J(DataRptrVec(&indep[nIndepIdgas]));
		LinearSolvable<DataGptrVec>::solve(RHS, mp);
		
		//Update indep:
		for(int k=0; k<3; k++) indep[nIndepIdgas+k] = I(epsMF[k]);
		
		//Cleanup:
		alphaN.clear();
	}
};

double FluidMixture::minimizeEpsMF()
{	if(!polarizable) return 0.;
	if(!epsMFsolver) epsMFsolver = std::make_shared<EpsMFsolver>(gInfo, component, Ceps);
	epsMFsolver->solve(state, rhoExternal);
	return 0.;
}
