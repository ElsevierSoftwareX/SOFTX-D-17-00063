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

struct EpsMFsolver : public LinearSolvable<DataGptr>
{
	const GridInfo& gInfo;
	const std::vector<const FluidComponent*>& component; //!< fluid components
	const double Cpol;
	RadialFunctionG Kkernel; //!< preconditioner
	DataRptrCollection CalphaN; //!< Cpol sum alpha_i N_i for each fluid component
	MinimizeParams mp;
	
	EpsMFsolver(const GridInfo& gInfo, const std::vector<const FluidComponent*>& component, const double Cpol)
	: gInfo(gInfo), component(component), Cpol(Cpol)
	{
		//Initialize preconditioner:
		unsigned nGradial = unsigned(ceil(gInfo.GmaxGrid/gInfo.dGradial))+5;
		std::vector<double> KkernelSamples(nGradial);
		for(unsigned i=0; i<nGradial; i++)
		{	double G = i*gInfo.dGradial, G2=G*G;
			//Compute diagonal part of the hessian (K^-1 - chi):
			double diagH = G2/(4*M_PI);
			for(const FluidComponent* c: component)
				diagH += G2 * c->idealGas->get_Nbulk() * c->molecule.getAlphaTot() * pow(c->molecule.mfKernel(G), 2);
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
	
	DataGptr chi(const DataGptr& phi) const
	{	DataGptr chiPhi;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const RadialFunctionG& mfKernel = component[ic]->molecule.mfKernel;
			chiPhi += mfKernel*divergence(J(CalphaN[ic] * I(mfKernel*gradient(phi))));
		}
		return chiPhi;
	}
	
	DataGptr hessian(const DataGptr& phi) const
	{	return (-1./(gInfo.detR*4*M_PI))*L(phi) - chi(phi);
	}
	
	DataGptr precondition(const DataGptr& v) const
	{	return Kkernel * v;
	}
	
	void solve(DataRptrCollection& indep, const DataGptr& rhoExternal, vector3<> Eexternal)
	{	static StopWatch watch("EpsMFsolve");
	
		//Compute site densities and RHS of equation:
		CalphaN.assign(component.size(), 0);
		DataGptr RHS; vector3<> Prot0; double epsInf=1.;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *(component[ic]);
			DataRptrCollection N(c.molecule.sites.size()); DataRptrVec P;
			c.idealGas->getDensities(&indep[c.offsetIndep], &N[0], P);
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.chargeKernel)
					RHS += s.chargeKernel(0) * (c.molecule.mfKernel * J(N[i]));
				if(s.polKernel)
				{	CalphaN[ic] += Cpol * N[i] * pow(s.polKernel(0),2);
					//Polarizability contributions to rhoMF and PMF:
					if(rhoExternal)
					{	RHS += Cpol * divergence( s.polKernel(0) * (c.molecule.mfKernel * J(N[i]
							* I(gradient(s.polKernel*coulomb(rhoExternal))) ) ) );
					}
				}
				epsInf += 4*M_PI * (integral(N[i])/gInfo.detR) * Cpol * pow(s.polKernel(0),2);
			}
			if(P) for(int k=0; k<3; k++) Prot0[k] += integral(P[k])/gInfo.detR;
		}
		watch.start();
		
		//Linear solve:
		unsigned nIndepIdgas = component.back()->offsetIndep + component.back()->idealGas->nIndep;
		nullToZero(indep, gInfo);
		state = -Linv(O(divergence(J(DataRptrVec(&indep[nIndepIdgas])))));
		LinearSolvable<DataGptr>::solve(RHS, mp);
		
		//Update indep:
		DataGptrVec epsMF = -gradient(state);
		vector3<> epsMF0 = (Eexternal - 4*M_PI*Prot0)/epsInf;
		for(int k=0; k<3; k++) indep[nIndepIdgas+k] = I(epsMF[k]) + epsMF0[k];
		
		//Cleanup:
		CalphaN.clear();
		watch.stop();
	}
};

double FluidMixture::minimizeEpsMF()
{	if(!polarizable) return 0.;
	if(!epsMFsolver) epsMFsolver = std::make_shared<EpsMFsolver>(gInfo, component, Cpol);
	epsMFsolver->solve(state, rhoExternal, Eexternal);
	return 0.;
}
