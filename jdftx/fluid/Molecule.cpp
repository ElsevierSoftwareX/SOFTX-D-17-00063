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

#include <fluid/ErfFMTweight.h>
#include <fluid/Molecule.h>
#include <core/Operators.h>


Molecule::Site::Site(string name, int atomicNumber) : name(name), Rhs(0), atomicNumber(atomicNumber), Znuc(0), sigmaNuc(0), Zelec(0), aElec(0), alpha(0), aPol(0), initialized(false)
{	
}

Molecule::Site::~Site()
{	if(initialized)
	{	elecKernel.free();
		chargeKernel.free();
		if(Rhs)
		{	w0.free();
			w1.free();
			w2.free();
			w3.free();
			w1v.free();
			w2m.free();
		}
	}
}

void Molecule::Site::setup(const GridInfo& gInfo)
{
	//TODO: initialize charge and elec kernels
	if(Rhs)
	{	//TODO: Initialize hard sphere kernels
		//ErfFMTweight erfFMTweight(sphereRadius, sphereSigma);
		//applyFuncGsq(gInfo, erfFMTweight, w0->data, w1->data, w2->data, w3->data, w1v->data, w2m->data);
	}
	initialized = true;
}


void Molecule::setup(const GridInfo& gInfo)
{	for(auto& site: sites) if(!*site) site->setup(gInfo);
}


bool Molecule::isMonoatomic() const
{	return (sites.size()==1) && (sites[0]->positions.size()==1);
}

double Molecule::getCharge() const
{	double Q = 0.0;
	for(const auto& site: sites)
		if(site->chargeKernel)
			Q += site->chargeKernel.Gzero * site->positions.size();
	return Q;
}

vector3<> Molecule::getDipole() const
{	vector3<> P;
	for(const auto& site: sites)
		if(site->chargeKernel)
			for(const vector3<>& r: site->positions)
				P += site->chargeKernel.Gzero * r;
	return P;
}

std::map<double,int> Molecule::getBonds() const
{	std::map<double,int> bond;
	for(const auto& site1: sites)
	{	double R1 = site1->Rhs;
		if(R1) for(vector3<> pos1: site1->positions)
		{	for(const auto& site2: sites)
			{	double R2 = site2->Rhs;
				if(R2) for(vector3<> pos2: site2->positions)
				{	if(fabs(R1+R2-(pos1-pos2).length()) < 1e-6*(R1+R2))
						bond[R1*R2/(R1+R2)]++;
				}
			}
		}
	}
	//correct for double counting:
	for(auto& bondEntry: bond)
	{	assert(bondEntry.second % 2 == 0);
		bondEntry.second /= 2;
	}
	return bond;
}
