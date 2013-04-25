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


SiteProperties::SiteProperties() : Rhs(0), atomicNumber(0), Znuc(0), sigmaNuc(0), Zelec(0), aElec(0), alpha(0), aPol(0)
{	
}

SiteProperties::~SiteProperties()
{	
	elecKernel.free();
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

void SiteProperties::setup(const GridInfo& gInfo)
{
	//TODO: initialize charge and elec kernels
	if(Rhs)
	{	//TODO: Initialize hard sphere kernels
		//ErfFMTweight erfFMTweight(sphereRadius, sphereSigma);
		//applyFuncGsq(gInfo, erfFMTweight, w0->data, w1->data, w2->data, w3->data, w1v->data, w2m->data);
	}
}



double Molecule::getCharge() const
{	double Q = 0.0;
	for(const SiteGroup& group: siteGroup)
		if(group.siteProp->chargeKernel)
			Q += group.siteProp->chargeKernel.Gzero * group.pos.size();
	return Q;
}

vector3<> Molecule::getDipole() const
{	vector3<> P;
	for(const SiteGroup& group: siteGroup)
		if(group.siteProp->chargeKernel)
			for(const vector3<>& r: group.pos)
				P += group.siteProp->chargeKernel.Gzero * r;
	return P;
}

std::map<double,int> Molecule::getBonds() const
{	std::map<double,int> bond;
	for(const SiteGroup& group1: siteGroup)
	{	double R1 = group1.siteProp->Rhs;
		if(R1) for(vector3<> pos1: group1.pos)
		{	for(const SiteGroup& group2: siteGroup)
			double R2 = group2.siteProp->Rhs;
			if(R2) for(vector3<> pos2: group2.pos)
			{	if(fabs(R1+R2-(pos1-pos2).length()) < 1e-6*(R1+R2))
					bond[R1*R2/(R1+R2)]++;
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
