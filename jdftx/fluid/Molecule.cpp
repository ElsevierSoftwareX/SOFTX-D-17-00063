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
		polKernel.free();
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


//Fourier transform of cuspless exponential
inline double cusplessExpTilde(double G, double norm, double a)
{	double aG = a*G;
	double den = 1./(1.+aG*aG);
	return norm * den*den*den;
}

//Fourier transform of gaussian
inline double gaussTilde(double G, double norm, double sigma)
{	double sigmaG = sigma*G;
	return norm * exp(-0.5*sigmaG*sigmaG);
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


/*
void ConvCoupling::setRadialKernel(SiteProperties& s)
{	
	ifstream ifs(s.kernelFilename.c_str());
	if(!ifs.is_open()) die("Can't open radial electron density file '%s' for reading.\n", s.kernelFilename.c_str());
	logPrintf("\nReading radial electron density model for %s site from '%s':\n",s.siteName.c_str(),s.kernelFilename.c_str());
	 
	std::vector<double> rVec,nVec; //to store values of radius and density
	double deltaRMin = 10.0; //minimum distance between radial grid points.
	
	while (!ifs.eof())
	{
		double r,n;
		string line;
		getline(ifs, line);
		istringstream(line) >> r >> n;
		logPrintf("r: %.12lf n: %.12lf",r,n);
		rVec.push_back(r);
		nVec.push_back(n);	
		
		int Nr = rVec.size();
		double deltaR = rVec[Nr]-rVec[Nr-1];
		
		if (deltaR <= 0.0)
		{
			die("ERROR reading electron density model for %s site from %s:\n"
				"Radial gridpoints must be in ascending order\n", s.siteName.c_str(), s.kernelFilename.c_str());
		}
		else if (deltaR < deltaRMin)
			deltaRMin = deltaR;
		
		if( ifs.fail() || ifs.bad())
			die("ERROR reading electron density model for %s site from %s\n", s.siteName.c_str(), s.kernelFilename.c_str());
		
	}
	
	RadialFunctionR KernelR(rVec.size());
	KernelR.r = rVec;
	KernelR.f = nVec;
	KernelR.initWeights();
	
	double dG = 0.02; //The uniform G-spacing for the radial function
	int ngridG = 2 * ceil(2.0*M_PI/(deltaRMin*dG)); //set the number of G grid points
	RadialFunctionG KernelG;
	
	//transform to fourier space
	KernelR.transform(0, dG, ngridG, KernelG);
	
	//allocate the kernel if not already allocated
	if(!s.couplingElecKernel)
	{
		s.couplingElecKernel = new RealKernel(gInfo);	
		//Interpolate fourier space function onto gInfo
		radialFunctionG(KernelG, *s.couplingElecKernel);
	}
	else //assumes radial function is an addition to an exponential (or other analytical function) 
	{
		RealKernel tmp(gInfo);	
		//Interpolate fourier space function onto gInfo
		radialFunctionG(KernelG, tmp);
		for(int i=0; i<gInfo.nG; i++)
			s.couplingElecKernel->data[i] += tmp.data[i];
	}
	
	double Kernel0 = s.couplingElecKernel->data[0];
	if (fabs(Kernel0-s.couplingZnuc-s.chargeZ*s.chargeKernel->data[0])>1e-12)
		logPrintf("Warning: classical site charges and electron kernel charges not balanced in convolution coupling.\n");

	logPrintf("Final electron density model for %s site has net charge %.12lf and site charge %.12lf\n",
				s.siteName.c_str(), Kernel0, Kernel0-s.couplingZnuc);
				
	s.couplingElecKernel->set();	
}
*/

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
			Q += site->chargeKernel(0) * site->positions.size();
	if(fabs(Q) < 1e-12) return 0.; //simplify neutrality testing
	else return Q;
}

vector3<> Molecule::getDipole() const
{	vector3<> P;
	for(const auto& site: sites)
		if(site->chargeKernel)
			for(const vector3<>& r: site->positions)
				P += site->chargeKernel(0) * r;
	if(P.length() < 1e-12) return vector3<>(); //simplify polarity testing
	else return P;
}

double Molecule::getVhs() const
{	double Vhs = 0;
	for(const auto& site: sites)
		if(site->Rhs > 0.)
			Vhs += (4.*M_PI/3)*pow(site->Rhs,3) * site->positions.size();
	return Vhs;
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
