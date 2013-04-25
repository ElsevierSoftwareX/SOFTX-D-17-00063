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

#ifndef JDFTX_FLUID_MOLECULE_H
#define JDFTX_FLUID_MOLECULE_H

#include <electronic/RadialFunction.h>
#include <core/vector3.h>
#include <core/string.h>
#include <core/Data.h>
#include <vector>
#include <map>

//! Properties of a site in a multi-site molecule model
struct SiteProperties
{
	string name; //!< site name
	double Rhs; //!< hard sphere radius
	int atomicNumber; //!< necessary for vdW parameters
	double Znuc, sigmaNuc; //!< magnitude of the nuclear charge (positive) and corresponding gaussian width
	double Zelec, aElec; string elecFilename; //!< magnitude of electron charge (positive) and corresponding cuspless-exponential width (or override with file)
	double alpha, aPol; //!< isotropic polarizability and corresponding cuspless-exponential width

	SiteProperties();
	~SiteProperties();
	void setup(const GridInfo& gInfo); //!< initialize the radial functions from the properties specified above

	RadialFunctionG w0, w1, w2, w3, w1v, w2m; //!< Hard sphere weight functions
	RadialFunctionG elecKernel, chargeKernel; //!< Electron density and net charge density kernels for the sites
};


//! Group of sites that are identical under molecule symmetries
struct SiteGroup
{	std::shared_ptr<SiteProperties> siteProp; //!< Site properties: can be shared between multiple symmetry classes within one molecule and same species in different molecules
	std::vector< vector3<> > pos; //!< Positions w.r.t molecular origin in the reference orientation
};

//! Molecule: a collection of sites
struct Molecule
{	string name; //!< An identifier for molecule (used for EnergyComponent labels)
	std::vector<SiteGroup> siteGroup; //!< list of symmetry-equivalent site groups
	
	double getCharge() const; //!< total charge on molecule
	vector3<> getDipole() const; //!< total dipole moment on molecule

	//! Get the harmonic sum of radii for spheres in contact, with the multiplicities for each such pair
	std::map<double,int> getBonds() const;
};

#endif // JDFTX_FLUID_MOLECULE_H
