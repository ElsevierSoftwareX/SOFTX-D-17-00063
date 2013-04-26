/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#include <electronic/FluidSolverParams.h>

FluidSolverParams::FluidSolverParams()
: verboseLog(false),
linearDielectric(false), linearScreening(false),
vdwScale(0.75)
{
}

void FluidSolverParams::setPCMparams()
{
	//Set PCM fit parameters:
	switch(pcmVariant)
	{	case PCM_SLSA13:
		{	nc = 1.2e-3;
			sigma = sqrt(0.5);
			cavityTension = 0.;
			assert(fluidType == FluidNonlocalPCM);
			initWarnings += "WARNING: Nonlocal PCM is highly experimental!\n";
			break;
		}
		case PCM_SGA13:
		{	nc = 7e-3;
			sigma = 0.6;
			cavityTension = 0.;
			initWarnings += "WARNING: PCM variant SGA13 is highly experimental!\n";
			break;
		}
		case PCM_GLSSA13:
		{	
			switch(solventName)
			{
				case CHCl3:
				{	nc = 2.4e-05;
					sigma = 0.6;
					cavityTension = -9.066e-6;
					break;
				}
				case CCl4:
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 1.15e-4;
							sigma = 0.6;
							cavityTension = -8.99e-06;
							break;
						case FluidNonlinearPCM:
							die("\nERROR: You can't use NonlinearPCM with CCl4 as it does not have a permanent dipole moment!\n");
						default: //Other fluids do not use these parameters
							break;
					}
					break;
				}
				default: // For water and unparametrized fluids
				{	switch(fluidType)
					{	case FluidLinearPCM:
							nc = 3.7e-4;
							sigma = 0.6;
							cavityTension = 5.4e-6;
							break;
						case FluidNonlinearPCM:
							nc = 1.0e-3;
							sigma = 0.6;
							cavityTension = 9.5e-6;
							break;
						default: //Other fluids do not use these parameters
							break;
					}
					if(solventName != H2O)
					{	initWarnings +=
							"WARNING: PCM variant GLSSA has not been parametrized for this solvent; using bulk\n"
							"   surface tension as effective cavity tension and water parameters for cavity size.\n";
						cavityTension = sigmaBulk;
					}
					break;
				}
			}
			break;
		}
		case PCM_LA12:
		case PCM_PRA05:
		{	nc = 7e-4;
			sigma = 0.6;
			cavityTension = 0.;
			if(fluidType == FluidNonlinearPCM)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has not been parametrized for NonlinearPCM; using LinearPCM fit parameters.\n";
			if( (fluidType==FluidLinearPCM || fluidType==FluidNonlinearPCM) && solventName != H2O)
				initWarnings += "WARNING: PCM variant LA12/PRA05 has been fit only for H2O; using nc and sigma from H2O fit.\n";
			break;
		}
	}
}

bool FluidSolverParams::needsVDW() const
{	switch(fluidType)
	{	case FluidNone:
			return false;
		case FluidLinearPCM:
		case FluidNonlinearPCM:
			return (pcmVariant == PCM_SGA13);
		case FluidNonlocalPCM:
		default: //All explicit fluid functionals
			return true;
	}
}

