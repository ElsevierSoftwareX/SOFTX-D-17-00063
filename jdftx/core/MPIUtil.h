/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_MPIUTIL_H
#define JDFTX_CORE_MPIUTIL_H

//! MPI wrapper class
class MPIUtil
{
public:
	int nProcesses, iProcess; //!< number of processes and rank
	bool isHead() { return iProcess==0; } //!< makes code more readable
	
	MPIUtil(int argc, char** argv);
	~MPIUtil();
};


#endif // JDFTX_CORE_MPIUTIL_H