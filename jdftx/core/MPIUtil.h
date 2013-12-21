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

#include <core/string.h>
#include <core/scalar.h>
#include <cstdlib>

//! MPI wrapper class
class MPIUtil
{
	int nProcs, iProc;
public:
	int iProcess() const { return iProc; } //!< rank of current process
	int nProcesses() const { return nProcs; }  //!< number of processes
	bool isHead() const { return iProc==0; } //!< whether this is the root process (makes code more readable)
	void exit(int errCode) const; //!< global exit (kill other MPI processes as well)
	
	//Broadcast functions:
	void bcast(int* data, size_t nData, int root=0) const;
	void bcast(double* data, size_t nData, int root=0) const;
	void bcast(complex* data, size_t nData, int root=0) const;
	void bcast(string& s, int root=0) const;
	
	//Reduce functions:
	enum ReduceOp { ReduceMin, ReduceMax, ReduceSum, ReduceProd };
	void allReduce(int* data, size_t nData, ReduceOp op) const;
	void allReduce(double* data, size_t nData, ReduceOp op) const;
	void allReduce(complex* data, size_t nData, ReduceOp op) const;
	
	MPIUtil(int argc, char** argv);
	~MPIUtil();
};


#endif // JDFTX_CORE_MPIUTIL_H