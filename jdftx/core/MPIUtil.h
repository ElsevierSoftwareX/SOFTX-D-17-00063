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
#include <cstdlib>
#include <cstdio>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

//! MPI wrapper class
class MPIUtil
{
	int nProcs, iProc;
public:
	int iProcess() const { return iProc; } //!< rank of current process
	int nProcesses() const { return nProcs; }  //!< number of processes
	bool isHead() const { return iProc==0; } //!< whether this is the root process (makes code more readable)

	MPIUtil(int argc, char** argv);
	~MPIUtil();
	void exit(int errCode) const; //!< global exit (kill other MPI processes as well)

	//Broadcast functions:
	void bcast(int* data, size_t nData, int root=0) const;
	void bcast(bool* data, size_t nData, int root=0) const;
	void bcast(double* data, size_t nData, int root=0) const;
	void bcast(string& s, int root=0) const;

	//Reduce functions (safe mode gaurantees identical results irrespective of round-off (but could be slower)):
	enum ReduceOp { ReduceMin, ReduceMax, ReduceSum, ReduceProd, ReduceLAnd, ReduceBAnd, ReduceLOr, ReduceBOr, ReduceLXor, ReduceBXor };
	void allReduce(int* data, size_t nData, ReduceOp op) const;
	void allReduce(bool* data, size_t nData, ReduceOp op) const;
	void allReduce(double* data, size_t nData, ReduceOp op, bool safeMode=false) const;

	//File access (tiny subset of MPI-IO, using byte offsets alone, and made to closely resemble stdio):
	#ifdef MPI_ENABLED
	typedef MPI_File File;
	#else
	typedef FILE* File;
	#endif
	void fopenRead(File& fp, const char* fname, size_t fsizeExpected=0, const char* fsizeErrMsg=0) const; //!< open file for reading and optionally check file size
	void fopenWrite(File& fp, const char* fname) const; //!< open file for writing
	void fclose(File& fp) const;
	void fseek(File fp, long offset, int whence) const; //!< syntax consistent with fseek from stdio
	void fread(void *ptr, size_t size, size_t nmemb, File fp) const;
	void fwrite(const void *ptr, size_t size, size_t nmemb, File fp) const;
};


#endif // JDFTX_CORE_MPIUTIL_H