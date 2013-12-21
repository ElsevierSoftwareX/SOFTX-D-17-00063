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

#include <core/MPIUtil.h>
#include <core/Util.h>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

MPIUtil::MPIUtil(int argc, char** argv)
{
	#ifdef MPI_ENABLED
	int rc = MPI_Init(&argc, &argv);
	if(rc != MPI_SUCCESS) { printf("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc); }
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
	#else
	//No MPI:
	nProcs = 1;
	iProc = 0;
	#endif
}

MPIUtil::~MPIUtil()
{
	#ifdef MPI_ENABLED
	MPI_Finalize();
	#endif
}

void MPIUtil::exit(int errCode) const
{
	#ifdef MPI_ENABLED
	MPI_Abort(MPI_COMM_WORLD, errCode);
	#else
	exit(errCode);
	#endif
}

void MPIUtil::bcast(int* data, size_t nData, int root) const
{	
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Bcast(data, nData, MPI_INT, root, MPI_COMM_WORLD);
	#endif
}

void MPIUtil::bcast(double* data, size_t nData, int root) const
{
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Bcast(data, nData, MPI_DOUBLE, root, MPI_COMM_WORLD);
	#endif
}

void MPIUtil::bcast(complex* data, size_t nData, int root) const
{	if(nProcs>1) bcast((double*)data, 2*nData, root);
}

void MPIUtil::bcast(string& s, int root) const
{	if(nProcs>1)
	{
		#ifdef MPI_ENABLED
		//Synchronize length of string:
		int len = s.length();
		bcast(&len, 1, root);
		if(iProc!=root) s.resize(len);
		//Bcast content:
		MPI_Bcast(&s[0], len, MPI_CHAR, root, MPI_COMM_WORLD);
		#endif
	}
}

#ifdef MPI_ENABLED
static MPI_Op mpiOp(MPIUtil::ReduceOp op)
{	switch(op)
	{	case MPIUtil::ReduceMax: return MPI_MAX;
		case MPIUtil::ReduceMin: return MPI_MIN;
		case MPIUtil::ReduceSum: return MPI_SUM;
		case MPIUtil::ReduceProd: return MPI_PROD;
		default: assert(!"Unknown reduction operation");
	}
	return 0;
}
#endif

void MPIUtil::allReduce(int* data, size_t nData, MPIUtil::ReduceOp op) const
{
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Allreduce(MPI_IN_PLACE, data, nData, MPI_INT, mpiOp(op), MPI_COMM_WORLD);
	#endif
}

void MPIUtil::allReduce(double* data, size_t nData, MPIUtil::ReduceOp op) const
{
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Allreduce(MPI_IN_PLACE, data, nData, MPI_DOUBLE, mpiOp(op), MPI_COMM_WORLD);
	#endif
}

void MPIUtil::allReduce(complex* data, size_t nData, MPIUtil::ReduceOp op) const
{	if(nProcs>1) allReduce((double*)data, 2*nData, op);
}
