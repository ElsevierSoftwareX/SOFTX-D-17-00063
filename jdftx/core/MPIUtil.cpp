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

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

MPIUtil::MPIUtil(int argc, char** argv)
{
	#ifdef MPI_ENABLED
	int rc = MPI_Init(&argc, &argv);
	if(rc != MPI_SUCCESS) { printf("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc); }
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &iProcess);
	#else
	//No MPI:
	nProcesses = 1;
	iProcess = 0;
	#endif
}

MPIUtil::~MPIUtil()
{
	#ifdef MPI_ENABLED
	MPI_Finalize();
	#endif
}
