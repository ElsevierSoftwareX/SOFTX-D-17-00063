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
#include <vector>
#include <algorithm>

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

//----------------------- Broadcast routines -------------------------------

void MPIUtil::bcast(int* data, size_t nData, int root) const
{	
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Bcast(data, nData, MPI_INT, root, MPI_COMM_WORLD);
	#endif
}

void MPIUtil::bcast(bool* data, size_t nData, int root) const
{	
	#ifdef MPI_ENABLED
	if(nProcs>1)
	{	std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		MPI_Bcast(&intCopy[0], nData, MPI_INT, root, MPI_COMM_WORLD);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
	#endif
}

void MPIUtil::bcast(double* data, size_t nData, int root) const
{
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Bcast(data, nData, MPI_DOUBLE, root, MPI_COMM_WORLD);
	#endif
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


//----------------------- Reduction routines -------------------------------

#ifdef MPI_ENABLED
static MPI_Op mpiOp(MPIUtil::ReduceOp op)
{	switch(op)
	{	case MPIUtil::ReduceMax: return MPI_MAX;
		case MPIUtil::ReduceMin: return MPI_MIN;
		case MPIUtil::ReduceSum: return MPI_SUM;
		case MPIUtil::ReduceProd: return MPI_PROD;
		case MPIUtil::ReduceLAnd: return MPI_LAND;
		case MPIUtil::ReduceBAnd: return MPI_BAND;
		case MPIUtil::ReduceLOr: return MPI_LOR;
		case MPIUtil::ReduceBOr: return MPI_BOR;
		case MPIUtil::ReduceLXor: return MPI_LXOR;
		case MPIUtil::ReduceBXor: return MPI_BXOR;
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

void MPIUtil::allReduce(bool* data, size_t nData, MPIUtil::ReduceOp op) const
{	
	#ifdef MPI_ENABLED
	if(nProcs>1)
	{	std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		MPI_Allreduce(MPI_IN_PLACE, &intCopy[0], nData, MPI_INT, mpiOp(op), MPI_COMM_WORLD);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
	#endif
}

void MPIUtil::allReduce(double* data, size_t nData, MPIUtil::ReduceOp op, bool safeMode) const
{
	#ifdef MPI_ENABLED
	if(nProcs>1)
	{	if(safeMode) //Reduce to root node and then broadcast result (to ensure identical values)
		{	MPI_Reduce(isHead()?MPI_IN_PLACE:data, data, nData, MPI_DOUBLE, mpiOp(op), 0, MPI_COMM_WORLD);
			bcast(data, nData, 0);
		}
		else //standard Allreduce
			MPI_Allreduce(MPI_IN_PLACE, data, nData, MPI_DOUBLE, mpiOp(op), MPI_COMM_WORLD);
	}
	#endif
}

//----------------------- File I/O routines -------------------------------

void MPIUtil::fopenRead(File& fp, const char* fname, size_t fsizeExpected, const char* fsizeErrMsg) const
{	if(fsizeExpected)
	{	off_t fsize = fileSize(fname);
		if(fsize != off_t(fsizeExpected))
			die("Length of '%s' was %ld instead of the expected %ld bytes.\n%s\n", fname, fsize, fsizeExpected, fsizeErrMsg ? fsizeErrMsg : "");
	}
	#ifdef MPI_ENABLED
	if(MPI_File_open(MPI_COMM_WORLD, (char*)fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "rb");
	if(!fp)
	#endif
		die("Error opening file '%s' for reading.\n", fname);
}

void MPIUtil::fopenWrite(File& fp, const char* fname) const
{
	#ifdef MPI_ENABLED
	if(MPI_File_open(MPI_COMM_WORLD, (char*)fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "wb");
	if(!fp)
	#endif
		 die("Error opening file '%s' for writing.\n", fname);
}

void MPIUtil::fclose(File& fp) const
{
	#ifdef MPI_ENABLED
	MPI_File_close(&fp);
	#else
	::fclose(fp);
	#endif
}

void MPIUtil::fseek(File fp, long offset, int whence) const
{
	#ifdef MPI_ENABLED
	int mpi_whence = 0;
	switch(whence)
	{	case SEEK_CUR: mpi_whence = MPI_SEEK_CUR; break;
		case SEEK_SET: mpi_whence = MPI_SEEK_SET; break;
		case SEEK_END: mpi_whence = MPI_SEEK_END; break;
		default: assert(!"Invalid seek offset mode");
	}
	if(MPI_File_seek(fp, offset, mpi_whence) != MPI_SUCCESS)
	#else
	if(::fseek(fp, offset, whence) != 0)
	#endif
		die("Error in file seek.\n");
}

void MPIUtil::fread(void *ptr, size_t size, size_t nmemb, File fp) const
{
	#ifdef MPI_ENABLED
	MPI_Status status;
	MPI_File_read(fp, ptr, size*nmemb, MPI_BYTE, &status);
	int count; MPI_Get_count(&status, MPI_BYTE, &count);
	if(size_t(count) != size*nmemb)
	#else
	if(::fread(ptr, size, nmemb, fp) != nmemb)
	#endif
		die("Error in file read.\n");
}

void MPIUtil::fwrite(const void *ptr, size_t size, size_t nmemb, File fp) const
{
	#ifdef MPI_ENABLED
	MPI_Status status;
	MPI_File_write(fp, (void*)ptr, size*nmemb, MPI_BYTE, &status);
	int count; MPI_Get_count(&status, MPI_BYTE, &count);
	if(size_t(count) != size*nmemb)
	#else
	if(::fwrite(ptr, size, nmemb, fp) != nmemb)
	#endif
		die("Error in file write.\n");
}
