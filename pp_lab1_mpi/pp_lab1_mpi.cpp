// pp_lab1_mpi.cpp: главный файл проекта.

#include "stdafx.h"
#include <iostream>
#include <string.h>
#include <intrin.h>
#include <math.h>
#include <mpi.h>

using namespace::std;
#pragma intrinsic(__rdtsc)

#define _USE_MATH_DEFINES

int proc_n;
int proc_count;
int **iA;
int **iB;
int **iC;
float **fA;
float **fB;
float **fC;
double **dA;
double **dB;
double **dC;
int listDim[] = { 10, 100, 200, 400, 800, 1000, 2000 };
double posl_res[3][7];
double par_res[3][7];


int main(int argc, char* argv[])
{
	unsigned __int64 t1, t2;
	double mpi_t1;
	double mpi_t2;
	cout.setf(std::ios::fixed);
	cout.precision(8);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_n);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

	// sequence program
	if (proc_n == 0)
	{
		iA = new int*[2000];
		for (int i = 0; i < 2000; i++)iA[i] = new int[2000];
		iB = new int*[2000];
		for (int i = 0; i < 2000; i++)iB[i] = new int[2000];
		iC = new int*[2000];
		for (int i = 0; i < 2000; i++)iC[i] = new int[2000];
		fA = new float*[2000];
		for (int i = 0; i < 2000; i++)fA[i] = new float[2000];
		fB = new float*[2000];
		for (int i = 0; i < 2000; i++)fB[i] = new float[2000];
		fC = new float*[2000];
		for (int i = 0; i < 2000; i++)fC[i] = new float[2000];
		dA = new double*[2000];
		for (int i = 0; i < 2000; i++)dA[i] = new double[2000];
		dB = new double*[2000];
		for (int i = 0; i < 2000; i++)dB[i] = new double[2000];
		dC = new double*[2000];
		for (int i = 0; i < 2000; i++)dC[i] = new double[2000];
		// инициализаци¤ int матрицы
		for (int i = 0; i < 2000; i++)
			for (int j = 0; j < 2000; j++)
			{
				iA[i][j] = rand();
				iB[i][j] = rand();
				iC[i][j] = 0;
			}
		// ѕоследовательный int
		for (int dim = 0; dim < 7; dim++)
		{
			t1 = __rdtsc();
			for (int i = 0; i < listDim[dim]; i++)
				for (int j = 0; j < listDim[dim]; j++)
				{
					iC[i][j] = 0;
					for (int k = 0; k < listDim[dim]; k++)
						iC[i][j] = iC[i][j] + iA[i][k] * iB[k][j];
				}
			t2 = __rdtsc();
			posl_res[0][dim] = (double)(t2 - t1) / 4000000000;
		}
		// инициализаци¤ float матрицы
		for (int i = 0; i < 2000; i++)
			for (int j = 0; j < 2000; j++)
			{
				fA[i][j] = (float)rand() / RAND_MAX + rand();
				fB[i][j] = (float)rand() / RAND_MAX + rand();
				fC[i][j] = 0;
			}
		// ѕоследовательный float
		for (int dim = 0; dim < 7; dim++)
		{
			t1 = __rdtsc();
			for (int i = 0; i < listDim[dim]; i++)
				for (int j = 0; j < listDim[dim]; j++)
				{
					fC[i][j] = 0;
					for (int k = 0; k < listDim[dim]; k++)
						fC[i][j] = fC[i][j] + fA[i][k] * fB[k][j];
				}
			t2 = __rdtsc();
			posl_res[1][dim] = (double)(t2 - t1) / 4000000000;
		}
		// инициализаци¤ double матрицы
		for (int i = 0; i < 2000; i++)
			for (int j = 0; j < 2000; j++)
			{
				dA[i][j] = (double)rand() / RAND_MAX + rand();
				dB[i][j] = (double)rand() / RAND_MAX + rand();
				dC[i][j] = 0;
			}
		// ѕоследовательный double
		for (int dim = 0; dim < 7; dim++)
		{
			t1 = __rdtsc();
			for (int i = 0; i < listDim[dim]; i++)
				for (int j = 0; j < listDim[dim]; j++)
				{
					dC[i][j] = 0;
					for (int k = 0; k < listDim[dim]; k++)
						dC[i][j] = dC[i][j] + dA[i][k] * dB[k][j];
				}
			t2 = __rdtsc();
			posl_res[2][dim] = (double)(t2 - t1) / 4000000000;
		}

		cout << "Sequence program\n";
		cout << "            10             100             200             400             800             1000            2000\n";
		cout << "int     ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[0][i] << "      ";
		cout << "\nfloat   ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[1][i] << "      ";
		cout << "\ndouble  ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[2][i] << "      ";
	}

	// MPI program
	if (proc_n == 0)
	{
		//int
		for (int dim = 0; dim < 7; dim++)
		{
			mpi_t1 = MPI_Wtime();
			int *a = new int[listDim[dim]];
			int *b = new int[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i += proc_count - 1)
			{
				int counter;
				if (listDim[dim] - i >= proc_count - 1) counter = proc_count;
				else counter = listDim[dim] - i + 1;
				for (int c = 1; c < counter; c++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						a[k] = iA[i + c - 1][k];
					}
					MPI_Ssend(a, listDim[dim], MPI_INT, c, 0, MPI_COMM_WORLD);
				}
				for (int j = 0; j < listDim[dim]; j++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						b[k] = iB[k][j];
					}
					if (i == 0)for (int c = 1; c < counter; c++)
					{
						MPI_Ssend(b, listDim[dim], MPI_INT, c, 0, MPI_COMM_WORLD);
					}
				}
				for (int c = 1; c < counter; c++)
				{
					MPI_Recv(b, listDim[dim], MPI_INT, c, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (int k = 0; k < listDim[dim]; k++) iC[i][k] = b[k];
				}
			}
			mpi_t2 = MPI_Wtime();
			par_res[0][dim] = mpi_t2 - mpi_t1;
		}
		for (int i = 0; i < 2000; i++) delete[] iA[i];
		delete[] iA;
		for (int i = 0; i < 2000; i++) delete[] iB[i];
		delete[] iB;
		for (int i = 0; i < 2000; i++) delete[] iC[i];
		delete[] iC;

		//float
		for (int dim = 0; dim < 7; dim++)
		{
			mpi_t1 = MPI_Wtime();
			float *a = new float[listDim[dim]];
			float *b = new float[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i += proc_count - 1)
			{
				int counter;
				if (listDim[dim] - i >= proc_count - 1) counter = proc_count;
				else counter = listDim[dim] - i + 1;
				for (int c = 1; c < counter; c++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						a[k] = fA[i + c - 1][k];
					}
					MPI_Ssend(a, listDim[dim], MPI_FLOAT, c, 0, MPI_COMM_WORLD);
				}
				for (int j = 0; j < listDim[dim]; j++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						b[k] = fB[k][j];
					}
					if (i == 0)for (int c = 1; c < counter; c++)
					{
						MPI_Ssend(b, listDim[dim], MPI_FLOAT, c, 0, MPI_COMM_WORLD);
					}
				}
				for (int c = 1; c < counter; c++)
				{
					MPI_Recv(b, listDim[dim], MPI_FLOAT, c, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (int k = 0; k < listDim[dim]; k++) fC[i][k] = b[k];
				}
			}
			mpi_t2 = MPI_Wtime();
			par_res[1][dim] = mpi_t2 - mpi_t1;
		}
		for (int i = 0; i < 2000; i++) delete[] fA[i];
		delete[] fA;
		for (int i = 0; i < 2000; i++) delete[] fB[i];
		delete[] fB;
		for (int i = 0; i < 2000; i++) delete[] fC[i];
		delete[] fC;

		//double
		for (int dim = 0; dim < 7; dim++)
		{
			mpi_t1 = MPI_Wtime();
			double *a = new double[listDim[dim]];
			double *b = new double[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i += proc_count - 1)
			{
				int counter;
				if (listDim[dim] - i >= proc_count - 1) counter = proc_count;
				else counter = listDim[dim] - i + 1;
				for (int c = 1; c < counter; c++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						a[k] = dA[i + c - 1][k];
					}
					MPI_Ssend(a, listDim[dim], MPI_DOUBLE, c, 0, MPI_COMM_WORLD);
				}
				for (int j = 0; j < listDim[dim]; j++)
				{
					for (int k = 0; k < listDim[dim]; k++)
					{
						b[k] = dB[k][j];
					}
					if (i == 0)for (int c = 1; c < counter; c++)
					{
						MPI_Ssend(b, listDim[dim], MPI_DOUBLE, c, 0, MPI_COMM_WORLD);
					}
				}
				for (int c = 1; c < counter; c++)
				{
					MPI_Recv(b, listDim[dim], MPI_DOUBLE, c, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for (int k = 0; k < listDim[dim]; k++) dC[i][k] = b[k];
				}
			}
			mpi_t2 = MPI_Wtime();
			par_res[2][dim] = mpi_t2 - mpi_t1;
		}
		for (int i = 0; i < 2000; i++) delete[] dA[i];
		delete[] dA;
		for (int i = 0; i < 2000; i++) delete[] dB[i];
		delete[] dB;
		for (int i = 0; i < 2000; i++) delete[] dC[i];
		delete[] dC;
	}

	if (proc_n != 0)
	{
		//int
		for (int dim = 0; dim < 7; dim++)
		{
			int *a = new int[listDim[dim]];
			int **b;
			b = new int*[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i++)b[i] = new int[listDim[dim]];
			int *c = new int[listDim[dim]];
			for (int i = proc_n - 1; i < listDim[dim]; i += proc_count - 1)
			{
				MPI_Recv(a, listDim[dim], MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (i == proc_n - 1)for (int j = 0; j < listDim[dim]; j++)
				{
					MPI_Recv(b[j], listDim[dim], MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				for (int k = 0; k < listDim[dim]; k++)
				{
					for (int j = 0; j < listDim[dim]; j++)
					{
						c[k] += a[j] * b[k][j];
					}
				}
				MPI_Ssend(c, listDim[dim], MPI_INT, 0, proc_n, MPI_COMM_WORLD);
			}
		}

		//float
		for (int dim = 0; dim < 7; dim++)
		{
			float *a = new float[listDim[dim]];
			float **b;
			b = new float*[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i++)b[i] = new float[listDim[dim]];
			float *c = new float[listDim[dim]];
			for (int i = proc_n - 1; i < listDim[dim]; i += proc_count - 1)
			{
				MPI_Recv(a, listDim[dim], MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (i == proc_n - 1)for (int j = 0; j < listDim[dim]; j++)
				{
					MPI_Recv(b[j], listDim[dim], MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				for (int k = 0; k < listDim[dim]; k++)
				{
					for (int j = 0; j < listDim[dim]; j++)
					{
						c[k] += a[j] * b[k][j];
					}
				}
				MPI_Ssend(c, listDim[dim], MPI_FLOAT, 0, proc_n, MPI_COMM_WORLD);
			}
		}

		//double
		for (int dim = 0; dim < 7; dim++)
		{
			double *a = new double[listDim[dim]];
			double **b;
			b = new double*[listDim[dim]];
			for (int i = 0; i < listDim[dim]; i++)b[i] = new double[listDim[dim]];
			double *c = new double[listDim[dim]];
			for (int i = proc_n - 1; i < listDim[dim]; i += proc_count - 1)
			{
				MPI_Recv(a, listDim[dim], MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (i == proc_n - 1)for (int j = 0; j < listDim[dim]; j++)
				{
					MPI_Recv(b[j], listDim[dim], MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				for (int k = 0; k < listDim[dim]; k++)
				{
					for (int j = 0; j < listDim[dim]; j++)
					{
						c[k] += a[j] * b[k][j];
					}
				}
				MPI_Ssend(c, listDim[dim], MPI_DOUBLE, 0, proc_n, MPI_COMM_WORLD);
			}
		}
	}

	if (proc_n == 0)
	{
		cout << "\n\nParallel program \n";
		cout << "                  10             100             200             400             800             1000            2000\n";
		cout << "int           ";
		for (int i = 0; i < 7; i++)
			cout << par_res[0][i] << "      ";
		cout << "\nAcceleration  ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[0][i] / par_res[0][i] << "      ";
		cout << "\nfloat         ";
		for (int i = 0; i < 7; i++)
			cout << par_res[1][i] << "      ";
		cout << "\nAcceleration  ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[1][i] / par_res[1][i] << "      ";
		cout << "\ndouble        ";
		for (int i = 0; i < 7; i++)
			cout << par_res[2][i] << "      ";
		cout << "\nAcceleration  ";
		for (int i = 0; i < 7; i++)
			cout << posl_res[2][i] / par_res[2][i] << "      ";
		cout << "\n\n";
		system("PAUSE");
	}
	MPI_Finalize();
	return 0;
}
