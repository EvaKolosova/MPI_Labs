// mpi2.cpp: определяет точку входа для консольного приложения.
//
#include "stdafx.h"
#include "mpi.h" 
#include <iostream> 
#include <ctime>
#include <math.h>>
#include <Windows.h>
using namespace std;

int* result_creation(int column)
{
	int* res = new int[column];
	for (int i = 0; i<column; i++)
		res[i] = 0;
	return res;
}
void fill_matrix(double* matr, int row, int column, int coeff)
{
	for (int i = 0; i<row; i++)
	{
		for (int j = 0; j<column; j++)
		{
			matr[i * column + j] = rand() % coeff;
		}
	}
}
void fill_b(double* matr, int row, int coeff)
{
	for (int j = 0; j<row; j++)
		matr[j] = rand() % coeff;
}
void output_matrix(double* matr, double* b, int row, int column)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j <column; j++)
		{
			cout << matr[i * column + j] << " ";
		}
		cout << "|" << b[i] << " ";
		cout << endl;
	}
}
void output_res(double* res, int column)
{
	for (int i = 0; i<column; i++)
		cout << res[i] << endl;
	cout << endl;
}
int get_part_of_matrix(int size, int kolproc, int rank) //вычисляем части для процессов
{
	int part = size / kolproc;
	if (size%kolproc)
		part++;
	int othpart = part*kolproc - size;
	if (rank >= kolproc - othpart)
		part = part - 1;
	return part;
}

void Gauss(double** A, double* b, double* x, int n)
{
	int i, j, k;
	int size = n;
	double alfa;
	for (j = 0; j < size; j++) {
		for (i = j + 1; i < size; i++) {
			alfa = A[i][j] / A[j][j];
			for (k = j; k < size; k++) {
				A[i][k] -= alfa * A[j][k];
			}
			b[i] -= alfa * b[j];
		}
	}
	x[size - 1] = b[size - 1] / A[size - 1][size - 1];
	for (i = size - 2; i >= 0; i--) {
		double sum = 0;
		for (j = i + 1; j < size; j++) {
			sum += A[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / A[i][i];
	}
	if (n < 16)
	{
		cout << "Result of Gauss method:" << endl;
		output_res(x, n);
	}
}

int main(int argv, char *argc[])
{
	int procNum, procRank;
	setlocale(LC_ALL, "Russian");
	int n = 1000;
	double* A = new double[n*n];
	double* b = new double[n];
	double t1 = 0;//время

	fill_matrix(A, n, n, 9);
	fill_b(b, n, 9);

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	double tGauss1 = 0, tGauss2 = 0;
	double* bGauss = new double[n];
	for (int i = 0; i < n; i++)
	{
		bGauss[i] = b[i];
	}
	//if(procRank == 0)
	//{
	//	//---последовательно:	
	//	double** Agauss = new double *[n];
	//	for (int i = 0; i < n; i++)
	//		Agauss[i] = new double[n];
	//	cout << "Lin Gauss->\n";
	//	for (int i = 0; i < n; i++)
	//	{
	//		for (int j = 0; j < n; j++)
	//		{
	//			Agauss[i][j] = A[i*n + j];
	//			if (n < 16)
	//				cout << Agauss[i][j] << " ";
	//		}
	//		if (n < 16)
	//			cout << "\n";
	//	}
	//	cout << "\n";
	//	double* x = new double[n];

	//	tGauss1 = clock()/1000.0;
	//	//tGauss1 = MPI_Wtime();
	//	Gauss(Agauss, bGauss, x, n);
	//	tGauss2 = clock()/1000.0 - tGauss1;
	//	cout << "Lin Gauss - " << tGauss2 << "(sec)" << endl << endl;
	//}
	//---
	//*параллельно
	if (n < 16 && procRank == 0)
	{
		output_matrix(A, b, n, n);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int part = get_part_of_matrix(n, procNum, procRank);//вычисляем размер пакета с учетом b
	if (procRank == 0)
		cout << "part " << part << endl;
	int* NumberOfString = new int[part];//номера строк
	double* A1 = new double[(n + 1)*part];//расширенная матрица пакета
	double* X = new double[n];//решение
	double* tmp = new double[n + 1];

	for (int i = 0; i < part; i++)
	{
		NumberOfString[i] = procRank + procNum*i;
		for (int j = 0; j < n; j++)
		{
			A1[i*(n + 1) + j] = A[j + n*NumberOfString[i]];
			//cout << "A1:" << i*(n + 1) + j << " " << A1[i*(n + 1) + j] << "\t";
			//cout << "Rank " << procRank << endl;
		}
		//cout << endl;
		A1[i*(n + 1) + n] = b[NumberOfString[i]];
		//cout << "b" << i*(n + 1) + n << " " << A1[i*(n + 1) + n] << "\t"; cout << "Rank " << procRank << endl;
		//cout << "NumberOfString" << i << " " << NumberOfString[i] << endl;
	}

	//прямой ход
	int row = 0;
	for (int i = 0; i < n - 1; i++)
	{
		t1 = MPI_Wtime();
		//исключаем хi
		if (i == NumberOfString[row])
		{
			//cout << "NumberOfString[row]" << row << " " << NumberOfString[row] << endl;
			//рассылаем строку i, находящуюся в памяти текущего процесса
			MPI_Bcast(A1+row*(n + 1), n + 1, MPI_DOUBLE, procRank, MPI_COMM_WORLD);
			for (int j = 0; j <= n; j++)
			{
				tmp[j] = A1[row*(n + 1) + j];
				//cout << "tmp" << j << " " << tmp[j] << "\t";
			}
			cout << endl;
			row++;
		}
		else
		{
			//cout << "will Bcast" << endl;
			MPI_Bcast(tmp, n + 1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);
			//cout << "did Bcast" << endl;
		}
		//Вычитаем принятую строку из уравнений, хранящихся в текущем процессе
		for (int j = row; j < part; j++)
		{
			double scalar = A1[j*(n + 1) + i] / tmp[i];

			//cout << "A1:" << j*(n + 1) + i << " " << A1[j*(n + 1) + i] << endl;
			//cout << "tmp" << i << " " << tmp[i] << endl;
			//cout << "scal " << scalar << endl;

			for (int k = i; k < n + 1; k++)
			{
				A1[j*(n + 1) + k] -= scalar*tmp[k];
				//cout << "A1:" << j*(n + 1) + k << " " << A1[j*(n + 1) + k] << endl;
			}
		}
	}

	row = 0;
	for (int i = 0; i < n; i++)
	{
		X[i] = 0;
		if (i == NumberOfString[row])
		{
			X[i] = A1[row*(n + 1) + n];
			//cout << "X*" << i << " " << X[i] << endl;
			row++;
			//cout << "row " << row << endl;
		}
	}

	//Обратный ход
	row = part - 1;
	for (int i = n - 1; i > 0; i--)
	{
		if (row >= 0)
		{
			if (i == NumberOfString[row])
			{
				X[i] /= A1[row*(n + 1) + i]; //передаем найденное xi
				//cout << "X" << i << " " << X[i] << endl;
				MPI_Bcast(&X[i], 1, MPI_DOUBLE, procRank, MPI_COMM_WORLD);
				row--;
			}
			else
				MPI_Bcast(&X[i], 1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);
		}
		else
			MPI_Bcast(&X[i], 1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);

		for (int j = 0; j <= row; j++)//корректировка локальных xi
		{
			//cout << "Xl^0*" << NumberOfString[j] << " " << X[NumberOfString[j]] << endl;
			X[NumberOfString[j]] -= A1[j*(n + 1) + i] * X[i];
			//cout << "A1:" << j*(n + 1) + i << " " << A1[j*(n + 1) + i] << endl;
			//cout << "X" << i << " " << X[i] << endl;
			//cout << "Xl^" << NumberOfString[j] << " " << X[NumberOfString[j]] << endl;
		}
	}
	if (procRank == 0)
	{
		X[0] /= A1[row*(n + 1)]; //корректировка x0
		//cout << "X0 " << X[0] << endl;
	}
	MPI_Bcast(X, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Каждый процесс содержит корректный вектор решений

	delete[] tmp, NumberOfString, A1;
	double t2 = MPI_Wtime() - t1;

	if (procRank == 0)
	{
		cout << endl << "X[" << n << "]:" << endl;
		if (n < 16)
			output_res(X, n);
		else { cout << "Very big result - dont want to print it" << endl; }
		cout << "Parall Gauss: n " << n << ", proc " << procNum << ",  time " << t2 << " (sec)" << endl;
	}
	delete[] X;
	delete[] A;
	MPI_Finalize();
	return 0;
}