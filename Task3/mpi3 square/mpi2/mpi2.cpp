// mpi3.cpp : main project file.

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <thread>
#include <windows.h>
using namespace std;

int ProcNum, Rank, Grid, GridComm, ColComm, RowComm;
int GridC[2];

// инициализация матриц
void InitializationMatrix(double* MatrixA, double* MatrixB, int Size)
{
	for (int i = 0; i < Size; i++)
	{
		for (int j = 0; j < Size; j++)
		{
			MatrixA[i*Size + j] = rand() % 9;
			MatrixB[i*Size + j] = rand() % 9;
		}
	}
}

// Для блочного разделения матрицы между процессами процессорной решетки 
void Scatter(double* Matrix, double* MatrixBlock, int Size, int BSize)
{
	double* MatrixRow = new double[BSize*Size];//буфер для временного хранения горизонтальной полосы матрицы 

	//функцию Scatter вызываем только в тех процессах которые расположены в нулевом столбце решетки 
	if (GridC[1] == 0)  MPI_Scatter(Matrix, BSize*Size, MPI_DOUBLE, MatrixRow, BSize*Size, MPI_DOUBLE, 0, ColComm);

	for (int i = 0; i < BSize; i++)
		MPI_Scatter(&MatrixRow[i*Size], BSize, MPI_DOUBLE, &(MatrixBlock[i*BSize]), BSize, MPI_DOUBLE, 0, RowComm);

	delete[] MatrixRow;
}

//поблочно разделяем матрицу А(блоки сохр-ся в переменной Block) и матрицу В(блоки сохр-ся в переменной BlockB) 
void Distribution(double* MatrixA, double* MatrixB, double* Block, double* BlockB, int Size, int BSize)
{
	Scatter(MatrixA, Block, Size, BSize);
	Scatter(MatrixB, BlockB, Size, BSize);
}

// Создание коммуникатора в виде двумерной квадратной решетки 
// и коммуникаторов для каждой строки и каждого столбца решетки
void CreateGird()
{
	int DimSize[2];// Количество процессов в каждом измерении решетки
	int Periodic[2];// = 1 для каждого измерения, являющегося периодическим
	int Subdims[2];// = 1 для каждого измерения, оставляемого в подрешетке

	DimSize[0] = Grid;
	DimSize[1] = Grid;
	Periodic[0] = 0;
	Periodic[1] = 0;

	// Создание коммуникатора в виде квадратной решетки 
	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);

	// Определение координат процесса в решетке 
	MPI_Cart_coords(GridComm, Rank, 2, GridC);

	// Создание коммуникаторов для строк процессной решетки
	Subdims[0] = 0;// Фиксация измерения
	Subdims[1] = 1;// Наличие данного измерения в подрешетке

	//разделение решетки на подрешетки меньшей размерности
	MPI_Cart_sub(GridComm, Subdims, &RowComm);

	// Создание коммуникаторов для столбцов процессной решетки
	Subdims[0] = 1;
	Subdims[1] = 0;

	MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

//Функция выполняет рассылку блоков матрицы A по строкам процессорной решетки. Для этого в каждой строке решетки определяется ведущий процесс Pivot,
//осуществляющий рассылку. Для рассылки используется блок BlockA, переданный в процесс в момент начального распределения данных.
//Выполнение операции рассылки блоков осуществляется при помощи функции MPI_Bcast. Данная операция является 
//коллективной и ее локализация пределами отдельных строк решетки обеспечивается за счет использования коммуникаторов RowComm, 
//определенных для набора процессов каждой строки решетки в отдельности.
// Рассылка блоков матрицы А по строкам решетки процессов 
void BlockAComm(int i, double *BlockA, double* Block, int BSize)
{  // Определение ведущего процесса в строке процессной решетки		
	int Pivot = (GridC[0] + i) % Grid;

	if (GridC[1] == Pivot)
	{  // Копирование передаваемого блока в отдельный буфер памяти
		for (int j = 0; j < BSize*BSize; j++) 
			BlockA[j] = Block[j];
	}
	// Рассылка блока на все процессы (широковещательная рассылка) 
	MPI_Bcast(BlockA, BSize*BSize, MPI_DOUBLE, Pivot, RowComm);
}

//Функция выполняет циклический сдвиг блоков матрицы B по столбцам процессорной решетки.
//Каждый процесс передает свой блок следующему процессу Next в столбце процессов и получает блок, 
//переданный из предыдущего процесса Prev в столбце решетки.Выполнение операций передачи данных осуществляется 
//при помощи функции MPI_SendRecv_replace, которая обеспечивает все необходимые пересылки блоков, используя при этом один и 
//тот же буфер памяти pBblock.Кроме того, эта функция гарантирует отсутствие возможных тупиков, когда операции передачи данных 
//начинают одновременно выполняться несколькими процессами при кольцевой топологии сети.
// Циклический сдвиг блоков матрицы В вдоль столбца процессной решетки
void BlockBComm(double *BlockB, int BSize)
{
	MPI_Status Status;
	int Next = GridC[0] + 1;
	if (GridC[0] == Grid - 1) Next = 0;
	int Prev = GridC[0] - 1;
	if (GridC[0] == 0) Prev = Grid - 1;

	MPI_Sendrecv_replace(BlockB, BSize*BSize, MPI_DOUBLE, Next, 0, Prev, 0, ColComm, &Status);
}

// последовательное перемножение матриц
void Mult(double* MatrixA, double* MatrixB, double* MatrixR, int Size)
{
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
			for (int k = 0; k < Size; k++)
				MatrixR[i*Size + j] += MatrixA[i*Size + k] * MatrixB[k*Size + j];
}

// метод Фокса - перемножение блоков А и В
void Method(double* BlockA, double* Block, double* BlockB, double* BlockR, int BSize)
{
	for (int i = 0; i < Grid; i++)
	{
		//Рассылка блоков матрицы А по всем процессам своей строки решетки 
		BlockAComm(i, BlockA, Block, BSize);
		//Умножение матричных блоков на каждом процессе 
		Mult(BlockA, BlockB, BlockR, BSize);
		//Циклический сдвиг блоков матрицы В вдоль столбцов процессорной решетки 
		BlockBComm(BlockB, BSize);
	}
}

// соединение данных
void CollectionResult(double *Result, double* BlockR, int Size, int BSize)
{
	double *Row = new double[Size*BSize];

	//сбор данных из всех процессов в один корневой
	for (int i = 0; i < BSize; i++) MPI_Gather(&BlockR[i*BSize], BSize, MPI_DOUBLE, &Row[i*Size], BSize, MPI_DOUBLE, 0, RowComm);

	if (GridC[1] == 0) MPI_Gather(Row, BSize*Size, MPI_DOUBLE, Result, BSize*Size, MPI_DOUBLE, 0, ColComm);
	delete[]Row;
}

void PrintMatrix(double *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << matrix[i * cols + j] << " ";
		}
		cout << endl;
	}
}

int main(int argc, char *argv[])
{
	double *MatrixA = NULL, *BlockA = NULL; // матрица А
	double *MatrixB = NULL, *BlockB = NULL; // матрица В
	double *Result = NULL, *BlockR = NULL; // результат (матрица С)
	double *Block = NULL;
	int Size, BSize; //размер матричного блока
	double t1, t2;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	Grid = sqrt((double)ProcNum);

	if (ProcNum != Grid*Grid) { if (Rank == 0) cout << "Wrong count of Processes! It should be a perfect square." << endl; }
	else
	{
		if (Rank == 0) cout << "Matrix multiplication: " << endl;

		CreateGird();

		Size = stoi(string(argv[1]));

		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		BSize = Size / Grid;

		BlockA = new double[BSize*BSize];
		BlockB = new double[BSize*BSize];
		BlockR = new double[BSize*BSize];
		Block = new double[BSize*BSize];

		for (int i = 0; i < BSize*BSize; i++) BlockR[i] = 0;

		if (Rank == 0)
		{
			MatrixA = new double[Size*Size];
			MatrixB = new double[Size*Size];
			Result = new double[Size*Size];
			InitializationMatrix(MatrixA, MatrixB, Size);
			if (Size < 11)
			{
				cout << "Matrix A" << endl;
				PrintMatrix(MatrixA, Size, Size);
				cout << endl << "Matrix B" << endl;
				PrintMatrix(MatrixB, Size, Size);
				cout << endl;
			}
		}

		t1 = MPI_Wtime();

		Distribution(MatrixA, MatrixB, Block, BlockB, Size, BSize);

		Method(BlockA, Block, BlockB, BlockR, BSize);

		CollectionResult(Result, BlockR, Size, BSize);

		t2 = MPI_Wtime();

		if (Rank == 0)
		{
			if (Size < 11)
			{
				cout << endl <<  "Matrix C" << endl;
				PrintMatrix(Result, Size, Size);
				cout << endl;
			}
			cout << "Parallel method: " << t2 - t1 << endl;
		}

		delete[] MatrixA;
		delete[] MatrixB;
		delete[] Result;
		delete[] BlockA;
		delete[] BlockB;
		delete[] BlockR;
	}
	MPI_Finalize();
	//_getch();//Будет ждать пока пользователь не нажмёт любую клавишу
	return 0;
}