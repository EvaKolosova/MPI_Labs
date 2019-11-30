// mpi_kolosova.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mpi.h" 
#include <stdio.h>
#include <iostream> 
#include <fstream>
#include "function.h"

using namespace std;

int main(int argc, char *argv[])
{
	int procNum, procRank;
	int row = 10, column = 10;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	if (procRank == 0)//������ � ������� ������
	{
		int* summs = new int[row];
		int** matrix = create_new_matrix(row, column);

		fill_matrix(matrix, row, column, 9);
		fill_summs_null(summs, row);

		output_matrix(matrix, row, column);

		double t1 = 0;

		if (procNum == 1) //��� ��������� ������� �������, ��� ��� �� ����
		{
			t1 = MPI_Wtime();
			//cout << "0: i am lonely\n";

			summs = summing(summs, matrix, row, column);
		}
		else //����� �������� ������ ������ ������� ��� ������������
		{
			t1 = MPI_Wtime();
			int numberOfProcess;
			for (int indexRow = 0; indexRow < row; indexRow++)//���� �������� ������
			{
				//MPI_Scatter(matrix, procNum, MPI_INT, , , , , );
				numberOfProcess = 1 + indexRow % (procNum - 1);//������� ������ �������� ����� ������ ������
				MPI_Send(matrix[indexRow], column, MPI_INT, numberOfProcess, 0, MPI_COMM_WORLD);
				cout << "0: sent " << indexRow << " row to proc# " << numberOfProcess << endl;
			}
			for (int indexRow = 0; indexRow < row; indexRow++)//���� ��� ��������� ������
			{
				MPI_Status status;

				numberOfProcess = 1 + indexRow % (procNum - 1);
				MPI_Recv(summs + indexRow, 1, MPI_INT, numberOfProcess, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//cout << "0: i received " << indexRow << " row from proc# " << numberOfProcess << endl;
			}
			cout << endl;
		}
		cout << "Sums:" << endl;//����� �������� ����� � ���� �������

		output_summs(summs, row);
		cout << endl;

		double t2 = MPI_Wtime();
		cout << "Time of work is: " << t2 - t1 << endl;

		delete_matrix(matrix, row);
		delete[] summs;

	}
	else//���� �� ������� �����, ������������ ������ ���� ��������� �������
	{
		MPI_Status status;

		int* �ountSum = new int[column];
		int countRowsForProcess = row / (procNum - 1);
		if (row % (procNum - 1) != 0)
		{//���� ���� "�����" � ����� "������"-������ ������, ��� ���� ��������, �� �������� ��� � ����� �����, 
			//������� ���� �������������� ����� ������
			if (procRank <= row % (procNum - 1))
			{
				countRowsForProcess++;
			}
		}
		//���� ������������ ������ �����(�)
		for (int numberCountRow = 0; numberCountRow < countRowsForProcess; numberCountRow++)
		{
			MPI_Recv(�ountSum, row, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//cout << "proc " << procRank << " : recv " << "COLUMN" << endl;
			int sum = 0;
			for (int i = 0; i < row; i++)
			{
				sum += �ountSum[i];
			}
			MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			//cout << "proc " << procRank << " : sent sum" << sum << endl;
		}
		delete[] �ountSum;
	}

	MPI_Finalize();
}

