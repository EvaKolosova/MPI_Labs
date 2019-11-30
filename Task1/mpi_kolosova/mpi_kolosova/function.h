#include <iostream>
#include <cstdlib>
#include <stdio.h>

using namespace std;

int** create_new_matrix(int row, int column)
{
	int** matrix = new int*[row];
	for (int i = 0; i < row; i++)
		matrix[i] = new int[column];
	return matrix;
}

int* create_result(int row)
{
	int* result = new int[row];
	for (int i = 0; i < row; i++)
		result[i] = 0;
	return result;
}

void delete_matrix(int** matrix, int row)
{
	for (int i = 0; i < row; i++)
		delete[] matrix[i];
	delete[] matrix;
}

void fill_matrix(int** matrix, int row, int column, int coeff)
{
	for (int i = 0; i < row; i++) 
	{
		for (int j = 0; j < column; j++)
		{
			matrix[i][j] = rand() % coeff;
		}
	}
}

void fill_summs_null(int* summs, int row) 
{
	for (int i = 0; i < row; i++)
	{
		summs[i] = 0;
	}
}

void output_matrix(int** matrix, int row, int column)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void output_summs(int* summs, int row)
{
	for (int i = 0; i < row; i++)
	{
		cout << summs[i] << " ";
		cout << endl;
	}
}

int* summing(int* summs, int** matrix, int row, int column)
{
	for (int i = 0; i < row; i++)
		for (int j = 0; j < column; j++)
			summs[i] += matrix[i][j];
	return summs;
}