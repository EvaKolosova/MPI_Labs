// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO: reference additional headers your program requires here

//не совсем тот Гаусс
//double* gauss(double **a, double *y, int n)
//{
//	double *x, max;
//	int k, index;
//	const double eps = 0.00001;  // точность
//	x = new double[n];
//	k = 0;
//	while (k < n)
//	{
//		// Поиск строки с максимальным a[i][k]
//		max = abs(a[k][k]);
//		index = k;
//		for (int i = k + 1; i < n; i++)
//		{
//			if (abs(a[i][k]) > max)
//			{
//				max = abs(a[i][k]);
//				index = i;
//			}
//		}
//		// Перестановка строк
//		if (max < eps)
//		{
//			// нет ненулевых диагональных элементов
//			cout << "Null column ";
//			cout << index << " matr A" << endl;
//			return 0;
//		}
//		for (int j = 0; j < n; j++)
//		{
//			double temp = a[k][j];
//			a[k][j] = a[index][j];
//			a[index][j] = temp;
//		}
//		double temp = y[k];
//		y[k] = y[index];
//		y[index] = temp;
//		// Нормализация уравнений
//		for (int i = k; i < n; i++)
//		{
//			double temp = a[i][k];
//			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
//			for (int j = 0; j < n; j++)
//				a[i][j] = a[i][j] / temp;
//			y[i] = y[i] / temp;
//			if (i == k)  continue; // уравнение не вычитать само из себя
//			for (int j = 0; j < n; j++)
//				a[i][j] = a[i][j] - a[k][j];
//			y[i] = y[i] - y[k];
//		}
//		k++;
//	}
//	// обратная подстановка
//	for (k = n - 1; k >= 0; k--)
//	{
//		x[k] = y[k];
//		for (int i = 0; i < k; i++)
//			y[i] = y[i] - a[i][k] * x[k];
//	}
//	return x;
//}