// MPI.cpp: определяет точку входа для консольного приложения.
//
//#include "pch.h" 
#include "stdafx.h"
#include "mpi.h" 
#include "function.h"
#include <iostream> 
#include <fstream>

using namespace std; 

int main(int argv, char *argc[]) 
{ 
	int procNum, procRank; 
	setlocale(LC_ALL,"Russian");
	FILE* file = fopen("SUM_Matr.txt","w");
	int n = 5; //размер матрицы 
	double* A =new double[n*n];
	double* b=new double[n]; 
	//int* X = NULL;
	int sum=0;
	int cel=0;
	int ost=0;

	filling_array(A, n, n, 9);
	filling_b(b, n, 9);
	out_display_matr(A, b, n, n);
	//out_file_matr(file, A, b, n, n);

	MPI_Init(&argv, &argc); 
	double t=MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank); 
	MPI_Comm_size(MPI_COMM_WORLD, &procNum); 

	int part = get_part_of_matr(n, procNum, procRank);		//вычисляем размер пакета с учетом b
	cout<<"part "<<part<<endl;
	int* numcol=new int[part];	//номера локальных столбцов
	double* A1= new double[(n+1)*part];	//расширенная матрица части
	double* X = new double[n];				//решение
	double* tmp = new double[n+1];

	for(int i=0; i<part; i++)
	{	
		numcol[i]=procRank+procNum*i;
		for(int j=0; j<n; j++)
		{
			A1[i*(n+1)+j]=A[j*n+numcol[i]];//i
			
			cout<<"A1:"<<i*(n+1)+j<<" "<<A1[i*(n+1)+j]<<"\t";
			cout<<"Rank "<<procRank<<endl;
		}
		cout<<endl;
		A1[i*(n+1)+n]=b[numcol[i]];//i
		cout<<"b"<<i*(n+1)+n<<" "<<A1[i*(n+1)+n]<<"\t"; cout<<"Rank "<<procRank<<endl;
		
		cout<<"numcol"<<i<<" "<<numcol[i]<<endl;
	}
	//прямой ход
	int row=0;
	for(int i=0; i<n-1; i++)
	{
		//исключаем хi
		if(i==numcol[row])
		{
			cout<<"numcol[row]"<<row<<" "<<numcol[row]<<endl;
			//рассылаем строку i, находящуюся в памяти текущего процесса
			MPI_Bcast(&A1[row*(n+1)], n+1,MPI_DOUBLE,procRank,MPI_COMM_WORLD);
			for(int j=0; j<=n; j++)
			{
				tmp[j]=A1[row*(n+1)+j];
				cout<<"tmp"<<j<<" "<<tmp[j]<<"\t";
			}
			cout<<endl;
			row++;
		}
		else
		{
			cout<<"will do Bcast"<<endl;
			MPI_Bcast(tmp, n+1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);
			cout<<"do Bcast"<<endl;
		}
		//Вычитаем принятую строку из уравнений, хранящихся в текущем процессе
		for(int j=row; j<part; j++)
		{
			double scaling=A1[j*(n+1)+i]/tmp[i];

			cout<<"A1:"<<j*(n+1)+i<<" "<<A1[j*(n+1)+i]<<endl;
			cout<<"tmp"<<i<<" "<<tmp[i]<<endl;
			cout<<"scal "<<scaling<<endl;
			
			for(int k=i; k<n+1; k++)
			{
				A1[j*(n+1)+k]-=scaling*tmp[k];
				
				cout<<"A1:"<<j*(n+1)+k<<" "<<A1[j*(n+1)+k]<<endl;
			}
		}
	}

	//Инициализация неизвестных
	row=0;
	for(int i=0; i<n; i++)
	{
		X[i]=0;
		if(i==numcol[row])
		{
			X[i]=A1[row*(n+1)+n];

			cout<<"X*"<<i<<" "<<X[i]<<endl;

			row++;

			cout<<"row "<<row<<endl;
		}
	}

	//Обратный ход
	row=part-1;
	for(int i=n-1; i>0; i--)
	{
		if(row>=0)
		{
			if(i==numcol[row])
			{
				X[i]/=A1[row*(n+1)+i]; //передаем найденное xi
				cout<<"X"<<i<<" "<<X[i]<<endl;
				MPI_Bcast(&X[i],1,MPI_DOUBLE, procRank, MPI_COMM_WORLD);
				row--;
			}
			else
				MPI_Bcast(&X[i], 1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);
		}
		else
			MPI_Bcast(&X[i], 1, MPI_DOUBLE, i%procNum, MPI_COMM_WORLD);

		for(int j=0; j<=row; j++)//корректировка локальных xi
		{
			cout<<"Xl^0*"<<numcol[j]<<" "<<X[numcol[j]]<<endl;
			X[numcol[j]]-=A1[j*(n+1)+i]*X[i];
			cout<<"A1:"<<j*(n+1)+i<<" "<<A1[j*(n+1)+i]<<endl;
			cout<<"X"<<i<<" "<<X[i]<<endl;
			cout<<"Xl^"<<numcol[j]<<" "<<X[numcol[j]]<<endl;
		}
	}
	if (procRank==0)
	{
		X[0]/=A1[row*(n+1)]; //корректировка x0
		cout<<"X0 "<<X[0]<<endl;
	}
	MPI_Bcast(X, 1,MPI_DOUBLE, 0, MPI_COMM_WORLD); //Каждый процесс содержит корректный вектор решений

	delete[] tmp, numcol,A1;
	t=MPI_Wtime()-t;

	if(procRank==0)
	{
		printf("Gaus: n %d, proc %d, time(sec) %6f\n",n, procNum, t);
		printf("X[%d] ", n);
		for(int i=0; i<n; i++)
			printf("%f ", X[i]);
		printf("\n");
	}
	delete[] X;
	MPI_Finalize();
	return 0;

/*int rank,size;
    int count = 3;
    double t1,t2;
    MPI_Status status;
    MPI_Init(&argv,&argc);
    
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 
    MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
    double temp;
    int maxCol;
    double maxValue;   
    int numb;
    int* ind = new int[count];
    for(int i =0; i < count; i++)
	{
        ind[i] = i;
		cout<<"ind"<<i<<" "<<ind[i]<<"\t";
	}
    int N = count/size; //пакет данных
	cout<<"N "<<N<<endl;
    double** A = new double*[count];
    for(int i = 0; i < count; i++)
        A[i] = new double[count];//N
    
    double*B = new double [count];
    for(int i = 0; i < count; i++) 
	{
		B[i] = rand()%10;
		cout<<"B"<<B[i]<<"\n";
	}
	
    //srand(time(NULL) + rank);   
    for(int i = 0; i < count; i++)    
	{
        for(int j = 0; j < count;j++)//N
		{
            A[i][j] = rand()%10;
			cout<<"A"<< i<<j<<" "<<A[i][j]<<"\t";
		}
		cout<<endl;
	}

    double* L = new double[count];
    double orB = B[0];
        //if(rank == 0) printf("%f",orB);

//-------------------------
    t1 = MPI_Wtime();
for(int p = 0; p < size; p++)
{
    if(p == rank)
    {       
        for(int k = 0;k<N;k++)
        {       
			maxValue = A[N*p+k][k];
            maxCol = N*p + k;
            //поиск главного элемента
            for(int i = N*p + k + 1; i < count; i ++)
                if(abs(A[i][k]) > abs(maxValue)) {maxValue = A[i][k]; maxCol = i;}
            for(int i = 0; i < N; i++)
            {
                temp = A[N*p+k][i];
                A[N*p+k][i] = A[maxCol][i];
                A[maxCol][i] = temp;
            }
                
            if(maxCol!=N*p+k)
			{
				temp = ind[N*p+k];
				ind[N*p+k] =ind[maxCol];
				ind[maxCol] = temp;
				temp = B[N*p+k];
				B[N*p+k] = B[maxCol];
				B[maxCol] = temp;
			}
                
            //
            for(int j = N*p + k; j < count - 1; j++)
			{
                L[j] = A[j+1][k]/A[p*N+k][k];
            }
            for(int proc = p + 1;proc < size; proc++)
            {            
                MPI_Send(L,count,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
                MPI_Send(&k,1,MPI_DOUBLE,proc,2,MPI_COMM_WORLD);
                MPI_Send(&maxCol,1,MPI_INT,proc,22,MPI_COMM_WORLD);
                MPI_Send(ind,count,MPI_INT,proc,33,MPI_COMM_WORLD);
            
            }
        for(int i = N*p + k + 1; i < count; i++)
        {
            for(int j = k; j < N; j++)
			{
                A[i][j] = A[i][j] - A[N*p + k][j] * L[i - 1];
				if(abs(A[i][j]) < 0.00001) A[i][j] = 0;
            }
            B[i] = B[i] - B[N*p + k]*L[i - 1];
        }
        
        if(p != size - 1){
            for(int proc = p + 1;proc < size; proc++)            
                MPI_Send(B,count,MPI_DOUBLE,proc,3,MPI_COMM_WORLD);
        }
        
            
		}
    }
    if(rank > p)
    {       for(int k = 0;k<N;k++) {     
                MPI_Recv(&numb,1,MPI_DOUBLE,p,2,MPI_COMM_WORLD,&status);
                MPI_Recv(L,count,MPI_DOUBLE,p,1,MPI_COMM_WORLD,&status);
                MPI_Recv(B,count,MPI_DOUBLE,p,3,MPI_COMM_WORLD,&status);
                MPI_Recv(ind,count,MPI_INT,p,33,MPI_COMM_WORLD,&status);
                MPI_Recv(&maxCol,1,MPI_INT,p,22,MPI_COMM_WORLD,&status);
                for(int i = 0; i < N; i++){
                    temp = A[maxCol][i];
                    A[maxCol][i] = A[N*p + numb][i];
                    A[N*p + numb][i] = temp;}
                    for(int i = numb + 1; i < count; i++)
                        {
                            for(int j = numb; j < N; j++)
                                A[i][j] = A[i][j] - A[k][j] * L[i - 1];//к потому что по диагонали
                            //B[i] = B[i] - B[k]*L[i - 1];
                        }
            }
    
    }
 
}
 
//-------так как последний свободный вектор, верный
MPI_Bcast(B,count,MPI_DOUBLE,size - 1,MPI_COMM_WORLD);
 
MPI_Bcast(ind,count,MPI_INT,size - 1, MPI_COMM_WORLD);
double* X = new double[count];
orB = B[0];
//-----обратный ход
for(int p = size - 1; p >=0; p--)
{
    if(p == rank) 
        {
 
            X[N*p + N - 1] = B[N*p + N - 1]/ A[N*p + N-1][N-1];
            
                for(int i = N - 2; i>=0; i--){
                    for(int j = i + 1; j < N; j++)
                        B[N*p + i] = B[N*p + i] - A[N*p + i][j] * X[N*p+j];
                    
                    X[N*p + i] = B[N*p + i]/A[N*p + i][i];
                }
                for(int i = 0; i < N*p; i++)
                    for(int j = 0; j < N;j++)
                        B[i] = B[i] - A[i][j] * X[N*p + j]; 
                if(p>0){
                MPI_Send(B,count,MPI_DOUBLE,p-1,4,MPI_COMM_WORLD);
                MPI_Send(X,count,MPI_DOUBLE,p-1,5,MPI_COMM_WORLD);}         
                
        }
    if(rank == p - 1)
    {       
            MPI_Recv(B,count,MPI_DOUBLE,p,4,MPI_COMM_WORLD,&status);
            MPI_Recv(X,count,MPI_DOUBLE,p,5,MPI_COMM_WORLD,&status);                
        
    }
}
t2 = MPI_Wtime();
printf("time = %f\n",t2 - t1);                  
double ch = 0;
for(int i=0; i<count; i++)
	cout<<"X"<<i<<" "<<X[i]<<endl;
MPI_Bcast(X,count,MPI_DOUBLE,0,MPI_COMM_WORLD);
 
for(int p = 0; p < size; p++){
    if(rank == p) {
        for(int i = 0; i < N; i++)
            ch = ch +  A[0][i] * X[N*p + i];
        if(p != size - 1) MPI_Send(&ch,1,MPI_DOUBLE,p + 1,223,MPI_COMM_WORLD);
        if(p == size - 1) printf("pogreshnost = %2.16f\n",orB - ch);
    }
    if(rank == p + 1) 
        MPI_Recv(&ch,1,MPI_DOUBLE,p,223,MPI_COMM_WORLD,&status);        
        
    }
MPI_Finalize();*/
}
