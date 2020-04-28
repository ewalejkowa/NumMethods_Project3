#include<iostream>
#include<math.h>
#define SIZE sizeof(s21) / sizeof(double)
using namespace std;
double x[20];

void gauss2(double v[20][21], int n) 
{
	//int n = v->size();
	for (int i = 0; i<n; i++)
	{
		//szukanie maximum dla danego wiersza; wybór elementu podstawowoego(wybór czêœciowy) ograniczony do k tej kolumny dla k tego kroku
		double maxEl = abs(v[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++)
		{
			if (abs(v[k][i]) > maxEl) 
			{
				maxEl = abs(v[k][i]);
				maxRow = k;
			}
		}
		//zamieñ maksymaln¹ wiersz z obecnym
		for (int k = i; k<n + 1; k++) 
		{
			double tmp = v[maxRow][k];
			v[maxRow][k] = v[i][k];
			v[i][k] = tmp;
		}
		// odjêcie wielokrotnoœci aktualnego wiersza od pozostalych wierszy 
		for (int k = i + 1; k<n; k++) 
		{
			double c = 0;
			if (v[k][i] != 0 && v[i][i] != 0) c = -v[k][i] / v[i][i];
			for (int j = i; j<n + 1; j++)
			{
				if (i == j)
				{
					v[k][j] = 0;
				}
				else
				{
					v[k][j] += c * v[i][j];
				}
			}
		}
	}
	//wynikowe równania Ax=b - zapisanie wyników w x
	
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = 0;
		if (v[i][i] != 0) x[i] = v[i][n] / v[i][i];
		for (int k = i - 1; k >= 0; k--)
		{
			v[k][n] -= v[k][i] * x[i];
		}
	}
}

int main(){
	double s21[20] = 
	{
		//indeksy 2 lub 7
		 0.1233957541, 0.1699148622, 0.2457313981, 0.3760970395, 0.5982299829, 0.8762715691, 0.9978801826,
		0.9883261673, 0.9870309473, 0.9993667568, 0.9925572200, 0.9849207469, 0.9999956457, 0.9915374521,
		0.1937142444, 0.0034410015, 0.0071060741, 0.0005756823, 0.0040972271, 0.0067222053 
	};
	double freq[20] = { 2.160913, 2.184642, 2.208656, 2.232956, 2.257543, 2.282417, 2.307579, 2.333029, 2.358767, 2.384794, 2.411110,
						2.437714, 2.464608, 2.491789, 2.519259, 2.547017, 2.575062, 2.603393, 2.632010, 2.660913		
	};
	
	///obliczanie wspó³czynników m
	double macierz_M[SIZE][SIZE+1];
	double macierz_D[SIZE];
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE+1; j++)
		{
			macierz_M[i][j] = 0;
		}
	}
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			if (i == j)
			{
				macierz_M[i][j] = 2;
				if (i == 0){		
					macierz_M[i][j + 1] = 1;  //warunek brzegowy 4b
				}
				else if (i == SIZE-1){
					macierz_M[i][j - 1] = 1;   //warynek brzegowy 4b
				}
				else
				{
					macierz_M[i][j - 1] = (freq[j] - freq[j - 1]) / ((freq[j] - freq[j - 1]) + (freq[j + 1] - freq[j]));
					macierz_M[i][j + 1] = (freq[j + 1] - freq[j]) / ((freq[j] - freq[j - 1]) + (freq[j + 1] - freq[j]));
				}
			}
		}
	}

	//wyliczenie delty
	for (int i = 0; i <SIZE; i++)
	{
		if (i == 0)
		{
			macierz_M[0][SIZE] = (6 / (freq[1] - freq[0]))* (((s21[1] - s21[0]) / (freq[1] - freq[0]))); //warunek brzegowy 4b
			macierz_M[0][0] += (-((freq[1] - freq[0]) / 3));
			macierz_M[0][1] += (-((freq[1] - freq[0]) / 6));
		}
		else if (i == SIZE-1)
		{
			macierz_M[i][i+1] = (6 / (freq[i] - freq[i-1])) *(-((s21[i] -s21[i-1]) / (freq[i] - freq[i-1]))); //warynek brzegowy 4b
			macierz_M[i][i] += (-((freq[i] - freq[i-1]) / 3));
			macierz_M[i][i-1] += (-((freq[i] - freq[i-1]) / 6));
		}
		else
		{
			macierz_M[i][SIZE] = (6 / ((freq[i] - freq[i - 1]) + (freq[i + 1] - freq[i]))*
							(((s21[i + 1] - s21[i]) / (freq[i + 1] - freq[i])) - ((s21[i] - s21[i - 1]) / (freq[i] - freq[i - 1]))));
		}
	}
	//rozwi¹zanie -wyzanczenie wyznaczników m
	 gauss2(macierz_M,20);

	double a[SIZE], b[SIZE], c[SIZE], d[SIZE];
	//wyznaczenie poszczególnych wspó³czynników równañ
	for (int i = 0; i < SIZE-1; i++)
	{
		a[i] = s21[i];
		b[i] = (s21[i + 1] - s21[i]) / (freq[i + 1] - freq[i]) - ((2 * x[i] + x[i + 1]) / 6 *(freq[i + 1] - freq[i]));
		c[i] = x[i] / 2;
		d[i] = (x[i + 1] - x[i]) / (6 * (freq[i + 1] - freq[i]));
	}
	cout << "\n";
	int k = 0;
	double tablica[400];
	//wypisanie obliczonych wartoœci na okreslonym zakresie
	for (int i = 0; i < SIZE-1; i++)
	{
		double j = freq[i];
		while( j<freq[i + 1])
		{
			double f = (a[i] + b[i] * (j - freq[i]) + c[i] * (j - freq[i])*(j - freq[i]) + d[i] * (j - freq[i])*(j - freq[i])*(j - freq[i]));
			if (f>0) printf("%.15f\n", 20 * log(f) / log(10));
			tablica[k] = j;
			j += 0.002;
			k++;
		}
	}
	printf("\n");
	for (int i = 0; i < k; i++)
	{
		printf("%.15f\n" ,tablica[i]);
	}
	printf("\n%d\n", k); //iloœæ wyznaczonych charakterystyk
	//system("pause");
	return 0;
}