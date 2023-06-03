#include<vector>
#include<array>
#include<iostream>
#include<iomanip>


using namespace std;


double solution(double x)
{
	return 3 / (2 * exp(3*x) + 1);
}

double f(double x, double y)
{
	return -3 * y + y * y;
}

// x | y | abs(error)
using Table = vector<array<double, 3>>;

void initTaylorCoeffs(vector<double>& coeffs)
{
	coeffs.resize(12, 0.0);
	coeffs[0] = 1.0;
	coeffs[1] = -2.0;
	coeffs[2] = 1.0;
	coeffs[3] = 1.0;
	coeffs[4] = -5.0 / 4.0;
	coeffs[5] = -7.0 / 20.0;
	coeffs[6] = 49.0 / 40.0;
	coeffs[7] = -53.0 / 280.0;
	coeffs[8] = -443.0 / 448.0;
	coeffs[9] = 1259.0 / 2240.0;
	coeffs[10] = 14201.0 / 22400.0;
	coeffs[11] = -26171.0 / 35200.0;
}

double getValueFromTaylor(double x, const vector<double>& memb)
{
	const int n = static_cast<int>(memb.size()) - 1;
	double result = memb[n];
	for (int i = n - 2; i >= 0; --i)
	{
		result *= x;
		result += memb[i];
	}

	return result;
}


void initPoints(vector<double>& pointVector, double x0, double h, int M, int N)
{
	pointVector.resize(abs(M - N) + 1);
	for (int i = M; i <= N; ++i)
	{
		pointVector[i - M] = x0 + h * i;
	}
};

template<typename Iterator>
void initTable(Table& table, double x, Iterator begin, Iterator end)
{
	int size = end - begin;
	table.resize(size);

	Iterator it = begin;
	for (int i = 0; i < size; ++i)
	{
		table[i][0] = *it;
		if (abs(table[i][0] - x) < 1e-12)
		{
			table[i][1] = solution(x);
		}
		++it;
	}
}

void printTable(const Table& table, bool withErrorRow = true)
{
	int size = table.size();

	for (int i = 0; i < size; ++i)
	{
		cout << setw(14) << table[i][0] << " | " << setw(14) << table[i][1];
		if (withErrorRow)
		{
			cout << " | " << setw(14) << table[i][2];
		}
		cout << "\n";
	}
}

void fillTable(Table& table)
{
	int size = table.size();

	for (size_t i = 0; i < size; ++i)
	{
		table[i][1] = solution(table[i][0]);
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}

}

void TaylorMethod(Table& table, const vector<double>& memb)
{
	const size_t size = table.size();

	for (size_t i = 0; i < size; ++i)
	{
		table[i][1] = getValueFromTaylor(table[i][0], memb);
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}

}

void EulerMethodDefault(Table& table)
{
	const size_t size = table.size();
	const double h = abs(table[1][0] - table[0][0]);

	for (size_t i = 1; i < size; ++i)
	{
		double xPrev = table[i - 1][0];
		double yPrev = table[i - 1][1];

		table[i][1] = yPrev + h * f(xPrev, yPrev);
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}

}

void EulerMethodModifiedFirst(Table& table)
{
	const size_t size = table.size();
	const double h = abs(table[1][0] - table[0][0]);

	for (size_t i = 1; i < size; ++i)
	{
		double xPrev = table[i - 1][0];
		double yPrev = table[i - 1][1];
		double yHalf = yPrev + h / 2 * f(xPrev, yPrev);

		table[i][1] = yPrev + h * f(xPrev + h / 2, yHalf);
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}

}

void EulerMethodModifiedSecond(Table& table)
{
	const size_t size = table.size();
	const double h = abs(table[1][0] - table[0][0]);

	for (size_t i = 1; i < size; ++i)
	{
		double xPrev = table[i - 1][0];
		double yPrev = table[i - 1][1];
		double yTemp = yPrev + h * f(xPrev, yPrev);

		table[i][1] = yPrev + h / 2 * (f(xPrev, yPrev) + f(table[i][0], yTemp));
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}

}

void AdamsMethod(Table& table, const vector<double>& memb)
{
	const size_t size = table.size();

	Table tableTaylor;
	vector<double> X;
	initPoints(X, 0, 0.1, -2, 2);
	initTable(tableTaylor, 0, X.begin(), X.end());
	TaylorMethod(tableTaylor, memb);

	for (size_t i = 0; i < 5; ++i)
	{
		table[i] = tableTaylor[i];
	}

	for (size_t k = 5; k < size; ++k)
	{
		vector<double> l(5, 1.0);

		for (size_t i = 0; i < 5; ++i)
		{
			for (size_t j = 0; j < 5; ++j)
			{
				if (i == j) continue;
				l[i] *= (table[k][0] - table[j][0]) / (table[i][0] - table[j][0]);
			}
		}

		for (size_t i = 0; i < 5; ++i)
		{
			table[k][1] += table[i][1] * l[i];
		}

		table[k][2] = abs(table[k][1] - solution(table[k][0]));
	}
}

void RungeKuttaMethod(Table& table)
{
	const size_t size = table.size();
	const double h = abs(table[1][0] - table[0][0]);


	for (size_t i = 1; i < size; ++i)
	{
		double xPrev = table[i - 1][0];
		double yPrev = table[i - 1][1];

		double k1 = h * f(xPrev, yPrev);
		double k2 = h * f(xPrev + h / 2, yPrev + k1 / 2);
		double k3 = h * f(xPrev + h / 2, yPrev + k2 / 2);
		double k4 = h * f(xPrev + h, yPrev + k3);

		table[i][1] = yPrev + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
		table[i][2] = abs(table[i][1] - solution(table[i][0]));
	}
}

int main()
{
	cout << fixed << setprecision(8);
	int N;
	double h;

	cout << "#6. Solution of 1st order ODE" << endl;
	cout << "@@ #ID: 10" << endl << endl;
	cout << "@@ Parameters:" << endl;
	cout << "   y' = -3y + y^2" << endl << endl;
	cout << "   y(0) = 1" << endl;
	cout << "   N = "; cin >> N;
	cout << "   h = "; cin >> h;

	vector<double> pointsVector;
	initPoints(pointsVector, 0, h, -2, N);

	cout << endl;
	cout << "@@ Table of values:" << endl;
	cout << "        x      |       y(x)    " << endl;
	Table table;
	initTable(table, 0, pointsVector.begin(), pointsVector.end());
	fillTable(table);
	printTable(table, false);

	vector<double> memb;
	initTaylorCoeffs(memb);

	cout << endl;
	cout << "@@ Taylor:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableTaylor;
	initTable(tableTaylor, 0, pointsVector.begin(), pointsVector.end());
	TaylorMethod(tableTaylor, memb);
	printTable(tableTaylor);

	cout << endl;
	cout << "@@ Adams:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableAdams;
	initTable(tableAdams, 0, pointsVector.begin() + 2, pointsVector.end());
	AdamsMethod(tableAdams, memb);
	printTable(tableAdams);

	cout << endl;
	cout << "@@ Runge-Kutta:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableRungeKutta;
	initTable(tableRungeKutta, 0, pointsVector.begin() + 2, pointsVector.end());
	RungeKuttaMethod(tableRungeKutta);
	printTable(tableRungeKutta);

	cout << endl;
	cout << "@@ Euler default:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableEuler;
	initTable(tableEuler, 0, pointsVector.begin() + 2, pointsVector.end());
	EulerMethodDefault(tableEuler);
	printTable(tableEuler);

	cout << endl;
	cout << "@@ Euler I:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableEulerI;
	initTable(tableEulerI, 0, pointsVector.begin() + 2, pointsVector.end());
	EulerMethodModifiedFirst(tableEulerI);
	printTable(tableEulerI);

	cout << endl;
	cout << "@@ Euler II:" << endl;
	cout << "        x      |       y(x)     |   |y(x) - Sol(x)|" << endl;
	Table tableEulerII;
	initTable(tableEulerII, 0, pointsVector.begin() + 2, pointsVector.end());
	EulerMethodModifiedSecond(tableEulerII);
	printTable(tableEulerII);

	return 0;
}