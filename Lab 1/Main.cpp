#include <iostream>
#include <vector>


struct Package
{
	double FirstAnswer,
		   FinalAnswer,
		   LastLenght;
	int NumberOfIterations;
};

double F (double x)
{
	return 1.2 * pow(x, 4) + 2 * pow(x, 3) - 13 * pow(x, 2) - 14.2 * x - 24.1;
}

double FDerivative (double x)
{
	return 4.8 * pow(x, 3) + 6 * pow(x, 2) - 26 * x - 14.2;
}

// Root separation
std::vector <std::pair <double, double>> RootSeparation (std::vector <std::pair <double, double>> vect, int N)
{ 
	std::vector <std::pair <double, double>> newVector;
	std::vector <std::pair <double, double>>::iterator Iter = vect.begin();

	while (Iter != vect.end())
	{
		double h = (Iter->second - Iter->first) / N;
		double x1 = Iter->first, 
			   x2 = x1 + h, 
			   y1 = F(x1),
			   y2;

		while (x2 <= Iter->second)
		{
			y2 = F(x2);
			if (y1 * y2 <= 0)
				newVector.push_back(std::make_pair(x1, x2));
			x1 = x2;
			x2 = x1 + h;
			y1 = y2;
		}
		Iter++;
	}

	return newVector;
}

// Bisector method
Package BisectorMethod(std::pair <double, double> p, double e)
{
	Package final;
	double a = p.first, 
		   b = p.second, 
		   c;
	int counter = 0;
	bool FirstRootIdent = false;

	while (b - a > 2 * e)
	{
		counter++;
		c = (a + b) / 2;

		if (F(a) * F(c) <= 0)
		{
			if (!FirstRootIdent)
			{
				FirstRootIdent = true;
				final.FirstAnswer = p.first;
			}
			else
				final.FirstAnswer = p.second;

			b = c;
		}
		else
		{
			if (!FirstRootIdent)
			{
				FirstRootIdent = true;
				final.FirstAnswer = p.second;
			}

			a = c;
		}
	}

	final.FinalAnswer = (a + b) / 2;
	final.LastLenght = (b - a) / 2;
	final.NumberOfIterations = counter;

	return final;
}

// Newton's method
Package NewtonMethod(std::pair <double, double> p, double e)
{
	Package final;
	double x1 = p.first, 
		   x2 = p.second;
	int counter = 0;

	while (abs(x1 - x2) > e)
	{
		x1 = x2;
		x2 = x1 - F(x1) / FDerivative(x1);
		counter++;
	}

	final.FirstAnswer = p.first;
	final.FinalAnswer = x2;
	final.NumberOfIterations = counter;
	final.LastLenght = abs(x1 - x2);

	return final;
}

// Modified Newton's method
Package ModifiedNewtonMethod(std::pair <double, double> p, double e)
{
	Package final;
	double x1 = p.first, 
		   x2 = p.second;
	int counter = 0;

	while (abs(x1 - x2) > e)
	{
		x1 = x2;
		x2 = x1 - F(x1) / FDerivative(p.first);
		counter++;
	}

	final.FirstAnswer = p.first;
	final.FinalAnswer = x2;
	final.NumberOfIterations = counter;
	final.LastLenght = abs(x1 - x2);

	return final;
}

// Secant Method
Package SecantMethod(std::pair <double, double> p, double e)
{
	Package final;
	double x1,
		   x2 = p.first,
		   x3 = p.second;
	int counter = 0;

	while (abs(x3 - x2) > e)
	{
		counter++;
		x1 = x2;
		x2 = x3;
		x3 = x2 - F(x2) * (x2 - x1) / (F(x2) - F(x1));
	}

	final.FirstAnswer = p.first;
	final.FinalAnswer = x2;
	final.NumberOfIterations = counter;
	final.LastLenght = abs(x3 - x2);

	return final;
}


int main()
{
	std::vector <std::pair <double, double>> IntervalVector;
	std::vector <std::pair <double, double>>::iterator Iter;
	Package results;
	double A, B, e;
	int N, i = 1;
	bool out = false, WasThereSeparation = false;
	std::string answer;

	// Beggining
	std::cout << "@@ Numerical methods for solving nonlinear equations \n"
			  << '\n'
			  << "@@ Parameters: \n"
			  << "    f(x) = 1.2x^4 + 2x^3 - 13x^2 - 14.2x - 24.1 \n"
			  << "    A = ";
	std::cin >> A;
	std::cout << "    B = ";
	std::cin >> B;
	std::cout << "    e = ";
	std::cin >> e;

	IntervalVector.push_back(std::make_pair(A, B));

	// Root separation
	std::cout << '\n'
			  << "@@ Root separation procedure: \n"
			  << "    Hint: Choose N greater than 1 \n";

	while (true)
	{
		std::cout << "    N = ";
		std::cin >> N;

		while (N <= 1)
		{
			std::cout << "    Oh no! N is wrong. Finish entering? \n"
					  << "    # ans: ";
			std::cin >> answer;

			if (answer == "yes")
			{
				out = true;
				break;
			}
			else if (answer == "no")
			{
				out = false;
				std::cout << "    N = ";
				std::cin >> N;
				continue;
			}
			else
				std::cout << "    Not recognized :( \n";
		}

		if (out)
			break;

		WasThereSeparation = true;

		IntervalVector = RootSeparation(IntervalVector, N);

		std::cout << '\n'
			<< "/!\\ Number of intervals: " << IntervalVector.size() << std::endl;

		for (Iter = IntervalVector.begin(); Iter != IntervalVector.end(); Iter++)
			std::cout << "    [" << Iter->first << ", " << Iter->second << "]; \n";

		while (true)
		{
			std::cout << '\n'
					  << "    OK, let's split the intervals again? \n"
					  << "    # ans: ";
			std::cin >> answer;

			if (answer == "yes")
			{
				out = false;
				break;
			}
			else if (answer == "no")
			{
				out = true;
				break;;
			}
			else
				std::cout << "    Not recognized :( \n";
		}

		if (out)
			break;
	}

	if (!WasThereSeparation)
	{
		std::cout << '\n'
				  << "Oops, something went wrong. Start again!";
		return 404;
	}

	std::cout << "\n"
			  << "@@ Root Refinement: \n";

	for (Iter = IntervalVector.begin(); Iter != IntervalVector.end(); Iter++, i++)
	{
		// Bisector
		std::cout << "/!\\  #" << i << " root: \n";
		std::cout << "  Bisector Method: \n";

		results = BisectorMethod(*Iter, e);

		std::cout << "    Initial Root: " << results.FirstAnswer << std::endl;
		std::cout << "    Number of iterations: " << results.NumberOfIterations << std::endl;
		std::cout << "    Final Root: " << results.FinalAnswer << std::endl;
		std::cout << "    |x_m - x_{m - 1}|: " << results.LastLenght << std::endl;
		std::cout << "    Discrepancy ( |f(x)| ): " << abs(F(results.FinalAnswer)) << std::endl << std::endl;

		// Newton
		std::cout << "  Newton's Method: \n";

		results = NewtonMethod(*Iter, e);

		std::cout << "    Initial Root: " << results.FirstAnswer << std::endl;
		std::cout << "    Number of iterations: " << results.NumberOfIterations << std::endl;
		std::cout << "    Final Root: " << results.FinalAnswer << std::endl;
		std::cout << "    |x_m - x_{m - 1}|: " << results.LastLenght << std::endl;
		std::cout << "    Discrepancy ( |f(x)| ): " << abs(F(results.FinalAnswer)) << std::endl << std::endl;

		// Modified Newton
		std::cout << "  Modified Newton's Method: \n";

		results = ModifiedNewtonMethod(*Iter, e);

		std::cout << "    Initial Root: " << results.FirstAnswer << std::endl;
		std::cout << "    Number of iterations: " << results.NumberOfIterations << std::endl;
		std::cout << "    Final Root: " << results.FinalAnswer << std::endl;
		std::cout << "    |x_m - x_{m - 1}|: " << results.LastLenght << std::endl;
		std::cout << "    Discrepancy ( |f(x)| ): " << abs(F(results.FinalAnswer)) << std::endl << std::endl;

		// Secant
		std::cout << "  The Secant Method: \n";

		results = SecantMethod(*Iter, e);

		std::cout << "    Initial Root: " << results.FirstAnswer << std::endl;
		std::cout << "    Number of iterations: " << results.NumberOfIterations << std::endl;
		std::cout << "    Final Root: " << results.FinalAnswer << std::endl;
		std::cout << "    |x_m - x_{m - 1}|: " << results.LastLenght << std::endl;
		std::cout << "    Discrepancy ( |f(x)| ): " << abs(F(results.FinalAnswer)) << std::endl << std::endl;
	}

	return 0;
}
