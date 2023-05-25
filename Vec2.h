#include <array>
#include <vector>
using namespace std;


double NormArrays(array<double,2> X)
{
    return sqrt(pow(X[0],2) + pow(X[1],2));
}

array<double,2> ArraySubstract(const array<double,2>& O, const array<double,2>& E)
{
    array<double,2> Out;
    for(int i = 0; i < 2; i++){
        Out[i] = E[i] - O[i];
    }
    return Out;
}

array<double,2> ArrayProduct(const array<double,2>& O, const array<double,2>& E)
{
    array<double,2> Out;
    for(int i = 0; i < 2; i++){
        Out[i] = E[i] * O[i];
    }
    return Out;
}

double InnerProduct(const array<double,2> O, const array<double,2> E)
{
    double Out = 0.;
    for(int i = 0; i < 2; i++){
        Out += E[i] * O[i];
    }
    return Out;
}

double Sign(const double& X)
{
    return (X > 0) ? 1 : ((X < 0) ? -1 : 0);
}

array<double,2> ConstArrayProduct(const double X, const array<double,2>& E)
{
    array<double,2> Out;
    for(int i = 0; i < 2; i++){
        Out[i] = E[i] * X;
    }
    return Out;
}

array<double,2> UnitaryVector(const array<double,2>& I, const double& norm)
{
    array<double,2> Out;
    for(int i = 0; i < 2; i++){
        Out[i] = I[i]/norm;
    }
    return Out;
}

array<double,2> PolarRpz(const complex<double> E)
{
    array<double,2> ModAndPhase;
    ModAndPhase[0] = sqrt(pow(real(E),2) + pow(imag(E),2));
    ModAndPhase[1] = atan(imag(E)/real(E));
    return ModAndPhase;
}

double ArraySum(vector<double> vect)
{
    double s = vect.size(), sum = 0.;

    for(int i = 0; i < s; i++){
        sum += vect[i];
    }
    return sum;
}


double sinc(double X){
    if(X == 0){
        return 1;
    }
    return sin(Pi*X)/(Pi*X);
}


