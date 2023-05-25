#include "Vec2.h"
using namespace std;

class Wall
{
public:
    array<double,2> Origin, End, u, n;
    double NormU;
    bool diff;
    
    Wall(array<double,2> O, array<double,2> E)
    {
        Origin = O;
        End = E;
        
        u = ArraySubstract(Origin, End);
        NormU = NormArrays(u);
        u = UnitaryVector(u, NormU);
        n[0] = u[1]; n[1] = - u[0];
    }
};
