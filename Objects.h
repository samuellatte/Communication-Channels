#include "Vec2.h"
using namespace std;

class Wall
{
public:
    array<double,2> Origin, End, u, n;
    double NormU;
    
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

class Corner
{
public:
    array<double,2> Point;
    
    Corner(array<double,2> pt)
    {
        Point = pt;
    }
};

class Rectangle
{
public:
    array<double,2> bottomLeft, topRight;
    
    Rectangle(array<double,2> bL, array<double,2> tR)
    {
        bottomLeft = bL;
        topRight = tR;
    }
};

class MiddleStreet
{
public:
    array<double,2> Point;
    
    MiddleStreet(array<double,2> pt)
    {
        Point = pt;
    }
};
