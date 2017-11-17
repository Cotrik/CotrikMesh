#ifndef CONSTANTCONSTANT_HEADER
#define CONSTANTCONSTANT_HEADER
#define NOMINMAX
//#include <windows.h>
#include<math.h>

const double EPS = 1.0e-6;
const double PI = 3.1415926535898;
const float MIN_RANGE = -1e+20;
const float MAX_RANGE = +1e+20;

#define DISTANCE(a,b,c)	{a=sqrt(pow((b)[0]-(c)[0],2)+pow((b)[1]-(c)[1],2)+pow((b)[2]-(c)[2],2));}
#define CROSSVECTOR3(a,b,c)       {(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1]; \
	(a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2]; \
	(a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}
const int hex_tetra_table[8][4] =
{
{ 0, 3, 4, 1 },
{ 1, 0, 5, 2 },
{ 2, 1, 6, 3 },
{ 3, 2, 7, 0 },
{ 4, 7, 5, 0 },
{ 5, 4, 6, 1 },
{ 6, 5, 7, 2 },
{ 7, 6, 4, 3 }, };

#endif
