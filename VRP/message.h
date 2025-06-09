#pragma once
#include<vector>
#include<string>
#include<algorithm>
using namespace std;
extern int Vehicle_max;//Number of vehicles
extern int Capacity_max;//Maximum capacity of the vehicle
extern int Station_max;//The number of dielivery stations

extern double **distance_station;//The distance between dielivery stations
extern double *Early;//The earliest service time of the dielivery station
extern double *Later;//The latest service time of the dielivery station
extern int* Ser_time;//The time required for the service of the dielivery station

extern int* ca_station;//The demand of the dielivery stations


