#include"VRP-GA.h"
#include"message.h"
#include"readin.h"
#include <random>
#include <iostream>
using namespace std;
int main()
{
	//srand((unsigned)time(NULL));
	srand(static_cast<unsigned int>(time(NULL)));
	int Particle ;
	int Maxcycle ;
	float PC;
	float PM;
	float P;
	string filePath = "XX/XX";
	vector<string>files;
	string format = ".txt";
	GetAllFormatFiles(filePath, files, format);
	int size = files.size();
	for (int s = 0; s < size; s++)
	{
		ifstream fileProcess;
		fileProcess.open(files[s], ios::in);
		if (!fileProcess.is_open())
		{
			cout << "Failed to read the file" << endl;
		}
		string sstr = "XX/XX.txt";
		ofstream fout(sstr, ios::app);
		fout << s + 1 << " ";
		fileProcess >> Vehicle_max;
		fileProcess >> Capacity_max;
		fileProcess >> Station_max;
		distance_station = new double* [Station_max+1];
		for (int i = 0;i < Station_max + 1;i++)
		{
			distance_station[i] = new double[Station_max + 1];
		}

		int aa;
		int** coord;//Coordinates of dielivery stations
		coord = new int* [Station_max + 1];
		for (int i = 0;i < Station_max + 1;i++)
		{
			coord[i] = new int [2];
		}
		ca_station = new int[Station_max + 1];
		Early = new double[Station_max + 1];
		Later = new double[Station_max + 1];
		Ser_time = new int[Station_max + 1];
		for (int i = 0; i < Station_max+1; i++)
		{
			fileProcess >> aa;//Station number
			fileProcess >> coord[i][0];
			fileProcess >> coord[i][1];
			fileProcess >> ca_station[i];
			fileProcess >> Early[i];
			fileProcess >> Later[i];
			fileProcess >> Ser_time[i];
		}
		//Calculate the distance between dielivery stations based on the coordinates
		for (int i = 0;i < Station_max + 1;i++)
		{
			for (int j = 0;j < Station_max + 1;j++)
			{
				distance_station[i][j] = sqrt(pow(coord[i][0] - coord[j][0], 2) + pow(coord[i][1] - coord[j][1], 2));
			}
		}
		VRP_GA a(200, 500, 0.7, 0.6, 0.1);
		a.VRP_GAPP();
		
	}
}