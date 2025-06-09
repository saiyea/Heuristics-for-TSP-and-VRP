#include"VRP-GA.h"
#define MY_H_FILE
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include <numeric>
#include<iostream>
#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace std;
VRP_GA::VRP_GA(int iPS, int iMaxcycle, float iPC, float iPM, float iP)
{
	PS = iPS;
	Maxcycle = iMaxcycle;
	PC = iPC;
	PM = iPM;
	P = iP;
}

VRP_GA::~VRP_GA()
{

	for (int i = 0; i < PS; i++)
	{
		delete[] Chrom[i].JSequence;
		delete[] Chrom[i].MSequence;
	}
	delete[] Chrom; Chrom = NULL;


	for (int i = 0; i < PS; i++)
	{
		delete[] pTemp[i].JSequence;
		delete[] pTemp[i].MSequence;
	}
	delete[] pTemp; pTemp = NULL;


	delete[] minChrom.JSequence; minChrom.JSequence = NULL;
	delete[] minChrom.MSequence; minChrom.MSequence = NULL;
	delete[] bestChrom.MSequence; bestChrom.MSequence = NULL;


	for (int i = 0; i < PS; i++)
	{
		delete[] newPop[i].JSequence;
		delete[] newPop[i].MSequence;
	}
	delete[] newPop; newPop = NULL;





	cout << "destructor" << endl;


}

void VRP_GA::VRP_GAPP()
{
	cout << "pp run" << endl;
	Tmax = 10000;
	aveBest = 0;
	finalBest = 200000;

	valuecycle = new int[10];
	for (int gph = 0; gph < 1; gph++)
	{
		run_times(Tmax);
		valuecycle[gph] = minChrom.tFitness;
		aveBest += minChrom.tFitness;
		finalBest = min(finalBest, minChrom.tFitness);
	}
	aveBest = aveBest / 1;
	STD = 0;
	for (int gph = 0;gph < 1;gph++)
	{
		STD += pow(valuecycle[gph] - aveBest, 2);
	}
	STD = STD / 1;
	STD = pow(STD, 0.5);
	//Write to document
	string sstr = "XX/XX.txt";
	ofstream fout(sstr, ios::app);
	fout << finalBest << " " << aveBest << " " << STD << endl;
}

void VRP_GA::run_times(double t) 
{
	double start = clock();
	gen = 0;
	int couut = 1;
	int n = 1;
	Initialize();
	Copy(bestChrom, minChrom);
	double end = clock();
	do {
		cout << gen << " " << minChrom.tFitness << endl;
		Generationpp();
		if (bestChrom.tFitness < minChrom.tFitness)
		{
			Copy(bestChrom, minChrom);
		}
		gen++;
		end = clock();
	} while (gen < Maxcycle);

	minChrom.tFitness = newbestdecode(minChrom.JSequence, minChrom.MSequence);
}

void VRP_GA::Initialize()
{
	newPop = new ppChrom[PS];
	oldPop = new ppChrom[PS];

	minChrom.JSequence = new int[Station_max];
	bestChrom.JSequence = new int[Station_max];

	minChrom.MSequence = new int[Vehicle_max];
	bestChrom.MSequence = new int[Vehicle_max];

	
	
	for (int i = 0; i < PS; i++)
	{
		oldPop[i].JSequence = new int[Station_max];
		oldPop[i].MSequence = new int[Vehicle_max];
		for (int j = 0;j < Station_max;j++)
		{
			oldPop[i].JSequence[j] = 0;
		}
		for (int j = 0;j < Vehicle_max;j++)
		{
			oldPop[i].MSequence[j] = 0;
		}
		newPop[i].JSequence = new int[Station_max];
		newPop[i].MSequence = new int[Vehicle_max];
	}


	Chrom = new ppChrom[PS];
	for (int i = 0; i < PS; i++)
	{
		Chrom[i].JSequence = new int [Station_max];
		Chrom[i].MSequence = new int[Vehicle_max];
		for (int j = 0;j < Station_max;j++)
		{
			Chrom[i].JSequence[j] = j+1;
		}
		for (int j = 0;j < Vehicle_max;j++)
		{
			Chrom[i].MSequence[j] = 0;
		}
	}


	int temp1;
	int temp2;
	for (int i = 0;i < PS;i++)
	{
		for (int j = 0;j < Station_max;j++)
		{
			temp1 = rand() % Station_max;
			temp2 = rand() % Station_max;
			swap(Chrom[i].JSequence[temp1], Chrom[i].JSequence[temp2]);
		}
		temp1 = rand() % Vehicle_max;
		temp1 += 1;
		for (int j = 0;j < temp1;j++)
		{
			Chrom[i].MSequence[j] = 1;
		}
	}
	for (int i = 0;i < PS;i++)
	{
		Chrom[i].tFitness = newdecode(Chrom[i].JSequence, Chrom[i].MSequence);
		//cout << "tFitness: " << Chrom[i].tFitness << endl;
	}

	for (int i = 0; i < PS; i++)
	{
		Copy(Chrom[i], oldPop[i]);
	}
	oldPop=sorting(Chrom);
	Copy(oldPop[0], bestChrom);
}


ppChrom* VRP_GA::sorting(ppChrom* Parent)
{
	//Perform bubble sort
	vector<int> temp(PS, 0);
	for (int i = 0;i < PS;i++)
	{
		temp[i] = i;
	}
	
	pTemp = new ppChrom[PS];
	for (int i = 0; i < PS; i++)
	{
		pTemp[i].JSequence = new int[Station_max];
		pTemp[i].MSequence = new int[Vehicle_max];
		for (int j = 0;j < Station_max;j++)
		{
			pTemp[i].JSequence[j] = 0;
		}
		for (int j = 0;j < Vehicle_max;j++)
		{
			pTemp[i].MSequence[j] = 0;
		}
	}
	for (int i = 0;i < PS;i++)
	{
		Copy(Parent[i], pTemp[i]);
	}
	for (int i = 0;i < PS-1;i++)
	{
		for (int j = i + 1;j < PS;j++)
		{
			if (pTemp[i].tFitness > pTemp[j].tFitness)
			{
				swap(temp[i], temp[j]);
				swap(pTemp[i].tFitness, pTemp[j].tFitness);
			}
		}
	}

	for (int i = 0;i < PS;i++)
	{
		Copy(Parent[temp[i]], pTemp[i]);
	}
	return pTemp;
}

void VRP_GA::Generationpp()
{
	int mate1, mate2, mate, i, j;
	int topOldSize = static_cast<int>(PS * P);
	if (topOldSize % 2 == 1)
		topOldSize = topOldSize + 1; 

	int* oldPopIndex = new int[PS];
	int* oldPopFitness = new int[PS];


	for (int i = 0; i < PS; i++)
	{
		oldPopFitness[i] = oldPop[i].tFitness;//
		oldPopIndex[i] = i;
	}

	for (i = 0; i < PS - 1; i++) 
	{
		for (j = i + 1; j < PS; j++)
		{
			if (oldPopFitness[i] > oldPopFitness[j])
			{
				swap(oldPopFitness[i], oldPopFitness[j]);
				swap(oldPopIndex[i], oldPopIndex[j]);
			}
		}
	}

	for (j = 0; j < topOldSize; j++) 
	{
		Copy(oldPop[oldPopIndex[j]], newPop[j]);
	}

	// Crossing and mutation will be performed below.
	for (j = topOldSize; j < PS - 1; j += 2)
	{
		do {
			mate1 = TourSelection();
			mate2 = TourSelection();
		} while (mate1 == mate2);

		JCross(oldPop[mate1], oldPop[mate2], newPop[j], newPop[j + 1]); 
	}	

	for (j = topOldSize; j < PS; j++)
	{
		mate = TourSelection();
		JMutation(oldPop[mate], newPop[j]);
	}  

    //Formation of new populations
	for (i = 0; i < PS; i++)
	{
		Copy(oldPop[i], pTemp[i]);
		Copy(newPop[i], oldPop[i]);
		Copy(pTemp[i], newPop[i]);
	}

	for (int i = 0; i < PS; i++)
	{
		oldPopFitness[i] = oldPop[i].tFitness;
	}
	for (i = 0; i < PS - 1; i++)
	{
		for (j = i + 1; j < PS; j++)
		{
			if (oldPopFitness[i] > oldPopFitness[j])
			{
				swap(oldPopFitness[i], oldPopFitness[j]);
			}
		}
	}
	Copy(oldPop[0], bestChrom);

	delete[] oldPopFitness; oldPopFitness = NULL;
	delete[] oldPopIndex; oldPopIndex = NULL;
}

bool VRP_GA::flip(float probability)
{
	float ppp;
	//ppp = (float)(rand() / RAND_MAX + 1);
	float random = static_cast<float>(rand()) / static_cast<float>(RAND_MAX + 1);//linux
	ppp = abs(random);
	//	cout << "ppp: " << ppp;
	if (ppp <= probability)
		return true;
	else
		return false;
}

int VRP_GA::TourSelection()  
{
	int rand1, rand2;
	do
	{
		//rand1= static_cast(rand()) / static_cast<float>(RAND_MAX + 1);//linux
		rand1 = rand() % PS;
		rand2 = rand() % PS;
		//	cout << rand1 << " " << rand2 << endl;
	} while (rand1 == rand2);

	if (oldPop[rand1].tFitness <= oldPop[rand2].tFitness && flip((float)0.8))
		return rand1;
	else
		return rand2;
}

void VRP_GA::Copy(ppChrom& Parent, ppChrom& Child)
{
	for (int i = 0; i < Station_max; i++)
	{
		Child.JSequence[i] = Parent.JSequence[i];

	}

	for (int i = 0; i < Vehicle_max; i++)
	{
		
		Child.MSequence[i] = Parent.MSequence[i];
	}
	Child.tFitness = Parent.tFitness;
}


void VRP_GA::poxCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2)
{
	int j1 = 0, j2 = 0, j3 = 0, k;
	int* temp1 = new int [Station_max];
	int* temp2 = new int [Station_max];//Intermediate variable, storing chromosomes during the crossover process


	for (int i = 0; i < Station_max; i++)
	{
		temp1[i] = 0;
		temp2[i] = 0;
	}

	int* jobset1 = new int [Station_max];
	int* jobset2 = new int [Station_max];//Intermediate variable, storing chromosomes during the crossover process

	for (int i = 0; i < Station_max; i++)
	{
		jobset1[i] = parent1.JSequence[i];
		jobset2[i] = parent2.JSequence[i];
	}

	do
	{
		j1 = rand() % Station_max;
		j2 = rand() % Station_max;
	} while ((j1 >= j2) || (j1 == 0) || (j2 == Station_max - 1));

	for (int i = 0; i < j1; i++)
	{
		temp1[i] = parent1.JSequence[i];

		temp2[i] = parent2.JSequence[i];

	}

	for (int i = j2; i < Station_max; i++)
	{
		temp1[i]= parent1.JSequence[i];

		temp2[i] = parent2.JSequence[i];

	}


	for (int i = 0; i < Station_max; i++)
	{
		for (int j = 0; j < j1; j++)
		{
			if (jobset2[i] == temp1[j] )
			{
				jobset2[i] = 0;

			}

			if (jobset1[i]== temp2[j])
			{
				jobset1[i] = 0;

			}
		}
		for (int j = j2; j < Station_max; j++)
		{
			if (jobset2[i] == temp1[j] )
			{
				jobset2[i] = 0;
			}
			if (jobset1[i] == temp2[j])
			{
				jobset1[i] = 0;

			}
		}

	}
	k = j1;

	for (int i = 0; i < Station_max; i++)
	{
		if (jobset2[i] != 0 )
		{
			temp1[k] = jobset2[i];

			k++;
		}


	}
	k = j1;

	for (int i = 0; i < Station_max; i++)
	{
		if (jobset1[i] != 0)
		{
			temp2[k] = jobset1[i];
			k++;
		}

	}
	for (int i = 0; i < Station_max; i++)
	{
		child1.JSequence[i]= temp1[i];
		child2.JSequence[i]= temp2[i];

	}
	
	for (int i = 0; i < Vehicle_max; i++)
	{
		child1.MSequence[i] = parent1.MSequence[i];
		if (flip(0.2))
		{
			child1.MSequence[i] = parent2.MSequence[i];
		}
	}

	for (int i = 0; i < Vehicle_max; i++)
	{
		child2.MSequence[i] = parent2.MSequence[i];
		if (flip(0.2))
		{
			child2.MSequence[i] = parent1.MSequence[i];
		}
	}

	child1.tFitness = newdecode(child1.JSequence, child1.MSequence);
	child2.tFitness = newdecode(child2.JSequence, child2.MSequence);

	delete[]jobset1; jobset1 = NULL;
	delete[]jobset2; jobset2 = NULL;
	delete[]temp1;  temp1 = NULL;
	delete[]temp2;  temp2 = NULL;

}


void VRP_GA::joxCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2)
{
	int j1 = 0, j2 = 0, j3 = 0, k;
	int* temp1;
	int* temp2;
	temp1 = new int [Station_max];
	temp2 = new int [Station_max];//Intermediate variable, storing chromosomes during the crossover process


	for (int i = 0; i < Station_max; i++)
	{
		temp1[i] = 0;
		temp2[i] = 0;
	}
	int* jobset1;
	int* jobset2;
	jobset1 = new int [Station_max];
	jobset2 = new int [Station_max];//Intermediate variable, storing chromosomes during the crossover process


	for (int i = 0; i < Station_max; i++)
	{
		jobset1[i] = parent1.JSequence[i];
		jobset2[i] = parent2.JSequence[i];
	}


	do
	{
		j1 = rand() % Station_max;
		j2 = rand() % Station_max;
	} while ((j1 >= j2) || (j1 == 0) || (j2 == Station_max - 1));

	for (int i = j1; i < j2; i++)
	{
		temp1[i] = parent1.JSequence[i];

		temp2[i]= parent2.JSequence[i];

	}

	for (int i = 0; i < Station_max; i++)
	{
		for (int j = j1; j < j2; j++)
		{
			if (jobset2[i] == temp1[j])
			{
				jobset2[i] = 0;
			}

			if (jobset1[i] == temp2[j])
			{
				jobset1[i] = 0;
			}
		}
	}


	k = 0;
	int indexmm = 0;
	for (int i = 0; i < Station_max; i++)
	{
		if (jobset2[i] != 0 )
		{
			temp1[k] = jobset2[i];
			k++;
		}
		if (k >= j1)
		{
			indexmm = i;
			break;
		}
	}
	k = j2;


	for (int i = indexmm + 1; i < Station_max; i++)
	{
		if (jobset2[i]!= 0)
		{
			temp1[k] = jobset2[i];
			k++;
		}
	}

	k = 0;
	for (int i = 0; i < Station_max; i++)
	{
		if (jobset1[i] != 0 )
		{
			temp2[k] = jobset1[i];
			k++;
		}
		if (k >= j1)
		{
			indexmm = i;
			break;
		}
	}
	k = j2;

	for (int i = indexmm + 1; i < Station_max; i++)
	{
		if (jobset1[i]!= 0)
		{
			temp2[k] = jobset1[i];
			k++;
		}
	}


	for (int i = 0; i < Station_max; i++)
	{
		child1.JSequence[i] = temp1[i];
		child2.JSequence[i] = temp2[i];
	}

	

	for (int i = 0; i < Vehicle_max; i++)
	{
		child1.MSequence[i] = parent1.MSequence[i];
		if (flip(0.2))
		{
			child1.MSequence[i] = parent2.MSequence[i];
		}
	}

	for (int i = 0; i < Vehicle_max; i++)
	{
		child2.MSequence[i] = parent2.MSequence[i];
		if (flip(0.2))
		{
			child2.MSequence[i] = parent1.MSequence[i];
		}
	}
	
	child1.tFitness = newdecode(child1.JSequence, child1.MSequence);
	child2.tFitness = newdecode(child2.JSequence, child2.MSequence);



	delete[]jobset1; jobset1 = NULL;
	delete[]jobset2; jobset2 = NULL;
	delete[]temp1;  temp1 = NULL;
	delete[]temp2;  temp2 = NULL;
}


void VRP_GA::JCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2)
{
	int random01;
	random01 = rand() % 2;
	if (flip(PC))	
	{
		if (random01 == 0)
		{
			poxCross(parent1, parent2, child1, child2);
		}
		if (random01 == 1)
		{
			joxCross(parent1, parent2, child1, child2);
		}
	}
	else
	{
		Copy(parent1, child1);
		Copy(parent2, child2);
	}
}

//mutation
void VRP_GA::swapMutation(ppChrom& parent, ppChrom& child)
{
	int point1, point2;
	int a1, a2, b1, b2;

	int* jobset = new int [Station_max];
	
	for (int i = 0; i < Station_max; i++)
	{
		jobset[i] = parent.JSequence[i];
	}
	do
	{
		a1 = rand() % Station_max;
		a2 = rand() % Station_max;
	} while (a1 >= a2);

	do
	{
		b1 = rand() % Station_max;
		b2 = rand() % Station_max;
	} while (b1 >= b2);

	do {
		point1 = rand() % (a2 - a1);
	} while (0);

	do {
		point2 = rand() % (b2 - b1);
	} while (0);

	swap(jobset[point1], jobset[point2]);

//	cout << "MuJS:";
	for (int i = 0; i < Station_max; i++)
	{
		child.JSequence[i] = jobset[i];
	//	cout << child.JSequence[i] << ",";

	}
	//cout << endl;


	int* mset = new int[Vehicle_max];
	point1 = rand() % Vehicle_max;
	for (int i = 0; i < point1+1; i++)
	{
		mset[i] = 1;
	}


	for (int i = 0; i < Vehicle_max; i++)
	{
		child.MSequence[i] = mset[i];
	}

	child.tFitness = newdecode(child.JSequence, child.MSequence);


	delete[] jobset; jobset = NULL;
	delete[] mset; mset = NULL;

}

void VRP_GA::JMutation(ppChrom& parent, ppChrom& child)
{

	if (flip(PM))
	{
		swapMutation(parent, child);
	}
	else
	{
		Copy(parent, child);
	}
}

int VRP_GA::JRouttleSelection() 
{
	double totalSum, partSum;
	double random;
	int j = 0;
	partSum = 0.0;


	random = rand() / 32767.;
	totalSum = random * sumFitness;
	do
	{
		partSum = partSum + oldPop[j].tFitness;
		j = j + 1;
	} while ((partSum < totalSum) && (j < PS));
	return j - 1;
}



double VRP_GA::newdecode(int* coSequence, int* cmSequence)
{

	int* JLtemp = new int[Station_max];


	for (int i = 0; i < Station_max; i++)
	{
		JLtemp[i] = coSequence[i];
	}
	//cout << endl;
	int vehicle_num = 0;//Number of vehicles in chromosomes

	for (int i = 0;i < Vehicle_max;i++)
	{
		if (cmSequence[i] == 1)
		{
			vehicle_num += 1;
		}
	}


	int** VeRoutemp = new int* [vehicle_num];//The sequence of dielivery stations for vehicle travel
	int* Vehicle_num_stat = new int[vehicle_num];//The number of  each dielivery stations needs to travel through
	int* VeCatemp = new int[vehicle_num];//The capacity occupied by the vehicle
	double** Vehicle_station_time = new double* [vehicle_num];//The time when each vehicle arrives at each dielivery stations
	double** Vehicle_station_time_end = new double* [vehicle_num];//The time when each vehicle leaves the dielivery stations


	for (int i = 0;i < vehicle_num;i++)
	{
		VeRoutemp[i] = new int[Station_max];
		Vehicle_station_time[i] = new double[Station_max];
		Vehicle_station_time_end[i] = new double[Station_max];
		for (int j = 0;j < Station_max;j++)
		{
			VeRoutemp[i][j] = 0;
			Vehicle_station_time[i][j] = 0;
			Vehicle_station_time_end[i][j] = 0;
		}
		VeCatemp[i] = 0;
		Vehicle_num_stat[i] = 0;
	}

	int k = 0;
	vector<int> index(vehicle_num,0);
	//According to the current sequence, allocate to the vehicles and update simultaneously the time and route taken by the vehicles
	//Generate the idle time block of each vehicle
	double*** block_vehicle = new double** [vehicle_num];
	for (int i = 0;i < vehicle_num;i++)
	{
		block_vehicle[i] = new double* [Station_max+1];
		for (int j = 0;j < Station_max+1;j++)
		{
			block_vehicle[i][j] = new double[2];
			block_vehicle[i][j][0] = 0;
			block_vehicle[i][j][1] = Later[0];
		}
	}
	k = 0;
	//The first dielivery station
	VeRoutemp[0][0] = JLtemp[0];
	VeCatemp[0] += ca_station[JLtemp[0]];//Update the capacity sequence of the vehicle
	Vehicle_station_time[0][JLtemp[0]-1] = max(Early[JLtemp[0]], distance_station[0][JLtemp[0]]);;//The start service time of the first dielivery station
	Vehicle_station_time_end[0][JLtemp[0] - 1] = Vehicle_station_time[0][JLtemp[0] - 1] + Ser_time[JLtemp[0]];//The service completion time of the first dielivery station
	Vehicle_num_stat[0] += 1;

	//Update the insertable time block: The first vehicle is divided into two idle time periods，[0,early]，[later,later[0]]
	block_vehicle[0][0][0] = 0;
	block_vehicle[0][0][1] = Vehicle_station_time[0][JLtemp[0] - 1];
	block_vehicle[0][1][0] = Vehicle_station_time_end[0][JLtemp[0] - 1];
	block_vehicle[0][1][1] = Later[0];

	//Traverse each free block of each vehicle
	k += 1;//dielivery station number
	index[0] += 1;
	int activity = 0;
	for (int j = 0;j < vehicle_num;j++)
	{
		if (VeCatemp[j] <= Capacity_max)
		{
			for (int g = 0;g < Vehicle_num_stat[j] + 1;g++)//Traverse the free blocks in sequence
			{
				if (k < Station_max)
				{
					if ((block_vehicle[j][g][0] <= Later[JLtemp[k]] && Later[JLtemp[k]] <= block_vehicle[j][g][1]) || (block_vehicle[j][g][0] <= Early[JLtemp[k]] && Early[JLtemp[k]] <= block_vehicle[j][g][1]))
					{
						//1. If the free block is the first free block, that is, it starts from the distribution center
						if (g == 0 && distance_station[0][JLtemp[k]] <= block_vehicle[j][g][1]&& distance_station[0][JLtemp[k]]<= Later[JLtemp[k]])
						{
							Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], distance_station[0][JLtemp[k]]);
							if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][VeRoutemp[j][g]] <= block_vehicle[j][g][1])
							{
								VeRoutemp[j][index[j]] = JLtemp[k];
								VeCatemp[j] += ca_station[JLtemp[k]];
								Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
								Vehicle_num_stat[j] += 1;
								activity = 1;
							}
						}
						//2. If the free block is the middle free block, that is, both sides are occupied by distribution stations
						if (g > 0 && g < Vehicle_num_stat[j])
						{
							if (block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= block_vehicle[j][g][1]&& block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]<= Later[JLtemp[k]])
							{
								Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]);
								if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][VeRoutemp[j][g]] <= block_vehicle[j][g][1])
								{
									VeRoutemp[j][index[j]] = JLtemp[k];
									VeCatemp[j] += ca_station[JLtemp[k]];
									Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
									Vehicle_num_stat[j] += 1;
									activity = 1;
								}
							}
						}
						// 3. If the free block is the last free block, that is, the end point is the distribution center
						if (g == Vehicle_num_stat[j] && block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= block_vehicle[j][g][1]&& block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]<= Later[JLtemp[k]])//或者是Later[0]一样的。
						{
							Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]);
							if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][0] <= block_vehicle[j][g][1])
							{
								VeRoutemp[j][index[j]] = JLtemp[k];
								VeCatemp[j] += ca_station[JLtemp[k]];
								Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
								Vehicle_num_stat[j] += 1;
								activity = 1;
							}
						}
					}
					//It represents a successful insertion in one attempt. 
					//Here, it is necessary to update the time blocks occupied by the delivery stations and the time sequence on each vehicle
					if (activity == 1)
					{
						//Update the time block sequence of this vehicle, from the front to the back row according to the size of its start service time
						vector<double> temp(Vehicle_num_stat[j], 0);
						for (int i = 0;i < Vehicle_num_stat[j];i++)
						{
							temp[i] = Vehicle_station_time[j][VeRoutemp[j][i] - 1];
						}
						for (int i = 0;i < Vehicle_num_stat[j] - 1;i++)
						{
							for (int i1 = i + 1;i1 < Vehicle_num_stat[j];i1++)
							{
								if (temp[i] > temp[i1])
								{
									swap(VeRoutemp[j][i], VeRoutemp[j][i1]);
									swap(temp[i], temp[i1]);
								}
							}
						}
						//Blank block update
						block_vehicle[j][0][0] = 0;
						block_vehicle[j][0][1] = Vehicle_station_time[j][VeRoutemp[j][0] - 1];
						for (int i = 1;i < Vehicle_num_stat[j];i++)
						{
							block_vehicle[j][i][0] = Vehicle_station_time_end[j][VeRoutemp[j][i - 1] - 1];
							block_vehicle[j][i][1] = Vehicle_station_time[j][VeRoutemp[j][i] - 1];
						}
						block_vehicle[j][Vehicle_num_stat[j]][0] = Vehicle_station_time_end[j][VeRoutemp[j][Vehicle_num_stat[j] - 1] - 1];;
						block_vehicle[j][Vehicle_num_stat[j]][1] = Later[0];
						k += 1;
						index[j] += 1;
						g = -1;
						//If j=0, it indicates that the task has been completed on vehicle 0 at this time. 
						// If j=1, it indicates that it doesn't work on vehicle 0 but is completed on vehicle 1. 
						// It's time to start the selection of the next delivery station.
						if (j > 0)
						{
							j = -1;
						}
						activity = 0;
					}
				}
			}
		}
	}

	int mark = 0;
	if (k < Station_max)
	{
		mark = 1;
	}

	//Calculate the objective function = the number of vehicles used *30+ the distance traveled
	double OBJ = 0;
	if (mark == 1)
	{
		OBJ = 50000;
	}
	if (mark == 0)
	{
		for (int i = 0;i < vehicle_num;i++)
		{
			if (Vehicle_num_stat[i] > 0)
			{
				OBJ += distance_station[0][VeRoutemp[i][0]];
				for (int j = 1;j < Vehicle_num_stat[i];j++)
				{
					OBJ += distance_station[VeRoutemp[i][j - 1]][VeRoutemp[i][j]];
				}
				OBJ += distance_station[VeRoutemp[i][Vehicle_num_stat[i] - 1]][0];
			}
		}
	}
	OBJ += vehicle_num * 30;

	//Release memory
	delete[] JLtemp; JLtemp= NULL;
	delete[] Vehicle_num_stat; Vehicle_num_stat = NULL;
	delete[] VeCatemp; VeCatemp = NULL;
	

	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] VeRoutemp[i];
	}
	delete[] VeRoutemp; VeRoutemp = NULL;

	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] Vehicle_station_time[i];
	}
	delete[] Vehicle_station_time; Vehicle_station_time = NULL;

	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] Vehicle_station_time_end[i];
	}
	delete[] Vehicle_station_time_end; Vehicle_station_time_end = NULL;


	for (int i = 0;i < vehicle_num;i++)
	{
		for (int j = 0;j < Station_max + 1;j++)
		{
			delete[] block_vehicle[i][j];
		}
		delete[] block_vehicle[i];
	}
	delete[] block_vehicle; block_vehicle = NULL;


	return OBJ;
}

//Similar to the above, only the detailed information of the optimal solution is output
double VRP_GA::newbestdecode(int* coSequence, int* cmSequence)
{

	int* JLtemp = new int[Station_max];


	for (int i = 0; i < Station_max; i++)
	{
		JLtemp[i] = coSequence[i];
		cout << JLtemp[i] << ",";
	}
	cout << endl;
	int vehicle_num = 0;

	for (int i = 0;i < Vehicle_max;i++)
	{
		if (cmSequence[i] == 1)
		{
			vehicle_num += 1;
		}
	}

	double* VeTtemp = new double[vehicle_num];
	int** VeRoutemp = new int* [vehicle_num];
	int* Vehicle_num_stat = new int[vehicle_num];
	int* VeCatemp = new int[vehicle_num];
	double** Vehicle_station_time = new double* [vehicle_num];
	double** Vehicle_station_time_end = new double* [vehicle_num];


	for (int i = 0;i < vehicle_num;i++)
	{
		VeRoutemp[i] = new int[Station_max];
		Vehicle_station_time[i] = new double[Station_max];
		Vehicle_station_time_end[i] = new double[Station_max];
		for (int j = 0;j < Station_max;j++)
		{
			VeRoutemp[i][j] = 0;
			Vehicle_station_time[i][j] = 0;
			Vehicle_station_time_end[i][j] = 0;
		}
		VeTtemp[i] = 0;
		VeCatemp[i] = 0;
		Vehicle_num_stat[i] = 0;
	}

	int k = 0;
	vector<int> index(vehicle_num, 0);
	double*** block_vehicle = new double** [vehicle_num];
	for (int i = 0;i < vehicle_num;i++)
	{
		block_vehicle[i] = new double* [Station_max + 1];
		for (int j = 0;j < Station_max + 1;j++)
		{
			block_vehicle[i][j] = new double[2];
			block_vehicle[i][j][0] = 0;
			block_vehicle[i][j][1] = Later[0];

		}
	}

	k = 0;
	VeRoutemp[0][0] = JLtemp[0];
	VeCatemp[0] += ca_station[JLtemp[0]];
	Vehicle_station_time[0][JLtemp[0] - 1] = max(Early[JLtemp[0]], distance_station[0][JLtemp[0]]);
	Vehicle_station_time_end[0][JLtemp[0] - 1] = Vehicle_station_time[0][JLtemp[0] - 1] + Ser_time[JLtemp[0]];
	Vehicle_num_stat[0] += 1;

	block_vehicle[0][0][0] = 0;
	block_vehicle[0][0][1] = Vehicle_station_time[0][JLtemp[0] - 1];
	block_vehicle[0][1][0] = Vehicle_station_time_end[0][JLtemp[0] - 1];
	block_vehicle[0][1][1] = Later[0];

	k += 1;
	index[0] += 1;
	int activity = 0;
	for (int j = 0;j < vehicle_num;j++)
	{
		if (VeCatemp[j] <= Capacity_max)
		{
			for (int g = 0;g < Vehicle_num_stat[j] + 1;g++)
			{
				if (k < Station_max)
				{
					
					if ((block_vehicle[j][g][0] <= Later[JLtemp[k]] && Later[JLtemp[k]] <= block_vehicle[j][g][1]) || (block_vehicle[j][g][0] <= Early[JLtemp[k]] && Early[JLtemp[k]] <= block_vehicle[j][g][1]))
					{
						
						if (g == 0 && distance_station[0][JLtemp[k]] <= block_vehicle[j][g][1] && distance_station[0][JLtemp[k]] <= Later[JLtemp[k]])
						{
							Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], distance_station[0][JLtemp[k]]);
							if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][VeRoutemp[j][g]] <= block_vehicle[j][g][1])
							{
								VeRoutemp[j][index[j]] = JLtemp[k];
								VeCatemp[j] += ca_station[JLtemp[k]];
								Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
								Vehicle_num_stat[j] += 1;
								activity = 1;
							}
						}
						
						if (g > 0 && g < Vehicle_num_stat[j])
						{
							if (block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= block_vehicle[j][g][1] && block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= Later[JLtemp[k]])
							{
								Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]);
								if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][VeRoutemp[j][g]] <= block_vehicle[j][g][1])
								{
									VeRoutemp[j][index[j]] = JLtemp[k];
									VeCatemp[j] += ca_station[JLtemp[k]];
									Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
									Vehicle_num_stat[j] += 1;
									activity = 1;
								}
							}
						}
						if (g == Vehicle_num_stat[j] && block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= block_vehicle[j][g][1] && block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]] <= Later[JLtemp[k]])//或者是Later[0]一样的。
						{
							Vehicle_station_time[j][JLtemp[k] - 1] = max(Early[JLtemp[k]], block_vehicle[j][g][0] + distance_station[VeRoutemp[j][g - 1]][JLtemp[k]]);
							if (Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]] + distance_station[JLtemp[k]][0] <= block_vehicle[j][g][1])
							{
								VeRoutemp[j][index[j]] = JLtemp[k];
								VeCatemp[j] += ca_station[JLtemp[k]];
								Vehicle_station_time_end[j][JLtemp[k] - 1] = Vehicle_station_time[j][JLtemp[k] - 1] + Ser_time[JLtemp[k]];
								Vehicle_num_stat[j] += 1;
								activity = 1;
							}
						}
					}
					if (activity == 1)
					{
						vector<double> temp(Vehicle_num_stat[j], 0);
						for (int i = 0;i < Vehicle_num_stat[j];i++)
						{
							temp[i] = Vehicle_station_time[j][VeRoutemp[j][i] - 1];
						}
						for (int i = 0;i < Vehicle_num_stat[j] - 1;i++)
						{
							for (int i1 = i + 1;i1 < Vehicle_num_stat[j];i1++)
							{
								if (temp[i] > temp[i1])
								{
									swap(VeRoutemp[j][i], VeRoutemp[j][i1]);
									swap(temp[i], temp[i1]);
								}
							}
						}
						block_vehicle[j][0][0] = 0;
						block_vehicle[j][0][1] = Vehicle_station_time[j][VeRoutemp[j][0] - 1];
						for (int i = 1;i < Vehicle_num_stat[j];i++)
						{
							block_vehicle[j][i][0] = Vehicle_station_time_end[j][VeRoutemp[j][i - 1] - 1];
							block_vehicle[j][i][1] = Vehicle_station_time[j][VeRoutemp[j][i] - 1];
						}
						block_vehicle[j][Vehicle_num_stat[j]][0] = Vehicle_station_time_end[j][VeRoutemp[j][Vehicle_num_stat[j] - 1] - 1];;
						block_vehicle[j][Vehicle_num_stat[j]][1] = Later[0];
						k += 1;
						index[j] += 1;
						g = -1;
						if (j > 0)
						{
							j = -1;
						}
						activity = 0;
					}
				}
			}
		}
	}
	
	cout << "#####finish" << endl;
	cout << "k: " << k << endl;
	for (int j = 0;j < vehicle_num;j++)
	{
		for (int i = 0;i < Vehicle_num_stat[j];i++)
		{
			cout << VeRoutemp[j][i] << ",";
		}
		cout << endl;
	}

	int mark = 0;
	if (k < Station_max)
	{
		mark = 1;
	}


	double OBJ = 0;
	if (mark == 1)
	{
		OBJ = 50000;
	}
	
	if (mark == 0)
	{
		for (int i = 0;i < vehicle_num;i++)
		{
			if (Vehicle_num_stat[i] > 0)
			{
				OBJ += distance_station[0][VeRoutemp[i][0]];
				for (int j = 1;j < Vehicle_num_stat[i];j++)
				{
					OBJ += distance_station[VeRoutemp[i][j - 1]][VeRoutemp[i][j]];
				}
				OBJ += distance_station[VeRoutemp[i][Vehicle_num_stat[i] - 1]][0];
			}
		}
	}
	OBJ += vehicle_num * 30;
	cout << "OBJ: " << OBJ << endl;

	int sev_cout = 0;
	int capa_cout = 0;
	int capa_num = 0;
	for (int i = 0;i < vehicle_num;i++)
	{
		capa_cout = 0;
		for (int j = 0;j < Vehicle_num_stat[i];j++)
		{
			cout << "text: " << Early[VeRoutemp[i][j]] << "," << Later[VeRoutemp[i][j]] << "," << Vehicle_station_time[i][VeRoutemp[i][j] - 1] << "," << Vehicle_station_time_end[i][VeRoutemp[i][j] - 1] << endl;
			if ((Vehicle_station_time[i][VeRoutemp[i][j] - 1] >= Early[VeRoutemp[i][j]] && Vehicle_station_time[i][VeRoutemp[i][j] - 1] <= Later[VeRoutemp[i][j]]) || (Vehicle_station_time_end[i][VeRoutemp[i][j] - 1] >= Early[VeRoutemp[i][j]] && Vehicle_station_time_end[i][VeRoutemp[i][j] - 1] <= Later[VeRoutemp[i][j]]))
			{
				sev_cout += 1;
			}
			capa_cout += ca_station[VeRoutemp[i][j]];
		}
		if (capa_cout <= 200)
		{
			capa_num += 1;
		}
	}
	cout << "sev_cout: " << sev_cout << endl;
	cout << "capa_num: " << capa_num << endl;

	delete[] JLtemp; JLtemp = NULL;
	delete[] VeTtemp; VeTtemp = NULL;
	delete[] Vehicle_num_stat; Vehicle_num_stat = NULL;
	delete[] VeCatemp; VeCatemp = NULL;


	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] VeRoutemp[i];
	}
	delete[] VeRoutemp; VeRoutemp = NULL;

	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] Vehicle_station_time[i];
	}
	delete[] Vehicle_station_time; Vehicle_station_time = NULL;

	for (int i = 0; i < vehicle_num; i++)
	{
		delete[] Vehicle_station_time_end[i];
	}
	delete[] Vehicle_station_time_end; Vehicle_station_time_end = NULL;


	for (int i = 0;i < vehicle_num;i++)
	{
		for (int j = 0;j < Station_max + 1;j++)
		{
			delete[] block_vehicle[i][j];
		}
		delete[] block_vehicle[i];
	}
	delete[] block_vehicle; block_vehicle = NULL;

	return OBJ;
}