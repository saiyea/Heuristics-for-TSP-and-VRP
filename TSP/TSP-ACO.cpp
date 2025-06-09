#include <unordered_set>
#include <fstream>
#include <numeric>
#include <ctime>
#include<iostream>
#include<cmath>
#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include<io.h>
#include<algorithm>
#include<tuple>
using namespace std;

const double ALPHA = 1.0;
const double BETA = 2.0;
const double RHO = 0.4;
const double Q = 100.0;



//Read instances from the folder
void GetAllFormatFiles(string path, vector<string>& files, string format)
{

	long long hFile = 0;
	struct _finddata_t fileinfo;
	string p;
	//The _findfirst function looks up a file. (const char*, struct _finddata_t*) The first parameter is the file name, 
	// and "*.*" can be used to look up the file. The second parameter is the pointer of the structure. If the search results are found, it returns the file handle. If it fails, it returns -1.
	if ((hFile = _findfirst(p.assign(path).append("\\*" + format).c_str(), &fileinfo)) != -1) 
	{
		do
		{
			if ((fileinfo.attrib & _A_SUBDIR)) 
				
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					//files.push_back(p.assign(path).append("\\").append(fileinfo.name) );
					GetAllFormatFiles(p.assign(path).append("\\").append(fileinfo.name), files, format);
				}
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);

		_findclose(hFile);
	}
}


//Ant route calculation
vector<double> distance_cout_ant(const vector<vector<int>>& ant_route1, const vector<vector<double>>& city_distan1)
{
	int rows = ant_route1.size();
	int cols = ant_route1[0].size();
	vector<double> result(rows,0);

	for (int i = 0;i < rows;i++)
	{
		for (int j = 0;j < cols-1;j++)
		{
			result[i] += city_distan1[ant_route1[i][j]][ant_route1[i][j+1]];
		}
		result[i] = result[i] + city_distan1[ant_route1[i][cols - 1]][ant_route1[i][0]];//Return to the starting city
	}
	return result;
}

//Pheromone update
vector<vector<double>> update_pheromones(const vector<vector<double>> & phero1,const vector<vector<int>> &ant_route1,const vector<double> &ant_distan1)
{
	int rows = ant_route1.size();
	int cols = ant_route1[0].size();

	vector<vector<double>> result(cols, vector<double>(cols, 0));

	//1.Pheromone volatilization
	for (int i = 0;i < cols;i++)
	{
		for (int j = 0;j < cols;j++)
		{
			result[i][j] = phero1[i][j] * (1 - RHO);
		}
	}

	//2.1 Calculate the last 20% of the added pheromones. If so, change it to the reduced ones
	vector<double> temp_ant_distan(rows, 0);
	vector<int> temp_num(rows, 0);
	for (int i = 0;i < rows;i++)
	{
		temp_ant_distan[i] = ant_distan1[i];
		temp_num[i] = i;
	}

	for (int i = 0;i < rows-1;i++)
	{
		for (int j = i+1;j < rows;j++)
		{
			if (temp_ant_distan[i] < temp_ant_distan[j])
			{
				swap(temp_ant_distan[i], temp_ant_distan[j]);
				swap(temp_num[i], temp_num[j]);
			}
		}
	}
	for (int i = 0;i < int(rows*0.2);i++)
	{
		for (int j = 0;j < cols - 1;j++)
		{
			int temp1 = ant_route1[temp_num[i]][j];
			int temp2 = ant_route1[temp_num[i]][j + 1];
			result[temp1][temp2] -= Q / temp_ant_distan[i];
			result[temp2][temp1] -= Q / temp_ant_distan[i];
		}
	}

	//2.2 The optimal pheromone has increased again
	for (int i = 0;i < 10;i++)
	{
		for (int j = 0;j < cols - 1;j++)
		{
			int temp1 = ant_route1[temp_num[rows-1-i]][j];
			int temp2 = ant_route1[temp_num[rows-1- i]][j + 1];
			result[temp1][temp2] += Q / ant_distan1[temp_num[rows-1 - i]];
			result[temp2][temp1] += Q / ant_distan1[temp_num[rows -1- i]];
		}
	}


	//2.3 All pheromones increase
	for (int i = 0;i < rows;i++)
	{
		for (int j = 0;j < cols-1;j++)
		{
			int temp1 = ant_route1[i][j];
			int temp2 = ant_route1[i][j + 1];
			result[temp1][temp2] += Q / ant_distan1[i];
			result[temp2][temp1] += Q / ant_distan1[i];
		}
		result[ant_route1[i][cols - 1]][ant_route1[i][0]] += Q / ant_distan1[i];
		result[ant_route1[i][0]][ant_route1[i][cols - 1]] += Q / ant_distan1[i];
	}

	return result;
}

//The next city that ants visit
int visit_city_update(const vector<int> visit1, const vector<vector<double>>Pselect1, const int allvisit1)
{
	int rows = Pselect1.size();
	//Count the unvisited cities
	vector<int>no_visit_city(rows - allvisit1, rows);
	int k = 0;
	int cout_no_equal;
	for (int i = 0;i < rows;i++)
	{
		cout_no_equal=0;
		for (int j = 0;j < allvisit1;j++)
		{
			if (i != visit1[j])
			{
				cout_no_equal += 1;
			}
		}
		if (cout_no_equal == allvisit1)
		{
			no_visit_city[k] = i;
			//cout << "k: " << k << "no_visit_city: " << no_visit_city[k];
			k++;
		}
	}
	//Calculate the probability of not visiting the city
	double SUM;//Denominator
	SUM = 0;
	for (int i = 0;i < rows - allvisit1;i++)
	{
		SUM += Pselect1[visit1[allvisit1 - 1]][no_visit_city[i]];
	}
	//cout << "allvisit1: " << allvisit1 << endl;
	vector<vector<double>> Pselect_actually(rows, vector<double>(rows, 0));//The probability of choosing city j for the current city i.
	for (int i = 0;i < rows - allvisit1;i++)
	{
		Pselect_actually[visit1[allvisit1 - 1]][no_visit_city[i]] = Pselect1[visit1[allvisit1 - 1]][no_visit_city[i]] / SUM;

	}
	int Novisit = rows - allvisit1;//There is no number of cities visited
	vector<int>numbervisit(Novisit, 0);
	//Bubble sort
	for (int j = 0;j < Novisit;j++)
	{
		numbervisit[j] = no_visit_city[j];
	}
	for (int j = 0;j < Novisit - 1;j++)
	{
		for (int e = j + 1;e < Novisit;e++)
		{
			if (Pselect_actually[visit1[allvisit1 - 1]][no_visit_city[j]] < Pselect_actually[visit1[allvisit1 - 1]][no_visit_city[e]])
			{
				swap(Pselect_actually[visit1[allvisit1 - 1]][no_visit_city[j]], Pselect_actually[visit1[allvisit1 - 1]][no_visit_city[e]]);
				swap(numbervisit[j], numbervisit[e]);
			}

		}
	}
	//Select the next city to visit, including both random methods and roulette methods
	int selectvisit_next;
	double ppss = static_cast<double>(std::rand()) / RAND_MAX;
	//cout << "ppss:" << ppss << endl;
	if (ppss >= 0.1)
	{
		//The weight assignment for roulette is Q
		if (Novisit >= 6)
		{
			vector<int> Q(6, 0);
			for (int j = 0;j < 6;j++)
			{
				Q[j] = 6 - j;
			}
			int proulette = rand() % 21;
			int temprou = 0;//The process parameters of roulette
			for (int j = 0;j <6;j++)
			{
				temprou += Q[j];
				if (proulette <= temprou)
				{
					selectvisit_next = numbervisit[j];
					break;
				}
			}
		}
		else
		{
			selectvisit_next = numbervisit[0];
		}
	}
	else
	{
		int aass = rand() % Novisit;
		selectvisit_next = numbervisit[aass];
	}
	return selectvisit_next;

}

//The solution construction process of ants
vector<vector<int>> update_route_ante(const int ant1,const vector<vector<double>> phero1, const vector<vector<double>>city_distan)
{
	int rows = phero1.size();
	int cols = ant1;

	vector<vector<int>> visit(cols, vector<int>(rows, rows));//Mark whether to visit the city
	vector<vector<double>> Pselect(rows,vector<double>(rows, 0));

	//1.Calculate the probability of the ant's next city choice£»1.1 Calculate the molecules of the probability formula first
	for (int j = 0;j < rows - 1;j++)
	{
		for (int r = j + 1;r < rows;r++)
		{
			Pselect[j][r] = (pow(phero1[j][r], ALPHA)) * (pow((1/city_distan[j][r]), BETA));//The larger the pheromone and the smaller the distance, the greater the probability of being selected
			Pselect[r][j] = (pow(phero1[r][j], ALPHA)) * (pow((1/city_distan[r][j]), BETA));
			//cout << "Pselect: " << Pselect[j][r] << ", " << Pselect[r][j] << endl;
		}
	}
	for (int j = 0;j < rows;j++)
	{
		Pselect[j][j] = 0;
	}
	for (int i = 0;i < ant1;i++)
	{
		//1.2 Randomly initialize the first departure city of ants
		visit[i][0] = rand() % rows;
		vector<int> visit_i(rows, 0);
		int allvisit=1;//The number of cities visited

		for (int j = 0;j < rows-1;j++)
		{
			//cout << "visit: " << visit[i][j] << endl;
			for (int e = 0;e < rows;e++)
			{
				visit_i[e] = visit[i][e];
			}
			visit[i][j+1] = visit_city_update(visit_i, Pselect, allvisit);
			allvisit += 1;
			//cout << "visit: " << visit[i][j] << ", " << visit[i][j + 1] << endl;
		}
	}
	return visit;
}

//Save the optimal solution of this round and the global optimal solution
tuple<double, double, vector<int>, vector<int>> solution_best(double gbest1,const vector<int> &gbest1_route, const vector<vector<int>>& ant_route1, const vector<double>& ant_distan1)
{
	double result1 = 0.0;
	double result2 = 0.0;
	int row = ant_route1.size();
	int cols = ant_route1[0].size();
	vector<int> result3(cols, 0);
	vector<int> result4(cols, 0);

	vector<int> nums;
	vector<double>target(row,0);
	for (int i = 0;i < row;i++)
	{
		nums.push_back(i);
	}
	for (int i = 0;i < row;i++)
	{
		target[i] = ant_distan1[i];
	}


	for (int i = 0;i < row;i++)
	{
		for (int j = i + 1;j < row-1;j++)
		{
			if (target[i] >= target[j])
			{
				swap(target[i], target[j]);
				swap(nums[i], nums[j]);
			}
		}
	}
	
	result1 = target[0];
	for (int j = 0;j < cols;j++)
	{
		result3[j] = ant_route1[nums[0]][j];
	}
	
	if (gbest1 > result1)
	{
		result2 = result1;
		for (int j = 0;j < cols;j++)
		{
			result4[j] = result3[j];
		}
	}
	else
	{
		result2 = gbest1;
		result4 = gbest1_route;
	}


	return { result1, result2, result3, result4 };


}

//Randomly initialize the ant sequence
vector<vector<int>> random_update_route(const vector<vector<int>> temp_ant_route,const vector<double> ant_distan1, const vector<vector<double>> city_distan1)
{
	//Because the optimal solution has not been updated for a long time, the ant is reinitialized
	int rows = temp_ant_route.size();
	int cols = temp_ant_route[0].size();
	vector<vector<int>> result(rows, vector<int>(cols, 1));
	for (int c = 0;c < rows;c++)
	{
		vector<int> nums;
		for (int i = 0;i < cols;i++)
		{
			nums.push_back(i);
		}
		for (int i = 0;i < cols;i++)
		{
			int j = rand() % (cols);
			swap(nums[i], nums[j]);
		}
		for (int j = 0;j < cols;j++)
		{
			result[c][j] = nums[j];
		}
	}

	//Sort the current ant paths and replace the worst ants
	vector<int> numer_ant(rows, 0);
	for (int i = 0;i < rows;i++)
	{
		numer_ant[i] = i;
	}
	vector<double> temp_ant_distan1(rows, 0);
	temp_ant_distan1 = ant_distan1;
	for (int i = 0;i < rows-1;i++)
	{
		for (int j = i + 1;j < rows;j++)
		{
			if (temp_ant_distan1[i] < temp_ant_distan1[j])
			{
				swap(temp_ant_distan1[i], temp_ant_distan1[j]);
				swap(numer_ant[i], numer_ant[j]);
			}
		}
	}
	for (int i = 40;i < rows;i++)
	{
		for (int j = 0;j < cols;j++)
		{
			result[i][j] = temp_ant_route[numer_ant[i]][j];
		}
	}
	return result;
}


//Neighborhood search near the optimal solution
vector<vector<int>> local_search_ant(const vector<vector<int>>temp_ant_route, const vector<vector<double>>city_distan1, const vector<double> ant_distan1)
{
	int rows = temp_ant_route.size();
	int cols= temp_ant_route[0].size();
	vector<vector<int>> result1(rows, vector<int>(cols, 0));
	result1 = temp_ant_route;
	vector<int> result(cols, 0);
	//1.1Find the longest adjacent distance of the optimal solution.
	vector<int> numer_ant(rows, 0);
	for (int i = 0;i < rows;i++)
	{
		numer_ant[i] = i;
	}
	vector<double> temp_ant_distan1(rows, 0);
	temp_ant_distan1 = ant_distan1;
	for (int i = 0;i < rows - 1;i++)
	{
		for (int j = i + 1;j < rows;j++)
		{
			if (temp_ant_distan1[i] > temp_ant_distan1[j])
			{
				swap(temp_ant_distan1[i], temp_ant_distan1[j]);
				swap(numer_ant[i], numer_ant[j]);
			}
		}
	}

	//Best:temp_ant_distan1[0],Best route: temp_ant_route[numer_ant[0]]
	vector<double> distance_in_route_city(cols, 0);
	for (int i = 0;i < cols-1;i++)
	{
		int temp1 = temp_ant_route[numer_ant[0]][i];
		int temp2 = temp_ant_route[numer_ant[0]][i + 1];
		distance_in_route_city[i] = city_distan1[temp1][temp2];
	}
	//Find the maximum value in distance_in_route_city
	vector<int>nnmm(rows, 0);
	for (int i = 0;i < rows;i++)
	{
		nnmm[i] = i;
	}
	for (int i = 0;i < rows - 1;i++)
	{
		for (int j = i + 1;j < rows;j++)
		{
			if (distance_in_route_city[i] < distance_in_route_city[j])
			{
				swap(distance_in_route_city[i], distance_in_route_city[j]);
				swap(nnmm[i], nnmm[j]);
			}
		}
	}


	//distance_in_route_city[0] is the maximum distance, and nnmm[0] corresponds to the starting city
	//Find the nearest city to it. If it is not its previous starting point, move to the next arrival point. For other sequences, move backward in sequence
	int temp1 = temp_ant_route[numer_ant[0]][nnmm[0]];
	double temp2=10000;
	int temp3;
	for (int i = 0;i < cols;i++)
	{
		if (i != temp1&& city_distan1[temp1][i]<temp2)
		{
			temp2 = city_distan1[temp1][i];
			temp3 = i;
		}
	}
	//If the previous city of this city was not temp3 and this city is not the last one
	if (nnmm[0] + 1!=cols&&nnmm[0]!=0&&temp_ant_route[numer_ant[0]][nnmm[0] - 1] != temp3)
	{
		int k = 0;
		for (int i = 0;i < nnmm[0]+1;i++)
		{
			if (i != temp3)
			{
				result[k] = temp_ant_route[numer_ant[0]][i];
				k++;
			}
		}
		result[k] = temp3;
		for (int i = nnmm[0] + 2;i < cols;i++)
		{
			if (i != temp3)
			{
				k++;
				result[k] = temp_ant_route[numer_ant[0]][i];
			}
		}
		//Recalculate the changes before and after the local search.
		double distance_local = 0.0;
		for (int j = 0;j < cols - 1;j++)
		{
			distance_local += city_distan1[result[j]][result[j+1]];
		}
		distance_local += city_distan1[result[cols-1]][result[0]];
		//cout << "distance_local: " << temp_ant_distan1[0] << ", " << distance_local << endl;
		if (distance_local < temp_ant_distan1[0])
		{
			for (int i = 0;i < cols;i++)
			{
				result1[numer_ant[0]][i] = result[i];
			}
		}
	}

	return result1;
}

int main()
{
	int maxiter = 1000;//Maximum number of iterations
	srand(static_cast<unsigned int>(time(NULL)));
	//1. Read data
	int city = 0;
	string filePath = "XX/XX";
	vector<string>files;
	string format = ".txt";
	GetAllFormatFiles(filePath, files, format);
	int size = files.size();
	string sstr= "XX/XX/XX.txt";
	for (int s = 0;s < size;s++)
	{
		ifstream fileProcess;
		fileProcess.open(files[s], ios::in);
		if (!fileProcess.is_open())
		{
			cout << "Failed to read the file" << endl;
		}
		fileProcess >> city;
		int *city_x;
		int* city_y;
		city_x = new int[city];
		city_y = new int[city];
		int a;
		for (int i = 0;i < city;i++)
		{
			fileProcess >> a;
			fileProcess >> city_x[i];
			fileProcess >> city_y[i];
		}
		vector<double> pbest(maxiter,5000);//The optimal solution of this round 
		vector<double> gbest(maxiter, 5000);//Global optimal solution
		vector<int> pbest_route(city, 0);//The corresponding route of the optimal solution in this round
		vector<int> gbest_route(city, 0);//The corresponding route of the global optimal solution


		//2.Initialize the algorithm parameters
		int ant = city;

		//Initialize the pheromone matrix: between cities
		vector<vector<double>> phero(city, vector<double>(city, 1.0));
		//Initialize the distance matrix
		vector<vector<double>> city_distan(city, vector<double>(city, 0.0));
		for (int i = 0;i < city;i++)
		{
			for (int j = 0;j < city;j++)
			{
				city_distan[i][j] = sqrt(pow(city_x[i] - city_x[j], 2) + pow(city_y[i] - city_y[j], 2));
			}
		}


		//Initialize the route that each ant needs to take
		vector<vector<int>> ant_route(ant, vector<int>(city, 1));
		for (int c = 0;c < ant;c++)
		{
			vector<int> nums;
			for (int i = 0;i < city;i++)
			{
				nums.push_back(i);
			}
			for (int i = 0;i < city;i++)
			{
				int j = rand() % (city);
				swap(nums[i], nums[j]);
			}
			for (int j = 0;j < city;j++)
			{
				ant_route[c][j] = nums[j];
			}
		}
		//Initialization: Select the city that is currently closest to it for each city. Start from City 0
		ant_route[0][0] = 0;
		int minIndex = 0;

		for (int i = 0;i < city;i++)
		{
			if (i != 0 && city_distan[0][i] < city_distan[0][minIndex])
			{
				minIndex = i;
			}
		}
		ant_route[0][1] = minIndex;

		int visit = 0;
		for (int i = 0;i < city - 1;i++)
		{
			minIndex = i + 1;
			for (int j = 0;j < city;j++)
			{
				for (int k = 0;k < visit;k++)
				{
					if (j != ant_route[0][k])
					{
						if (city_distan[i][j] < city_distan[i][minIndex])
						{
							minIndex = i;
						}
					}
				}
			}
			ant_route[0][i + 1] = minIndex;
			visit += 1;
		}
		//Calculate the route length of each ant
		vector<double> ant_distan(ant, 0);
		ant_distan=distance_cout_ant(ant_route, city_distan);



		double gbest1=5000;
		vector<vector<int>> temp_ant_route(ant, vector<int>(city, 1));
		vector<int> gbest1_route(city, 0);
		for (int epoch = 0;epoch < maxiter;epoch++)
		{
			//Save the optimal solution of this round and the global optimal solution
			gbest1_route = gbest_route;
			tie(pbest[epoch], gbest[epoch], pbest_route, gbest_route) = solution_best(gbest1, gbest1_route, ant_route, ant_distan);
			//cout << "pbest: " << pbest << ",gbest: " << gbest << endl;
			cout << "gbest " << epoch<<": "<<gbest[epoch] << endl;
			gbest1 = gbest[epoch];

			//Update pheromone
			phero = update_pheromones(phero, ant_route, ant_distan);
			//cout <<"phero: "<< phero[0][0] << phero[1][0] << endl;

			//The solution construction of each ant
			ant_route = update_route_ante(ant,phero, city_distan);

			//Calculate the route of the ant
			ant_distan = distance_cout_ant(ant_route, city_distan);


			//Neighborhood optimization of the optimal solution
			temp_ant_route = ant_route;
			ant_route = local_search_ant(temp_ant_route, city_distan, ant_distan);

			
			//If there is no improvement in the optimal solution after 100 rounds, reinitialize some of the solutions
			if (epoch > 100 && gbest[epoch] == gbest[epoch - 100]) 
			{
				temp_ant_route = ant_route;
				ant_route = random_update_route(temp_ant_route, ant_distan, city_distan);
			}
			//Calculate the route of the ant
			ant_distan = distance_cout_ant(ant_route, city_distan);
		}

		cout << "finalist best: " << gbest[maxiter-1] << endl;
		for (int i = 0;i < city;i++)
		{
			cout << gbest_route[i] << " ";
		}
		cout << endl;
	}
	return 0;
}