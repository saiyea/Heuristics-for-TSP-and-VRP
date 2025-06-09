#pragma once
#include"message.h"
#include<iostream>
#include<ctime>
#include<vector>
#include<cmath>
#include<fstream>
#include<cstdlib>


using namespace std;
struct ppChrom//Chromosome
{
	int* JSequence;//The sequence string of the docking station
	int* MSequence;//Whether the vehicle uses a sequence
	double tFitness;
};

class VRP_GA
{
public:
	//friend class ppSolution;
	VRP_GA() :x(0), y(0) {                          //Default constructor
		cout << "Default Constructor called." << endl;
	}
	int x, y;
	VRP_GA(int iPS, int iMaxcycle, float iPC, float iPM, float iP);
	~VRP_GA();

	//variable
	void VRP_GAPP();
	int gen;
	double aveBest;
	double ARPD;
	double STD;
	double finalBest;
	int PS;//Population number
	int Maxcycle;//Maximum number of iterations
	float PC;//Cross probability
	float PM;//Mutation probability
	float P;//Select proportion
	int* valuecycle;

	ppChrom* Chrom;//°üº¬¸öÌå
	ppChrom* newPop;
	ppChrom* oldPop;
	ppChrom* pTemp;
	ppChrom minChrom;
	ppChrom bestChrom;
	int bestFitness;        
	float avgFitness;
	int sumFitness;

	double Tmax;//The time limit of the algorithm

	//function
	void run_times(double t);

	void Initialize(); //Individual initialization
	ppChrom* sorting(ppChrom* parent);
	void Generationpp();    

	//Calculate the fitness value
	double newdecode(int* coSequence, int* cmSequence);                    
	double newbestdecode(int* coSequence, int* cmSequence);              

	//Replicate chromosome
	void Copy(ppChrom& Parent, ppChrom& Child);

	//Competitive choice
	int TourSelection();      
	//Roulette choice
	int JRouttleSelection();

	//Cross-mutation operation
	void poxCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2);
	void joxCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2);
	void JCross(ppChrom& parent1, ppChrom& parent2, ppChrom& child1, ppChrom& child2);
	void swapMutation(ppChrom& parent, ppChrom& child);
	void JMutation(ppChrom& parent, ppChrom& child);

	//Used to generate crossover, mutation probabilities or selection probabilities
	bool flip(float probability);   

	ofstream dataFile;





};


