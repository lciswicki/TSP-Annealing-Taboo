#pragma once
#include "Util.h"
#include <vector>
#include <string>


/**
 * @brief Calculates objective function for a permutation
 * O(n)
 * @param perm starting permutation
 * @return Objective value of the permutation
 */
int Objective_Function(std::vector<int> perm, std::vector<std::vector<int>> weights);

/**
 * @brief Inverts a part of a vector
 * O(n)
 * @param perm starting permutation
 * @param i start of inversion
 * @param j end of inversion
 * @return std::vector<int> inverted vector
 */
std::vector<int> invert(std::vector<int> perm, int i, int j);

double NewLength(std::vector<int> perm_old, double length_old, int i, int j, std::vector<std::vector<int>> weights);

class Annealing {
private:
	int size;
	std::vector<std::vector<int>> coords;
	std::vector<std::vector<int>> weights;
	std::vector<int> RndNeighbour(std::vector<int>);
	double ProbAnn(int new_value, int old_value, double temp);
	std::vector<int> BestRndNeighbour(std::vector<int> perm, int neighb_count);

public:
	int iterations = 0;
	Annealing(std::string data_path);
	~Annealing();
	int CalcAnnealing(double start_temp, int max_epoch_count,
		int max_epoch_steps, double temp_multip, int neighb_count);

};

class TabuSearch {
private:
	int size;
	std::vector<std::vector<int>> coords;
	std::vector<std::vector<int>> weights;
	std::vector<int> BestRndNeighbour(std::vector<int> perm, int neighb_count);
public:
	int iterations = 0;
	TabuSearch(std::string input_path);
	~TabuSearch();
	int CalcTabu(int list_length, int max_stale_iter, int neigbh_count);
};