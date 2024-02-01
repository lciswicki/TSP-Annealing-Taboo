#include <iostream>
#include "Util.h"
#include "TSP.h"
#include <string>
#include <vector>
#include <fstream>

int main() {

	int repetitions = 100;
	std::string input_path = "./data/icw1483.tsp";
	int n = 1083;


	TabuSearch tab(input_path);
	int list_length = n * n;	// n^2
	int max_stale = 100;			// n
	int neighbour_count = 10000;	// n

	std::vector<int> tabu_res(repetitions);

	std::ofstream output_file;
	output_file.open("./results/tab.txt");
	output_file << "i result\n";


	for (int i = 0; i < repetitions; i++) {
		tabu_res[i] = tab.CalcTabu(list_length, max_stale, neighbour_count);
		output_file << i << " " << tabu_res[i] << std::endl;
	}
	output_file.close();


	std::vector<int> ann_res(repetitions);
	output_file.open("./results/ann.txt");
	output_file << "i result\n";

	Annealing ann(input_path);
	double temp_mult = 0.90;	// 0.9
	double start_temp = n;	// n	
	int epoch = n;			// 10 * n
	int samples = n;			// n

	for (int i = 0; i < repetitions; i++) {
		ann_res[i] = ann.CalcAnnealing(start_temp, epoch, samples, temp_mult, 1000);
		output_file << i << " " << ann_res[i] << std::endl;
	}
	output_file.close();





}