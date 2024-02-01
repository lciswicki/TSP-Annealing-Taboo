#include "Util.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <chrono>


void PrintVectorInt(std::vector<int> v) {
	for (int i = 0; i < v.size(); i++) {
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

void PrintVectorDouble(std::vector<double> v) {
	for (int i = 0; i < v.size(); i++) {
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

void Print2DVectorInt(std::vector<std::vector<int>> v) {
	for (int i = 0; i < v.size(); i++) {
		for (int j = 0; j < v[0].size(); j++) {
			std::cout << v[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void Print2DVectorDouble(std::vector<std::vector<double>> v) {
	for (int i = 0; i < v.size(); i++) {
		for (int j = 0; j < v[0].size(); j++) {
			std::cout << v[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

std::vector<int> String_Split(std::string s) {
    std::vector<int> words(3);
    std::string delimiter = " ";
    int pos = 0;
    for (int i = 0; i < 3; i++) {
        pos = s.find(delimiter);
        words[i] = stoi(s.substr(0, pos));
        s.erase(0, pos + delimiter.length());
    }
    return words;
}

std::vector<std::vector<int>> Parser(std::string file_name) {
    std::string line;
    std::string word;
    std::ifstream file_reader(file_name);
    int problem_size = 0;
    // iterate over entry lines
    for (int i = 0; i < 5; i++)
        getline(file_reader, line);
    // get problem size
    for (int i = 0; i < 3; i++)
        file_reader >> word;
    problem_size = stoi(word);
    // vector of coordinates of cities
    std::vector<std::vector<int>> coordinates(problem_size, std::vector<int>(2));
    // iterate over lines left
    for (int i = 0; i < 3; i++) {
        getline(file_reader, line);
    }
    // now read the data
    std::vector<int> split_line(3);
    for (int i = 0; i < problem_size; i++) {
        getline(file_reader, line);
        split_line = String_Split(line);
        coordinates[i][0] = split_line[1];
        coordinates[i][1] = split_line[2];
    }
    return coordinates;
}

std::vector<std::vector<int>> Weights(std::vector<std::vector<int>> coords) {
    int c_size = coords.size();
    std::vector<std::vector<int>> weights(c_size, std::vector<int>(c_size));
    for (int i = 0; i < c_size; i++) {
        for (int j = 0; j < c_size; j++) {
            weights[i][j] = sqrt(pow((coords[i][0] - coords[j][0]), 2) + pow((coords[i][1] - coords[j][1]), 2));
        }
    }
    return weights;
}

std::vector<int> Perm_Gen(int n) {
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::vector<int> permutation(n);
    for (int i = 0; i < n; i++)
        permutation[i] = i;
    shuffle(permutation.begin(), permutation.end(), rng);
    return permutation;
}

void To_File(std::string path, std::vector<int> cycle) {
    std::ofstream cycle_file(path);
    std::ostream_iterator<int> output_iterator(cycle_file, "\n");
    std::copy(std::begin(cycle), std::end(cycle), output_iterator);
    cycle_file.close();
}