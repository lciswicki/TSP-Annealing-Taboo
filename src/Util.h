#pragma once
#include <vector>
#include <string>
#include <chrono>
#include <random>


void PrintVectorInt(std::vector<int> v);
void PrintVectorDouble(std::vector<double> v);
void Print2DVectorInt(std::vector<std::vector<int>> v);
void Print2DVectorDouble(std::vector<std::vector<double>> v);

/**
 * @brief Splits input string using delimiters
 * Util function for Parser
 * @param s string to be split
 * @return std::vector<int> strings after splitting
 */
std::vector<int> String_Split(std::string s);
/**
 * @brief Parses a test file into a coordinates vector
 * @param file_name path to the test file
 * @return std::vector<std::vector<int>> coordinates vector
 */
std::vector<std::vector<int>> Parser(std::string file_name);

/**
 * @brief Calculates euclidean distance between vertices
 * and presents it as weights of edges
 * O((coords.size()) ^ 2)
 * @param coords coordinates of vertices
 * @return Distances between vertices as weights
 */
std::vector<std::vector<int>> Weights(std::vector<std::vector<int>> coords);

/**
 * @brief Generates a random permutation of a vector
 * @param n length of the vector
 * @return std::vector<int> a permutated vector
 */
std::vector<int> Perm_Gen(int n);

void To_File(std::string path, std::vector<int> cycle);