#include "TSP.h"
#include "Util.h"
#include <vector>
#include <chrono>
#include <random>
#include <iostream>


int Objective_Function(std::vector<int> perm, std::vector<std::vector<int>> weights) {
    int p_length = perm.size(); // permutation length
    int c_weight = 0;   // cycle weight
    for (int i = 0; i < p_length - 1; i++) {
        c_weight += weights[perm[i]][perm[i + 1]];
    }
    // adding an edge from last to first to close the cycle
    c_weight += weights[perm[p_length - 1]][perm[0]];
    return c_weight;
}

std::vector<int> invert(std::vector<int> perm, int i, int j) {
    int p_size = perm.size();
    if (i < 0 || i >= p_size) {
        std::cout << "Error: " << i << " out of bound" << std::endl;
        return perm;
    }
    if (j < 0 || j >= p_size) {
        std::cout << "Error: " << j << " out of bound" << std::endl;
        return perm;
    }
    if (j < i) {
        std::cout << "Error: " << j << " smaller than " << i << std::endl;
        return perm;
    }
    int temp = 0;
    // division with an upper ceiling
    int half_range = ((j - i) + 2 - 1) / 2;
    for (int x = 0; x < half_range; x++) {
        temp = perm[i + x];
        perm[i + x] = perm[j - x];
        perm[j - x] = temp;
    }
    return perm;
}

double NewLength(std::vector<int> perm_old, double length_old, int i, int j, std::vector<std::vector<int>> weights) {
    int p_size = perm_old.size();
    double length_new = length_old;
    if (i == 0 && j == p_size - 1) {
        // cycle length doesnt change
    }
    else if (i == 0) {
        length_new -= weights[perm_old[p_size - 1]][perm_old[i]];
        length_new -= weights[perm_old[j]][perm_old[j + 1]];
        length_new += weights[perm_old[p_size - 1]][perm_old[j]];
        length_new += weights[perm_old[i]][perm_old[j + 1]];
    }
    else if (j == p_size - 1) {
        length_new -= weights[perm_old[i - 1]][perm_old[i]];
        length_new -= weights[perm_old[j]][perm_old[0]];
        length_new += weights[perm_old[i - 1]][perm_old[j]];
        length_new += weights[perm_old[i]][perm_old[0]];
    }
    else {
        length_new -= weights[perm_old[i - 1]][perm_old[i]];
        length_new -= weights[perm_old[j]][perm_old[j + 1]];
        length_new += weights[perm_old[i - 1]][perm_old[j]];
        length_new += weights[perm_old[i]][perm_old[j + 1]];
    }
    return length_new;
}

Annealing::Annealing(std::string input_path) {
    coords = Parser(input_path);
    weights = Weights(coords);
    size = coords.size();
}

Annealing::~Annealing() {

}

std::vector<int> Annealing::RndNeighbour(std::vector<int> previous_perm) {
    int p_size = previous_perm.size();
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<std::mt19937::result_type> distrib_i(0, p_size - 2);
    int i = distrib_i(rng);
    std::uniform_int_distribution<std::mt19937::result_type> distrib_j(i, p_size - 1);
    int j = distrib_j(rng);
    std::vector<int> inverted_perm = invert(previous_perm, i, j);
    return inverted_perm;
}

std::vector<int> Annealing::BestRndNeighbour(std::vector<int> perm, int neighb_count) {
    int length_start = Objective_Function(perm, weights);
    int size_p = perm.size();
    int length_new = 0;
    int length_best = length_start;
    int i_best = 0;
    int j_best = 1;

    std::vector<int> temp_v(size_p);
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    int i = 0;
    int j = 0;
    for (int x = 0; x < neighb_count; x++) {
        std::uniform_int_distribution<std::mt19937::result_type> distrib_i(0, size_p - 2);
        i = distrib_i(rng);
        std::uniform_int_distribution<std::mt19937::result_type> distrib_j(i, size_p - 1);
        j = distrib_j(rng);

        length_new = length_start;
        if (i == 0 && j == size_p - 1) {
            // cycle length doesnt change
        }
        else if (i == 0) {
            length_new -= weights[perm[size_p - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[j + 1]];
            length_new += weights[perm[size_p - 1]][perm[j]];
            length_new += weights[perm[i]][perm[j + 1]];
        }
        else if (j == size_p - 1) {
            length_new -= weights[perm[i - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[0]];
            length_new += weights[perm[i - 1]][perm[j]];
            length_new += weights[perm[i]][perm[0]];
        }
        else {
            length_new -= weights[perm[i - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[j + 1]];
            length_new += weights[perm[i - 1]][perm[j]];
            length_new += weights[perm[i]][perm[j + 1]];
        }

        if (length_new < length_best) {
            length_best = length_new;
            i_best = i;
            j_best = j;
        }
    }

    if (length_best >= length_start) {
        return perm;
    }

    else {
        return invert(perm, i_best, j_best);
    }
}


double Annealing::ProbAnn(int new_value, int old_value, double temp) {
    return exp((old_value - new_value) / temp);
}

int Annealing::CalcAnnealing(double start_temp, int max_epoch_count,
    int max_epoch_steps, double temp_multip, int neighb_count) {
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    std::vector<int> current_perm = Perm_Gen(coords.size());
    std::vector<int> new_perm;
    int current_value = Objective_Function(current_perm, weights);
    int new_value;
    int epoch_counter = 0;
    int steps_in_epoch = 0;
    double temp = start_temp;
    double prob;
    double rnd_number;
    int p_size = current_perm.size();

    iterations = 0;

    std::cout << "Starting perm length = " << current_value << std::endl;
    std::cout << "Epoch Length" << std::endl;
    while (epoch_counter < max_epoch_count) {
        steps_in_epoch = 0;
        std::cout << epoch_counter << " " << current_value << std::endl;
        while (steps_in_epoch < max_epoch_steps) {
            //std::uniform_int_distribution<std::mt19937::result_type> distrib_i(0, p_size - 2);
            //int i = distrib_i(rng);
            //std::uniform_int_distribution<std::mt19937::result_type> distrib_j(i, p_size - 1);
            //int j = distrib_j(rng);
            new_perm = BestRndNeighbour(current_perm, neighb_count);
            new_value = Objective_Function(new_perm, weights);

            if (new_value < current_value) {
                current_value = new_value;
                current_perm = new_perm;
            }
            else {
                prob = ProbAnn(new_value, current_value, temp);
                rnd_number = distrib(rng);
                if (rnd_number < prob) {
                    current_value = new_value;
                    current_perm = new_perm;
                }
            }
            steps_in_epoch++;
            iterations++;
        }
        temp *= temp_multip;
        epoch_counter++;
    }


    std::cout << "Final perm length = " << current_value << std::endl;
    return current_value;
}


TabuSearch::TabuSearch(std::string input_path) {
    coords = Parser(input_path);
    weights = Weights(coords);
    size = coords.size();
}

TabuSearch::~TabuSearch() {

}

std::vector<int> TabuSearch::BestRndNeighbour(std::vector<int> perm, int neighb_count) {
    int length_start = Objective_Function(perm, weights);
    int size_p = perm.size();
    int length_new = 0;
    int length_best = length_start;
    int i_best = 0;
    int j_best = 1;

    std::vector<int> temp_v(size_p);
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    int i = 0;
    int j = 0;
    for (int x = 0; x < neighb_count; x++) {
        std::uniform_int_distribution<std::mt19937::result_type> distrib_i(0, size_p - 2);
        i = distrib_i(rng);
        std::uniform_int_distribution<std::mt19937::result_type> distrib_j(i, size_p - 1);
        j = distrib_j(rng);

        length_new = length_start;
        if (i == 0 && j == size_p - 1) {
            // cycle length doesnt change
        }
        else if (i == 0) {
            length_new -= weights[perm[size_p - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[j + 1]];
            length_new += weights[perm[size_p - 1]][perm[j]];
            length_new += weights[perm[i]][perm[j + 1]];
        }
        else if (j == size_p - 1) {
            length_new -= weights[perm[i - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[0]];
            length_new += weights[perm[i - 1]][perm[j]];
            length_new += weights[perm[i]][perm[0]];
        }
        else {
            length_new -= weights[perm[i - 1]][perm[i]];
            length_new -= weights[perm[j]][perm[j + 1]];
            length_new += weights[perm[i - 1]][perm[j]];
            length_new += weights[perm[i]][perm[j + 1]];
        }

        if (length_new < length_best) {
            length_best = length_new;
            i_best = i;
            j_best = j;
        }
    }

    if (length_best >= length_start) {
        return perm;
    }

    else {
        return invert(perm, i_best, j_best);
    }
}

int TabuSearch::CalcTabu(int list_length, int max_stale_iter, int neigbh_count) {
    std::vector<std::vector<int>> tabu_list;
    std::vector<int> current_perm = Perm_Gen(size);
    std::vector<int> new_perm;
    int current_length = Objective_Function(current_perm, weights);
    int new_length;
    int stale_iter = 0;

    iterations = 0;

    std::cout << "Starting perm length = " << current_length << std::endl;
    std::cout << "List_size Length" << std::endl;

    while (tabu_list.size() < list_length && stale_iter < max_stale_iter) {
        new_perm = BestRndNeighbour(current_perm, neigbh_count);
        new_length = Objective_Function(new_perm, weights);
        if (new_length < current_length &&
            (std::find(tabu_list.begin(), tabu_list.end(), new_perm) == tabu_list.end())){
            current_perm = new_perm;
            current_length = new_length;
            tabu_list.push_back(current_perm);
            stale_iter = 0;
        }
        if (tabu_list.size() % 100 == 0) {
            std::cout << " " << tabu_list.size() << " " << current_length << std::endl;
        }
        
        stale_iter++;
    }
    std::cout << "Final perm length = " << current_length << std::endl;
    return current_length;
}

