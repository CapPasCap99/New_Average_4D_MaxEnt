#ifndef INIT_H
#define INIT_H

#include "global.h"
#include "moves.h"
using namespace Eigen;

std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations(number_of_threads);

std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_inter(number_of_threads);
std::vector<std::unordered_map<std::pair<int, int>, int, pair_hash>> contacts_lin(number_of_threads);
std::vector<std::unordered_map<std::vector<int>, std::vector<int>, vec_hash>> locations_lin(number_of_threads);


void initialize(int thread_num, int step) {
//    int stage = thread_num % number_of_stages;
    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    // Set the values of is_replicated//
    for (int i = 0; i < pol_length; i++) {
        if(lin_length[thread_num]==0){
            is_replicated[thread_num][i] = false;
            continue;
        }
        if (oriC + lin_length[stage] >= pol_length) {
            if (i < oriC - lin_length[stage] && i > (oriC + lin_length[stage]) % pol_length) {
                is_replicated[thread_num][i] = false;
            }
        }
        else {
            if (i > oriC + lin_length[stage]) {
                is_replicated[thread_num][i] = false;
            }
            else if (i < oriC - lin_length[stage]) {
                is_replicated[thread_num][i] = false;
            }
        }
    }

    // Construct initial configurations if these are not taken from saved configurations//
    if(!initConfig) {
        int xi{0}, yi{0}, zi{0};
        Vector3i monomer;
        for (int i = 0; i < pol_length; i++) {
            if (i % 2 == 0) {
                zi -= 1;
            } else {
                zi += 1;
            }
            monomer = {xi, yi, zi};
            polymer[thread_num].push_back(monomer);
            lin_polymer[thread_num].push_back(monomer);

            //JM: The values of the linear polymer are set to zero if they are past the replication fork
            if (oriC + lin_length[stage] >= pol_length) {
                if (i < oriC - lin_length[stage] && i > (oriC + lin_length[stage]) % pol_length) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                }
            } else {
                if (i > oriC + lin_length[stage]) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                } else if (i < oriC - lin_length[stage]) {
                    lin_polymer[thread_num][i] = {0, 0, 0};
                }
            }
        }
    }

    // Construct initial configurations from saved configurations//
    if (initConfig) {
        check_input_configuration_compatibility(configuration_data_folder);
        read_configuration("ring",  thread_num, step);
        read_configuration("strand", thread_num, step);
    }


    // set contacts and locations //
    for (int i = 0; i < pol_length; i++) {
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(i);
    }
    if (lin_length[stage] != 0) {
        if (lin_length[stage] < pol_length/2) {
            if (oriC + lin_length[stage] >= pol_length) {
                for (int i = 0; i < pol_length; i++) {
                    if (i > oriC - lin_length[stage] || i < (oriC + lin_length[stage]) % pol_length) {
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                            }
                        }
                        if (locations[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_inter[thread_num][{elem, i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}].push_back(i);
                    }
                }
            } else {
                for (int i = 0; i < pol_length; i++) {
                    if (i > oriC - lin_length[stage] && i < oriC + lin_length[stage]) {
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                            }
                        }
                        if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                                contacts_inter[thread_num][{elem, i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(i);
                    }
                }
            }
        }
        else {
            for (int i = 0; i < pol_length; i++) {
                if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                    for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                        contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                    }
                }
                if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                    for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                        contacts_inter[thread_num][{elem, i}] = 0;
                    }
                }
                locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(i);
            }
        }
    }
}


void burn_in(int thread_num, int n_steps) {
    for (int m = 0; m < n_steps ; m++) {   //burn-in
        move(thread_num, m);
    }
//    int stage = thread_num % number_of_stages;
    int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)

    std::vector<double> zeroVec(pol_length, 0);
    std::fill(total_contacts[thread_num].begin(), total_contacts[thread_num].end(), zeroVec);
    z_close[thread_num] = 0;
    z_far[thread_num] = 0;
    z_close_squared[thread_num] = 0;
    z_far_squared[thread_num] = 0;
    for(int i=0; i<sites_constrained_mean.size(); i++){
        z_mean_data[thread_num][i] = 0;
    }
    for(int i=0; i<sites_constrained_separation.size(); i++){
        z_separation_data[thread_num][i] = 0;
    }
    //z[thread_num] = 0;
    //z_lin[thread_num] = 0;
    locations[thread_num].clear();
    locations_lin[thread_num].clear();
    contacts[thread_num].clear();
    contacts_lin[thread_num].clear();
    contacts_inter[thread_num].clear();
    for (int i = 0; i < pol_length; i++) {
        if (locations[thread_num].find({ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2] }]) {
                contacts[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
            }
        }
        locations[thread_num][{polymer[thread_num][i][0],polymer[thread_num][i][1],polymer[thread_num][i][2]}].push_back(i);
    }
    if (lin_length[stage] != 0) {
        if (lin_length[stage] < pol_length/2) {
            if (oriC + lin_length[stage] >= pol_length) {
                for (int i = 0; i < pol_length; i++) {
                    if (i > oriC - lin_length[stage] || i < (oriC + lin_length[stage]) % pol_length) {
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                            }
                        }
                        if (locations[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}]) {
                                contacts_inter[thread_num][{elem, i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1], lin_polymer[thread_num][i][2]}].push_back(i);
                    }
                }
            } else {
                for (int i = 0; i < pol_length; i++) {
                    if (i > oriC - lin_length[stage] && i < oriC + lin_length[stage]) {
                        if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                            for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                                contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                            }
                        }
                        if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                            for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                                contacts_inter[thread_num][{elem, i}] = 0;
                            }
                        }
                        locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(i);
                    }
                }
            }
        }
        else {
            for (int i = 0; i < pol_length; i++) {
                if (locations_lin[thread_num].find({lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}) !=locations_lin[thread_num].end()) {
                    for (auto elem : locations_lin[thread_num][{lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}]) {
                        contacts_lin[thread_num][{std::min(elem, i), std::max(elem, i)}] = 0;
                    }
                }
                if (locations[thread_num].find({ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }) != locations[thread_num].end()) {
                    for (auto elem : locations[thread_num][{ lin_polymer[thread_num][i][0],lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2] }]) {
                        contacts_inter[thread_num][{elem, i}] = 0;
                    }
                }
                locations_lin[thread_num][{lin_polymer[thread_num][i][0], lin_polymer[thread_num][i][1],lin_polymer[thread_num][i][2]}].push_back(i);
            }
        }
    }
}

#endif
