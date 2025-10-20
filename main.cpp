//HALLOOO
// Code accompanying the paper 'Learning the dynamic organization of a replicating bacterial chromosome from time-course Hi-C data' by Harju et al. (2025)
// Here, the energies of the MaxEnt model are obtained via iterative Monte Carlo simulations.

// The code is structured as follows:

// global.h
// - Initiation of global variables

// files.h
// - Reading of input files, with checks on correct data structures. 
// - Writing of output files.

// Initialization.h
// - Define starting set of polymer configurations, either via reading in saved configurations, or by default construction.
// - Configure polymer replication stage.
// - Burn in polymer configurations.

// RandomGenerator.h
// - Definition of the RandomGenerator class, used for random number generation during Monte Carlo simulations.

// moves.h
// - All polymer moves used during the Monte Carlo simulation, with supporting functions.

// energy_changes.h
// - Computation of energy changes associated with a Monte Carlo move.

// main.cpp
// - Setting of simulation parameters
// - Orchestration of iterative Monte Carlo simulation

// CHANGES TO THE 4D CODE:
//  -Make lin_length and length number of stages long.
//  -Previously, number_of_threads and number_of_stages were used interchangeably, now we make sure we use one or the other and convert from one to the other using replicates_per_stage.
//  -Include number_of_threads>>number_of_threads_in_parallel
//  -Instead of having :
//         total_contacts(nb_threads,bin,bin)  --averaged over all threads-->  final_contacts(bin,bin)
//         change to:
//         total_contacts(nb_threads,bin,bin)  --averaged over replicates per stage-->  final_contacts(stage,bin,bin) --averaged over stages-->  final_contacts_averaged(bin,bin)
//  -Average over the energies/contacts and compare with averaged out Hi-C map to update energies
//  -Actually simulates all stages at once



#include <thread>
#include <chrono>
#include <iomanip>
#include "RandomGenerator.h"
#include "global.h"
#include "files.h"
#include "initialization.h"

////// Simulation properties //////
const int number_of_threads_in_parallel = 48; 
//const int number_of_threads = 2880; // Total number of simulations (must be divisible by number_of_stages)
const int number_of_threads = 96; // Total number of simulations (must be divisible by number_of_stages)
const std::vector<int> stages = {0, 10, 30, 45, 60, 75}; //C.Crescentus
const int number_of_stages = stages.size(); //number of replication stages
//const std::vector<int> lin_length = {0, 72, 292, 456, 624, 788}; //C.Crescentus number of replicated sites (physical ones)
//const std::vector<double> length = {18, 19, 22, 24, 28, 31}; //C.Crescentus cell length
const int replicates_per_stage = number_of_threads / number_of_stages;


const std::string bacteria_name = "crescentus_separations_3"; //Determines which input files will be used
const int bin_num = 405; //Crescentus
const int reduction_factor = 4; //Crescentus
//CB: It is weirdly hard coded in that the oriC is not 0, for other bacteria check if the logic makes sense with other oriC or just translate them when anaylising the results
const int oriC = 1220; // Crescentus   //position of the origin of replication. Determines where replication is initiated. 
const int Ter = 410;
double radius{ 3.2 }; // Crescentus


const int pol_length = reduction_factor*bin_num;
const int mc_moves = 4e6;
//const int mc_moves = 3.8e7;
const int burn_in_steps = 3e7;
//const int burn_in_steps = 3e3;
//const int it_steps = 60; //JM: number of steps inverse algorithm
const int it_steps =80;
bool boundary_cond = true;
bool constrain_pol = true; // constrains one origin to always be the closer one to a cell pole

///////// initial configurations ////////////
bool initConfig = false;

//int init_config_number = 60; //iteration number of the initial configuration used
int init_config_number = 15; //iteration number of the initial configuration used
const int res{1}; //JM: how often the mean positions are sampled

///////// initial energies //////////////
bool initEnerg = false; //load initial energies?
bool burn_in_done = false;
//std::string energy_data_folder = "home/Documents/test/4D/HiC_data/2025-08-24_1604_averaged/Energies";
std::string old_simulation_folder="2025-08-25_2249_averaged/";
std::string energy_data_iteration = std::to_string(init_config_number);

///////// input Hi-C data ////////////
std::string dir = "/home/capucine/Documents/test/4D/HiC_data/";
std::string HiC_file = dir + "t_averaged_normalized_rotated.txt";


// ////// Learning rates //////
// // Why so many learning rates?? //
// double learning_rate{ 0.5 }; // contacts updating (used) the 1st
// //double learning_rate{ 0.2 }; // contacts updating (used) the 2nd
// //double learning_rate{ 0.1 }; // contacts updating (used) the 3rd
// //double learning_rate{ 0.05 }; // contacts updating (used) the 4th
// //double learning_rate{ 0.02 }; // contacts updating (used) the 5th

// double learning_rate_close { 0.5 };
// double learning_rate_far { 0 };
// double learning_rate_close_var { 0 };
// double learning_rate_far_var { 0 };
// double learning_rate_means { 0 };
// double learning_rate_separations { 0.5 }; //for ori-ori distance the 1st(used)
// //double learning_rate_separations { 0.2 }; //for ori-ori distance the 2nd
// //double learning_rate_separations { 0.05 }; //for ori-ori distance the 3rd

double learning_rate = 0.1;
double learning_rate_close = 0.1;
double learning_rate_far = 0.1;
double learning_rate_close_var = 0.0;
double learning_rate_far_var = 0.0;
double learning_rate_means = 0.1;
double learning_rate_separations = 0.1;
double update_cap_factor {0.1};
int update_cap_onset {5};

/////// positional constraints ///////
// The only constrained sites are the Oris right? //
std::vector<int> sites_constrained_mean = {Ter}; //sites for which constraints on mean are imposed. mean of what ?
std::vector<int> sites_constrained_separation = {oriC}; //sites for which constraints on separation are imposed. what separation ?
const int n_constrained_mean = sites_constrained_mean.size();
const int n_constrained_separation = sites_constrained_separation.size();
std::unordered_map<int, int> sites_constrained_mean_map; // GG: keeps index of a site in the "sites_constrained_mean" list
std::unordered_map<int, int> sites_constrained_separation_map;
bool include_replicated_in_mean = false; //GG: include (z1+z2)/2 as contribution to mean z for replicated sites?

////// to know if we have data for a given site at a given stage (to avoid updating corresponding energies if we don't) ///////
std::vector<std::vector<bool>> is_constrained_ori(number_of_stages, std::vector<bool>(4, false));
std::vector<std::vector<bool>> is_constrained_mean(number_of_stages, std::vector<bool>(n_constrained_mean, false));
std::vector<std::vector<bool>> is_constrained_separation(number_of_stages, std::vector<bool>(n_constrained_separation, false));

/////// cellular properties ///////
std::vector<int> lin_length(number_of_stages, 0); //  progress of the replication forks at either of the two chromosome arms. If the value is above the bin number, replication is at 100% //JM: not clear what values larger than the number of bins imply, i.e. does it matter how much above the bin number it is?

std::vector<double> length(number_of_stages, 0); // JM: cell lengths. Ecoli: 8-21

std::vector<double> offset {0.5,0.5}; //JM: offsets in x and y directions, determines the center of the cylinder


/////// Set properties of randomly generated numbers ////////
std::vector<RandomGenerator> generators;

/////// variables to be used later ///////
//std::string output_folder;

std::vector<int> pole (number_of_threads,0);  // JM: stores the z value of the close cell pole, redundant in the current implementation because always the same

/////// Experimental constraints ////////
std::vector<double> xp_z_close(number_of_stages);// Experimental constraints [units of cell length]
std::vector<double> xp_z_far(number_of_stages);
std::vector<double> xp_z_close_var(number_of_stages);
std::vector<double> xp_z_far_var(number_of_stages);

std::vector<double> xp_z_close_simunits(number_of_stages); // Constraints in the units of simulation lattice (needed for calculating variance energy)
std::vector<double> xp_z_far_simunits(number_of_stages);

std::vector<std::vector<double>> target_means(number_of_stages,std::vector<double>(n_constrained_mean,0.)); //[stage][index of site in (sites_constrained_mean)]
std::vector<std::vector<double>> target_separations(number_of_stages,std::vector<double>(n_constrained_separation,0.)); //[stage][]

///////// Storing ori z-coordinate statistics  ///////////////
std::vector<double> z_close(number_of_threads,0);
std::vector<double> z_far(number_of_threads,0);
std::vector<double> z_close_squared(number_of_threads,0);
std::vector<double> z_far_squared(number_of_threads,0);
std::vector<double> z_close_var(number_of_threads,0);
std::vector<double> z_far_var(number_of_threads,0);


std::vector<double> z_close_tot(number_of_stages,0);
std::vector<double> z_far_tot(number_of_stages,0);
std::vector<double> z_close_squared_tot(number_of_stages,0);
std::vector<double> z_far_squared_tot(number_of_stages,0);

////////// Interaction, position, separation ENERGIES ///////////
std::vector<std::vector<double>> Interaction_E(pol_length, std::vector<double>(bin_num, 0));

std::vector<double> alpha(number_of_stages,0); // Energy close origin, for each stage of replication
std::vector<double> beta(number_of_stages,0); // Energy far origin
std::vector<double> alpha2(number_of_stages,0);// Energy close origin (variance), for each stage of replication
std::vector<double> beta2(number_of_stages,0);// Energy far origin (variance)

std::vector<std::vector<double>> energ_coeff_mean(number_of_stages, std::vector<double>(sites_constrained_mean.size(),0));// = {{}}; //[stage] equivalent of Lucas's alpha/beta, multiplies (z_close + z_far)
std::vector<std::vector<double>> energ_coeff_separation(number_of_stages, std::vector<double>(sites_constrained_separation.size(),0));// = {{}}; //[stage] separation energy - multiplies (|z_far - z_close|)
std::vector<std::unordered_map<int, double>> energy_mean_map(number_of_stages, std::unordered_map<int, double>()); //[stage] site_index -> energy coefficient
std::vector<std::unordered_map<int, double>> energy_separation_map(number_of_stages, std::unordered_map<int, double>());

/////////// Storing z-coordinate statistics of constrained sites /////////////////
std::vector<std::vector<double>> z_mean_data(number_of_threads, std::vector<double>(sites_constrained_mean.size(),0));
std::vector<std::vector<double>> z_separation_data(number_of_threads, std::vector<double>(sites_constrained_separation.size(),0));
std::vector<std::vector<double>> z_mean_data_tot(number_of_stages, std::vector<double>(sites_constrained_mean.size(),0));
std::vector<std::vector<double>> z_separation_data_tot(number_of_stages, std::vector<double>(sites_constrained_separation.size(),0));
std::vector<std::vector<unsigned long long>> z_separation_samples;
std::vector<std::vector<unsigned long long>> z_mean_samples;

std::vector<std::vector<Eigen::Vector3i>> polymer(number_of_threads);
std::vector<std::vector<Eigen::Vector3i>> lin_polymer(number_of_threads);

std::vector<std::vector<bool>> is_replicated(number_of_threads, std::vector<bool>(pol_length, true)); // to have an easy check for which sites are replicated. Unreplicated are set to "false" in "initialization.h"


std::vector<std::vector<double>> final_contacts_averaged(pol_length, std::vector<double>(bin_num, 0));
std::vector<std::vector<std::vector<double>>> final_contacts(number_of_stages, std::vector< std::vector<double>>(bin_num, std::vector<double>(bin_num, 0))); // This stores the contact frequencies by stages during the simulation
std::vector<std::vector<std::vector<double>>> total_contacts(number_of_threads, std::vector< std::vector<double>>(bin_num, std::vector<double>(bin_num, 0))); // This stores the contact frequencies during the simulation
std::vector<std::vector<double>> xp_contacts(bin_num, std::vector<double>(bin_num, 0));

/////// functions to be used later ///////
void run(int thread_num, int move_num);
void normalize();
void update_E(int step);
void normalize_measured_z();
void update_alpha_beta(int step);
void update_energ_coeff_mean (int step);
void update_energ_coeff_separation (int step);
void clean_up();
void set_unconstrained_energies_to_zero();

/////// script ///////

int main() {
    // Resize contact maps
    total_contacts.resize(number_of_threads, std::vector<std::vector<double>>(bin_num, std::vector<double>(bin_num, 0.0)));
    contacts.resize(number_of_threads);
    contacts_lin.resize(number_of_threads);
    contacts_inter.resize(number_of_threads);

    // Resize location maps
    locations.resize(number_of_threads);
    locations_lin.resize(number_of_threads);

    // Resize polymers
    polymer.resize(number_of_threads);
    lin_polymer.resize(number_of_threads);

    // Resize scalar vectors
    z_close.resize(number_of_threads, 0.0);
    z_far.resize(number_of_threads, 0.0);
    z_close_squared.resize(number_of_threads, 0.0);
    z_far_squared.resize(number_of_threads, 0.0);

    is_replicated.resize(number_of_threads, std::vector<bool>(pol_length, true));


    // Resize z-mean and z-separation data
    z_mean_data.resize(number_of_threads, std::vector<double>(sites_constrained_mean.size(), 0.0));
    z_separation_data.resize(number_of_threads, std::vector<double>(sites_constrained_separation.size(), 0.0));

    z_separation_samples.resize(number_of_threads,std::vector<unsigned long long>(n_constrained_separation, 0ULL));
    z_mean_samples.resize(number_of_threads,std::vector<unsigned long long>(n_constrained_mean, 0ULL));

    is_constrained_mean.assign(number_of_stages, std::vector<bool>(n_constrained_mean, false));
    is_constrained_separation.assign(number_of_stages, std::vector<bool>(n_constrained_separation, false));
    energ_coeff_mean.assign(number_of_stages, std::vector<double>(n_constrained_mean, 0.0));
    energ_coeff_separation.assign(number_of_stages, std::vector<double>(n_constrained_separation, 0.0));
    std::cout << "bin number: " << bin_num << '\n' << "polymer length: " << pol_length << '\n' << "mc moves: " << mc_moves <<'\n' << "boundaries: " << boundary_cond << '\n' << "learning rate: " << learning_rate << '\n' << "saved in " << dir <<'\n';

    //Fill site index maps (needed for reading in input data)
    for(int i=0; i<sites_constrained_mean.size(); i++){
        sites_constrained_mean_map[ sites_constrained_mean[i] ] = i;
    }
    for(int i=0; i<sites_constrained_separation.size(); i++){
        sites_constrained_separation_map[ sites_constrained_separation[i] ] = i;
    }

    //Read in constraints//
    read_input_data(); 
    for (int s=0; s<number_of_stages; ++s) {
        std::cout << "target ori separation for stage s=" << s << " is: ";
        for (double v1 : target_separations[s]) std::cout << v1 << ' ';
        std::cout << '\n';

        std::cout << "target Ter position for stage s=" << s << " is: ";
        for (double v2 : target_means[s]) std::cout << v2 << ' ';
        std::cout << '\n';
    }
    read_file(xp_contacts, HiC_file);

    //Read in energies//
    if(initEnerg) {
        check_input_energies_compatibility(old_simulation_folder, energy_data_iteration);
        // GG: throws an error if the stages in the input file are not the same as in the simulation. Could be modified to automatically find needed stages in the input folder...

        read_interaction_energies(old_simulation_folder, energy_data_iteration);
        //read_ori_energies(old_simulation_folder, energy_data_iteration);
        //read_position_energies(old_simulation_folder, energy_data_iteration);
        read_separation_energies(old_simulation_folder, energy_data_iteration);

        set_unconstrained_energies_to_zero();
    }


    //Create output directory//
    time_t time_now = time(0);   // get time now
    struct tm * now = localtime( & time_now );
    char buffer [20];
    strftime(buffer, 20, "%Y-%m-%d_%H%M", now);
    std::string buffer_str = buffer; // converts to string

    std::string general_output_folder = dir + buffer_str + "_averaged";
    std::system(("mkdir -p " + general_output_folder + "/Contacts").c_str());
    std::system(("mkdir -p " + general_output_folder + "/Energies").c_str());
    std::system(("mkdir -p " + general_output_folder + "/Positions").c_str());
    
    for (int stage = 0; stage < number_of_stages; ++stage) {
        std::string stage_folder = buffer_str + "_stage_" + std::to_string(stages[stage]);
        std::string full_path = dir + stage_folder;

        std::system(("mkdir -p " + full_path).c_str());
        std::system(("mkdir -p " + full_path + "/Configurations").c_str());
        std::system(("mkdir -p " + full_path + "/Contacts").c_str());
        std::system(("mkdir -p " + full_path + "/Energies").c_str());
        std::system(("mkdir -p " + full_path + "/Positions").c_str());

    }

    // create and seed random generators //

    // for (int i = 0; i < number_of_threads; i++) {
    //     int stage = i / replicates_per_stage; // Get the stage for this simulation
    //     generators.push_back(RandomGenerator(int(time_now) + i, pol_length, lin_length[stage], oriC));
    // }

    generators.reserve(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        int stage = i / replicates_per_stage;
        generators.emplace_back(int(time_now) + i, pol_length, lin_length[stage], oriC);
    }


    for (int s=0; s<number_of_stages; ++s) {
        double offset_stage = (int(length[s]) % 2) / 2.0; // stage-specific
        xp_z_close_simunits[s] = xp_z_close[s]*(length[s] + 2*radius)- (length[s]/2 + radius - offset_stage);
        xp_z_far_simunits[s]   = xp_z_far[s]  *(length[s] + 2*radius) - (length[s]/2 + radius - offset_stage);
    }


    //Save the input parameters of the simulation in "sim_params.txt"//
    get_sim_params(general_output_folder);
    std::cout<<"Initialising..."<<std::endl;
    //Initialize configurations//
    for (int i = 0; i < number_of_threads; i += number_of_threads_in_parallel) {
        // Launch a batch of threads
        int batch_size = std::min(number_of_threads_in_parallel, number_of_threads - i);
        std::vector<std::thread> iniThreads(batch_size);
        for (int j = 0; j < batch_size; ++j) {
            iniThreads[j] = std::thread(initialize, i + j, init_config_number);
        }
        // Join all threads in this batch
        for (auto& t : iniThreads) {
            t.join();
        }
        // Debug check for threads 48–58 after initialization
        for (int t = i; t < i + batch_size; ++t) {
            if (t >= 48 && t <= 58) {
                std::cout << "[Check] Thread " << t
                        << ": polymer size = " << polymer[t].size()
                        << ", lin_polymer size = " << lin_polymer[t].size()
                        << ", generators valid = " << (t < generators.size())
                        << ", is_replicated size = " << is_replicated[t].size()
                        << std::endl;
            }
        }

        for (int t = i; t < i + batch_size; ++t) {
            if (t >= 48 && t <= 58) {
                std::cout << "[Check] Thread " << t
                        << ": contacts size = " << contacts[t].size()
                        << ", contacts_lin size = " << contacts_lin[t].size()
                        << ", contacts_inter size = " << contacts_inter[t].size()
                        << ", locations size = " << locations[t].size()
                        << ", locations_lin size = " << locations_lin[t].size()
                        << std::endl;
            }
        }

    }

    auto start = std::chrono::high_resolution_clock::now();
    std::cout<<"Burning in..."<<std::endl;

    // Burn in configurations in parallel batches
    for (int i = 0; i < number_of_threads; i += number_of_threads_in_parallel) {
        int batch_size = std::min(number_of_threads_in_parallel, number_of_threads - i);
        std::vector<std::thread> burnThreads(batch_size);
        for (int j = 0; j < batch_size; ++j) {
            burnThreads[j] = std::thread(burn_in, i + j, burn_in_steps);
        }
        for (auto& t : burnThreads) {
            t.join();
        }
    }
    // After burn-in steps
    burn_in_done = true;


    std::cout << "Running MC moves..." << std::endl;

    std::cout<<"Start gradient descent..."<<std::endl;

    // Perform inverse algorithm //
    for (auto step = 0; step <= it_steps; step++) {

        if (step < 20) {
            learning_rate = learning_rate_separations = learning_rate_close =learning_rate_means = 0.5;
        } else if (step < 40) {
            learning_rate = learning_rate_separations = learning_rate_close =learning_rate_means = 0.2;
        } else if (step < 60) {
            learning_rate = learning_rate_separations = learning_rate_close =learning_rate_means = 0.1;
        }
        else if (step < 70) {
            learning_rate = learning_rate_separations = learning_rate_close=learning_rate_means  = 0.05;
        } else { // step < 60 (and also step==60 due to <= it_steps)
            learning_rate =  0.05;
            learning_rate_separations = learning_rate_close =learning_rate_means = 0.02;
        }
        // also zero per-thread accumulators and counters
        for (int l = 0; l < number_of_threads; ++l) {
            for (int i = 0; i < n_constrained_separation; ++i) {
                z_separation_data[l][i] = 0.0;
                z_separation_samples[l][i] = 0ULL;
            }
        }
        for (int l = 0; l < number_of_threads; ++l) {
            for (int i = 0; i < n_constrained_mean; ++i) {
                z_mean_data[l][i] = 0.0;
                z_mean_samples[l][i] = 0ULL;
            }
        }


        //Forward run//

        //copying positioning energies from vectors to unordered maps (for convenience when calculating delta_E)
        for(int stage=0; stage<number_of_stages; stage++) {
            for (int i = 0; i < sites_constrained_mean.size(); i++) {
                energy_mean_map[stage][sites_constrained_mean[i]] = energ_coeff_mean[stage][i];
            }
            for (int i = 0; i < sites_constrained_separation.size(); i++) {
                energy_separation_map[stage][sites_constrained_separation[i]] = energ_coeff_separation[stage][i];
            }
        }


        std::cout << "Running MC moves..." << std::endl;

        for (int i = 0; i < number_of_threads; i += number_of_threads_in_parallel) {
            int batch_size = std::min(number_of_threads_in_parallel, number_of_threads - i);
            std::vector<std::thread> runThreads(batch_size);
            for (int j = 0; j < batch_size; ++j) {
                runThreads[j] = std::thread(run, i + j, std::round(mc_moves * sqrt(step + 1)));
            }
            for (auto& t : runThreads) {
                t.join();
            }
        }


        
        for (int s = 0; s < number_of_stages; ++s) {
            int i = replicates_per_stage * s;
            double cell = length[s] + 2*radius;
            // z at oriC (ring always exists; lin only if replicated at this stage)
            double zc = static_cast<double>(polymer[i][oriC][2]);
            double zf = (lin_length[s] > 0)
                        ? static_cast<double>(lin_polymer[i][oriC][2])
                        : zc;
            
            double z_pos_ter = static_cast<double>(polymer[i][Ter][2]);
            double z_pos_ter_norm = z_pos_ter / cell;
            double sep_sim  = std::abs(zf - zc);
            double sep_norm = sep_sim / cell;

            std::cout << "step=" << step
                    << " t=" << i
                    << " s=" << s
                    << " cell=" << cell
                    << " lin_len=" << lin_length[s]
                    << " z_ring=" << zc
                    << " z_lin="  << zf
                    << " sep_sim=" << sep_sim
                    << " sep_norm="<< sep_norm
                    << " z_pos_ter_sim=" << z_pos_ter
                    << " z_pos_ter_norm="<< z_pos_ter_norm
                    << '\n';
        }

        // Step 1: average over replicates per stage
        for (int s = 0; s < number_of_stages; s++) {
            for (int r = 0; r < replicates_per_stage; r++) {
                int idx = s * replicates_per_stage + r;
                for (int i = 0; i < bin_num; i++) {
                    for (int j = 0; j < bin_num; j++) {
                        final_contacts[s][i][j] += total_contacts[idx][i][j];
                    }
                }
            }
            // Divide by number of replicates to get average per stage
            for (int i = 0; i < bin_num; i++) {
                for (int j = 0; j < bin_num; j++) {
                    final_contacts[s][i][j] /= replicates_per_stage;
                }
            }
        }

        //std::vector<double> w = {1., 1., 1., 1., 0.9270073, 0.35766423}; // biological weights (prob not right)
        std::vector<double> w = {1., 1., 1., 1., 1., 1.}; // uniform weights
        double w_sum = std::accumulate(w.begin(), w.end(), 0.0);
        for (auto& wi : w) wi /= w_sum; // Normalize to sum to 1

        // Step 2: weighted average over stages
        for (int s = 0; s < number_of_stages; s++) {
            for (int i = 0; i < bin_num; i++) {
                for (int j = 0; j < bin_num; j++) {
                    final_contacts_averaged[i][j] += w[s] * final_contacts[s][i][j];
                }
            }
        }

        //Ensures final_contacts_averaged is symmetric//
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                final_contacts_averaged[j][i] = final_contacts_averaged[i][j];
            }
        }

        // Ensure final_contacts[stage] is symmetric for all stages
        for (int s = 0; s < number_of_stages; s++) {
            for (int i = 0; i < bin_num; i++) {
                for (int j = 0; j < bin_num; j++) {
                    final_contacts[s][j][i] = final_contacts[s][i][j];
                }
            }
        }

        //Saves configuration of the strand and ring for each thread with a separate file for each stage at the end of each step
        for (int l = 0; l < number_of_threads; l++) {
            int s = l / replicates_per_stage;
            std::string output_folder = buffer_str + "_stage_" + std::to_string(stages[s]);

            get_configuration(step, "strand", l, output_folder);
            get_configuration(step, "ring", l, output_folder);
        }

        // std::cout << "[CHK] sites_sep.size()=" << sites_constrained_separation.size() << '\n';
        // for (int l = 0; l < std::min(4, number_of_threads); ++l) {
        //     double maxv = 0.0;
        //     for (double v : z_separation_data[l]) maxv = std::max(maxv, v);
        //     std::cout << "[CHK] thread " << l
        //             << " z_sep_len=" << z_separation_data[l].size()
        //             << " max=" << maxv << '\n';
        // }


        for (int l = 0; l < std::min(3, number_of_threads); ++l) {
            int s = l / replicates_per_stage;
            double cell = length[s] + 2*radius;

            double zc = polymer[l][oriC][2];
            double zf = (lin_length[s] > 0 ? lin_polymer[l][oriC][2] : polymer[l][oriC][2]);

            double sep_sim  = std::abs(zf - zc);
            double sep_norm = sep_sim / cell;

            // std::cout << std::fixed << std::setprecision(4)
            //         << "[UNIT] s="<< s
            //         << " sep_sim=" << sep_sim
            //         << " cell="    << cell
            //         << " sep_norm="<< sep_norm
            //         << " target="  << (target_separations[s].empty()? -1 : target_separations[s][0])
            //         << '\n';
        }
        //calculation of positional constraint values after forward run//


        // inside the step loop, before using them
        z_mean_data_tot.assign(number_of_stages,std::vector<double>(n_constrained_mean, 0.0));
        z_separation_data_tot.assign(number_of_stages, std::vector<double>(n_constrained_separation, 0.0));

        std::vector<std::vector<double>> stage_sum_sep(
            number_of_stages, std::vector<double>(n_constrained_separation, 0.0));
        std::vector<std::vector<unsigned long long>> stage_count_sep(
            number_of_stages, std::vector<unsigned long long>(n_constrained_separation, 0ULL));

        std::vector<std::vector<double>> stage_sum_val(
            number_of_stages, std::vector<double>(n_constrained_mean, 0.0)); // Update to track mean data
        std::vector<std::vector<unsigned long long>> stage_count_val(
            number_of_stages, std::vector<unsigned long long>(n_constrained_mean, 0ULL)); // Count for mean data
        


        for (int l=0; l < number_of_threads; l++) { //JM: calculation of average z coordinates & squares of both oris
            int stage = l / replicates_per_stage;
            if (lin_length[stage]==0) {
                z_close_tot[stage] += z_close[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
                z_close_squared_tot[stage] += z_close_squared[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
            }
            else {
                z_close_tot[stage] += z_close[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
                z_far_tot[stage] += z_far[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
                z_close_squared_tot[stage] += z_close_squared[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
                z_far_squared_tot[stage] += z_far_squared[l]/number_of_threads*number_of_stages/std::round(mc_moves * sqrt(step+1))*res;
            }
            // Accumulate the mean z positions for each site across threads
            for (int i = 0; i < n_constrained_mean; ++i) {
                stage_sum_val[stage][i]   += z_mean_data[l][i];  // Accumulate z_mean_data
                stage_count_val[stage][i] += z_mean_samples[l][i]; // Track the number of samples for each site
                //std::cout <<"stage_count_val of stage "<< stage << " is " << stage_count_val[stage][i] <<" and stage_sum_val  "<< stage_sum_val[stage][i] ;
            }

            // Accumulate the separation data 
            for (int i = 0; i < n_constrained_separation; ++i) {
                stage_sum_sep[stage][i]   += z_separation_data[l][i];
                stage_count_sep[stage][i] += z_separation_samples[l][i];
            }   
        }       
        // Convert to per-stage averages for mean data
        for (int s = 0; s < number_of_stages; ++s) {
            for (int i = 0; i < n_constrained_mean; ++i) {
                
                // Debug print to check the values of stage_count_val and stage_sum_val
                // std::cout << "Checking stage_count_val[" << s << "][" << i << "]: " << stage_count_val[s][i] << std::endl;
                // std::cout << "Checking stage_sum_val[" << s << "][" << i << "]: " << stage_sum_val[s][i] << std::endl;

                // Debugging the mean value calculation process
                double mean_val = stage_count_val[s][i] ? stage_sum_val[s][i] / double(stage_count_val[s][i]) : 0.0;
                // std::cout << "Calculated mean_val[" << s << "][" << i << "] = " << mean_val << std::endl;

                // std::cout <<"stage_count_val of stage "<< s << " is " << stage_count_val[s][i] <<" and stage_sum_val  "<< stage_sum_val[s][i] << std::endl;

                // Optional unit fix if sampling was in sim units:
                // mean_val /= (length[s] + 2*radius);

                z_mean_data_tot[s][i] = mean_val;  // Store the averaged value for each stage
            }
        }

        // Convert to per-stage averages for separation data.
        for (int s = 0; s < number_of_stages; ++s) {
            for (int i = 0; i < n_constrained_separation; ++i) {
                
                // Debug print to check the values of stage_count_sep and stage_sum_sep
                // std::cout << "Checking stage_count_sep[" << s << "][" << i << "]: " << stage_count_sep[s][i] << std::endl;
                // std::cout << "Checking stage_sum_sep[" << s << "][" << i << "]: " << stage_sum_sep[s][i] << std::endl;

                // std::cout <<"stage_count_sep of stage "<< s << " is " << stage_count_sep[s][i] <<" and stage_sum_sep  "<< stage_sum_sep[s][i] << std::endl;

                double mean_sep = stage_count_sep[s][i] ? stage_sum_sep[s][i] / double(stage_count_sep[s][i]) : 0.0;

                // Optional unit fix if sampling was in sim units:
                // mean_sep /= (length[s] + 2*radius);

                z_separation_data_tot[s][i] = mean_sep;
            }
        }


        for (int s=0; s<number_of_stages;s++) {
            z_close_var[s] = z_close_squared_tot[s] - pow(z_close_tot[s],2);
            z_far_var[s] = z_far_squared_tot[s] - pow(z_far_tot[s],2);
        }

        for (int s = 0; s < number_of_stages; ++s) {
            // std::cout << "[SEP] stage " << s
            //             << " mean_sep_sim=" << z_separation_data_tot[s][0]
            //             << " n=" << stage_count_sep[s][0] << '\n';
            // std::cout << "[POS] stage " << s
            //             << " mean_pos_sim=" << z_mean_data_tot[s][0]
            //             << " n=" << stage_count_val[s][0] << '\n';
        }



        //update simulation parameters//
        normalize();
        update_E(step); //JM: Updates pairwise interaction energies
        // before normalize_measured_z()
        std::cout << "[PRE] sep s0=" << z_separation_data_tot[0][0] << '\n';
        std::cout << "[PRE] pos s0=" << z_mean_data_tot[0][0] << '\n';        
        if (n_constrained_separation > 0) {
            for (int s = 0; s < number_of_stages; ++s) {
                double pre_sep = z_separation_data_tot[s][0];
                double cell = length[s] + 2*radius;
                std::cout << "[UNITSCHK] s="<<s<<" pre sep="<<pre_sep
                        << " post≈"<< pre_sep / cell
                        << " cell="<<cell << '\n';
            }
        }
        if (n_constrained_mean > 0) {
            for (int s = 0; s < number_of_stages; ++s) {
                double pre_mean = z_mean_data_tot[s][0];
                double cell = length[s] + 2*radius;
                double offset_stage = (int(length[s]) % 2) / 2.0;
                double left_pole_z  = -(length[s]/2 + radius - offset_stage);

                std::cout << "[UNITSCHK] s="<<s<<" pre mean="<<pre_mean
                        << " post≈"<< (pre_mean - left_pole_z) / cell
                        << "different post (is it better)"<< pre_mean / cell
                        << " cell="<<cell << '\n';
            }
        }
        normalize_measured_z();//GG: express measured values of z in the units of cell length, and shift such that z=0 at the left pole
        std::cout << "[POST] sep s0=" << z_separation_data_tot[0][0]
                << " (target " << target_separations[0][0] << ")\n";
        std::cout << "[POST] pos s0=" << z_mean_data_tot[0][0]
                << " (target " << target_means[0][0] << ")\n";                

        update_energ_coeff_mean(step);
        update_energ_coeff_separation(step);
        update_alpha_beta(step);

        //save simulation parameters//
        get_final_contacts(step, buffer_str);
        get_final_contacts_averaged(step, general_output_folder);
        get_energ_coeff_mean(step,general_output_folder);
        get_energ_coeff_separation(step,general_output_folder);
        get_alpha_beta(step,general_output_folder); //JM: Saves new values alpha and beta
        get_z_lin_far_close(step,general_output_folder); //JM: Saves new values ori positions
        get_energies_plot(step,general_output_folder); //JM: Saves pairwise interaction energies
        for (int i=0;i< n_constrained_mean;i++){
            get_z_mean_rest(i, step,general_output_folder);
        }
        for (int i=0;i< n_constrained_separation;i++){
            get_z_separation_rest(i, step,general_output_folder);
        }

        clean_up(); //JM: sets everything to zero again before next round
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}

void run(int thread_num, int move_num) {
    int m = 0;
    while (m < move_num) {    //perform polymer simulation
        move(thread_num, m); //GG: "m" is passed by reference now
    }
    for (const auto& elem : contacts[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
    for (const auto& elem : contacts_lin[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
    for (const auto& elem : contacts_inter[thread_num]) { //add the contacts remaining at the end of the simulation (during a the simulation, a contact is only added to 'total_contacts' when a contact is removed)
        total_contacts[thread_num][std::min(elem.first.first,elem.first.second)][std::max(elem.first.first,elem.first.second)] += double(move_num) - elem.second;
    }
}

void normalize() {
    double sum = 0;
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                sum += final_contacts_averaged[i][j];
            }
        }
    }
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            final_contacts_averaged[i][j] *= float(bin_num) / sum;
        }
    }

    // Normalize final_contacts per stage
    for (int s = 0; s < number_of_stages; s++) {
        double sum_stage = 0;
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                if (std::abs(i - j) > 1 && !(i == 0 && j == bin_num - 1) && !(i == bin_num - 1 && j == 0)) {
                    sum_stage += final_contacts[s][i][j];
                }
            }
        }
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j < bin_num; j++) {
                final_contacts[s][i][j] *= float(bin_num) / sum_stage;
            }
        }
    }

}

void update_E(int step) {
    for (int i = 0; i < bin_num; i++) {                 //update energy map
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                auto comp = xp_contacts[i][j];
                if (comp != 0) {
                    if (step>update_cap_onset) { //JM we impose an upper limit on the update after the first few stages, to facilitate convergence
                        double prop_update = learning_rate * (final_contacts_averaged[i][j] - comp) / pow(std::abs(comp), 0.5);
                        if (abs(prop_update)<update_cap_factor * learning_rate){
                            Interaction_E[i][j] +=prop_update;
                        }
                        else{
                            if (prop_update<0){
                                Interaction_E[i][j] -=update_cap_factor * learning_rate;
                            }
                            else{
                                Interaction_E[i][j] +=update_cap_factor * learning_rate;
                            }
                        }
                    }
                    else{
                        Interaction_E[i][j] +=
                                learning_rate * (final_contacts_averaged[i][j] - comp) / pow(std::abs(comp), 0.5);
                    }
                }
            }
        }
    }
    double meanE = 0;                                   //mean energy
    double sum = 0;
    for (int i = 0; i < bin_num; i++) {
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                auto comp = xp_contacts[i][j];
                meanE += Interaction_E[i][j] * comp;
                sum += comp;
            }
        }
    }
    if (sum != 0) {
        meanE /= sum;
    }
    // std::cout << std::setw(4) << std::left << step << std::setw(10)  << std::setprecision(3) << std::right << meanE << "\n";
    for (int i = 0; i < bin_num; i++) {                 //set average energy to zero
        for (int j = 0; j < bin_num; j++) {
            if (std::abs(i-j) > 1 && !(i==0 && j==bin_num - 1) && !(i==bin_num - 1 && j==0)) {
                Interaction_E[i][j] -= meanE;
            }
        }
    }
}

void normalize_measured_z(){ // express measured values of z in the units of cell length, also shifts such that z=0 at the left pole
    for (int s=0; s<number_of_stages; s++){
        double total_cell_length = length[s] + 2*radius;
        double offset_stage = (int(length[s]) % 2) / 2.0;
        double left_pole_z  = -(length[s]/2 + radius - offset_stage);
        z_close_tot[s] = (z_close_tot[s] - left_pole_z)/total_cell_length;
        z_far_tot[s] = (z_far_tot[s] - left_pole_z)/total_cell_length;
        z_close_var[s] = z_close_var[s] / pow(total_cell_length, 2);
        z_far_var[s] = z_far_var[s] / pow(total_cell_length, 2);
        for(int i=0; i<sites_constrained_mean.size(); i++){

            // z_mean_data_tot[s][i] = (z_mean_data_tot[s][i] - left_pole_z) / total_cell_length;
            z_mean_data_tot[s][i] = z_mean_data_tot[s][i]  / total_cell_length;
        }
        for(int i=0; i<sites_constrained_separation.size(); i++){
            z_separation_data_tot[s][i] = z_separation_data_tot[s][i] / total_cell_length;
        }
    }
}

// void update_alpha_beta (int step) {
//     for (int s = 0; s < number_of_stages; s++) {
//         if(is_constrained_ori[s][0]){
//             alpha[s] += learning_rate_close*(z_close_tot[s]-xp_z_close[s]);
//         }
//         if(is_constrained_ori[s][2]){
//             alpha2[s] += learning_rate_close_var * (z_close_var[s] - xp_z_close_var[s]);
//         }
//         if(is_replicated[s][oriC]){
//             if(is_constrained_ori[s][1]){
//                 beta[s] += learning_rate_far * (z_far_tot[s]-xp_z_far[s]);
//             }
//             if(is_constrained_ori[s][3]){
//                 beta2[s] += learning_rate_far_var * (z_far_var[s]- xp_z_far_var[s]);
//             }
//         }
//     }
// }

void update_alpha_beta(int step) {
    for (int s = 0; s < number_of_stages; ++s) {
        // "close" ori terms (exist regardless of replication)
        if (is_constrained_ori[s][0]) {
            alpha[s] += learning_rate_close * (z_close_tot[s] - xp_z_close[s]);
        }
        if (is_constrained_ori[s][2]) {
            alpha2[s] += learning_rate_close_var * (z_close_var[s] - xp_z_close_var[s]);
        }

        // Check if ANY replicate in this stage has the ori replicated
        bool any_rep = false;
        for (int r = 0; r < replicates_per_stage && !any_rep; ++r) {
            const int idx = s * replicates_per_stage + r;  // thread index
            any_rep |= is_replicated[idx][oriC];
        }

        // "far" ori terms only make sense once replication occurred
        if (any_rep) {
            if (is_constrained_ori[s][1]) {
                beta[s]  += learning_rate_far * (z_far_tot[s] - xp_z_far[s]);
            }
            if (is_constrained_ori[s][3]) {
                beta2[s] += learning_rate_far_var * (z_far_var[s] - xp_z_far_var[s]);
            }
        }
    }
}



void update_energ_coeff_mean (int step) {
    for(int stage=0; stage<number_of_stages; stage++) {
        for (int i = 0; i < sites_constrained_mean.size(); i++) {
            if (step % 10 == 0) {
                // std::cout << " The z_mean_position at stage "<< stage << " and iteration step= "<< step <<" is :" << z_mean_data_tot[stage][i] << std::endl;
            }
            if(is_constrained_mean[stage][i]){
                // std::cout<<"Energy is being updated!"<<std::endl;
                energ_coeff_mean[stage][i] += learning_rate_means * (z_mean_data_tot[stage][i] - target_means[stage][i]);
            }
            // for (double v : target_means[stage]) std::cout << v << ' ';
            // std::cout << '\n';

            // std::cout << "energ_coeff_mean for stage s=" << stage << " is: ";
            // for (double v : energ_coeff_mean[stage]) std::cout << v << ' ';
            // std::cout << '\n';

            // std::cout << "z_mean_data_tot for stage s=" << stage << " is: ";
            // for (double v : z_mean_data_tot[stage]) std::cout << v << ' ';
            // std::cout << '\n';
            // std::cout << "n_sites=" << sites_constrained_mean.size() << '\n';  
        }
    }
}

// void update_energ_coeff_separation (int step) {
//     for(int stage=0; stage<number_of_stages; stage++) {
//         for (int i = 0; i < sites_constrained_separation.size(); i++) {
//             // if(is_replicated[stage][sites_constrained_separation[i]] && is_constrained_separation[stage][i]) {
//             //     energ_coeff_separation[stage][i] += learning_rate_separations * (z_separation_data_tot[stage][i] - target_separations[stage][i]);
//             //}
//             bool any_rep=false;
//             for (int r=0; r<replicates_per_stage; ++r){
//                 int idx = stage*replicates_per_stage + r;
//                 any_rep |= is_replicated[idx][ sites_constrained_separation[i] ];
//                 }
//                 if (any_rep && is_constrained_separation[stage][i]) {
//                     energ_coeff_separation[stage][i] += learning_rate_separations * (z_separation_data_tot[stage][i] - target_separations[stage][i]);
//                 }
//         }
//     }
// }

void update_energ_coeff_separation (int step) {
    for (int stage = 0; stage < number_of_stages; ++stage) {
        for (int i = 0; i < (int)sites_constrained_separation.size(); ++i) {
            if (!is_constrained_separation[stage][i]) continue;   // early skip

            int site = sites_constrained_separation[i];
            bool any_rep = false;
            for (int r = 0; r < replicates_per_stage && !any_rep; ++r) {
                int idx = stage * replicates_per_stage + r;
                any_rep |= is_replicated[idx][site];
            }
            if (!any_rep) continue;

            energ_coeff_separation[stage][i] +=
                learning_rate_separations *
                (z_separation_data_tot[stage][i] - target_separations[stage][i]);
        }
    //     for (double v : target_separations[stage]) std::cout << v << ' ';
    //     std::cout << '\n';

    //     std::cout << "energ_coeff_separation for stage s=" << stage << " is: ";
    //     for (double v : energ_coeff_separation[stage]) std::cout << v << ' ';
    //     std::cout << '\n';

    //     std::cout << "z_separation_data_tot for stage s=" << stage << " is: ";
    //     for (double v : z_separation_data_tot[stage]) std::cout << v << ' ';
    //     std::cout << '\n';
    // std::cout << "n_sites=" << sites_constrained_separation.size() << '\n';    
    }
}

void set_unconstrained_energies_to_zero() { //for safety: makes sure that the initial energies with no corresponding experimental constraints are 0
    for (int s=0; s<number_of_stages; s++){
        if(!is_constrained_ori[s][0]){
            alpha[s] = 0;
        }
        if(!is_constrained_ori[s][1]){
            beta[s] = 0;
        }
        if(!is_constrained_ori[s][2]){
            alpha2[s] = 0;
        }
        if(!is_constrained_ori[s][3]){
            beta2[s] = 0;
        }
        for(int i=0; i<n_constrained_mean; i++) {
            if (!is_constrained_mean[s][i]) {
                energ_coeff_mean[s][i] = 0;
            }
        }
        for(int i=0; i<n_constrained_separation; i++) {
            if(!is_constrained_separation[s][i]){
                energ_coeff_separation[s][i] = 0;
            }
        }
    }
}

void clean_up() {
    std::vector<double> zeroVec(bin_num, 0);
    std::fill(final_contacts_averaged.begin(), final_contacts_averaged.end(), zeroVec);
    for (int s = 0; s < number_of_stages; s++) {
        std::fill(final_contacts[s].begin(), final_contacts[s].end(), zeroVec);
        z_close_var[s] = 0;
        z_far_var[s] = 0;
        z_close_squared_tot[s]=0;
        z_far_squared_tot[s]=0;
        z_close_tot[s]=0;
        z_far_tot[s]=0;
        for(int i=0; i<sites_constrained_mean.size(); i++){
            z_mean_data_tot[s][i] = 0;
        }
        for(int i=0; i<sites_constrained_separation.size(); i++){
            z_separation_data_tot[s][i] = 0;
        }
    }

    for (int l = 0; l < number_of_threads; l++) {
        int stage = l / replicates_per_stage;
        std::fill(total_contacts[l].begin(), total_contacts[l].end(), zeroVec);
        z_close[l] = 0;
        z_far[l] = 0;
        z_close_squared[l] = 0;
        z_far_squared[l] = 0;

        for(int i=0; i<sites_constrained_mean.size(); i++){
            z_mean_data[l][i] = 0;
        }
        for(int i=0; i<sites_constrained_separation.size(); i++){
            z_separation_data[l][i] = 0;
        }

        contacts[l].clear();
        contacts_lin[l].clear();
        contacts_inter[l].clear();
        for (int i = 0; i < pol_length; i+=reduction_factor) {
            int red_i=i/reduction_factor;
            if (locations[l].find({polymer[l][i][0], polymer[l][i][1], polymer[l][i][2]}) != locations[l].end()) {
                for (auto elem : locations[l][{polymer[l][i][0], polymer[l][i][1],polymer[l][i][2]}]) {
                    contacts[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                }
            }
        }
        for (int i = 0; i < pol_length; i+=reduction_factor) {
            int red_i=i/reduction_factor;
            if (oriC + lin_length[stage] >= pol_length) {
                if (i > oriC - lin_length[stage] || i < (oriC + lin_length[stage]) % pol_length) {
                    if (locations_lin[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) != locations_lin[l].end()) {
                        for (auto elem : locations_lin[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_lin[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                        }
                    }
                    if (locations[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) !=
                        locations[l].end()) {
                        for (auto elem : locations[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_inter[l][{elem, red_i}] = 0;
                        }
                    }
                }
            } else {
                if (i > oriC - lin_length[stage] && i < oriC + lin_length[stage]) {
                    if (locations_lin[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) != locations_lin[l].end()) {
                        for (auto elem : locations_lin[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_lin[l][{std::min(elem, red_i), std::max(elem, red_i)}] = 0;
                        }
                    }
                    if (locations[l].find({lin_polymer[l][i][0], lin_polymer[l][i][1],lin_polymer[l][i][2]}) !=
                        locations[l].end()) {
                        for (auto elem : locations[l][{lin_polymer[l][i][0],lin_polymer[l][i][1],lin_polymer[l][i][2]}]) {
                            contacts_inter[l][{elem, red_i}] = 0;
                        }
                    }
                }
            }
        }
    }
}
