#include<vector>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <functional>
#include <algorithm>
#include <tuple>

#include "stats.hpp"
#include "trajectory.hpp" 

double calculate_pvalue(Random & random, Trajectory const & trajectory, const int max_num_bootstraps=10000, const int min_numerator_counts=100);

const char delim=',';

int main(int argc, char * argv[]){

    bool DISABLED=false;
    if(argc>1){
        DISABLED = true;
        std::cerr << "Warning: skipping pvalue calculation & setting p=0!" << std::endl; 
    }
    
    // random number generator
    // deterministic seed to ensure reproducibility
    // once pipeline is completed
    auto random = create_random(42); 
    
    // used for reporting purposes
    int num_processed = 0;
    int num_passed = 0;
    int num_surprising = 0;
    
    // iterate over all trajectory records in file
    std::string line;
    
    // first print header
    std::getline(std::cin,line);
    std::cout << line << std::endl;
    
    while(std::getline(std::cin,line)){
        
        // parse trajectory record
        std::stringstream line_stream(line);
        std::string site_id;
        std::string item;
        std::string subitem;
        
        std::vector<double> alts;
        std::vector<double> depths;
        
        // read site_id
        std::getline(line_stream, site_id, '\t');
        
        // read alts and depths trajectory
        auto trajectory = Trajectory();
        while(std::getline(line_stream, item, '\t')){
            // loop over samples
            double alt, depth;
            std::stringstream item_stream(item);
            std::getline(item_stream, subitem, ',');
            alt = std::stof(subitem);
            std::getline(item_stream, subitem, ',');
            depth = std::stof(subitem);
            trajectory.push_back( Sample{ alt, depth } );
            //std::cerr << alt << "," << depth << " ";
        }
        //std::cerr << std::endl;
        
        double pvalue;
            
        num_processed+=1;
                
        if(!passes_filter(trajectory)){
            pvalue = 2; // shouldn't happen!    
        }
        else{
            // actually calculate autocorrelation score
            num_passed+=1;
            if(!DISABLED){
                pvalue = calculate_pvalue(random, trajectory); 
            }
            else{
                pvalue = 0;
            }
        }
            
        // increment number of "surprising" trajectories
        // used only for display purposes during run
        if(pvalue <= 5e-02){
            num_surprising += 1; 
        }
        
        if(num_processed % 10000 == 0){
                std::cerr << num_processed/10000 << "0k trajectories processed, " << num_passed << " passed, " << num_surprising << " surprising!\n";
        }
        
        // print annotated line
        std::cout << site_id << "|" << pvalue;
        for(auto const & sample : trajectory){
            std::cout << "\t" << sample.alt << "," << sample.depth;
        }
        std::cout << std::endl;  
    }
    
    std::cerr << "Finished: " << num_processed << " trajectories processed, " << num_passed << " passed, " << num_surprising << " surprising!\n";
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculates pvalue for trajectory based on likelihood ratio test
//
// Returns: pvalue
//
///////////////////////////////////////////////////////////////////////////////


double calculate_pvalue(Random & random, Trajectory const & observed_trajectory, const int max_num_bootstraps, const int min_numerator_counts){
     
    if(!passes_filter(observed_trajectory)){
        std::cerr << "Trajectory did not pass filter, should not get here!" << std::endl;
        return 2;
    }
    
    // null model is based on resampling of observed trajectory
    auto trajectory_generator = BinomialTrajectoryGenerator(observed_trajectory);
    
    double observed_LRT =  calculate_LRT(trajectory_generator.unmasked_observed_trajectory);
    
    //std::cerr << observed_LRT << std::endl;
    
    // calculate # bootstrapped Ts > T
    int num_greater_LRTs = 0;
    int current_num_bootstraps = 0;    
    
    //std::cerr << "Starting bootstraps!" << std::endl;
    while(current_num_bootstraps < max_num_bootstraps){
        
        current_num_bootstraps+=1;
        
        //std::cerr << "Bootstrap " << current_num_bootstraps << std::endl;
            
        trajectory_generator.generate_bootstrapped_trajectory(random);
        
        // calculate bootstrapped test statistics
        double bootstrapped_LRT = calculate_LRT(trajectory_generator.unmasked_bootstrapped_trajectory);
        
        //std::cerr << bootstrapped_LRT << std::endl;
        
        // compare to observed values
        if( bootstrapped_LRT >= observed_LRT ){
            // bootstrapped trajectory is at least as extreme
            num_greater_LRTs+=1;
        }
        
        if( num_greater_LRTs > min_numerator_counts ){ 
            break;
        }
    }   
    
    // calculate pvalues
    double pvalue = (num_greater_LRTs+1.0)/(current_num_bootstraps+1.0);
    return pvalue;
}