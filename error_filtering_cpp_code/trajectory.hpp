#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include<vector>
#include <iostream>
#include <random>
#include <functional>
#include <algorithm>
#include <tuple>

#include "stats.hpp"

const double max_error_rate = 0.01;
const double prevalence_threshold = 0.5;

// A sample could represent a technical replicate, 
// the sum of all technical replicates for a timepoint
// or the sum of all techncial replicates for a host
class Sample{
    public:
        double alt;
        double depth;
};

// a Trajectory object is just a collection of samples.
// (order is not important; vector used for performance)
typedef std::vector<Sample> Trajectory; 

// calculates the average frequency using the binomial MLE
inline double calculate_average_frequency(Trajectory const & trajectory){

    // first calculate average frequency
    double total_alts = 0;
    double total_depths = 0;
    for(auto const & sample : trajectory){
        total_alts += sample.alt;
        total_depths += sample.depth;
    }
    double avg_frequency = total_alts*1.0/total_depths;
    return avg_frequency;
}

inline double calculate_error_rate(Trajectory const & trajectory){
    return std::min( calculate_average_frequency(trajectory), max_error_rate);
}

// We assume that the trajectory has been pre-processed 
// such that low-coverage samples have depth=0
inline bool filter_sample(Sample const & sample){
    if(sample.depth==0){
        return true;
    }
    else{
        return false;
    } 
}

/* Filters out trajectories that do not have at least minimum fraction
   of covered samples */
inline bool passes_filter(Trajectory const & trajectory){  

    // we assume that sites have been pre-filtered by the pipe_snps() function
    // that feeds data to this program.
    return true;

    int num_covered_samples = 0;
    int num_total_samples = 0;
    for(auto const & sample : trajectory){
        num_total_samples += 1;
        if(sample.depth > 0.5){
            num_covered_samples += 1;
        }    
    }
    
    if(num_covered_samples < prevalence_threshold*num_total_samples){
        return false;
    }
    else{
        return true;
    }
}

class BinomialAltGenerator{
    public:
        int D;
        double p;
        BinomialAltGenerator(): D(0), p(0), sample_A(D, p) {};
        BinomialAltGenerator(int D, double p): D(D), p(p), sample_A(D, p) {};
        std::binomial_distribution<int> sample_A;
};

class BinomialTrajectoryGenerator{
    public:
    
        Trajectory bootstrapped_trajectory; // a full bootstrapped replica of the observed trajectory
        Trajectory unmasked_observed_trajectory; // only non-filtered points in the observed trajectory
        Trajectory unmasked_bootstrapped_trajectory; // only non-filtered points in the bootstrapped trajectory 

        double p;
        
        BinomialTrajectoryGenerator(Trajectory const & trajectory){
    
            p = calculate_error_rate(trajectory);
            
            //std::cerr << "Error rate = " << p << std::endl;
            
            for(int i=0,imax=trajectory.size();i<imax;++i){
                Sample const & sample = trajectory[i];
                bootstrapped_trajectory.push_back(sample);
                if(!filter_sample(sample)){
                    
                    unmasked_indices.push_back(i);
                    unmasked_observed_trajectory.push_back(sample);
                    unmasked_bootstrapped_trajectory.push_back(sample);
                    unmasked_alt_generators.push_back(BinomialAltGenerator((int) sample.depth, p));
                }
            }
        };
        
        Trajectory & generate_bootstrapped_trajectory(Random & random){
            int num_attempts = 0;
            do{
                num_attempts+=1;
                //std::cerr << num_attempts << std::endl;
                // resample alts 
                for(int i=0,imax=unmasked_bootstrapped_trajectory.size(); i<imax; ++i){
                    //std::cerr << "Before sampling: " << unmasked_alt_generators[i].D << " " << unmasked_alt_generators[i].p << std::endl;
                    int new_alt = unmasked_alt_generators[i].sample_A(random);
                    //std::cerr << new_alt << std::endl;
                    unmasked_bootstrapped_trajectory[i].alt = new_alt;
                    bootstrapped_trajectory[unmasked_indices[i]].alt = new_alt;    
                }
            }while(!passes_filter(bootstrapped_trajectory));
            return bootstrapped_trajectory;
        };
      
    private:
        std::vector<int> unmasked_indices;
        std::vector<BinomialAltGenerator> unmasked_alt_generators; // used for generating bootstrapped trajectories
        
};

// Calculates a (log) likelihood ratio test statistic comparing the null model
// (all samples are binomial variables with the same avg) vs the alternative
// (all samples are binomial variables with a sample-specific avg)
inline double calculate_LRT(Trajectory const & trajectory){
    
    
    // cap the avg frequency at avg_frequency threshold
    double p = calculate_error_rate(trajectory);
    
    double log_LRT = 0;
    
    if(p<=0 || p>= 1){
        return log_LRT;
    }
    
    for(auto const & sample : trajectory){
    
        // calculate upper and lower alt thresholds
        // (so that we don't waste statistical power on
        //  changes that we don't value scientifically)
        double expected_alt = sample.depth*p;
        double stddev_alt = std::sqrt(sample.depth*p*(1-p));
        
        double lower_alt_discreteness_threshold = std::floor(expected_alt)-0.5;
        double lower_alt_stddev_threshold = expected_alt-stddev_alt;
        double lower_alt_threshold = std::min(lower_alt_discreteness_threshold, lower_alt_stddev_threshold);
        lower_alt_threshold = std::max(lower_alt_threshold, 0.0);
    
        double upper_alt_discreteness_threshold = std::ceil(expected_alt)+0.5;
        double upper_alt_stddev_threshold = expected_alt+stddev_alt;
        double upper_alt_threshold = std::max(upper_alt_discreteness_threshold, upper_alt_stddev_threshold);
        upper_alt_threshold = std::min(upper_alt_threshold, sample.depth);
    
        if( (sample.alt < lower_alt_threshold) || (sample.alt > upper_alt_threshold) ){
        
            double pi = sample.alt/sample.depth;
            
            if(pi>0){
                log_LRT += sample.depth*pi*std::log( pi/p ); 
            }
            
            if(pi<1){
                log_LRT += sample.depth*(1-pi)*std::log( (1-pi)/(1-p) );
            }
        }
    }
    
    return log_LRT;
}

#endif