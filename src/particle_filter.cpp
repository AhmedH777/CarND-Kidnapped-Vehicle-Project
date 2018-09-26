/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	// Set Number of particles
	num_particles = 50;
	// Set standard deviations for x, y, and theta.
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];


	// These lines creates a normal (Gaussian) distribution.
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	for (int i = 0; i < num_particles; ++i)
	 {
		struct Particle current_particle;

		current_particle.id = i;
		current_particle.x = x + dist_x(gen);
		current_particle.y = y + dist_y(gen);
		current_particle.theta = theta + dist_theta(gen);
		current_particle.weight = 1.0;

		particles.push_back(current_particle);
		weights.push_back(current_particle.weight);
	 }

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    	
	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	double dist, curve, change_rate;

	//buffer values
	dist = velocity * delta_t;
	curve = yaw_rate * delta_t;
	change_rate = velocity/yaw_rate;

	// Set standard deviations for x, y, and theta.
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];


	// These lines creates a normal (Gaussian) distribution.
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	

	for(int i = 0; i < num_particles; i++)
	{
		if (fabs(yaw_rate) < 0.00001)
		{
			particles[i].x 		+= ( dist * cos(particles[i].theta) ) + dist_x(gen);
			particles[i].y 		+= ( dist * sin(particles[i].theta) ) + dist_y(gen);
			particles[i].theta 	+= dist_theta(gen);
		}
		else
		{
			particles[i].x += ( change_rate * ( sin(particles[i].theta + curve) - sin(particles[i].theta) ) ) + dist_x(gen);
			particles[i].y += ( change_rate * ( cos(particles[i].theta) - cos(particles[i].theta + curve) ) ) + dist_y(gen);
			particles[i].theta += curve + dist_theta(gen);
		}
	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	
	
	for(unsigned i = 0; i < observations.size(); i++)
	{
		double min_dist = numeric_limits<double>::max();
		int id = -1;

		for(unsigned j = 0; j < predicted.size(); j++)
		{
			double cur_distance;
			cur_distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			
			if(cur_distance < min_dist)
			{
				min_dist = cur_distance;
				id = predicted[j].id;
			}
		}
		observations[i].id = id;
	}


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// For Each Particle


	for(unsigned index = 0; index < particles.size(); index++)
	{
		struct Particle cur_particle;
		cur_particle = particles[index];

		// First: Transform all Observed Landmarks To Particle Coordinates
		std::vector<LandmarkObs> pred_landmarks;
		std::vector<LandmarkObs> obs_trans_landmarks;

		for(unsigned i = 0; i < observations.size(); i++)
		{
		    LandmarkObs trans_obs;
		    trans_obs.id = observations[i].id;
		    trans_obs.x = observations[i].x * cos(cur_particle.theta) - observations[i].y * sin(cur_particle.theta) + cur_particle.x;
		    trans_obs.y = observations[i].x * sin(cur_particle.theta) + observations[i].y * cos(cur_particle.theta) + cur_particle.y;

		    obs_trans_landmarks.push_back(trans_obs);
		}

		// Second : Find map Landmarks within sensor FOV
		for(unsigned i = 0; i < map_landmarks.landmark_list.size();i++)
		{
			
			if(dist(cur_particle.x,
				cur_particle.y,
				map_landmarks.landmark_list[i].x_f,
				map_landmarks.landmark_list[i].y_f)
			   < sensor_range) 
			
			{
				LandmarkObs temp;
				temp.id = map_landmarks.landmark_list[i].id_i;
				temp.x  = map_landmarks.landmark_list[i].x_f;
				temp.y  = map_landmarks.landmark_list[i].y_f;
				
				pred_landmarks.push_back(temp);

			}
		}


		// Third : Find Associations between Predicted and Observed Landmarks

		dataAssociation(pred_landmarks,obs_trans_landmarks);

		// Fourth : Update Particle weights

		particles[index].weight = 1.0;

		for(unsigned i =0; i < obs_trans_landmarks.size(); i++)
		{
		    double predicted_x, predicted_y, observed_x, observed_y, std_x, std_y;

		    observed_x = obs_trans_landmarks[i].x;
		    observed_y = obs_trans_landmarks[i].y;
		    std_x = std_landmark[0];
		    std_y = std_landmark[1];


		    //get the predicted landmark that is associated with the observed landmark
		    for(unsigned j =0; j < pred_landmarks.size(); j++)
		    {
		        if (pred_landmarks[j].id == obs_trans_landmarks[i].id)
		        {
		            predicted_x = pred_landmarks[j].x;
		            predicted_y = pred_landmarks[j].y;
		        }
		    }

		    //Update the particle weights by multiplying the gaussian multivariate distributions
		    particles[index].weight *= 1.0/(2* M_PI * std_x * std_y) *
		            exp(-((pow(observed_x - predicted_x, 2) / (2 * pow(std_x,2))) +
		                (pow(observed_y - predicted_y, 2) / (2 * pow(std_y,2)))));

        	}

		weights[index] = particles[index].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distributioon

	    default_random_engine gen;

	    //This vector contains the rsampled particles
	    std::vector<Particle> resampled_particles;

	    //choose random index
	    uniform_int_distribution<int> int_dist(0,num_particles-1);

	    int index = int_dist(gen);

	    double beta = 0.0;

	    // extract the maximum weight
	    double max_weight = *max_element(weights.begin(), weights.end());

	    uniform_real_distribution<double> double_dist(0,max_weight);


	    //run the resampling wheel
	    for (int i = 0; i< num_particles; i++)
	    {

		beta = beta + 2.0 * double_dist(gen);

		while (weights[index] < beta)
		{
		    beta = beta - weights[index];
		    index = (index+1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	    }

	    particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
