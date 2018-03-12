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
#include <map>
#include <assert.h>
#include "particle_filter.h"
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// Gaussian noise generator
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i<num_particles;i++)
	{
		// Sample from Gaussian distributions
		double sample_x = dist_x(gen);
        double sample_y = dist_y(gen);
        double sample_theta = dist_theta(gen);
		double weight = 1.0;		
		Particle new_particle{i,sample_x, sample_y, sample_theta, weight};
		particles.push_back(new_particle);
	}

	//for (auto particle:particles)
	//	cout<<particle.id<<','<<particle.x<<','<<particle.y<<','<<particle.theta<<','<<particle.weight<<endl;
	
	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
		// Gaussian noise generator
	default_random_engine gen;
	normal_distribution<double> noise_x(0, std_pos[0]);
	normal_distribution<double> noise_y(0, std_pos[1]);
	normal_distribution<double> noise_theta(0, std_pos[2]);

	double delta_x, delta_y, delta_theta;
	for (auto &p:particles)
	{
		
		// Calculate prediction
		if (fabs(yaw_rate)<0.000001)
		{
            delta_x=velocity*cos(p.theta)*delta_t;
            delta_y=velocity*sin(p.theta)*delta_t;
          
		}
		else
		{
			delta_x=(velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
			delta_y=(velocity/yaw_rate)*(-cos(p.theta+yaw_rate*delta_t)+cos(p.theta));
		}

		delta_theta=yaw_rate*delta_t;

		p.x = p.x + delta_x + noise_x(gen);
		p.y = p.y + delta_y + noise_y(gen);
		p.theta = p.theta + delta_theta + noise_theta(gen);

	}	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	

	weights.clear();
				
	for (auto &particle:particles)
	{	
		// convert observations vehicle coordinate to map coordinate	
		vector<LandmarkObs> observations_remapped; //observation in map coorindates
		for (auto obs:observations)
		{
			LandmarkObs obs_remapped;
			obs_remapped.id=-1;
			obs_remapped.x = particle.x + obs.x * cos(particle.theta) - obs.y * sin(particle.theta);
			obs_remapped.y = particle.y + obs.x * sin(particle.theta) + obs.y * cos(particle.theta);
			observations_remapped.push_back(obs_remapped);
		}

		// use only the landmark that is within sensor range
		vector<Map::single_landmark_s> landmark_in_range;
		for (auto &landmark:map_landmarks.landmark_list)
		{
			if (dist(landmark.x_f, landmark.y_f, particle.x, particle.y)<=sensor_range)
			{
				landmark_in_range.push_back(landmark);
			}
		} 

		// Data association. Match each observation to a landmark. This involves 2 rounds of nearest-neighbor algorithm
		// 1. shortlist potential landmark
		vector<int> maybe_ids;
		int maybe_id;
		for (auto &obs:observations_remapped)
		{
			double minDist=10000;
			double distance;
			for (auto &landmark:landmark_in_range)
			{
				distance=dist(obs.x, obs.y, landmark.x_f, landmark.y_f);
				if (distance<minDist)
				{
					minDist=distance;
					maybe_id=landmark.id_i;
				}
			}
			maybe_ids.push_back(maybe_id);
		}

		// remove duplicate in potential landmark id
		sort(maybe_ids.begin(), maybe_ids.end());
		maybe_ids.erase( unique( maybe_ids.begin(), maybe_ids.end() ), maybe_ids.end());

		// 2. match landmark candidate to only one observation
		for (auto id_i:maybe_ids)
		{
			// assume landmark id increment from 1, otherwise will use index to landmark_list instead of id_i
			Map::single_landmark_s landmark=map_landmarks.landmark_list[id_i-1];
			double minDist=10000;
			double distance;			
			int obs_id;
			for (int i=0; i<observations_remapped.size(); i++)
			{
				distance=dist(observations_remapped[i].x, observations_remapped[i].y, landmark.x_f, landmark.y_f);
				if (distance<minDist)
				{
					minDist=distance;	
					obs_id=i;
				}
			}	
			observations_remapped[obs_id].id=id_i;		
		}

		// Calculate probabilities
		particle.weight=1.0;
		for (auto &obs:observations_remapped)
		{
			if (obs.id==-1) continue; //ignore if this observation has not been matched with a landmark

			float x=obs.x;
			float y=obs.y;
			float ux=map_landmarks.landmark_list[obs.id-1].x_f;
			float uy=map_landmarks.landmark_list[obs.id-1].y_f;

			// multivariate Gaussian pdf
			double xterm=(x-ux)*(x-ux)/(2*std_landmark[0]*std_landmark[0]);
			double yterm=(y-uy)*(y-uy)/(2*std_landmark[1]*std_landmark[1]);
			double pdf=(1/(2*M_PI*std_landmark[0]*std_landmark[1]))*exp(-(xterm+yterm));
			particle.weight*=pdf;		

		}
		weights.push_back(particle.weight);

	}
	
	// normalize weights
	double sum_weights=0.0;
	for (auto w:weights)
	{
		sum_weights+=w;
	}

	for (auto &w:weights)
	{
		w/=sum_weights;
	}	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//cout<<"resample"<<endl;
	// set up the random bit
    std::random_device rd;
    std::mt19937 gen(rd());

    std::discrete_distribution<> d(weights.begin(), weights.end());
	std::vector<Particle> resampled_particles;
    for(int n=0; n<num_particles; ++n) {
		int sampled_idx=d(gen);
		resampled_particles.push_back(particles[sampled_idx]);
    }

	particles=move(resampled_particles);

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].id <<" " << particles[i].x << " " << particles[i].y << " " << particles[i].theta <<"\n";
	}
	dataFile.close();
}
