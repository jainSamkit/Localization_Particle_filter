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
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;
	normal_distribution <double> dist_x(x,std[0]);
	normal_distribution <double> dist_y(y,std[1]);
	normal_distribution <double> dist_theta(theta,std[2]);
	for(int i=0;i<num_particles;i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		weights.push_back(1.0);
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	normal_distribution <double> dist_x(0.0,std_pos[0]);
	normal_distribution <double> dist_y(0.0,std_pos[1]);
	normal_distribution <double> dist_theta(0.0,std_pos[2]);
	for(int i=0;i<num_particles;i++)
	{
		cout<<"samkit"<<endl;
		particles[i].x +=velocity*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta))/yaw_rate;
		particles[i].y +=velocity*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t))/yaw_rate;
		particles[i].theta+=yaw_rate*delta_t;

		particles[i].x +=dist_x(gen);
		particles[i].y +=dist_y(gen);
		particles[i].theta +=dist_theta(gen);
		
	}
	for(int i=0;i<num_particles;i++)
	{
		cout<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].theta<<endl;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(int i = 0;i<observations.size();i++)
	{
		LandmarkObs o = observations[i];
		double dist_ = INT_MAX;
		for(int j=0;j<predicted.size();j++)
		{
			if(dist(o.x,o.y,predicted[j].x,predicted[j].y) < dist_)
			{
				o.id = predicted[j].id;
				dist_ = dist(o.x,o.y,predicted[j].x,predicted[j].y);
			}
		}
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
	vector<LandmarkObs>predictions;
	vector<Map::single_landmark_s>landmark_list_ = map_landmarks.landmark_list;
	double weight_sum = 0;
	for(int i=0;i<num_particles;i++)
	{
		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{

			if(fabs(landmark_list_[j].x_f-particles[i].x)<=sensor_range && fabs(landmark_list_[j].y_f-particles[i].y))
			{
				predictions.push_back(LandmarkObs{landmark_list_[j].id_i,landmark_list_[j].x_f,landmark_list_[j].y_f});
			}
		}	
		double theta_=particles[i].theta;
		double x_ = particles[i].x;
		double y_ = particles[i].y;

		vector<LandmarkObs> transform_obs = observations;
		for(int j=0;j<transform_obs.size();j++)
		{
			transform_obs[j].x = x_ + cos(theta_)*transform_obs[j].x - sin(theta_)*transform_obs[j].y;
			transform_obs[j].y = y_ + sin(theta_)*transform_obs[j].x + cos(theta_)*transform_obs[j].y;
		}
		dataAssociation(predictions,transform_obs);
		for(int j=0;j<transform_obs.size();j++)
		{
			particles[i].associations.push_back(transform_obs[j].id);
			particles[i].sense_x.push_back(transform_obs[j].x);
			particles[i].sense_y.push_back(transform_obs[j].y);
		}
		double weight_prob = 1.0;
		for(int j=0;j<transform_obs.size();j++)
		{
			int count=0;
			for(int k=0;k<predictions.size();j++)
			{
				if(transform_obs[j].id == predictions[k].id)
				{
					double std_x = std_landmark[0];
					double std_y = std_landmark[1];
					double p_x = predictions[k].x;
					double p_y = predictions[k].y;
					double o_x = transform_obs[j].x;
					double o_y = transform_obs[j].y;

					// weight_prob * = get_gauss_prob(std_landmark[0],std_landmark[1],transform_obs[j].x,predictions[k].x,transform_obs[j].y,predictions[k].y);
					weight_prob*=( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(p_x-o_x,2)/(2*pow(std_x, 2)) + (pow(p_y-o_y,2)/(2*pow(std_y, 2))) ) );
					count+=1;
				}
				if(count==1)
					break;
			}
		}
		weights[i] = weight_prob;
		weight_sum +=weight_prob;
		
	}
	

	for(int i=0;i<num_particles;i++)
	{
		weights[i] /=weight_sum;
		particles[i].weight = weights[i];
	}
	for(int i=0;i<weights.size();i++)
	{
		cout<<weights[i]<<" ";
	}
	cout<<endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	uniform_int_distribution<int> uniintdist(0,num_particles-1);
	auto index = uniintdist(gen);

	double max_weight = *max_element(weights.begin(),weights.end());
	uniform_real_distribution<double>uniintdist_(0.0,max_weight);

	double beta = 0.0;
	vector<Particle>rmsp_particles;
	for(int i=0;i<num_particles;i++)
	{
		beta+=uniintdist_(gen);
		while(beta>weights[index])
		{
			beta-=weights[index];
			index=(index+1)%num_particles;
		}
		rmsp_particles.push_back(particles[index]);
	}
	particles = rmsp_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
