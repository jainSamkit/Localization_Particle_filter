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
	num_particles = 10;
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
	
	default_random_engine gen;
	normal_distribution <double> dist_x(0.0,std_pos[0]);
	normal_distribution <double> dist_y(0.0,std_pos[1]);
	normal_distribution <double> dist_theta(0.0,std_pos[2]);
	for(int i=0;i<num_particles;i++)
	{
		// default_random_engine gen;
		// normal_distribution <double> dist_x(0.0,std_pos[0]);
		// normal_distribution <double> dist_y(0.0,std_pos[1]);
		// normal_distribution <double> dist_theta(0.0,std_pos[2]);
		double thresh = 0.001;
		double pred_x=0.0,pred_y=0.0,pred_theta=0.0;
		Particle p_ = particles[i];
		if(fabs(yaw_rate)<=thresh)
		{
			// particles[i].x +=velocity*cos(particles[i].theta)*delta_t;
			// particles[i].y +=velocity*sin(particles[i].theta)*delta_t;
			pred_x = p_.x + velocity*cos(p_.theta)*delta_t;
			pred_y = p_.y + velocity*sin(p_.theta)*delta_t;
			pred_theta = p_.theta;
		}
		else
		{
			pred_x = p_.x + (velocity/yaw_rate)*(sin(p_.theta+yaw_rate*delta_t) - sin(p_.theta));
			pred_y = p_.y + (velocity/yaw_rate)*(cos(p_.theta) - cos(p_.theta+yaw_rate*delta_t));
			pred_theta = p_.theta + yaw_rate * delta_t;
		}
		
		particles[i].x = pred_x + dist_x(gen);
		particles[i].y = pred_y + dist_y(gen);
		particles[i].theta = pred_theta + dist_theta(gen);
		
	}
	// for(int i=0;i<num_particles;i++)
	// {
	// 	cout<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].theta<<endl;
	// }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(int i = 0;i<observations.size();i++)
	{
		LandmarkObs o = observations[i];
		double dist_ = 99999.0;
		for(int j=0;j<predicted.size();j++)
		{
			if(dist(o.x,o.y,predicted[j].x,predicted[j].y) < dist_)
			{
				observations[i].id = j;
				dist_ = dist(o.x,o.y,predicted[j].x,predicted[j].y);
			}
		}
		// cout<<"nearby to"<<" "<<i<<"st"<<" "<<"observation"<<" "<< "is"<<" "<<predicted[observations[i].id].x<<" "<<predicted[observations[i].id].y<<endl;
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
	double weight_sum = 0.0;
	for(int i=0;i<num_particles;i++)
	{
		vector<LandmarkObs>predictions;
		vector<Map::single_landmark_s>landmark_list_ = map_landmarks.landmark_list;

		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{
			double err_x = landmark_list_[j].x_f-particles[i].x;
			double err_y = landmark_list_[j].y_f-particles[i].y;
			double error_ = sqrt(err_x*err_x + err_y*err_y);
			if(fabs(error_)<=sensor_range)
			{
				predictions.push_back(LandmarkObs{landmark_list_[j].id_i,landmark_list_[j].x_f,landmark_list_[j].y_f});
			}
		}	
		// cout<<"Size of prediction"<<" "<<predictions.size()<<endl;
		// cout<<"Size of landmarks list"<<" "<<landmark_list_.size()<<endl;
		// cout<<map_landmarks.landmark_list.size()-prediction_count<<endl;
		double theta_=particles[i].theta;
		double x_ = particles[i].x;
		double y_ = particles[i].y;

		vector<LandmarkObs> transform_obs = observations;
		for(int j=0;j<transform_obs.size();j++)
		{
			// cout<<"samkit"<<endl;
			transform_obs[j].x = x_ + cos(theta_)*observations[j].x - sin(theta_)*observations[j].y;
			transform_obs[j].y = y_ + sin(theta_)*observations[j].x + cos(theta_)*observations[j].y;

		}

		// for(int j=0;j<observations.size();j++)
		// {
		// 	cout<<"local"<<" "<<observations[j].x<<" "<<observations[j].y<<endl;
		// 	cout<<"global"<<" "<<transform_obs[j].x<<" "<<transform_obs[j].y<<endl;
		// }
		// cout<<"Tranformed"<<endl;
		// for(int j=0;j<transform_obs.size();j++)
		// {
		// 	cout<<transform_obs[j].x<<" "<<transform_obs[j].y<<endl;
		// }
		
		// cout<<"Before dataAssociation"<<endl;
		dataAssociation(predictions,transform_obs);
		// cout<<"After dataAssociation"<<endl;
		for(int j=0;j<transform_obs.size();j++)
		{
			particles[i].associations.push_back(predictions[transform_obs[j].id].id);
			// particles[i].associations.push_back(transform_obs[j].id);
			particles[i].sense_x.push_back(transform_obs[j].x);
			particles[i].sense_y.push_back(transform_obs[j].y);
		}
		double weight_prob = 1.0;
		// cout<<"Before weight_prob"<<endl;
		// cout<<transform_obs.size()<<endl;
		for(int j=0;j<transform_obs.size();j++)
		{
			// cout<<"Inside"<<endl;
			int id_ = transform_obs[j].id;
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];
			double p_x = predictions[id_].x;
			double p_y = predictions[id_].y;
			double o_x = transform_obs[j].x;
			double o_y = transform_obs[j].y;
			// default_random_engine gen;
			// normal_distribution <int> prob(0.0,1.0);

			weight_prob*=( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(p_x-o_x,2)/(2*pow(std_x, 2)) + (pow(p_y-o_y,2)/(2*pow(std_y, 2))) ) );
			// weight_prob*=prob(gen);
			// int count=0;
		}
		weights[i] = weight_prob;
		particles[i].weight = weights[i];
		weight_sum +=weight_prob;
	}
	// cout<<"After weight_prob"<<endl;

	for(int i=0;i<num_particles;i++)
	{
		weights[i] /=weight_sum;
		particles[i].weight = weights[i];
	}
	// cout<<"Printing weights"<<endl;
	// for(int i=0;i<weights.size();i++)
	// {
	// 	cout<<weights[i]<<" ";
	// }
	// cout<<endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	uniform_int_distribution<int> uniintdist(0,num_particles-1);
	auto index = uniintdist(gen);

	double max_weight = *max_element(weights.begin(),weights.end());
	uniform_real_distribution<double>uniintdist_(0.0,2*max_weight);

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
	cout<<"Printing resampled weights"<<endl;
	// for(int i=0;i<num_particles;i++)
	// {
	// 	cout<<particles[i].weight<<" ";
	// }
	// cout<<endl;
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
