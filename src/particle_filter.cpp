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
#include <cfloat>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  cout << "x: " << x << ", y = " << y << "theta = " << theta << endl;

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	default_random_engine gen;

	particles.clear();
	num_particles = 100;

	for(int i = 0; i < num_particles; i++) {
		Particle p;
    p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
    weights.push_back(1.0);
	}

	is_initialized = true;

	cout << "Created # of particles: " << particles.size() << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	this->predictionCount++;
	cout << "Prediction: " << predictionCount << endl;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

  default_random_engine gen;

	double change_factor = velocity / yaw_rate;

	if (yaw_rate == 0) {

		for (int index = 0; index < particles.size(); index++) {

			Particle &p = particles[index];
			p.x = p.x + velocity * delta_t * cos(p.theta);
			p.y = p.y + velocity * delta_t * sin(p.theta);

		}
	} else {
		for (int index = 0; index < particles.size(); index++) {
			Particle &p = particles[index];
			p.x = p.x + change_factor * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)) + dist_x(gen);
			p.y = p.y + change_factor * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t)) + dist_y(gen);
			p.theta = p.theta + yaw_rate * delta_t + dist_theta(gen);
		}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//Implement nearest neighbour association
	//iterate over landmarks and find closest measurement


	for(int observationsIndex = 0; observationsIndex < observations.size(); observationsIndex++) {
		auto &observation = observations[observationsIndex];
		double closest = DBL_MAX;

		for(int index = 0; index < predicted.size(); index ++) {
		  auto landmark = predicted[index];

			double observation_dist = dist(landmark.x, landmark.y, observation.x, observation.y);

			if(observation_dist < closest) {
				closest = observation_dist;
				observation.id = landmark.id; //Is this right or maybe just create map
				//observations[observationsIndex] = observation;
			}

		}

	}

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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for(int index = 0; index < particles.size(); index++) {

		Particle &p = particles[index];
		//auto landmarks_in_range = landmarksInRange(p, sensor_range, map_landmarks);

		std::vector<LandmarkObs> mappedObservations;
		mappedObservations = mapObservationsToWorldCoordinates(observations, p);


		dataAssociation(predictObservations(map_landmarks,p, mappedObservations.size()), mappedObservations);
		p.weight = 1.0;
		for(int observationIndex = 0; observationIndex< mappedObservations.size(); observationIndex++) {
			p.weight *= measurementProbability(mappedObservations[observationIndex], map_landmarks, std_landmark);
		}
		this->weights[index] = p.weight;

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());

  vector<Particle> resampled;
  for(int index = 0; index < this->particles.size(); index++) {
    resampled.push_back(particles[distribution(gen)]);
  }

  particles = resampled;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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

double ParticleFilter::measurementProbability(LandmarkObs observation, Map map, double stdLandmark[]) {

	double max_probability = 0.0;
	for(int mapLandmarkIndex = 0; mapLandmarkIndex < map.landmark_list.size(); mapLandmarkIndex++) {
		auto landMark = map.landmark_list[mapLandmarkIndex];
		if (landMark.id_i == observation.id) {
			auto probability = gaussian(observation, landMark, stdLandmark);
			//cout << "Found matching landmark for observation " << observation.id << "with prob: " << probability << endl;
			return gaussian(observation, landMark, stdLandmark);
		}
	}
	//cout << "No matching landmark for observation " << observation.id << endl;
	return 0.0; //Should not come here?
}

double
ParticleFilter::gaussian(LandmarkObs observation, Map::single_landmark_s landmark, double std_dev[]) {
	double exponent1 = pow(observation.x-landmark.x_f, 2) / (2*std_dev[0]*std_dev[0]);
	double exponent2 = pow(observation.y-landmark.y_f, 2) / (2*std_dev[1]*std_dev[1]);
	return exp(-(exponent1 + exponent2)) / (2 * M_PI * std_dev[0] * std_dev[1]);
}

LandmarkObs ParticleFilter::mapToWorldCoordinates(LandmarkObs &obs, Particle &particle) {
	auto mapped = LandmarkObs();
	mapped.x = particle.x + (obs.x * cos(particle.theta) - obs.y * sin(particle.theta));
	mapped.y = particle.y + (obs.x * sin(particle.theta) +  obs.y * cos(particle.theta));
	mapped.id = obs.id;
	return mapped;
}

std::vector<LandmarkObs>
ParticleFilter::mapObservationsToWorldCoordinates(std::vector<LandmarkObs> observations,
																									Particle particle) {
	auto mapped = std::vector<LandmarkObs>();
	for(int index = 0; index < observations.size(); index ++) {
		mapped.push_back(mapToWorldCoordinates(observations[index], particle));
	}

	return mapped;
}

std::vector<LandmarkObs, std::allocator<LandmarkObs>>
ParticleFilter::landmarksInRange(Particle particle, double range, Map map) {
	auto inRange =  vector<LandmarkObs, allocator<LandmarkObs>>();
	for(int index = 0; index < map.landmark_list.size(); index ++) {
		auto landmark = map.landmark_list[index];
		if(dist(particle.x, particle.y, landmark.x_f, landmark.y_f) <= range || true) {
			LandmarkObs landmarkInRange;
			landmarkInRange.x = landmark.x_f;
			landmarkInRange.y = landmark.y_f;
			landmarkInRange.id = landmark.id_i;
			inRange.push_back(landmarkInRange);
		}
	}
	return inRange;
}

std::vector<LandmarkObs, std::allocator<LandmarkObs>>
ParticleFilter::predictObservations(Map map, Particle &particle, unsigned long size) {


	auto predicted = vector<LandmarkObs, allocator<LandmarkObs>>() ;
	for(int index = 0; index < map.landmark_list.size(); index++) {
		LandmarkObs landmarkObs = LandmarkObs();
		landmarkObs.x = map.landmark_list[index].x_f;
		landmarkObs.y = map.landmark_list[index].y_f;
		landmarkObs.id = map.landmark_list[index].id_i;
		predicted.push_back(landmarkObs);
	}

	/*
	Particle p = particle;
	struct sorter
	{
			inline bool operator() (const LandmarkObs& landmarkObs1, const LandmarkObs& landmarkObs2)
			{
				return (dist(p.x, p.y, landmarkObs1.x, landmarkObs1.y) < dist(p.x, p.y, landmarkObs2.x, landmarkObs2.y));
			}
	};

	sort(predicted.begin(), predicted.end(), sorter());
	*/
	sort( predicted.begin( ), predicted.end( ), [ particle]( const LandmarkObs& landmarkObs1, const LandmarkObs& landmarkObs2 )
	{
			return (dist(particle.x, particle.y, landmarkObs1.x, landmarkObs1.y) < dist(particle.x, particle.y, landmarkObs2.x, landmarkObs2.y));
	});

	predicted.resize(size);
	return predicted;

}
