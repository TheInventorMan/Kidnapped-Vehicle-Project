/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 500;  // TODO: Set the number of particles

  std::default_random_engine gen;
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

 // This line creates a normal (Gaussian) distribution for x
 normal_distribution<double> dist_x(x, std_x);
 normal_distribution<double> dist_y(y, std_y);
 normal_distribution<double> dist_theta(theta, std_theta);

 for (int i = 0; i < num_particles; i++) {
   Particle new_particle;

   new_particle.id = i;
   new_particle.x = x + dist_x(gen);
   new_particle.y = y + dist_y(gen);
   new_particle.theta = theta + dist_theta(gen);
   new_particle.weight = 1.0;

   particles.push_back(new_particle);
 }
 is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   std::default_random_engine gen;

   for (int i = 0; i < num_particles; i++){
     double x0 = particles[i].x;
     double y0 = particles[i].y;
     double t0 = particles[i].theta;
     double x1, y1, t1;

     //check if yaw_rate close to zero
     if(fabs(yaw_rate) < 0.0001){
       x1 = x0 + velocity*delta_t*cos(t0);
       y1 = y0 + velocity*delta_t*sin(t0);
       t1 = t0;
     }else{
       x1 = x0 + (velocity/yaw_rate)*(sin(t0 + yaw_rate*delta_t)-sin(t0));
       y1 = y0 + (velocity/yaw_rate)*(cos(t0)-cos(t0 + yaw_rate*delta_t));
       t1 = t0 + yaw_rate*delta_t;
     }
     normal_distribution<double> dist_x(x1, std_pos[0]);
     normal_distribution<double> dist_y(y1, std_pos[1]);
     normal_distribution<double> dist_theta(t1, std_pos[2]);
     
     particles[i].x = dist_x(gen);
     particles[i].y = dist_y(gen);
     particles[i].theta = dist_theta(gen);
     
   }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
   double lowest_dist;
   double obs_x, obs_y;
   double pred_x, pred_y;
   double distance;

   for(int o = 0; o < observations.size(); o++){
     lowest_dist = 10000;
     obs_x = observations[o].x;
     obs_y = observations[o].y;

      for(int p = 0; p < predicted.size(); p++){
          pred_x = predicted[p].x;
          pred_y = predicted[p].y;
          distance = dist(obs_x,obs_y,pred_x,pred_y);
          if (distance < lowest_dist){
            lowest_dist = distance;
            observations[o].id = predicted[p].id;
          }
      }
   }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   double std_x = std_landmark[0];
   double std_y = std_landmark[1];

   double part_x, part_y, part_t;

   for(int i = 0; i < num_particles; i++){
   //for each particle, store particle coordinates
      part_x = particles[i].x;
      part_y = particles[i].y;
      part_t = particles[i].theta;

      particles[i].weight = 1.0;

      //for each landmark, check distance with sensor_range, add to new predictions vector
      vector<LandmarkObs> predictions;
      for(int l = 0; l < map_landmarks.landmark_list.size(); l++){
          double land_x = map_landmarks.landmark_list[l].x_f;
          double land_y = map_landmarks.landmark_list[l].y_f;

          if (dist(part_x, part_y, land_x, land_y) <= sensor_range){
            LandmarkObs new_pred;
            new_pred.id = map_landmarks.landmark_list[l].id_i;
            new_pred.x = land_x;
            new_pred.y = land_y;
            predictions.push_back(new_pred);
          }
      }

      //transform obs into vehicle frame
     double obs_x, obs_y;
     vector<LandmarkObs> transformed;
      for(int o = 0; o < observations.size(); o++){
          obs_x = observations[o].x;
          obs_y = observations[o].y;
		  LandmarkObs new_obs;
          new_obs.id = o;
          new_obs.x = part_x + cos(part_t)*obs_x - sin(part_t)*obs_y;
          new_obs.y = part_y + sin(part_t)*obs_x + cos(part_t)*obs_y;
        
          transformed.push_back(new_obs);
      }
    //for each observation, set to transformed observation
    dataAssociation(predictions, transformed);
    //for each observation, store x,y,id
    double mu_x, mu_y, temp_w;
    for(int o = 0; o < observations.size(); o++){
      for(int p = 0; p < predictions.size(); p++){
        //loop through predictions, get prediction x,y if obs.id matches
        if(observations[o].id == predictions[p].id){
          mu_x = predictions[p].x;
          mu_y = predictions[p].y;
          obs_x = observations[o].x;
          obs_y = observations[o].y;

          temp_w = 1/(2 * M_PI * std_x * std_y);
          temp_w *= exp(-pow(obs_x-mu_x, 2)/(2 * pow(std_x, 2)) - (pow(obs_y-mu_y, 2)/(2 * pow(std_y, 2))));

          particles[i].weight *= temp_w;
        }
      }
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   double max_w = 0.0;
   for(int i = 0; i < particles.size(); i++){
     weights.push_back(particles[i].weight);
     if (particles[i].weight > max_w){
       max_w = particles[i].weight;
     }
   }
   vector<Particle> resampled;
   std::default_random_engine gen;
   std::discrete_distribution<int> w_distro(weights.begin(), weights.end());
   int idx;
   //std::discrete_distribution<int> idx_distro(0, num_particles-1);

   for(int i = 0; i < num_particles; i++){
     idx = w_distro(gen);
     resampled.push_back(resampled[idx]);
   }
  particles = resampled;

}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
