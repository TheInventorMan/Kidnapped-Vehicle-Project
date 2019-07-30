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
#include <map>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::cout;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
        /**
         * TODO: Set the number of particles. Initialize all particles to
         *   first position (based on estimates of x, y, theta and their uncertainties
         *   from GPS) and all weights to 1.
         * TODO: Add random Gaussian noise to each particle.
         * NOTE: Consult particle_filter.h for more information about this method
         *   (and others in this file).
         */
        //cout << "begin init" << std::endl;
        num_particles = 100; // TODO: Set the number of particles
        particles.resize(num_particles);
        weights.resize(num_particles);

        std::default_random_engine gen;
        double std_x = std[0];
        double std_y = std[1];
        double std_theta = std[2];

        // This line creates a normal (Gaussian) distribution for x
        normal_distribution<double> dist_x(x, std_x);
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);

        for (int i = 0; i < particles.size(); i++) {
                Particle new_particle;

                new_particle.id = i;
                new_particle.x = x + dist_x(gen);
                new_particle.y = y + dist_y(gen);
                new_particle.theta = theta + dist_theta(gen);
                new_particle.weight = 1.0;

                particles[i] = new_particle;
                weights[i] = 1.0;
        }
        is_initialized = true;
        //cout << "end init" << std::endl;
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
        //cout << "begin prediction" << std::endl;
        std::default_random_engine gen;
        double x0, y0, t0, x1, y1, t1;

        normal_distribution<double> dist_x(0, std_pos[0]);
        normal_distribution<double> dist_y(0, std_pos[1]);
        normal_distribution<double> dist_theta(0, std_pos[2]);

        for (int i = 0; i < particles.size(); i++) {
                x0 = particles[i].x;
                y0 = particles[i].y;
                t0 = particles[i].theta;

                //check if yaw_rate not close to zero
                if(fabs(yaw_rate) > 0.0001) {
                        x1 = x0 + (velocity/yaw_rate)*(sin(t0 + yaw_rate*delta_t)-sin(t0));
                        y1 = y0 + (velocity/yaw_rate)*(cos(t0)-cos(t0 + yaw_rate*delta_t));
                        t1 = t0 + yaw_rate*delta_t;
                } else {
                        x1 = x0 + velocity*delta_t*cos(t0);
                        y1 = y0 + velocity*delta_t*sin(t0);
                        t1 = t0;
                }

                particles[i].x = x1 + dist_x(gen);
                particles[i].y = y1 + dist_y(gen);
                particles[i].theta = t1 + dist_theta(gen);

        }
        //cout << "end prediction" << std::endl;
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
        //cout << "begin updateWeights" << std::endl;
        double std_x = std_landmark[0];
        double std_y = std_landmark[1];

        unsigned int i;
        double part_x, part_y, part_th;

        unsigned int o;
        unsigned int t;
        double obs_x, obs_y, new_obs_x, new_obs_y, lowest_dist;
        double mu_x, mu_y, temp_w, mvx, mvy;
        int nearest_assoc;

        unsigned int l;
        double land_x, land_y, dist_to_land, nearest_land_x, nearest_land_y;

        for(i = 0; i < particles.size(); i++) {
                //for each particle, store particle coordinates
                //cout << i << std::endl;
                vector<int> new_associations;
                vector<double> new_sense_x;
                vector<double> new_sense_y;
                vector<LandmarkObs> transformed;

                part_x = particles[i].x;
                part_y = particles[i].y;
                part_th = particles[i].theta;
                particles[i].weight = 1.0;
                weights[i] = 1.0;
                //cout << i << std::endl;
                for(o = 0; o < observations.size(); o++) {
                        LandmarkObs new_obs;
                        obs_x = observations[o].x;
                        obs_y = observations[o].y;
                        new_obs.x = part_x + cos(part_th)*obs_x - sin(part_th)*obs_y;
                        new_obs.y = part_y + sin(part_th)*obs_x + cos(part_th)*obs_y;
                        transformed.push_back(new_obs);
                }

                lowest_dist = 1000000.0;
                for(t = 0; t < transformed.size(); t++) {
                        LandmarkObs t_obs = transformed[t];
                        //cout << map_landmarks.landmark_list.size() << std::endl;
                        for(l = 0; l < map_landmarks.landmark_list.size(); l++) {
                                land_x = map_landmarks.landmark_list[l].x_f;
                                land_y = map_landmarks.landmark_list[l].y_f;
                                //dist_to_land = dist(new_obs_x, new_obs_y, land_x, land_y);
                                dist_to_land = sqrt(pow((land_y-t_obs.y),2.0) + pow((land_x-t_obs.x),2.0));
                                //cout << dist_to_land << std::endl;
                                //if (dist_to_land <= sensor_range) {
                                if(dist_to_land < lowest_dist) {
                                        lowest_dist = dist_to_land;
                                        nearest_land_x = land_x;
                                        nearest_land_y = land_y;
                                        nearest_assoc = map_landmarks.landmark_list[l].id_i;
                                        //cout << particles[i].x << " " << land_x << " " << particles[i].y << " " << land_y << std::endl;
                                }
                                //}
                        }
                        //cout << obs_x << " " << obs_y << " " << lowest_dist << std::endl;
                        mu_x = nearest_land_x;
                        mu_y = nearest_land_y;
                        //cout << new_obs_x << " " << mu_x << " " << new_obs_y << " " << mu_y << std::endl;
                        //cout << new_obs_x << " " << nearest_land_x << " " << new_obs_y << " " << nearest_land_y << std::endl;

                        temp_w = 1/(2 * M_PI * std_x * std_y);
                        mvx = pow(t_obs.x-mu_x, 2)/(2 * pow(std_x, 2));
                        mvy = pow(t_obs.y-mu_y, 2)/(2 * pow(std_y, 2));
                        temp_w *= exp(-(mvx+mvy));

                        //cout << temp_w << std::endl;
                        if(temp_w > 0) {
                                particles[i].weight *= temp_w;
                                weights[i] *= temp_w;
                        }
                        new_associations.push_back(nearest_assoc);
                        new_sense_x.push_back(new_obs_x);
                        new_sense_y.push_back(new_obs_y);
                        //if(lowest_dist == 1000000.0){
                        //	particles[i].weight = 0.0;
                        //    weights[i] = 0.0;
                        //}

                }
                particles[i].associations = new_associations;  //push_back(nearest_assoc);
                particles[i].sense_x = new_sense_x; //push_back(new_obs_x);
                particles[i].sense_y = new_sense_y; //push_back(new_obs_y);
                //cout << "particle " << i << " weight updated" << std::endl;
        }
        //cout << "end updateWeights" << std::endl;
}

void ParticleFilter::resample() { //revisit
        /**
         * TODO: Resample particles with replacement with probability proportional
         *   to their weight.
         * NOTE: You may find std::discrete_distribution helpful here.
         *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
         */
        //cout << "begin resample" << std::endl;
        //get maximum weight
        double max_w = 0.0;
        for(unsigned int i = 0; i < particles.size(); i++) {
                weights[i] = particles[i].weight;
                //cout << weights[i] << std::endl;
                if (particles[i].weight > max_w) {
                        max_w = particles[i].weight;
                }
        }

        //cout << "check1" << std::endl;
        vector<Particle> resampled;

        std::default_random_engine gen;
        std::uniform_int_distribution<int> idx_distro(0, num_particles-1);
        std::uniform_real_distribution<double> w_distro(0.0, max_w);

        int idx = idx_distro(gen);
        double beta = 0.0;

        for(int i = 0; i < num_particles; i++) {
                beta += 2.0*w_distro(gen);
                while(beta > weights[idx]) {
                        beta -= weights[idx];
                        idx = (idx + 1) % num_particles;
                }
                resampled.push_back(particles[idx]);
        }
        particles = resampled;
        //cout << "end resample" << std::endl;
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
        particle.associations = associations;
        particle.sense_x = sense_x;
        particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
        vector<int> v = best.associations;
        std::stringstream ss;
        copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
        string s = ss.str();
        s = s.substr(0, s.length()-1); // get rid of the trailing space
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
        s = s.substr(0, s.length()-1); // get rid of the trailing space
        return s;
}
