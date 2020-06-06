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
#include<stdlib.h>
#include <time.h>
#include "helper_functions.h"
 
using namespace std;
// declare a random engine to be used across multiple and various method calls
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if (is_initialized) {
    return;
}
  num_particles = 100; // TODO: Set the number of particles
  
    // normal distribution x with std_x 
  normal_distribution<double> dist_x(x, std[0]);
   // normal distribution y with std_y
  normal_distribution<double> dist_y(y, std[1]);
    //normal distribution theta with std_theta
  normal_distribution<double> dist_theta(theta, std[2]);
   // =========================================================
  //  Creating particles and the normally distributed state 
  // ==========================================================
  for(int i = 0; i < num_particles; i++){
    Particle particle;
    particle.id = i;
    particle.x=dist_x(gen);
    particle.y=dist_y(gen);
    particle.theta=dist_theta(gen);
    // assign weight to particles 
    particle.weight=1;
    particles.push_back(particle);
    }
   /* set is_initialized to True */
  is_initialized=true;

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
  // normal distribution x with zero mean and each std.
  normal_distribution<double> dist_x(0, std_pos[0]);
  // normal distribution y with std_y
  normal_distribution<double> dist_y(0, std_pos[1]);
  //normal distribution theta with std_theta
  normal_distribution<double> dist_theta(0, std_pos[2]);
   // =========================================
  //  calculate new state
  // =========================================
  for(int i=0;i < num_particles;i++){
    // if yaw_rate is greater then 0
    if(fabs(yaw_rate)>0.00001){
    
      particles[i].x=particles[i].x+(velocity/yaw_rate)*(sin(particles[i].theta+(yaw_rate*delta_t))-sin(particles[i].theta))+dist_x(gen);
      particles[i].y=particles[i].y+(velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+(yaw_rate*delta_t)))+dist_y(gen);      
      particles[i].theta = particles[i].theta + (yaw_rate * delta_t);
  
    }
    // if yaw_rate is zero
    else{

      particles[i].x=particles[i].x+velocity*delta_t*cos(particles[i].theta);
      particles[i].y=particles[i].y+velocity*delta_t*sin(particles[i].theta);
       }
    }

}
// 
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
  for(unsigned int i=0;i < observations.size();i++){
   
     //for each observation, initialize the min distance to a big number
    double mindistance = numeric_limits<double>::max();
    
  int id_in_map=-1;
    for(unsigned int j=0;j<predicted.size();j++){
     //distance calculation with helper function
    double distance=dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
      
      // if distance is smaller than the distance, then save the id and iterate all the predicted value
      //finally find the nearest precited to GT value. 
      if(distance<mindistance){
      mindistance=distance;
      id_in_map= predicted[j].id; 
      }
    }
      //assign observed measurement to particular landmark.
    observations[i].id=id_in_map;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
 // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
 //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
 // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
 //   according to the MAP'S coordinate system. You will need to transform between the two systems.
 //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
 //   and the following is a good resource for the actual equation to implement (look at equation 
 //   http://planning.cs.uiuc.edu/node99.html

  for (int i = 0; i < num_particles; i++) {

    // get particle x, y coordinates
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;

    // create a vector to hold map landmark locations 
    vector<LandmarkObs> predictions;

    // for each map landmark...
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      // get id and x,y coordinates
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;
      
      double disinrange=dist(lm_x,lm_y,p_x,p_y);
      
      //find landmarks in particle's range
      if(disinrange<sensor_range){
 
        LandmarkObs temp;
        temp.x=lm_x;
        temp.y=lm_y;
        temp.id=lm_id;
        // add prediction to vector
        predictions.push_back(temp);
      }
    }

    // transform the list of observations from vehicle coordinates to map coordinates
    vector<LandmarkObs> transformed_os;
    for (unsigned int j = 0; j < observations.size(); j++) {
      double transx = p_x + cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y;
      double transy = p_y + sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y;
      LandmarkObs temp;
      temp.x=transx;
      temp.y=transy;
      temp.id=observations[j].id;
      transformed_os.push_back(temp);
      
    }

    
     dataAssociation(predictions, transformed_os);

    // init weight
    particles[i].weight = 1.0;

    for (unsigned int j = 0; j < transformed_os.size(); j++) {
      
      //declaring variables 
      double o_x, o_y, pr_x, pr_y;
      o_x = transformed_os[j].x;
      o_y = transformed_os[j].y;

      int associated_prediction = transformed_os[j].id;

      // get the x,y coordinates of prediction associated with current observation
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (predictions[k].id == associated_prediction) {
          pr_x = predictions[k].x;
          pr_y = predictions[k].y;
        }
      }

      
      double stdLandmarkRange = std_landmark[0];
      double stdLandmarkBearing = std_landmark[1];
      double dX = pr_x-o_x;
      double dY = pr_y-o_y;
      // calculate weight for observation with multivariate Gaussian
      double obs_w = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp( -( pow(dX,2)/(2*pow(stdLandmarkRange, 2)) + (pow(dY,2)/(2*pow(stdLandmarkBearing, 2))) ) );

      // product of obersvation weight with total observations weight
      particles[i].weight *= obs_w;
    }
  }
} 
  
  

void ParticleFilter::resample() {
  
  vector<Particle> new_particles;
  vector<double> weights;
  
  double maxweight=0;
  double beta=0;

  for(int i=0;i < num_particles;i++){
    weights.push_back(particles[i].weight);
    if(particles[i].weight>maxweight){
    maxweight=particles[i].weight;
    }  
  }
  
    // Creating distributions.
  uniform_real_distribution<float> unirealdist(0.0, maxweight);
  uniform_int_distribution<int> uniintdist(0, num_particles - 1);
  
  auto index = uniintdist(gen);
  
  for(int j=0;j < num_particles;j++){
    beta += unirealdist(gen) * 2.0;
    while(beta > weights[index]){
    beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
    
  }
  particles =new_particles;
  

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
    //Clear the previous associations

  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

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