//Joel Jaffe February 2024
//Based on https://www.youtube.com/watch?v=xiUpAeos168&list=PLZ1w5M-dmhlGWtqzaC2aSLfQFtp0Dz-F_&index=3
//Programming Chaos on YouTube

#include "al/app/al_App.hpp"
#include "al/system/al_Time.hpp"
#include "al/math/al_Random.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/app/al_DistributedApp.hpp"
#include "al/app/al_OmniRendererDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"
//#include "al/scene/al_DistributedScene.hpp"
//#include "al_ext/statedistribution/al_Serialize.hpp"
using namespace al;

#include <iostream>
#include <string>
using namespace std;

float fMap(float value, float in_min, float in_max, float out_min, float out_max) { // custom mapping function
  return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

Vec3f randomVec3f(float scale) { // <- Function that returns a Vec2f containing random coords
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} 

struct Particle { // Particle struct
  int type;
  Vec3f position; 
  Vec3f velocity;
};

struct State {
  // state member variables
  float pointSize;
  Pose camera;
  Parameter simScale{"/simScale", "", 0.5f, 0.f, 1.f}; // <- creates GUI parameter
  Parameter springConstant{"/springConstant", "", 0.4, 0.0, 1.0}; // <- creates GUI parameter

  static const int numTypes = 6; // numTypes
  int numParticles = 1000; // numParticles (1000 seems to be the limit for my M2 Max)
  float colorStep = 1.f / numTypes; // colorStrep
  float K = 0.05; // make smaller to slow sim (0.05 looks good for simScale around 1)
  float friction = 0.6; // make smaller to slow sim (0.6 looks good for a simScale around 1)
  float forces[numTypes][numTypes]; // forces table
  float minDistances[numTypes][numTypes]; // minDistances table
  float radii[numTypes][numTypes]; // radii table

  vector<Particle> swarm; // swarm vector

  // state methods
  void seed() {
    for (int i = 0; i < numParticles; i++) {  // for each iter...
      Particle particle; // initialzie a particle
      particle.type = rnd::uniformi(0, numTypes - 1); // give random type
      particle.position = randomVec3f(simScale); // give random pos within simScale 
      particle.velocity = 0; // give initial velocity of 0
      swarm.push_back(particle); // append to swarm vector
    }
  }

  void setParameters(int numTypes) { // define setParams function (seems to be working)
    for (int i = 0; i < numTypes; i++) {
      for (int j = 0; j < numTypes; j++) {
        forces[i][j] = rnd::uniform<float>(.01, .003); // .01, .003 for simScale of 1
        if (rnd::uniformi(1, 100) < 50) {
          forces[i][j] *= -1;
        }
        minDistances[i][j] = rnd::uniform<float>(.1, .05); // .1, .05 for simScale of 1
        radii[i][j] = rnd::uniform<float>(.5, .15); // .5, .15 for simScale of 1
        cout << "forces[" << i << "][" << j << "]: " << forces[i][j] << endl;
        cout << "minDistances[" << i << "][" << j << "]: " << minDistances[i][j] << endl;
        cout << "radii[" << i << "][" << j << "]: " << radii[i][j] << endl;
      }
    }
  }

  void update() {
    for (int i = 0; i < numParticles; i++) { // for each particle, ~60fps...
      float dis = 0; // initialize dis
      Vec3f direction = 0; // initialize direction
      Vec3f totalForce = 0; // initialize totalForce
      Vec3f acceleration = 0; // initialize acceleration

      float springForceMag = springConstant * (swarm[i].position.mag() - simScale); // create sphere 
      Vec3f normalizedNegative = -swarm[i].position / swarm[i].position.mag(); // create sphere
      acceleration += springForceMag * normalizedNegative; // create sphere
      
      for (int j = 0; j < numParticles; j++) {
        if (i == j) {continue;} // don't have particles calculate forces on themselves

        direction *= 0; // set direction to 0
        direction += swarm[j].position; // direction = j particle 
        direction -= swarm[i].position; // direction -= current particle
        dis = direction.mag(); // Euclidian distance between particles
        direction.normalize(); // normalize to unit vector

        if (dis < minDistances[swarm[i].type][swarm[j].type]) { // separation forces (personal space)
          Vec3f force = direction; 
          force *= abs(forces[swarm[i].type][swarm[j].type]) * -3; // calculate repulsion force based on type
          force *= fMap(dis, 0, minDistances[swarm[i].type][swarm[j].type], 1, 0); // map based on distance
          force *= K; // scale by K
          totalForce += force; // add to totalForce
        } 
        
        if (dis < radii[swarm[i].type][swarm[j].type]) { // love/hate forces
          Vec3f force = direction; 
          force *= forces[swarm[i].type][swarm[j].type]; // calculate force based on type
          force *= fMap(dis, 0, radii[swarm[i].type][swarm[j].type], 1, 0); // map based on distance (flip last two arguments?)
          force *= K; // scale by K
          totalForce += force; // add to totalForce
        } 
        
      } 

      acceleration += totalForce; // integrate totalForce
      swarm[i].velocity += acceleration; // integrate acceleration
      swarm[i].velocity *= friction; // apply friction here?
      swarm[i].position += swarm[i].velocity; // integrate velocity
      //swarm[i].velocity *= friction; // or apply friction here?
    } 
  }
};

class swarmOrb : public DistributedAppWithState<State> {
public:

  float channelLeft = 0;
  float channelRight = 0;

  void onInit() override {
    auto cuttleboneDomain =
        CuttleboneStateSimulationDomain<State>::enableCuttlebone(this);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      quit();
    }
  }

  void onCreate() {
  if (isPrimary()) {
    state().seed();
    state().setParameters(state().numTypes);
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertex(state().swarm[i].position);
      verts.color(HSV(state().swarm[i].type * state().colorStep, 1.f, 1.f));
    }
    verts.primitive(Mesh::POINTS);
  } else { // if != primary...
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertex(state().swarm[i].position);
      verts.color(HSV(state().swarm[i].type * state().colorStep, 1.f, 1.f));
    }
    verts.primitive(Mesh::POINTS);
  }
  }

  bool freeze = false; // <- for pausing sim
  float phase = 0;
  void onAnimate(double dt) {
  if (isPrimary()) {
    if (freeze) return; // <- if freeze is true, then pause sim
    phase += dt;
    if (phase > 3) {
      state().setParameters(state().numTypes);
      phase = 0;
    }
    state().update();
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertices()[i] = state().swarm[i].position;
    }
  } else { // if not primary... 
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertices()[i] = state().swarm[i].position;
    }
  }
  } 
  
  bool onKeyDown(const Keyboard &k) override {
  if (isPrimary()) {
    if (k.key() == ' ') { // <- on spacebar, freeze or unfreeze simulation
      freeze = !freeze; // <- invert state of freeze
    }
    if (k.key() == '1') { // <- on 1, setParams
      state().setParameters(state().numTypes);
    }
    return true;
  }
  }

  void onSound(AudioIOData& io) override{
  if (isPrimary()) {
    float maxSamp = 0;
    while(io()) { 
      channelLeft = io.in(0);
      channelRight = io.in(1);
      io.out(2) = channelLeft;
      io.out(3) = channelRight;
      float mixDown = abs(channelLeft + channelRight);
      if (mixDown > maxSamp) {
        maxSamp = mixDown;
      }
      state().pointSize = 10 * maxSamp;
    }
  }
  }

  void onDraw(Graphics &g) { 
    g.clear(0); // black background
    g.pointSize(state().pointSize); // set pointSize
    g.meshColor(); // color vertices based on type
    g.draw(verts); // draw verts
  }

private:
  Mesh verts;
};

int main() {
  swarmOrb app;
  AudioDevice alloAudio = AudioDevice("AlloAudio");
  alloAudio.print();
  app.configureAudio(alloAudio, 48000, 128, 4, 2);
  app.start();
}