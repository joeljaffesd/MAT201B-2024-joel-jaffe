// Joel Jaffe February 2024
//Based on https://www.youtube.com/watch?v=xiUpAeos168&list=PLZ1w5M-dmhlGWtqzaC2aSLfQFtp0Dz-F_&index=3
//Programming Chaos on YouTube

#include "al/app/al_App.hpp"
#include "al/system/al_Time.hpp"
#include "al/math/al_Random.hpp"
#include "al/app/al_GUIDomain.hpp"
using namespace al;

#include <iostream>
using namespace std;

struct Particle {
  int type;
  Vec2f position;
  Vec2f velocity;
};

float fMap(float value, float in_min, float in_max, float out_min, float out_max) {
  return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

float fWrap(float input, float scope) {
    float lowerBound = -1 * scope;
    float range = 2 * scope;
    return lowerBound + fmodf((input - lowerBound + range), range);
}

Vec2f randomVec2f(float scale) { // <- Function that returns a Vec3f containing random coords
  return Vec2f(rnd::uniformS(), rnd::uniformS()) * scale;
} 

struct MyApp : public App {
  Parameter timeStep{"/timeStep", "", 0.444, 0.01, 0.6}; // <- creates GUI parameter
  Parameter simScale{"/simScale", "", 0.9, 1.f, 100.f}; // <- creates GUI parameter
  Mesh verts;

  static const int numTypes = 6;
  int numParticles = 1000;
  float colorStep = 1.f / numTypes;
  float K = 0.0000005; // make smaller to slow sim
  float friction = 0.0000085; // make smaller to slow sim
  float forces[numTypes][numTypes];
  float minDistances[numTypes][numTypes];
  float radii[numTypes][numTypes];
  vector<Particle> swarm;

  void setParameters(int numTypes) {
    for (int i = 0; i < numTypes; i++) {
      for (int j = 0; j < numTypes; j++) {
        forces[i][j] = rnd::uniform<float>(1, 3);
        if (rnd::uniformi(1, 100) < 50) {
          forces[i][j] *= -1;
        }
        minDistances[i][j] = rnd::uniform<float>(50, 30);
        radii[i][j] = rnd::uniform<float>(250, 70);
      }
    }
  }

  void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(timeStep); // add parameter to GUI
    gui.add(simScale); // add parameter to GUI
  }

  void onCreate() {
    verts.primitive(Mesh::POINTS); // skin mesh as points
    for (int i = 0; i < numParticles; i++) {
      Particle particle;
      particle.type = rnd::uniformi(1, 6);
      particle.position = randomVec2f(simScale);
      particle.velocity = 0;
      swarm.push_back(particle);

      verts.vertex(particle.position);
      verts.color(HSV(particle.type * colorStep, 1.f, 1.f));
    }
    setParameters(numTypes);
  }


  bool freeze = false; // <- for pausing sim
  void onAnimate(double dt) {
    if (freeze) return; // <- if freeze is true, then pause sim
    dt = timeStep; // override dt

    for (int i = 0; i < numParticles; i++) {

      float dis = 0;
      Vec2f direction = 0;
      Vec2f totalForce = 0;
      Vec2f acceleration = 0;
       
      for (int j = 0; j < numParticles; j++) {
        if (i == j) {continue;}
        direction *= 0;
        direction += swarm[j].position;
        direction -= swarm[i].position;
        
        if (direction[0] > simScale) {
          direction[0] -= 2 *simScale;
        }
        if (direction[0] < -1 * simScale) {
          direction[0] += 2 * simScale;
        }
        if (direction[1] > simScale) {
          direction[1] -= 2 * simScale;
        }
        if (direction[1] < -1 * simScale) {
          direction[1] += 2 * simScale;
        }
        
        dis = direction.mag();
        direction.normalize();

        if (dis < minDistances[swarm[i].type][swarm[j].type]) {
          Vec2f force = direction;
          force *= abs(forces[swarm[i].type][swarm[j].type]) * -3;
          force *= fMap(dis, 0, minDistances[swarm[i].type][swarm[j].type], 1, 0);
          force *= K;
          totalForce += force;
        } 

        if (dis < radii[swarm[i].type][swarm[j].type]) {
          Vec2f force = direction;
          force *= forces[swarm[i].type][swarm[j].type];
          force += fMap(dis, 0, radii[swarm[i].type][swarm[j].type], 1, 0);
          force *= K;
          totalForce += force;
        } 

      } 

      acceleration += totalForce; // integrate totalForce
      swarm[i].velocity += acceleration; // integrate acceleration
      swarm[i].position += swarm[i].velocity; // integrate velocity
      
      swarm[i].position[0] = fWrap(swarm[i].position[0], simScale); // wrapping x-dim
      swarm[i].position[1] = fWrap(swarm[i].position[1], simScale); // wrapping y-dim
      
      swarm[i].velocity *= friction; // apply friction
      verts.vertices()[i] = swarm[i].position; // skin mesh
      
    } 

  } 
  
  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') { // <- on spacebar, freeze or unfreeze simulation
      freeze = !freeze; // <- invert state of freeze
    }
    if (k.key() == '1') { // <- on 1
      setParameters(numTypes);
    }
    return true;
  }

  void onDraw(Graphics &g) {
    g.clear(0);
    g.camera(Viewpoint::IDENTITY); // Ortho [-1:1] x [-1:1]
    g.pointSize(10);
    g.meshColor();
    g.draw(verts);
    g.camera(Viewpoint::ORTHO_FOR_2D); // Ortho [0:width] x [0:height]
  }
};

int main() {
  MyApp app;
  app.configureAudio(48000, 512, 2, 0);
  app.start();
}