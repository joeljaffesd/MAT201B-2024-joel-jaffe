// Oraculum: Omnispherical Audio-Reactive Visuals for Electric Instruments
// Joel Jaffe February 2024

// Particle life based on
// https://www.youtube.com/watch?v=xiUpAeos168&list=PLZ1w5M-dmhlGWtqzaC2aSLfQFtp0Dz-F_&index=3
// Programming Chaos on YouTube

// This implementation is optimized for use in UCSB's Allosphere

// TO DO
// -disable stereo rendering
// -tune background color for sphere (hide projector borders)

#include "al/app/al_App.hpp"
#include "al/system/al_Time.hpp"
#include "al/math/al_Random.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/app/al_DistributedApp.hpp"
#include "al/app/al_OmniRendererDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"
using namespace al;

#include "Gamma/SamplePlayer.h"

#include <iostream>
using namespace std;

float fMap(float value, float in_min, float in_max, float out_min, float out_max) { // custom mapping function
  return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

float dBtoA (float dBVal) {
  float ampVal = pow(10.f, dBVal / 20.f);
  return ampVal;
}

float ampTodB (float ampVal) {
  float dBVal = 20.f * log10f(abs(ampVal));
  return dBVal;
}

Vec3f randomVec3f(float scale) { // <- Function that returns a Vec2f containing random coords
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} 

class EnvFollower { // implementation based on https://www.musicdsp.org/en/latest/Analysis/97-envelope-detector.html
public:
  EnvFollower (float sampleRate, float attackMs, float releaseMs) {
    attackCoef = exp(-1.0f / (sampleRate * (attackMs / 1000)));
    releaseCoef = exp(-1.0f / (sampleRate * (releaseMs / 1000)));
    envOut = 0.f;
  }

  float processBuffer (float samples[], int size) {
    for (int i = 0; i < size; i++) {
      float envIn = fabs(samples[i]);
      if (envOut < envIn) {
        envOut = envIn + attackCoef * (envOut - envIn);
      } else {
        envOut = envIn + releaseCoef * (envOut - envIn);
      }
    }
    return envOut;
  }

private:
  float attackCoef;
  float releaseCoef;
  float envOut;
};

struct Particle { // Particle struct
  int type; // <- this is the problem.
  Vec3f position; 
  Vec3f velocity;
};

struct SimulationState {
  // state() member variables
  float pointSize;
  float simScale = 0.5f;
  float springConstant = 0.4f;

  static const int numTypes = 6; // numTypes
  static const int numParticles = 1000; // numParticles (1000 seems to be the limit for my M2 Max)
  float colorStep = 1.f / numTypes; // colorStep
  float K = 0.05; // make smaller to slow sim (0.05 looks good for simScale around 1)
  float friction = 0.6; // make smaller to slow sim (0.6 looks good for a simScale around 1)
  float forces[numTypes][numTypes]; // forces table
  float minDistances[numTypes][numTypes]; // minDistances table
  float radii[numTypes][numTypes]; // radii table

  Particle swarm[numParticles];

  // state() methods
  void seed() {
    for (int i = 0; i < numParticles; i++) {  // for each iter...
      swarm[i].type = rnd::uniformi(0, numTypes - 1); // give random type
      swarm[i].position = randomVec3f(simScale); // give random pos within simScale 
      swarm[i].velocity = 0; // give initial velocity of 0
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
        //cout << "forces[" << i << "][" << j << "]: " << forces[i][j] << endl;
        //cout << "minDistances[" << i << "][" << j << "]: " << minDistances[i][j] << endl;
        //cout << "radii[" << i << "][" << j << "]: " << radii[i][j] << endl;
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

class swarmOrb : public DistributedAppWithState<SimulationState> {
public:
  Parameter volControl{"volControl", "", 0.f, -96.f, 6.f};
  Parameter rmsMeter{"/rmsMeter", "", -96.f, -96.f, 0.f};
  Parameter dBThresh{"/dBThresh", "", -4.f, -96.f, 0.f}; // -17 great for Tom Sawyer
  Parameter envAttack{"/envAttack", "", 32.f, 1.f, 50.f}; // 7.2 looks good!, 11.8!
  Parameter envRelease{"/envRelease", "", 399.f, 10.f, 500.f}; // 270 looks good! 399!
  Parameter pointSizeOffset{"/pointSizeOffset", "", 2.f, 0.f, 10.f};
  Parameter clearValue{"/clearValue", "", 0.f, 0.f, 1.f};
  ParameterBool audioOutput{"audioOutput", "", false, 0.f, 1.f};
  ParameterBool filePlayback{"filePlayback", "", false, 0.f, 1.f};
  gam::SamplePlayer<float, gam::ipl::Linear, gam::phsInc::Loop> player;

  void onInit() override {
    auto cuttleboneDomain =
      CuttleboneStateSimulationDomain<SimulationState>::enableCuttlebone(this);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      quit();
    }
  if (isPrimary()) {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(volControl); // add parameter to GUI
    gui.add(rmsMeter); // add parameter to GUI
    gui.add(dBThresh); // add parameter to GUI
    gui.add(envAttack); // add parameter to GUI
    gui.add(envRelease); // add parameter to GUI
    gui.add(pointSizeOffset); // add parameter to GUI
    gui.add(clearValue); // add parameter to GUI
    gui.add(audioOutput); // add parameter to GUI
    gui.add(filePlayback); // add parameter to GUI
    
    //load file to player
    player.load("../Resources/HuckFinn.wav");
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
      verts.vertex(state().swarm[i].position); // initializes wrong, overridden by onAnimate
      verts.color(HSV(state().swarm[i].type * state().colorStep, 1.f, 1.f)); // initializes wrong, overridden by onAnimate
    }
    verts.primitive(Mesh::POINTS);
  }
  displayMode(Window::DisplayMode::DOUBLE_BUF | Window::DisplayMode::ALPHA_BUF | Window::DisplayMode::DEPTH_BUF);
  cout << "Display Mode: " << displayMode() << endl;
  }

  bool freeze = false; // <- for pausing sim
  float phase = 0;
  void onAnimate(double dt) {
  if (isPrimary()) {
    if (freeze) return; // <- if freeze is true, then pause sim
    state().update(); // <- simulation step
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertices()[i] = state().swarm[i].position; // update mesh
    }
  } else { // if not primary... 
    for (int i = 0; i < state().numParticles; i++) {
      verts.vertices()[i] = state().swarm[i].position; // update mesh
      verts.colors()[i] = HSV(state().swarm[i].type * state().colorStep, 1.f, 1.f); // update mesh
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
    if (k.key() == 'm') { // <- on m, muteToggle
      audioOutput = !audioOutput;
      cout << "Mute Status: " << audioOutput << endl;
    }
    if (k.key() == 'p') { // <- on p, playTrack
      filePlayback = !filePlayback;
      player.reset();
      cout << "File Playback: " << filePlayback << endl; 
    }
    return true;
  }
  }

  bool hold = false;
  int boomCounter = 0;
  void onSound(AudioIOData& io) override{
  if (isPrimary()) {
    EnvFollower envFollow (io.framesPerSecond(), envAttack, envRelease); // instance of EnvFollower
    
    // audio throughput
    while(io()) { 
    if (filePlayback) {
      for (int i = 0; i < io.channelsOut(); i++) {
        if (i % 2 == 0) {
          io.out(i) = player(0) * dBtoA(volControl) * audioOutput;
        } else {
          io.out(i) = player(1) * dBtoA(volControl) * audioOutput;
        }
      }
    } else {
      for (int i = 0; i < io.channelsOut(); i++) {
        io.out(i) = io.in(0) * dBtoA(volControl) * audioOutput; // <- feedback risk!
      }
    }
    }

    // audio analysis, including:
    // -amplitude thresholding to trigger state().setParamaters
    // -enveloper follower to drive pointSize
    // -RMS calculation for rmsMeter
    float myBuffer [io.framesPerBuffer()]; // initialize myBuffer
    for (int i = 0; i < io.framesPerBuffer(); i++) {
      if (filePlayback) {
        myBuffer[i] = io.out(0, i) + io.out(1, i); // populate myBuffer with samples
      } else {
        myBuffer[i] = io.in(0, i); // populate myBuffer with samples
      }
    }

    float threshAmp = dBtoA(dBThresh); // set threshAmp
    float bufferPower = 0; // initialize bufferPower
    for (int i = 0; i < io.framesPerBuffer(); i++) {
      myBuffer[i] *= dBtoA(volControl); // scale by volControl
      bufferPower += myBuffer[i] * myBuffer[i];
    }

    // Envelope following
    state().pointSize = 200 * log(1 + envFollow.processBuffer(myBuffer, io.framesPerBuffer())) + pointSizeOffset;

    // RMS calculation + amplitude thresholding
    bufferPower /= io.framesPerBuffer();
    rmsMeter = ampTodB(bufferPower);
    if (bufferPower > threshAmp && !hold) {
      boomCounter += 1;
      cout << "BOOM!" << boomCounter << endl;
      state().setParameters(state().numTypes);
      hold = true;
    } else if (bufferPower < threshAmp && hold) {
      hold = false;
    }
  }
  }

  void onDraw(Graphics &g) { 
    g.clear(clearValue); // black background
    g.pointSize(state().pointSize); // set pointSize
    g.meshColor(); // color vertices based on type
    g.draw(verts); // draw verts
  }

private:
  Mesh verts;
};

int main() {
  swarmOrb app;

  if (al::Socket::hostName() == "ar01") { // if in AlloSphere...
    AudioDevice alloAudio = AudioDevice("ECHO X5");
    alloAudio.print();
    app.player.rate(1.0 / alloAudio.channelsOutMax());
    app.configureAudio(alloAudio, 44100, 128, alloAudio.channelsOutMax(), 2);
  } else { // if not... 
    AudioDevice alloAudio = AudioDevice("AlloAudio");
    alloAudio.print();
    app.player.rate(1.0 / alloAudio.channelsOutMax());
    app.configureAudio(alloAudio, 44100, 128, alloAudio.channelsOutMax(), 2);
  }
  
  app.start();
}