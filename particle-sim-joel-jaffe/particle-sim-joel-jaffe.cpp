// Karl Yerkes 
// 2022-01-20

// Modified by Joel Jaffe
// 2024-01-30

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

Vec3f randomVec3f(float scale) { // <- Function that returns a Vec3f containing random coords
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} 

string slurp(string fileName);  // loads glsl files ig??

struct AlloApp : App {
  Parameter pointSize{"/pointSize", "", 1.47, 0.0, 2.0}; // <- creates GUI parameter (OSC name, OSCGroup, default, min, max)
  Parameter timeStep{"/timeStep", "", 0.444, 0.01, 0.6}; // <- creates GUI parameter
  Parameter dragValue{"/dragValue", "", 0.33, 0.0, 1.0}; // <- creates GUI parameter
  Parameter springConstant{"/springConstant", "", 0.4, 0.0, 1.0}; // <- creates GUI parameter
  Parameter coulombConstant{"/coulombConstant", "", 8.92, 0, 10}; // <- creates GUI parameter
  Parameter loveConstant{"/loveConstant", "", 27, 0, 50}; // <- creates GUI parameter
  Parameter radius{"/radius", "", 2.5, 0.0, 5.0}; // <- creates GUI parameter

  ShaderProgram pointShader; // <- Defines shader for points

  //  simulation state
  Mesh mesh;  // position *is inside the mesh* mesh.vertices() are the positions
  vector<Vec3f> velocity; // <- creates Vec3f to hold velocity values
  vector<Vec3f> acceleration; // <- creates Vec3f to hold acceleration values
  vector<float> mass; // <- creates vector of floats to hold mass values

  void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(pointSize); // <- add parameter to GUI
    gui.add(timeStep);   // <- add parameter to GUI
    gui.add(dragValue); // <- add parameter to GUI
    gui.add(springConstant); // <- add parameter to GUI
    gui.add(coulombConstant); // <- add parameter to GUI
    gui.add(loveConstant); // <- add parameter to GUI
    gui.add(radius); // <- add parameter to GUI
  }

  void onCreate() override {
    // compile shaders // <- skins the points
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // set initial conditions of the simulation

    // c++11 "lambda" function
    auto randomColor = []() { return HSV(rnd::uniform(), 1.0f, 1.0f); }; // <- generates random colors

    mesh.primitive(Mesh::POINTS); // <- displays vertices as points
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?
    for (int _ = 0; _ < 1000; _++) { // <- adds 1000 vertcies to mesh
      mesh.vertex(randomVec3f(5)); // <- for each iter, generate vertex of scale 5
      mesh.color(randomColor()); // <- for each vertex, generate random color

      // float m = rnd::uniform(3.0, 0.5);
      float m = 3 + rnd::normal() / 2; // <- for each vertex, assign random mass
      if (m < 0.5) m = 0.5; // <- for each iter, sets minimum for mass
      mass.push_back(m); // <- for each iter, append mass value to mass vector

      // using a simplified volume/size relationship
      mesh.texCoord(pow(m, 1.0f / 3), 0);  // s, t // <- for each iter, ??

      // separate state arrays
      velocity.push_back(randomVec3f(0.1)); // <- for each vertex, set initial velocity with random Vec3f
      acceleration.push_back(randomVec3f(1)); // <- for each vertex, set initial acceleration with random Vec3f
    }

    nav().pos(0, 0, 10); // <- sets camera pos
  }

  bool freeze = false; // <- for pausing sim

  void onAnimate(double dt) override {
    if (freeze) return; // <- if freeze is true, then pause sim

    // ignore the real dt and set the time step;
    dt = timeStep; // <- sets dt equal to timestep, a variable controlled in gui

    // Calculate forces

    // XXX you put code here that calculates forces and sets accelerations 
    // These are pair-wise. Each unique pairing of two particles
    // These are equal but opposite: A exerts a force on B while B exerts that
    // same amount of force on A (but in the opposite direction!) Use a nested
    // for loop to visit each pair once The time complexity is O(n^2)
    //
    // Vec3f has lots of operations you might use...
    // • +=
    // • -=
    // • +
    // • -
    // • .normalize() ~ Vec3f points in the direction as it did, but has length
    // • .normalize(float scale) ~ same but length `scale` 
    // • .mag() ~ length of the Vec3f 
    // • .magSqr() ~ squared length of the Vec3f 
    // • .dot(Vec3f f) 
    // • .cross(Vec3f f)

    vector<Vec3f> &position(mesh.vertices()); // <- create pointer to mesh for position of each vertex
    for (int i = 0; i < mesh.vertices().size(); i++) {
        // for each vertex, Calculate the Euclidean distance from the origin
        float euclidianDistanceFromOrigin = sqrt(pow(position[i][0], 2) + pow(position[i][1], 2) + pow(position[i][2], 2));
        // for each vertex, Calculate the spring force
        float springForceMagnitude = springConstant * (euclidianDistanceFromOrigin - radius);
        // for each dimension, Calculate acceleration
        for (int dim = 0; dim < 3; dim++) {
          float direction = -position[i][dim] / euclidianDistanceFromOrigin;
          acceleration[i][dim] += springForceMagnitude * direction;
        }
      
      for (int j = 0; j < mesh.vertices().size(); j++) {
      // for each pair, calculate Euclidean distance
      float euclideanDistance = sqrt(pow(position[i][0] - position[j][0], 2) + pow(position[i][1] - position[j][1], 2) + pow(position[i][2] - position[j][2], 2));
      // Avoid division by zero
        if (euclideanDistance > 1e-6) {
        // for each pair, calculate the electroForce
        float electroForceMagnitude = ((coulombConstant * 1e-5) / pow(euclideanDistance, 2));
        // for each dimension of i and j, calculate acceleration
        for (int dim = 0; dim < 3; dim++) {
            float direction = (position[j][dim] - position[i][dim]) / euclideanDistance;
            acceleration[i][dim] -= electroForceMagnitude * direction;
            acceleration[j][dim] += electroForceMagnitude * direction;  // Opposite force on particle j
          }
        }

        // Asymmetical force (love)
        Color colorI = mesh.colors()[i]; // <- retrieve color of vertex i
        Color colorJ = mesh.colors()[j]; // <- retrieve color of vertex j
        // If colorI is "red" and colorJ is "blue", apply attraction force
        if (colorI.r > 0.5 && colorI.g < 0.5 && colorI.b < 0.5 &&  // <- Red color condition
        colorJ.r < 0.5 && colorJ.g < 0.5 && colorJ.b > 0.5) {  // <- Blue color condition
        for (int dim = 0; dim < 3; dim++) {
          float direction = (position[j][dim] - position[i][dim]) / euclideanDistance;
          float loveForce = ((loveConstant * 1e-5) / pow(euclideanDistance, 2));
          acceleration[i][dim] += loveForce * direction; // <- red particles chase blue particles
          }
        }

        else if (colorI.r < 0.5 && colorI.g < 0.5 && colorI.b > 0.5 &&  // <- Blue color condition
        colorJ.r < 0.5 && colorJ.g > 0.5 && colorJ.b < 0.5) { // <- Green color condition
        for (int dim = 0; dim < 3; dim++) {
          float direction = (position[j][dim] - position[i][dim]) / euclideanDistance;
          float loveForce = ((loveConstant * 1e-5) / pow(euclideanDistance, 2));
          acceleration[i][dim] += loveForce * direction; // <- blue particles chase green particles
          }
        }
      }
    }
    
    // drag
    for (int i = 0; i < velocity.size(); i++) { // <- for size of velocity vector
      acceleration[i] -= velocity[i] * dragValue; // <- subtract velocity [i] from acceleration [i] (heat death), make variable
    }

    // Integration
    //vector<Vec3f> &position(mesh.vertices());
    for (int i = 0; i < velocity.size(); i++) {
      // "semi-implicit" Euler integration
      velocity[i] += acceleration[i] / mass[i] * dt;
      position[i] += velocity[i] * dt;
      // Explicit (or "forward") Euler integration would look like this:
      // position[i] += velocity[i] * dt;
      // velocity[i] += acceleration[i] / mass[i] * dt;
    }
    // clear all accelerations (IMPORTANT!!)
    for (auto &a : acceleration) a.set(0);
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') { // <- on spacebar, freeze or unfreeze simulation
      freeze = !freeze; // <- invert state of freeze
    }
    if (k.key() == '1') { // <- on 1
      // introduce some "random" forces
      for (int i = 0; i < velocity.size(); i++) { // <- for size of velocity vector
        // F = ma
        acceleration[i] = randomVec3f(1) / mass[i]; // <- set acceleration value at [i] equal to random / mass at [i]
      }
    }
    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear(0.3); 
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.blending(true);
    g.blendTrans();
    g.depthTesting(true);
    g.draw(mesh);
  }
};

int main() { // <- basic allolib main function
  AlloApp app;
  app.configureAudio(48000, 512, 2, 0);
  app.start();
}

string slurp(string fileName) { // <- file management function
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}
