#include "al/app/al_App.hpp"
#include "al/math/al_Random.hpp"
#include "al/graphics/al_Shapes.hpp" 
#include "al/math/al_Random.hpp"
#include "al/app/al_GUIDomain.hpp"
using namespace al;

#include <vector>
using namespace std;

/*
When a predator is too close, a prey _turns_ away. 
Keep creatures within a spherical boundary by _turning_ them toward the center when they stray too far.
Spawn food that attracts prey agents and disapears when consumed.
Spawn a small number of predators.

If you get the simulation working, animate the camera to alternate slowly between wide establishing shots and up close action shots. 
I.e., follow a predator for a bit from a "third person" point of view or watch some food when it spawns and turn away once it is consumed
*/

Vec3f randomVec3f(float scale) { // <- Function that returns a Vec3f containing random coords
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} 

struct MyApp : public al::App {
    Parameter personalSpace{"/personalSpace", "", 3.0, 0.0, 5.0}; // <- creates GUI parameter (OSC name, OSCGroup, default, min, max)
    Parameter cageSize{"/cageSize", "", 35, 0, 50}; // <- creates GUI parameter
    Parameter foodSense{"/foodSense", "", 20.0, 0, 20.0}; // <- creates GUI parameter

    Mesh mesh; // <- Create mesh
    Nav predator; // <- create nav for predator
    vector<Nav> preyFlock; // <- Create vector of navs for preyFlock
    vector<Nav> predatorFlock; // <- Create vector of navs for predatorFlock
    vector<Vec3f> flockColors; // <- Create vector to store flock colors
    vector<Vec3f> predColors; // <- Create vector to store predator colors
    vector<Vec3d> food; // <- Create vector to store food

    void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(personalSpace); // <- add parameter to GUI
    gui.add(cageSize); // <- add parameter to GUI
    gui.add(foodSense); // <- add parameter to GUI
  }

    void onCreate() {
        addCone(mesh); // <- Create cone mesh 
        mesh.generateNormals(); // <- Generate normals

        // Generate food
        int numFood = 3;
        food.resize(numFood);
        for (int i = 0; i < numFood; i++) {
          food[i] = al::rnd::ball<al::Vec3d>() * cageSize;
        }
        
        // Generate preyFlock
        int numBoids = 30; // <- Specify number of boids
        preyFlock.resize(numBoids);
        for (int i = 0; i < numBoids; i++) { 
            preyFlock[i].pos(randomVec3f(cageSize)); // <- Set initial position of prey
            flockColors.push_back(Vec3f(0, rnd::uniform(), rnd::uniform())); // <- generate random colors for each prey
        }

        // Generate predatorFlock
        int numPredators = 3; // <- Specify number of predators
        predatorFlock.resize(numPredators);
        for (int i = 0; i < numPredators; i++) { 
            predatorFlock[i].pos(randomVec3f(cageSize)); // <- Set initial position of predator
            predColors.push_back(Vec3f(rnd::uniform(), 0, 0)); // <- generate random colors for each predator
        }
        
        nav().pos(0, 0, 70); // <- Set initial position of camera
        nav().faceToward(0,0,0); // <- Make camera face origin
    }

    double phase = 0; // <- Initialize variable to accumulate dt
    double camPhase = 0;
    void onAnimate(double dt) {

      camPhase += dt;
      nav().pos(0, 0, 70 - camPhase); // <- Set initial position of camera
      nav().faceToward(0, 0, 0); // <- Make camera face origin

      phase += dt;
      if (phase > 4) {
        phase = 0;
      }

      for (int i = 0; i < food.size(); i++) {
        if (food[i].mag() > cageSize) { // if food is outside cage, 
            food[i] = al::rnd::ball<al::Vec3d>() * cageSize; // spawn new food
        }
      }
      

      double turn_rate = 0.05; // <- Set turn rate
      double move_rate = 0.04; // <- Set move rate

      
      float maxDistFromFood = 0;
      for (int i = 0; i < preyFlock.size(); i++) {
        Vec3f targetVec (0, 0, 0);

        // Manage Relationship to Food
        float closestFood = cageSize;
        for (int j = 0; j < food.size(); j++) {
          float distFromFood = (food[j] - preyFlock[i]).mag();
          if (distFromFood > maxDistFromFood) {
            maxDistFromFood = distFromFood;
          }
          if (maxDistFromFood > foodSense && phase > 3) {
            food[j] = al::rnd::ball<al::Vec3d>() * cageSize; // spawn new food
            }
          else if (distFromFood < foodSense) {
            if (distFromFood < closestFood) {
              closestFood = distFromFood;
              targetVec += food[j]; // <- prey faces food
            }
            if ((food[j] - preyFlock[i]).mag() < 0.1) { // if prey gets close enough to food
            food[j] = al::rnd::ball<al::Vec3d>() * cageSize; // eat it, spawn new food
            }
          }
        }

        // Manage neighbor groups
        Vec3f neighborAverage;
        Vec3f sum (0, 0, 0);
        float divisor = 0;
        for (int j = i + 1; j < preyFlock.size(); j++) {
          Vec3f toNeighbor = preyFlock[j].pos() - preyFlock[i].pos();
          float distance = toNeighbor.mag(); // calculate distance from neighbors
          if (distance < personalSpace) { // if in neighborhood
            targetVec -= toNeighbor; // keep away from neighbor
            sum += preyFlock[j].pos();
            divisor += 1;
            }
          }
          if (divisor > 0) { // pray faces (targetVec + neighborAverage) / 2
            neighborAverage = sum / divisor;
            targetVec += neighborAverage; 
            targetVec /= 2; 
          } 

          // Manage relationship to predators
          float fleeFactor = 1;
          for (int j = 0; j <predatorFlock.size(); j++) {
            Vec3f toPredator = predatorFlock[j].pos() - preyFlock[i].pos();
            float distance = toPredator.mag(); // calculate distance from predators
            if (distance < personalSpace) { // if danger close
              targetVec -= toPredator; // offset target
              fleeFactor += personalSpace - distance; // increase speed
              }
          }
            
          // Keep birds in cage
          if (targetVec.mag() > cageSize) { // if target is oustide cage
            targetVec /= 2; // set target inside cage
          }
          
          // "Integration"
          preyFlock[i].faceToward(targetVec, turn_rate); // <- pray faces target
          preyFlock[i].moveF(move_rate * fleeFactor); // <- set move rate
          preyFlock[i].step(); // <- pray moves towards target
        }

        for (int i = 0; i < predatorFlock.size(); i++) {
          Vec3f targetVec = preyFlock[0].pos();
          for (int j = 0; j < preyFlock.size(); j++) {
            float eucDist = (preyFlock[j].pos() - predatorFlock[i].pos()).mag();
            if (eucDist < (targetVec - predatorFlock[i].pos()).mag()) {
              targetVec = preyFlock[j].pos(); // set target to nearest prey
              }
          }
          for (int j = i + 1; j < predatorFlock.size(); j++) {
            Vec3f toNeighbor = predatorFlock[j].pos() - predatorFlock[i].pos();
            float distance = toNeighbor.mag();
            if (distance < personalSpace) {
              targetVec -= toNeighbor; // avoid crowding
              }
          }
          predatorFlock[i].faceToward(targetVec, turn_rate); // <- predator chases prey
          predatorFlock[i].moveF(move_rate / 2); // <- Predator moves half as fast as prey
          predatorFlock[i].step(); // <- advance predator
        }
    }

    void onDraw(al::Graphics& g) {
        g.depthTesting(true);
        g.lighting(true);
        g.clear(0);

        for (int i = 0; i < preyFlock.size(); i++) {
          g.color(flockColors[i][0], flockColors[i][1], flockColors[i][2]);
          g.pushMatrix();
          g.translate(preyFlock[i].pos());
          g.rotate(preyFlock[i].quat());
          g.scale(0.1);
          g.draw(mesh);
          g.popMatrix();
        }
       
        for (int i = 0; i < predatorFlock.size(); i++) {
          g.color(predColors[i][0], predColors[i][1], predColors[i][2]);
          g.pushMatrix();
          g.translate(predatorFlock[i].pos());
          g.rotate(predatorFlock[i].quat());
          g.scale(0.2);
          g.draw(mesh);
          g.popMatrix();
        }
        
        for (int i = 0; i < food.size(); i++) {
          g.color(0, 0, 0);
          g.pushMatrix();
          g.translate(food[i]);
          g.scale(0.1);
          g.scale(0.1);
          g.draw(mesh);
          g.popMatrix();
        }
    }
};

int main() {
    MyApp app;
    app.configureAudio(48000, 512, 2, 0);
    app.start();
}
