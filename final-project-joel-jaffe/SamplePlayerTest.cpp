// Karl Yerkes | MAT201B
// 2022-02-22  | Audio-reactive visualizer
//

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Functions.hpp"  // al::abs

#include "Gamma/SamplePlayer.h"  // XXX
#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"

using namespace al;

struct MyApp : App {
  Parameter value{"value", 0, 0, 1};

  gam::SamplePlayer<float, gam::ipl::Linear, gam::phsInc::Loop> player;
  gam::EnvFollow<> follow;

  Mesh sphere;

  void onCreate() override { addSphere(sphere); 
  
    follow.lag(0.5); 
    player.rate(0.5);

    std::cout << player.frameRate() << std::endl;
  
  }

  void onDraw(Graphics& g) override {
    g.clear(0.2);
    g.scale(0.9 * value.get());
    g.draw(sphere);
  }

  void onSound(AudioIOData& io) override {
    while (io()) {
      float l = player(0);
      float r = player(1);
      io.out(2) = l;
      io.out(3) = r;

      value.set(follow(l + r));
    }
  }

  void onInit() override {
    player.load("../Resources/HuckFinn.wav");

    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto& gui = GUIdomain->newGUI();
    gui.add(value);
    parameterServer() << value;
  }

  void onMessage(osc::Message& m) override { m.print(); }
  bool onKeyDown(const Keyboard& k) override { return false; }
  void onAnimate(double dt) override {}
};

int main() {
  MyApp app;
  AudioDevice alloAudio = AudioDevice ("AlloAudio");
  app.configureAudio(alloAudio, 44100, 128, alloAudio.channelsOutMax(), 2);
  app.start();
}
