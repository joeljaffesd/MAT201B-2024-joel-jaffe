
#include <iostream>
using namespace std;

#include "al/app/al_DistributedApp.hpp"
using namespace al;

struct foo : DistributedApp {

  float channelLeft = 0;
  float channelRight = 0;
  float scaleVal = 1;
  Mesh mesh;

  void onCreate() override {
    addSphere(mesh);
  }

  void onAnimate(double dt) override {
  }

  void onDraw(Graphics &g) override {
    g.clear(0, 0, 1);
    g.scale(scaleVal);
    g.draw(mesh);
  }

  bool onKeyDown(const Keyboard &k) override {
  }

  void onSound(AudioIOData& io) override{
    while(io()) { 
      channelLeft = io.in(0);
      channelRight = io.in(1);
      io.out(2) = channelLeft;
      io.out(3) = channelRight;
      //float monitor = static_cast <float> (io.out(2));
      //cout << monitor << endl;
      scaleVal = abs(channelLeft + channelRight);
      /*
      Generating control signals for natural-feeling motion

      For simScale, summing channels and sending to 
        peakamp(10ms)
        slide(1, 50)
      works pretty well. This scheme works less well with highly compressed audio 

      For setParams, sum channels and send to
        delta~
        if (abs($1) > $f2) {setParams}
      works decent, $f2 must be tuned

      In general, these amplitude-domain techniques have trouble 
      with high compressed audio. 
      It may be best to use spectral techniques instead
      */
    }
  }

};

int main() {
  foo app;
  AudioDevice alloAudio = AudioDevice("AlloAudio");
  alloAudio.print();
  app.configureAudio(alloAudio, 48000, 128, 4, 2);
  app.start();
}
