
#include <iostream>
using namespace std;

#include "al/app/al_DistributedApp.hpp"
using namespace al;

struct foo : DistributedApp {

  float channelLeft = 0;
  float channelRight = 0;
  float scaleVal = 0;
  Mesh mesh;

  float ampToDec (float in) {
    return 20.f * log10(in);
  }

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
    float maxSamp = 0;
    float lastSamp = 0;
    float sampleTime = 1 / 48000;
    float slideIn = exp(-sampleTime / 1);
    float slideOut = exp(-sampleTime / 50);

    while(io()) { 
      channelLeft = io.in(0);
      channelRight = io.in(1);
      io.out(2) = channelLeft;
      io.out(3) = channelRight;
      float mixDown = abs(channelLeft + channelRight);
      float mSlide = (0.5 * mixDown) + (0.5 * lastSamp);
      //float dbVal = abs(ampToDec(mSlide));
      //float mSlide = lastSamp + ((mixDown - lastSamp) * 0.5);
      //float mDifference = mixDown - lastSamp;
      //float scaledDiff = mDifference * ((slideOut + lastSamp) / (slideIn + lastSamp));
      //float mSlide = lastSamp + scaledDiff;
      if (mSlide > maxSamp) {
        maxSamp = mSlide;
      }
      lastSamp = mSlide;
      //scaleVal = maxSamp * 0.1;
    }
    scaleVal = abs(maxSamp) * 0.1f;
  }

};

int main() {
  foo app;
  AudioDevice alloAudio = AudioDevice("AlloAudio");
  alloAudio.print();
  app.configureAudio(alloAudio, 48000, 128, 4, 2);
  app.start();
}
