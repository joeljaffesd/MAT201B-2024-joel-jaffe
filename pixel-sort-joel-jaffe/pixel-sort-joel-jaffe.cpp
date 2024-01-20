#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/graphics/al_Image.hpp"
#include "al/io/al_File.hpp"
#include "al/graphics/al_Mesh.hpp"

using namespace std;
using namespace al;

struct MyApp : public App {
    Mesh targetMesh;
    Mesh currentMesh;
    double myClock = 0.0;
    double animationDuration = 2.0; // Animation duration in seconds
    double animationStartTime = 0.0; // Start time of the animation

    File file = File::currentPath() + "../image-joel-jaffe.png";
    Image image;

    float myLerp(float a, float b, float alpha) {
        return (1.0f - alpha) * a + alpha * b;
}

    void onCreate() {
        image = Image(file.path());
        targetMesh.primitive(Mesh::POINTS);
        currentMesh.primitive(Mesh::POINTS);
        for (int j = 0; j < image.height(); j++) {
            for (int i = 0; i< image.width(); i++) {
                auto pixel = image.at(i, j);
                currentMesh.vertex((1.0 * i) / image.width(), (1.0 * j) / image.height(), 0); //0-255
                currentMesh.color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);
            }
        }
        nav().pos(0.5, 0.5, 5);
    }

    void onAnimate(double dt) {
      myClock += dt;
      //change currentMesh to targetMesh over a set duration of time
        if (animationStartTime > 0.0) {
            double elapsed = myClock - animationStartTime;
            double alpha = elapsed / animationDuration;

            // Interpolate between the current and target configurations
            if (alpha < 1.0) {
                for (int i = 0; i < currentMesh.vertices().size(); i++) {
                        float targetX = targetMesh.vertices()[i][0];
                        float targetY = targetMesh.vertices()[i][1];
                        float targetZ = targetMesh.vertices()[i][2];

                        float currentX = currentMesh.vertices()[i][0];
                        float currentY = currentMesh.vertices()[i][1];
                        float currentZ = currentMesh.vertices()[i][2];

                        float newX = myLerp(currentX, targetX, alpha);
                        float newY = myLerp(currentY, targetY, alpha);
                        float newZ = myLerp(currentZ, targetZ, alpha);

                        currentMesh.vertices()[i] = {newX, newY, newZ};
                }
            } else {
                animationStartTime = 0.0; // Reset animation start time
            }
        }
    }

    void onDraw(al::Graphics& g) {
        g.clear(0.0);
        g.meshColor();
        g.pointSize(10);
        g.draw(currentMesh);
    }

      // This is called whenever a key is pressed.
    bool onKeyDown(const Keyboard& k) override {
    // Use a switch to do something when a particular key is pressed
    switch (k.key()) {
      case '1':
        //reset targetMesh
        targetMesh.reset();    
        //access pixels, display image
        for (int j = 0; j < image.height(); j++) {
            for (int i = 0; i< image.width(); i++) {
                auto pixel = image.at(i, j);
                targetMesh.vertex((1.0 * i) / image.width(), (1.0 * j) / image.height(), 0); //0-255
                targetMesh.color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0); //0.0-255.0
            }
        }
        animationStartTime = myClock;
        //cout status
        std::cout << "Displaying Input Image." << std::endl;
        break;

      case '2':
        //reset targetMesh
        targetMesh.reset(); 
        //access pixels, construct RGB Cube
        for (int j = 0; j < image.height(); j++) {
            for (int i = 0; i< image.width(); i++) {
                auto pixel = image.at(i, j);
                targetMesh.vertex(pixel.r / 255.0 , pixel.g / 255.0, pixel.b / 255.0);
                targetMesh.color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);
            }
        }
        //cout status
        animationStartTime = myClock;
        std::cout << "Reconfigured Pixels into an RGB Cube." << std::endl;
        break;

      case '3':
       //reset targetMesh
        targetMesh.reset(); 
        //access pixels, construct HSV Cylinder
        for (int j = 0; j < image.height(); j++) {
            for (int i = 0; i< image.width(); i++) {
                auto pixel = image.at(i, j);
                RGB rgb = RGB(pixel.r / 255.0 , pixel.g / 255.0 , pixel.b / 255.0);
                HSV hsv = HSV(rgb);

                float angle = hsv.h * M_2PI;
                float radius = hsv.s;
                float x = radius * cos(angle);
                float z = radius * sin(angle);

                targetMesh.vertex(x, hsv.v , z);
                targetMesh.color(rgb.r, rgb.g, rgb.b); 
            }
        }
        //cout status
        animationStartTime = myClock;
        std::cout << "Reconfigured Pixels into an HSV Cylinder." << std::endl;
        break;

      case '4':
        //reset targetMesh
        targetMesh.reset(); 
        //access pixels, construct Spiral Galaxy
        for (int j = 0; j < image.height(); j++) {
            for (int i = 0; i< image.width(); i++) {
                auto pixel = image.at(i, j);
                RGB rgb = RGB(pixel.r / 255.0 , pixel.g / 255.0 , pixel.b / 255.0);
                HSV hsv = HSV(rgb);

                float angle = hsv.h * M_2PI + j * 0.1;
                float radius = hsv.s * j * 0.01;
                float x = radius * cos(angle);
                float z = radius * sin(angle);

                targetMesh.vertex(x, hsv.v , z);
                targetMesh.color(rgb.r, rgb.g, rgb.b); 
            }
        }
        //cout status
        animationStartTime = myClock;
        std::cout << "Reconfigured Pixels into a Sprial Galaxy." << std::endl;
        break;
    }
    return true;
  }
};

int main() {
    MyApp app;
    app.configureAudio(48000, 512, 2, 0);
    app.start();
}