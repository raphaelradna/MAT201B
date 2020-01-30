// Raphael Radna
// MAT-201B W20
// Assignment 2

/* 1. Copy the starter code to a file named _particles-p1.cpp_. Implement a
 * stable[^stability] particle system based on _[Newton's law of universal
 * gravitation][]_ where all particles have a mass of 1. Balance max-force,
 * drag-factor, time-step, gravitational-constant */

#include "al/app/al_App.hpp"
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp" // gui.draw(g)

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

Vec3f rv(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

string slurp(string fileName); // forward declaration
float constrain(float f, float min, float max);

struct AlloApp : App {
  // add more GUI here
  Parameter pointSize{"/pointSize", "", 2.0, "", 0.0, 4.0};
  Parameter timeStep{"/timeStep", "", 0.1, "", 0.01, 2.0};
  ControlGUI gui;

  ShaderProgram pointShader;
  Mesh mesh; // vector<Vec3f> position is inside mesh

  // typedef al::Vec<float, 3> Vec3f;
  // typedef std::vector<Vec3f> Vertices;

  //  simulation state
  vector<Vec3f> velocity;
  vector<Vec3f> acceleration;
  vector<float> mass;

  const float g = 6.674 * pow(10, -11);

  void onCreate() override {
    // add more GUI here
    gui << pointSize << timeStep;
    gui.init();
    navControl().useMouse(false);

    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // set initial conditions of the simulation
    //

    // c++11 "lambda" function
    auto rc = []() { return HSV(rnd::uniform(), 1.0f, 1.0f); };

    mesh.primitive(Mesh::POINTS);
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?
    for (int _ = 0; _ < 1000; _++) {
      mesh.vertex(rv(1));
      mesh.color(rc());

      // float m = rnd::uniform(3.0, 0.5);
      float m = 1;
      mass.push_back(m);

      // using a simplified volume/size relationship
      mesh.texCoord(pow(m, 1.0f / 3), 0); // s, t

      // separate state arrays
      velocity.push_back(Vec3f(rv(1)));
      acceleration.push_back(Vec3f(0));
    }

    nav().pos(0, 0, 10);
  }

  bool freeze = false;
  void onAnimate(double dt) override {
    if (freeze)
      return;

    // ignore the real dt and set the time step;
    dt = timeStep;

    // Calculate forces

    // pair-wise and equal but opposite
    // nested for loop to visit each pair once
    // O(n*n)
    //
    for (int i = 0; i < mesh.vertices().size(); i++) {
      for (int j = i + 1; j < mesh.vertices().size(); j++) {
        // F = g*(m1m2/r^2)
        Vec3f r = mesh.vertices()[j] - mesh.vertices()[i];
        Vec3f f = r.normalize() * 0.1f * mass[i] * mass[j] / r.magSqr();
        // F = ma
        acceleration[i] += f / mass[i];
        acceleration[j] -= f / mass[j];
      }
      // drag
      acceleration[i] -= velocity[i] * 0.1f;
    }

    // Vec3f has
    // • +=
    // • -=
    // • .normalize()
    // • .normalize(float scale)
    // • .mag()
    // • .magSqr()
    // • .dot(Vec3f f)
    // • .cross(Vec3f f)

    // Integration
    //
    vector<Vec3f> &position(mesh.vertices());
    for (int i = 0; i < velocity.size(); i++) {
      // "backward" Euler integration
      velocity[i] += acceleration[i] / mass[i] * dt;
      position[i] += velocity[i] * dt;

      // Explicit (or "forward") Euler integration would look like this:
      // position[i] += velocity[i] * dt;
      // velocity[i] += acceleration[i] / mass[i] * dt;
    }

    // clear all accelerations (IMPORTANT!!)
    for (auto &a : acceleration)
      a.zero();
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '1') {
      // introduce some "random" forces
      for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        acceleration[i] = rv(1) / mass[i];
      }
    }

    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear();
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.blending(true);
    g.blendModeTrans();
    g.depthTesting(true);
    g.draw(mesh);
    gui.draw(g);
  }
};

int main() { AlloApp().start(); }

string slurp(string fileName) {
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}

float constrain(float f, float min, float max) {
  if (f < min) {
    f = min;
  }
  if (f > max) {
    f = max;
  }
  return f;
};
