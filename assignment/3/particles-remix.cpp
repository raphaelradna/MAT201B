// Raphael Radna
// MAT-201B W20
// Assignment 3

/* 2. Choose code from another student's repository---Make a copy of their
 * `particles-p?.cpp` and add it to your of repository, clearly marking the copy
 * with the name of their repo and the commit you copied from. Attribute the
 * copy to them by name. Make some changes. Sumbit this as `particles-remix.cpp`
 */

/* Original code 'particles-p4.cpp' by Stejara Dinulescu
 * https://github.com/sdinulescu/MAT201B
 * 220d57f60d4f7c08ea9892454106fd75a75268dd
 */

/* Changes I made:
 *  – Point masses now 'swell', modulated by a sine function; I edited the
 *  geometry shader slightly to represent this visually
 *  – Added 'frequency' and 'depth' parameters for the modulation of the mass
 *  – Altered parameter ranges and defaults to better demonstrate new
 *  functionality
 *  – Initialized particles with a range of masses in order to show a variety of
 * interactions
 *  – Limited range of particle colors
 */

#include "al/app/al_App.hpp"
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp" // gui.draw(g)

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

int partNum = 1000;

Vec3f rv(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} // assigns randomness in three dimensions, multiplies it by scale

string slurp(string fileName); // forward declaration

struct AlloApp : App {
  Parameter pointSize{"/pointSize", "", 0.5, "", 0.0, 5.0};
  Parameter timeStep{
      "/timeStep", "", 0.05, "",
      0.01,        0.6}; // simplest way to not get NANs, keep timestep small
  Parameter gravConst{"/gravConst", "", 0.75, "", 0, 1};
  Parameter dragFactor{"/dragFactor", "", 0.075, "", 0.01, 0.99};
  // Parameter gravityBound{"/gravityBound", "", 0.8, "", 0.01, 0.99};
  Parameter maxAccel{"/maxAccel", "", 10, "", 0, 20};
  Parameter scaleVal{"/scaleVal", "", 0.8, "", 0, 2};
  // add GUI params here
  Parameter frequency{"/frequency", "", 1, "", 0, 5};
  Parameter depth{"/depth", "", 1, "", 0, 5};

  ControlGUI gui;

  ShaderProgram pointShader;
  Mesh mesh; // simulation state position is located in the mesh (positions are
             // the direct simulation states that we use to draw)

  // simulation state
  vector<Vec3f> velocity;
  vector<Vec3f> acceleration;
  vector<float> mass, massMod;

  void reset() { // empty all containers
    mesh.reset();
    velocity.clear();
    acceleration.clear();

    // c++11 "lambda" function
    // seed random number generators to maintain determinism
    rnd::Random<> rng;
    rng.seed(42);
    auto rc = [&]() { return HSV(rng.uniform() / 2 + 0.5f, 1.0f, 1.0f); };
    auto rv = [&](float scale) -> Vec3f {
      return Vec3f(rng.uniformS(), rng.uniformS(), rng.uniformS()) * scale;
    };

    mesh.primitive(Mesh::POINTS);
    for (int r = 0; r < partNum; r++) { // create 1000 points, put it into mesh
      mesh.vertex(rv(5));
      mesh.color(rc());

      float m = 5 + rnd::uniformS() * 10; // mass = 1
      mass.push_back(m);
      massMod.push_back(m);

      // set texture coordinate to be the size of the point (related to the
      // mass) using a simplified volume size relationship -> V = 4/3 * pi * r^3
      // pow is power -> m^(1/3)
      mesh.texCoord((4 / 3) * 3.14 * pow(m, 1.0f / 3), 0);
      // pass in an s, t (like x, y)-> where on an image do we want to grab the
      // color from for this pixel (2D texture) normalized between 0 and 1

      // separate state arrays for v and a

      // here, instead of a random velocity, use the cross product to make an
      // initial velocity that "tends to spin" around the yaxis this is so that
      // they don't all crash into each other
      velocity.push_back(rv(0.1));   // start with some small velocity
      acceleration.push_back(rv(0)); // start with no acceleration
    }
  }

  void onCreate() override {
    // add more GUI here
    gui << pointSize << timeStep << gravConst << dragFactor << maxAccel
        << scaleVal << frequency << depth; // stream operator
    gui.init();
    navControl().useMouse(false);

    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // set initial conditions

    reset();

    nav().pos(0, 0, 50); // push camera back
  }

  bool freeze =
      false; // state that freezes simulation -> doesn't run onAnimate if frozen
  int t = 0;
  void onAnimate(double dt) override {
    if (freeze)
      return;

    // ignore the real dt (real number of seconds that have passed since the
    // last call) and set the time step numerical simulation stability is due to
    // small, regular timesteps
    dt = timeStep; // simulation time, not wall time

    for (int i = 0; i < partNum; i++) {
      // float m = sinf();
      // mass.at(i) = m;

      massMod.at(i) =
          mass.at(i) * (float)depth +
          sinf(((float)t * (float)frequency / 60.0f) * 6.283185308) *
              mass.at(i) * (float)depth +
          0.01;
      // mesh.texCoord1s()[i] = ((4 / 3) * 3.14 * pow(m, 1.0f / 3), 0);
      // above line segfaults, I think due to attempt to change mesh after it's
      // compiled?
    }
    if (t >= 1 / (float)frequency * 60)
      t = 0;
    t++;

    // *********** Calculate forces ***********

    // gravity
    for (int i = 0; i < partNum;
         i++) { // nested for loops (for each particle, calculate force with all
                // other particles but itself one at a time)
      for (int j = 1 + i; j < partNum; j++) {

        rnd::Random<> rng;
        auto rv = [&](float scale) -> Vec3f {
          return Vec3f(rng.uniformS(), rng.uniformS(), rng.uniformS()) * scale;
        };

        Vec3f distance(
            mesh.vertices()[j] -
            mesh.vertices()[i]); // calculate distances between particles
        Vec3f gravityVal = gravConst * massMod[i] * massMod[j] *
                           distance.normalize() /
                           pow(distance.mag(), 2); // F = G * m1 * m2 / r^2

        // using a random multiplier here in order to show a more naturalistic,
        // "swimming" motion of the particles
        acceleration[i] += gravityVal * rv(scaleVal) / massMod[i];
        acceleration[j] -= gravityVal * rv(scaleVal) / massMod[j];
      }
    }

    // drag -> stabilizes simulation
    for (int i = 0; i < partNum; i++) {
      // force of drag is proportional to the opposite of velocity * small
      // amount normally, it is v^2 -> can change the statement to
      // vector[i].mag() * velocity[i] * 0.07; take a bit of acceleration away
      // proportional to what the velocity is
      acceleration[i] -= velocity[i] * dragFactor;
    }

    // limit acceleration
    for (int i = 0; i < acceleration.size(); i++) {
      float m = acceleration[i].mag();
      if (m > maxAccel) {
        acceleration[i].normalize(maxAccel);
        // cout << "Limiting Acceleration: " << m << " -> " << (float)maxAccel
        //      << endl;
      }
    }

    // Integration -> don't mess with this
    vector<Vec3f> &position(
        mesh.vertices()); // reference (alias) to mesh.vertices()
    for (int i = 0; i < partNum; i++) {
      // "backward" Euler integration (semi-implicit, more stable)
      velocity[i] += acceleration[i] / massMod[i] * dt;
      position[i] +=
          velocity[i] *
          dt; // this is actually changing mesh.vertices() -> interchangeable

      // Explicit (or "forward") Euler integration would look like this:
      // position[i] += velocity[i] * dt;
      // velocity[i] += acceleration[i] / mass[i] * dt;
    }

    // clear all accelerations (IMPORTANT!!) -> accelerations have been used and
    // counted already
    for (auto &a : acceleration)
      a.zero();
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '1') {
      // introduce some "random" forces
      for (int i = 0; i < partNum; i++) {
        acceleration[i] = rv(1) / massMod[i]; // rv gives a force -> divide by
                                              // mass to get acceleration
        // cout << acceleration[i] << endl;
      }
    }

    if (k.key() == 'r') {
      reset();
    }

    return true;
  }

  void onDraw(Graphics &g) override {
    float scaledSize = pointSize / 50;
    g.clear();
    g.shader(pointShader);
    g.shader().uniform("pointSize", scaledSize);
    g.shader().uniform(
        "pointSizeAgain",
        scaledSize * (float)depth +
            sinf(((float)t * (float)frequency / 60.0f) * 6.283185308) *
                scaledSize * (float)depth);
    g.blending(true);
    g.blendModeTrans();
    g.depthTesting(true);
    g.draw(mesh);
    gui.draw(g);
  }
};

int main() { AlloApp().start(); }

// slurp definition
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