// Raphael Radna
// MAT-201B W20
// Assignment 2

/* 4. Copy the file _particles-p3.cpp_ to _particles-p4.cpp_. Start with one of
 * the systems above. Create your own twist on the system. Make something
 * different and interesting.
 */

#include "al/app/al_App.hpp"
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp" // gui.draw(g)

using namespace al;

#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

Vec3f rv(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

string slurp(string fileName); // forward declaration
float constrain(float f, float min, float max);

struct celestialBody {
  float d; // distance from origin
  float m; // mass
  float r; // radius
  Color c;

  celestialBody(float d_, float m_, float r_, Color c_) {
    d = d_;
    m = m_;
    r = r_;
    c = c_;
  }
};

struct AlloApp : App {
  // add more GUI here
  Parameter pointSize{"/pointSize", "", 2.5, "", 0.0, 4.0};
  Parameter timeStep{"/timeStep", "", 0.01, "", 0.01, 2.0};
  Parameter symmetry{"/symmetry", "", 0.5, "", 0.0, 1.0};
  ControlGUI gui;

  ShaderProgram pointShader;
  Mesh mesh; // vector<Vec3f> position is inside mesh

  // typedef al::Vec<float, 3> Vec3f;
  // typedef std::vector<Vec3f> Vertices;

  //  simulation state
  vector<Vec3f> velocity;
  vector<Vec3f> acceleration;
  vector<float> mass, massMod;
  vector<celestialBody> bodies;

  const float g = 1.190588106E-5; // true value is E-19

  void onCreate() override {
    // add more GUI here
    gui << pointSize << timeStep << symmetry;
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

    // values taken from
    // https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size and
    // https://www.universetoday.com/15462/how-far-are-the-planets-from-the-sun
    //
    celestialBody Sun(0, 333000, 109.3, Color(0.992157, 0.427451, 0.070588));
    bodies.push_back(Sun);
    celestialBody Jupiter(5.2, 317.83, 10.97, rc());
    bodies.push_back(Jupiter);
    celestialBody Saturn(9.58, 95.162, 9.14, rc());
    bodies.push_back(Saturn);
    celestialBody Neptune(30.1, 17.147, 3.865, rc());
    bodies.push_back(Neptune);
    celestialBody Uranus(19.2, 14.536, 3.981, rc());
    bodies.push_back(Uranus);
    celestialBody Earth(1, 1, 1, rc());
    bodies.push_back(Earth);
    celestialBody Venus(0.772, 0.815, 0.9499, rc());
    bodies.push_back(Venus);
    celestialBody Mars(1.52, 0.107, 0.532, rc());
    bodies.push_back(Mars);
    celestialBody Mercury(0.387, 0.0553, 0.3829, rc());
    bodies.push_back(Mercury);
    celestialBody Pluto(39.5, 0.0022, 0.186, rc());
    bodies.push_back(Pluto);

    mesh.primitive(Mesh::POINTS);
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?

    for (celestialBody cb : bodies) { // add planets
      float dist = sqrt(cb.d);        // reel things in a lil
      float size = sqrt(cb.r);        // scale things down a lil

      mesh.vertex(Vec3f(dist, 0, 0));
      mesh.color(cb.c);

      mass.push_back(cb.m);

      mesh.texCoord(size, 0); // modify relative sizes

      velocity.push_back(Vec3f(0, sqrt(g * Sun.m / dist), 0));
      acceleration.push_back(Vec3f(0));
    }

    for (int i = 0; i < 190; i++) { // add moons
      mesh.vertex(rv(5));
      mesh.color(HSV(0.0, 0.0, 1.0));

      float m = 10.0 + rnd::uniform(0.0, 10.0);
      mass.push_back(m);

      mesh.texCoord(pow(m, 1.0f / 2), 0); // s, t
      velocity.push_back(rv(3));
      acceleration.push_back(Vec3f(0));
    }

    for (int i = 0; i < 800; i++) { // add comets
      mesh.vertex(rv(5));
      mesh.color(HSV(0.0, 0.0, 0.25));

      float m = 5.0 + rnd::uniform(0.0, 5.0);
      mass.push_back(m);

      mesh.texCoord(pow(m, 1.0f / 2), 0); // s, t

      velocity.push_back(rv(5));
      acceleration.push_back(Vec3f(0));
    }

    velocity.at(0) = Vec3f(0); // set sun vel to 0

    for (int i = 0; i < mass.size(); i++) {
      massMod.push_back(mass[i]);
    }

    nav().pos(0, 0, 20);
  }

  bool freeze = false;
  int t = 0;
  void onAnimate(double dt) override {
    if (freeze)
      return;

    // ignore the real dt and set the time step;
    dt = timeStep; // a realistic g value makes things move very slowly
    // Calculate forces

    for (int i = 0; i < mass.size(); i++) {
      massMod.at(i) = mass.at(i) * abs(sinf((t / 60.0f) * 6.283185308));
    }
    t++;
    if (t >= 60) {
      t = 0;
    }

    // pair-wise and equal but opposite
    // nested for loop to visit each pair once
    // O(n*n)
    //

    for (int i = 0; i < mesh.vertices().size(); i++) {
      for (int j = i + 1; j < mesh.vertices().size(); j++) {
        Vec3f r = mesh.vertices()[j] - mesh.vertices()[i];
        Vec3f f = r * g * massMod[i] * massMod[j] / pow(r.mag(), 3);
        // idk why the orbits look better w/o dividing by the mass
        acceleration[i] += 2 * (float)symmetry * f;       // / mass[i];
        acceleration[j] -= 2 * (1 - (float)symmetry) * f; // / mass[j];
      }
      // drag
      // acceleration[i] -= velocity[i] * 0.0001f;
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
