// Raphael Radna
// MAT-201B W20
// Assignment 4

/* 1. Choose your best work so far from the homeworks. Port this work to use
 * `al::DistributedAppWithState<SharedState>` rather than `al::App`. Submit this
 * by pushing `distributed-work.cpp`  to your repository.
 */

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

#include "al/app/al_DistributedApp.hpp" // #include "al/app/al_App.hpp"
#include "al/graphics/al_Font.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"

using namespace al;

#include <fstream>
#include <vector>

using namespace std;

// Distributed App provides a simple way to share states and parameters between
// a simulator and renderer apps using a single source file.
//
// You must define a state to be broadcast to all listeners. This state is
// synchronous but unreliable information i.e. missed data should not affect the
// overall state of the application.
//
struct Particle {
  Vec3f position, velocity, acceleration;
  float mass, massMod;
};
const int N = 1000;

struct SharedState {
  uint16_t frameCount{0};
  Particle particle[N]; // copy from a vector<Particle>
};

Vec3f rv(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
} // assigns randomness in three dimensions, multiplies it by scale

string slurp(string fileName); // forward declaration

// must do this instead of including a constructor in Particle
//
void initParticle(Particle &p, Vec3f p_, Vec3f v_, Vec3f a_, float m_,
                  float mM_);

// Inherit from DistributedApp and template it on the shared state data struct
//
class MyApp : public DistributedAppWithState<SharedState> {
  Parameter pointSize{"/pointSize", "", 0.5, "", 0.0, 5.0};
  Parameter timeStep{"/timeStep", "", 0.05, "", 0.01, 0.6};
  Parameter gravConst{"/gravConst", "", 0.75, "", 0, 1};
  Parameter dragFactor{"/dragFactor", "", 0.075, "", 0.01, 0.99};
  Parameter maxAccel{"/maxAccel", "", 10, "", 0, 20};
  Parameter scaleVal{"/scaleVal", "", 0.8, "", 0, 2};
  Parameter frequency{"/frequency", "", 1, "", 0, 5};
  Parameter depth{"/depth", "", 1, "", 0, 5};

  /* DistributedApp provides a parameter server. In fact it will
   * crash if you have a parameter server with the same port,
   * as it will raise an exception when it is unable to acquire
   * the port
   */

  ControlGUI gui;

  ShaderProgram pointShader;
  Mesh mesh; // simulation state position is located in the mesh (positions are
             // the direct simulation states that we use to draw)

  vector<Particle> particle;

  void reset() { // empty all containers
    mesh.reset();
    particle.clear();

    // c++11 "lambda" function
    // seed random number generators to maintain determinism
    rnd::Random<> rng;
    rng.seed(42);
    auto rc = [&]() { return HSV(rng.uniform() / 2 + 0.5f, 1.0f, 1.0f); };
    auto rv = [&](float scale) -> Vec3f {
      return Vec3f(rng.uniformS(), rng.uniformS(), rng.uniformS()) * scale;
    };

    mesh.primitive(Mesh::POINTS);

    for (int r = 0; r < N; r++) { // create 1000 points, put it into mesh
      Particle p;
      float m = 5 + rnd::uniformS() * 10;
      initParticle(p, rv(5), rv(0.1), Vec3f(0), m, m);

      mesh.vertex(p.position);
      mesh.color(rc());
      // set texture coordinate to be the size of the point (related to the
      // mass) using a simplified volume size relationship -> V = 4/3 * pi * r^3
      // pow is power -> m^(1/3)
      mesh.texCoord((4 / 3) * 3.14 * pow(p.mass, 1.0f / 3), 0);
      // pass in an s, t (like x, y)-> where on an image do we want to grab the
      // color from for this pixel (2D texture) normalized between 0 and 1

      particle.push_back(p);
    }
  }

  // You can keep a pointer to the cuttlebone domain
  // This can be useful to ask the domain if it is a sender or receiver
  std::shared_ptr<CuttleboneStateSimulationDomain<SharedState>>
      cuttleboneDomain;

  void onCreate() override {
    cuttleboneDomain =
        CuttleboneStateSimulationDomain<SharedState>::enableCuttlebone(this);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      quit();
    }

    gui << pointSize << timeStep << gravConst << dragFactor << maxAccel
        << scaleVal << frequency << depth;
    gui.init();

    // DistributedApp provides a parameter server.
    // This links the parameters between "simulator" and "renderers"
    // automatically
    parameterServer() << pointSize << timeStep << gravConst << dragFactor
                      << maxAccel << scaleVal << frequency << depth;

    navControl().useMouse(false);

    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    reset();
    nav().pos(0, 0, 10); // push camera back
  }

  bool freeze =
      false; // state that freezes simulation -> doesn't run onAnimate if frozen
  int t = 0;
  void onAnimate(double dt) override {
    if (cuttleboneDomain->isSender()) {
      state().frameCount++;
      navControl().active(!isImguiUsingInput());
    } else {
    }

    if (freeze)
      return;

    // ignore the real dt (real number of seconds that have passed since the
    // last call) and set the time step numerical simulation stability is due to
    // small, regular timesteps
    dt = timeStep; // simulation time, not wall time

    for (Particle p : particle) {
      p.massMod = p.mass * (float)depth +
                  sinf(((float)t * (float)frequency / 60.0f) * 6.283185308) *
                      p.mass * (float)depth +
                  0.01;
    }

    if (t >= 1 / (float)frequency * 60)
      t = 0;
    t++;

    // *********** Calculate forces ***********

    // gravity
    for (int i = 0; i < N; i++) {
      for (int j = 1 + i; j < N; j++) {

        rnd::Random<> rng;
        auto rv = [&](float scale) -> Vec3f {
          return Vec3f(rng.uniformS(), rng.uniformS(), rng.uniformS()) * scale;
        };

        Particle &p1(particle.at(i));
        Particle &p2(particle.at(j));
        // Vec3f distance(mesh.vertices()[j] - mesh.vertices()[i]);
        Vec3f distance(p2.position - p1.position);
        Vec3f gravityVal = gravConst * p1.massMod * p2.massMod *
                           distance.normalize() /
                           pow(distance.mag(), 2); // F = G * m1 * m2 / r^2

        // using a random multiplier here in order to show a more naturalistic,
        // "swimming" motion of the particles

        // acceleration[i] += gravityVal * rv(scaleVal) / massMod[i];
        // acceleration[j] -= gravityVal * rv(scaleVal) / massMod[j];
        p1.acceleration += gravityVal * rv(scaleVal) / p1.massMod;
        p2.acceleration -= gravityVal * rv(scaleVal) / p2.massMod;

        // drag -> stabilizes simulation

        // force of drag is proportional to the opposite of velocity * small
        // amount normally, it is v^2 -> can change the statement to
        // vector[i].mag() * velocity[i] * 0.07; take a bit of acceleration away
        // proportional to what the velocity is
        p1.acceleration -= p1.velocity * dragFactor;

        // limit acceleration
        //
        if (p1.acceleration.mag() > maxAccel)
          p1.acceleration.normalize(maxAccel);
      }
    }

    // // drag -> stabilizes simulation
    // for (int i = 0; i < N; i++) {
    //   // force of drag is proportional to the opposite of velocity * small
    //   // amount normally, it is v^2 -> can change the statement to
    //   // vector[i].mag() * velocity[i] * 0.07; take a bit of acceleration
    //   away
    //   // proportional to what the velocity is
    // acceleration[i] -= velocity[i] * dragFactor;
    // }

    // // limit acceleration
    // for (int i = 0; i < acceleration.size(); i++) {
    //   float m = acceleration[i].mag();
    //   if (m > maxAccel) {
    //     acceleration[i].normalize(maxAccel);
    //     // cout << "Limiting Acceleration: " << m << " -> " <<
    //     (float)maxAccel
    //     //      << endl;
    //   }
    // }

    // Integration -> don't mess with this

    // reference (alias) to mesh.vertices()
    //
    vector<Vec3f> &position(mesh.vertices());

    for (int i = 0; i < N; i++) {
      Particle &p(particle.at(i));

      p.velocity += p.acceleration / p.massMod * dt;
      p.position += p.velocity * dt;

      mesh.vertices()[i] = p.position;

      // clear all accelerations (IMPORTANT!!) -> accelerations have been used
      // and counted already
      //
      p.acceleration.zero();
    }

    // // clear all accelerations (IMPORTANT!!) -> accelerations have been used
    // and
    // // counted already
    // for (auto &a : acceleration)
    //   a.zero();

    // for (int i = 0; i < N; i++) {
    //   // "backward" Euler integration (semi-implicit, more stable)
    //   velocity[i] += acceleration[i] / massMod[i] * dt;
    //   position[i] +=
    //       velocity[i] *
    //       dt; // this is actually changing mesh.vertices() -> interchangeable

    //   // Explicit (or "forward") Euler integration would look like this:
    //   // position[i] += velocity[i] * dt;
    //   // velocity[i] += acceleration[i] / mass[i] * dt;
    // }

    // // clear all accelerations (IMPORTANT!!) -> accelerations have been used
    // and
    // // counted already
    // for (auto &a : acceleration)
    //   a.zero();
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '1') { // introduce some "random" forces
      for (int i = 0; i < N; i++)
        particle.at(i).acceleration = rv(1) / particle.at(i).massMod;
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

    // Draw th GUI on the simulator only
    if (isPrimary()) {
      gui.draw(g);
    }
  }
};

int main() {
  MyApp app;
  app.start();
}

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

void initParticle(Particle &p, Vec3f p_, Vec3f v_, Vec3f a_, float m_,
                  float mM_) {
  p.position = p_;
  p.velocity = v_;
  p.acceleration = a_;
  p.mass = m_;
  p.massMod = mM_;
}