// Raphael Radna
// MAT-201B W20
// Assignment 3

/* 1. Using your code from the previous assignment as a guide, implement
 * [boids][] -style flocking. Submit by pushing `boids.cpp` to your github
 * repository.
 */

#include "al/app/al_App.hpp"
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp" // gui.draw(g)

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

Vec3f rv(float scale = 1.0f) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

string slurp(string fileName); // forward declaration

// suggested structure for Agent
//
struct Agent : Pose {
  Vec3f heading, center, separation; // of the local flock! averages
  unsigned flockCount{1};
};

struct AlloApp : App {
  // add more GUI here
  Parameter moveRate{"/moveRate", "", 1.0, "", 0.0, 2.0};
  Parameter turnRate{"/turnRate", "", 1.0, "", 0.0, 2.0};
  Parameter localRadius{"/localRadius", "", 0.4, "", 0.01, 0.9};
  Parameter size{"/size", "", 1.0, "", 0.0, 2.0};
  Parameter ratio{"/ratio", "", 1.0, "", 0.0, 2.0};
  Parameter align{"/align", "", 0.5, "", 0.0, 1.0};
  Parameter cohere{"/cohere", "", 0.5, "", 0.0, 1.0};
  Parameter separate{"/separate", "", 0.5, "", 0.0, 1.0};

  ControlGUI gui;

  ShaderProgram shader;
  Mesh mesh;

  vector<Agent> agent;

  void onCreate() override {
    // add more GUI here
    gui << moveRate << turnRate << localRadius << size << ratio << align
        << cohere << separate;
    gui.init();
    navControl().useMouse(false);

    // compile shaders
    shader.compile(slurp("../tetrahedron-vertex.glsl"),
                   slurp("../tetrahedron-fragment.glsl"),
                   slurp("../tetrahedron-geometry.glsl"));

    mesh.primitive(Mesh::POINTS);

    for (int _ = 0; _ < 1000; _++) {
      Agent a;
      a.pos(rv());
      a.faceToward(rv());
      agent.push_back(a);
      //
      mesh.vertex(a.pos());
      mesh.normal(a.uf());
      const Vec3f &up(a.uu());
      mesh.color(up.x, up.y, up.z);
    }

    nav().pos(0, 0, 10);
  }

  void onAnimate(double dt) override {
    int N = agent.size();
    for (unsigned i = 0; i < N; i++) {
      agent[i].flockCount = 1;
      for (unsigned j = 1 + i; j < N; j++) {
        float distance = (agent[j].pos() - agent[i].pos()).mag();
        if (distance < localRadius) {
          agent[i].heading += agent[j].uf();
          agent[i].center += agent[j].pos();
          agent[i].separation += agent[j].pos() - agent[i].pos();
          agent[i].flockCount++;
        }
      }
      agent[i].heading /= agent[i].flockCount;
      agent[i].center /= agent[i].flockCount;
      agent[i].separation /= agent[i].flockCount;
    }

    // only once the above loop is done do we have good data on average headings
    // and centers

    //
    // put code here
    //

    // move the agents along (KEEP THIS CODE)
    //
    for (unsigned i = 0; i < N; i++) {
      Vec3f a = agent[i].heading * (float)align;
      Vec3f c = (agent[i].center - agent[i].pos()).normalize() * (float)cohere;
      Vec3f s = (agent[i].separation * -1).normalize() * (float)separate;

      agent[i].faceToward(a + c + s, 0.03 * turnRate);
      agent[i].pos() += agent[i].uf() * moveRate * 0.002;
    }

    // respawn agents if they go too far (MAYBE KEEP)
    //
    for (unsigned i = 0; i < N; i++) {
      if (agent[i].pos().mag() > 1.1) {
        agent[i].pos(rv());
        agent[i].faceToward(rv());
      }
    }

    // visualize the agents
    //
    vector<Vec3f> &v(mesh.vertices());
    vector<Vec3f> &n(mesh.normals());
    vector<Color> &c(mesh.colors());
    for (unsigned i = 0; i < N; i++) {
      v[i] = agent[i].pos();
      n[i] = agent[i].uf();
      const Vec3d &up(agent[i].uu());
      c[i].set(up.x, up.y, up.z);
    }
  }

  void onDraw(Graphics &g) override {
    g.clear(0.0, 0.0, 0.0);
    // gl::depthTesting(true); // or g.depthTesting(true);
    // gl::blending(true); // or g.blending(true);
    // gl::blendTrans(); // or g.blendModeTrans();
    g.shader(shader);
    g.shader().uniform("size", size * 0.03);
    g.shader().uniform("ratio", ratio * 0.2);
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