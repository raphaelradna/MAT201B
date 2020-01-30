//Raphael Radna
//MAT-201B W20
//Assignment 1

#include "al/app/al_App.hpp"
#include "al/graphics/al_Image.hpp"  // al::Image
#include "al/math/al_Random.hpp"
#include "al/ui/al_ControlGUI.hpp"  // gui.draw(g)

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

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

float scale(float value=0, float inLow=0, float inHigh=1, float outLow=0, float outHigh=1) {
  float normValue = (value - inLow) / (inHigh - inLow);
  return (normValue * (outHigh - outLow)) + outLow;
}

float lerp(float a, float b, float t) {
  return a + t * (b - a);
}

Vec3f lerpVec3f(Vec3f start, Vec3f end, float p) {
  return Vec3f(lerp(start[0], end[0], p), lerp(start[1], end[1], p), lerp(start[2], end[2], p));
}

Vec3f rgbToHSV(Image::RGBAPix px) { //after https://www.rapidtables.com/convert/color/rgb-to-hsv.html
  float rgb[] = {px.r / 255.0f, px.g / 255.0f, px.b / 255.0f};
  float h = 0;
  float s = 0;
  float v = 0;
  float *cMax = max_element(rgb, rgb + 3);
  float *cMin = min_element(rgb, rgb + 3);
  float delta = *cMax - *cMin;

  if (cMax == &rgb[0])
    h = 60 * fmod(((rgb[1] - rgb[2]) / delta), 6);
  if (cMax == &rgb[1])
    h = 60 * (((rgb[2] - rgb[0]) / delta) + 2);
  if (cMax == &rgb[2])
    h = 60 * (((rgb[0] - rgb[1]) / delta) + 4);

  if (*cMax != 0) {
    s = delta / *cMax;
  } else {
    s = 0;
  }

  v = *cMax;

  return Vec3f(h, s, v);
}

float degToRad(float d) {
  return d * 3.141593 / 180;
}

struct AlloApp : App {
  Parameter pointSize{"/pointSize", "", 0.268, "", 0.0, 1.0};
  ControlGUI gui;
  ShaderProgram pointShader;
  Mesh mesh;
  vector<Vec3f> img, rgb, hsv, rr, *shape;

  void onCreate() override {
    gui << pointSize;
    gui.init();
    navControl().useMouse(false);

    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    mesh.primitive(Mesh::POINTS);

    const char* filename = "../color_spaces.png";
    auto imageData = Image(filename);
    if (imageData.array().size() == 0) {
      std::cout << "failed to load image" << std::endl;
      exit(1);
    }

    Image::RGBAPix pixel;
    const int W = imageData.width();
    const int H = imageData.height();
    for (int c = 0; c < W; c++)
      for (int r = 0; r < H; r++) {
        imageData.read(pixel, c, r);
        Vec3f imgPos((c - W / 2) * 0.005, (r - H / 2) * 0.005, 0);
        img.push_back(imgPos);
        Vec3f rgbPos(scale(pixel.r, 0, 255, -2.5, 2.5), scale(pixel.g, 0, 255, -2.5, 2.5), scale(pixel.b, 0, 255, -2.5, 2.5));
        rgb.push_back(rgbPos);
        Vec3f hsv_ = rgbToHSV(pixel);
        Vec3f hsvPos(hsv_[1] * sinf(degToRad(hsv_[0])) * 2.5, scale(hsv_[2], 0, 1, -2.5, 2.5), hsv_[1] * cosf(degToRad(hsv_[0])) * 2.5); //cylindrical to cartesian coordinates, after http://tutorial.math.lamar.edu/Classes/CalcIII/CylindricalCoords.aspx
        hsv.push_back(hsvPos);
        Vec3f rrPos(rgbPos[0] * hsvPos[0] / 2.5, rgbPos[1] * hsvPos[1] / 2.5, rgbPos[2] * hsvPos[2] / 2.5);
        rr.push_back(rrPos);
        mesh.vertex(imgPos);
        mesh.color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);
      }
    nav().pos(0, 0, 15);
  }

  bool doLerp = false;
  int t = 0;
  void onAnimate(double dt) override {
    if (doLerp) {
      if (t < 60) {
        for (int i = 0; i < mesh.vertices().size(); i++) {
          mesh.vertices()[i] = lerpVec3f(mesh.vertices()[i], shape->at(i), t / 60.0f);
        }
        t++;
      } else {
        doLerp = false;
      }
    }
  }

  // I tried to encapsulate the 'keyDown' behavior in a function; not sure why doesn't it work right.
  //
  void changeShape(vector<Vec3f> _shape) {
    shape = &_shape;
    t = 0;
    doLerp = true;
  }

  bool onKeyDown(const Keyboard& k) override {
    if (k.key() == '1') {
      // changeShape(img);
      shape = &img;
      t = 0;
      doLerp = true;
    }
    if (k.key() == '2') {
        shape = &rgb;
        t = 0;
        doLerp = true;
    }
    if (k.key() == '3') {
        shape = &hsv;
        t = 0;
        doLerp = true;
    }
    if (k.key() == '4') {
        shape = &rr;
        t = 0;
        doLerp = true;
    }
    return true;
  }

  void onDraw(Graphics& g) override {
    g.clear(0.01);
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.depthTesting(true);
    g.draw(mesh);
    gui.draw(g);
  }
};

int main() { AlloApp().start(); }
