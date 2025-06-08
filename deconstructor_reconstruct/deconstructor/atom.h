#pragma once
#include <string>

struct coordinate {
  double x, y, z;
  coordinate() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
  coordinate(double x, double y, double z) : x(x), y(y), z(z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }
};

struct atom {
  unsigned resID;
  std::string element;
  coordinate pos;

  atom() {
    resID = 0;
    element = "";
    pos = {0, 0, 0};
  }

  atom(int resID, std::string element, double x, double y, double z) {
    this->resID = resID;
    this->element = element;
    this->pos.x = x;
    this->pos.y = y;
    this->pos.z = z;
  }
};
