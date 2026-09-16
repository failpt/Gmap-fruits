#ifndef GMAP_FRUITS_TOMATO_H
#define GMAP_FRUITS_TOMATO_H

#include "../lmap.h"

#include <cmath>
#include <string>

namespace grammar_1 {
  struct Parameters {
    double height = 50;
    double width = 50;
    double angletop = 30;
    double anglebase = 10;
    int precision = 10;

    unsigned int trunk_order = 10;
    double trunk_size = 5;
    double pe_size = 10;
  };

  inline const CGAL::IO::Color trunk_color{60, 140, 55};
  inline const CGAL::IO::Color pericarp_color{235, 120, 35};

  inline double raw_trunk_height(const Parameters& p) {
    const double pi = std::acos(-1.0);
    const double vhei = pi / p.precision;
    double s = 0;
    for (int k = 1; k < p.precision; ++k) s += std::sin(k * vhei) * p.height;
    return s;
  }

  inline int steps(const Parameters& p) { return p.precision + 1; }

  inline lmap::LMap build_round_tomato(const Parameters& p = {}, int step = 0) {
    using namespace lmap;

    const double pi = std::acos(-1.0);

    Grammar g;

    g.define("height", p.height);
    g.define("width", p.width);
    g.define("precision", p.precision);
    g.define("angletop", p.angletop);
    g.define("anglebase", p.anglebase);
    g.define("vhei", pi / p.precision);
    g.define("vang", (p.angletop + p.anglebase) / p.precision);
    g.define("vsid", 2.0 / p.precision);

    const double height_sampling = (pi / p.precision) / 2.0;
    g.define("hsampling", height_sampling);

    g.define_volume("PE",
                    Volume_attributes(4, 0, 0, 0, p.pe_size, p.pe_size, p.pe_size),
                    pericarp_color);
    g.define_volume("TR",
                    Volume_attributes(p.trunk_order, 0, 0, 0, p.trunk_size, p.trunk_size, p.trunk_size),
                    trunk_color);

    // #axiom: TR(10, 0, 0, 0, 5, 5, 5)
    Frame axiom_frame = Frame::from_ez_ex(Point(0, 0, 0), Vector(0, 0, -1), Vector(1, 0, 0));
    g.set_axiom("TR",
                Volume_attributes(p.trunk_order, 0, 0, 0, p.trunk_size, p.trunk_size, p.trunk_size),
                axiom_frame);

    Rule p01;
    p01.name = "p01";
    p01.predecessor = "TR";
    p01.cond = [](const Context& c) { return c.stage < static_cast<int>(c.get("precision")); };
    p01.block2 = [](Context& c) {
      const double h = std::sin(c.stage * c.get("vhei")) * c.get("height") * c.get("hsampling");
      c.set("h", h);
    };
    p01.successor = [](Context& c) {
      c.add(Face_selector::end(), "TR", Attribute_overload{}.with_H(c.get("h")));
    };
    g.add_rule(p01);

    // The rule fires once, on every trunk volume at once, so it recovers each
    // ring's angle and length from the stage that created that volume: the
    // volume made at stage k carries ring k+1. The last trunk volume carries no
    // ring, which is what the old stage < precision guard used to express.
    Rule p02;
    p02.name = "p02";
    p02.predecessor = "TR";
    p02.cond = [](const Context& c) { return c.stage == static_cast<int>(c.get("precision")); };
    p02.successor = [](Context& c) {
      const int last = static_cast<int>(c.get("precision")) - 2;
      const int k = c.map->volume(c.self).stage;
      if (k > last) return;
      const double ring = k + 1;
      const double a = c.get("angletop") - ring * c.get("vang");
      const double x = ring * c.get("vsid") - 1.0;
      const double l = std::sqrt(std::max(0.0, 1.0 - x * x)) * c.get("width");
      c.add(Face_selector::all_sides(), "PE", Attribute_overload{}.with_aty(a).with_H(l));
    };
    g.add_rule(p02);

    Rule p03;
    p03.name = "p03";
    p03.predecessor = "PE";
    p03.cond = [](const Context& c) { return c.stage == static_cast<int>(c.get("precision")) + 1; };
    p03.successor = [](Context& c) {
      c.adjacency(Face_selector::all_sides(), "PE", Face_selector::all_sides());
    };
    g.add_rule(p03);

    LMap map;
    g.derive(map, step < 1 ? steps(p) : step);
    return map;
  }
}

#endif //GMAP_FRUITS_TOMATO_H