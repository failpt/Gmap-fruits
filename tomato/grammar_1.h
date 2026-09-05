#pragma once

#include "../lmap.h"

#include <cmath>
#include <string>

namespace tomato_g_1 {
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

  inline const CGAL::IO::Color kTrunkColor{60, 140, 55};
  inline const CGAL::IO::Color kPericarpColor{235, 120, 35};

  inline double raw_trunk_height(const Parameters& p) {
    const double pi = std::acos(-1.0);
    const double vhei = pi / p.precision;
    double s = 0;
    for (int k = 1; k < p.precision; ++k) s += std::sin(k * vhei) * p.height;
    return s;
  }

  enum class Step {
    Trunk = 1,
    Fan = 2,
    Fruit = 3
  };

  inline lmap::LMap build_round_tomato(const Parameters& p = {}, Step step = Step::Fruit) {
    const bool apply_p02 = (step != Step::Trunk);
    const bool apply_p03 = (step == Step::Fruit);
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

    const double kHeightSampling = (pi / p.precision) / 2.0;
    g.define("hsampling", kHeightSampling);

    g.define_volume("PE",
                    Volume_attributes(4, 0, 0, 0, p.pe_size, p.pe_size, p.pe_size),
                    kPericarpColor);
    g.define_volume("TR",
                    Volume_attributes(p.trunk_order, 0, 0, 0, p.trunk_size, p.trunk_size, p.trunk_size),
                    kTrunkColor);

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

    Rule p02;
    p02.name = "p02";
    p02.predecessor = "TR";
    p02.cond = [](const Context& c) { return c.stage < static_cast<int>(c.get("precision")); };
    p02.block2 = [](Context& c) {
      const double a = c.get("angletop") - c.stage * c.get("vang");
      const double x = c.stage * c.get("vsid") - 1.0;
      const double l = std::sqrt(std::max(0.0, 1.0 - x * x)) * c.get("width");
      c.set("a", a);
      c.set("l", l);
    };
    p02.successor = [](Context& c) {
      c.add(Face_selector::all_sides(), "PE", Attribute_overload{}.with_aty(c.get("a")).with_H(c.get("l")));
    };
    if (apply_p02) g.add_rule(p02);

    if (apply_p03) {
      Rule p03;
      p03.name = "p03";
      p03.predecessor = "PE";
      p03.cond = [](const Context& c) { return c.stage == static_cast<int>(c.get("precision")); };
      p03.successor = [](Context& c) {
        c.adjacency(Face_selector::all_sides(), "PE", Face_selector::all_sides());
      };
      g.add_rule(p03);
    }

    const int nb_stages = p.precision;
    LMap map;
    g.derive(map, nb_stages);
    return map;
  }
}