#ifndef GMAP_FRUITS_CHERRY_H
#define GMAP_FRUITS_CHERRY_H

#include "../lmap.h"

#include <string>

namespace grammar_2 {
  struct Parameters {
    double stone_size = 50;
    int subdiv = 5;
    double layer_H = 50;
    double layer_L = 20;
    double layer_W = 50;
  };

  inline const CGAL::IO::Color kStoneColor{205, 175, 130};
  inline const CGAL::IO::Color kEndocarpColor{240, 205, 195};
  inline const CGAL::IO::Color kMesocarpColor{198, 48, 56};
  inline const CGAL::IO::Color kExocarpColor{138, 18, 30};

  enum class Step {
    Cube = 0,
    Sphere = 1,
    Spikes = 2,
    Shell = 5,
    Fruit = 7
  };

  inline lmap::LMap build_cherry(const Parameters& p = {}, Step step = Step::Fruit) {
    using namespace lmap;

    Grammar g;

    g.define("subdiv", p.subdiv);
    g.define("H", p.layer_H);

    const Volume_attributes layer(4, 0, 0, 0, p.layer_H, p.layer_L, p.layer_W);
    const Volume_attributes stone(4, 0, 0, 0, p.stone_size, p.stone_size, p.stone_size);

    g.define_volume("EN", layer, kEndocarpColor);
    g.define_volume("ME", layer, kMesocarpColor);
    g.define_volume("EX", layer, kExocarpColor);
    g.define_volume("STONE", stone, kStoneColor);

    const double h = p.stone_size / 2;
    g.set_axiom("STONE", stone,
                Frame::from_ez_ex(Point(-h, -h, -h), Vector(0, 0, 1), Vector(1, 0, 0)));

    Rule p01;
    p01.name = "p01";
    p01.predecessor = "STONE";
    p01.cond = [](const Context& c) { return c.stage == 1; };
    p01.successor = [](Context& c) {
      c.subdivide(static_cast<int>(c.get("subdiv")));
    };
    g.add_rule(p01);

    Rule p02;
    p02.name = "p02";
    p02.predecessor = "STONE";
    p02.cond = [](const Context& c) { return c.stage == 2; };
    p02.successor = [](Context& c) {
      c.add(Face_selector::all(), "EN", Attribute_overload{}.with_H(c.get("H")));
    };
    g.add_rule(p02);

    const char* from[2] = { "EN", "ME" };
    const char* to[2] = { "ME", "EX" };
    const int when_add[2] = { 3, 4 };
    const char* name_add[2] = { "p03", "p04" };
    for (int i = 0; i < 2; ++i) {
      const std::string next = to[i];
      const int when = when_add[i];
      Rule r;
      r.name = name_add[i];
      r.predecessor = from[i];
      r.cond = [when](const Context& c) { return c.stage == when; };
      r.successor = [next](Context& c) {
        c.add(Face_selector::end(), next, Attribute_overload{}.with_H(c.get("H")));
      };
      g.add_rule(r);
    }

    const char* glue[3] = { "EN", "ME", "EX" };
    const int when_glue[3] = { 5, 6, 7 };
    const char* name_glue[3] = { "p05", "p06", "p07" };
    for (int i = 0; i < 3; ++i) {
      const std::string label = glue[i];
      const int when = when_glue[i];
      Rule r;
      r.name = name_glue[i];
      r.predecessor = label;
      r.cond = [when](const Context& c) { return c.stage == when; };
      r.successor = [label](Context& c) {
        c.adjacency(Face_selector::all_sides(), label, Face_selector::all_sides());
      };
      g.add_rule(r);
    }

    LMap map;
    g.derive(map, static_cast<int>(step));
    return map;
  }
}

#endif //GMAP_FRUITS_CHERRY_H