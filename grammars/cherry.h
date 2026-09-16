#ifndef GMAP_FRUITS_CHERRY_H
#define GMAP_FRUITS_CHERRY_H

#include "../lmap.h"

#include <string>

namespace grammar_2 {
  struct Parameters {
    double stone_size = 130;
    int subdiv = 5;

    double exo_size = 5;
    double endo_size = 50;
    double meso_size = 50;
  };

  inline const CGAL::IO::Color stone_color{205, 175, 130};
  inline const CGAL::IO::Color endocarp_color{198, 30, 46};
  inline const CGAL::IO::Color mesocarp_color{198, 48, 56};
  inline const CGAL::IO::Color exocarp_color{138, 18, 30};

  // One stage subdivides the stone, then each of the three layers takes one
  // stage to be added and one to be glued.
  inline int steps(const Parameters&) { return 7; }

  inline lmap::LMap build_cherry(const Parameters& p = {}, int step = 0) {
    using namespace lmap;

    Grammar g;

    g.define("subdiv", p.subdiv);

    const Volume_attributes meso(4, 0, 0, 0, p.meso_size, p.meso_size, p.meso_size);
    const Volume_attributes endo(4, 0, 0, 0, p.endo_size, p.endo_size, p.endo_size);
    const Volume_attributes exo(4, 0, 0, 0, p.exo_size, p.exo_size, p.exo_size);
    const Volume_attributes stone(4, 0, 0, 0, p.stone_size, p.stone_size, p.stone_size);

    g.define_volume("EN", endo, endocarp_color);
    g.define_volume("ME", meso, mesocarp_color);
    g.define_volume("EX", exo, exocarp_color);
    g.define_volume("STONE", stone, stone_color);

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
      c.add(Face_selector::all(), "EN");
    };
    g.add_rule(p02);

    const char* from[2] = { "EN", "ME" };
    const char* to[2] = { "ME", "EX" };
    const int when_add[2] = { 4, 6 };
    const char* name_add[2] = { "p03", "p04" };
    for (int i = 0; i < 2; ++i) {
      const std::string next = to[i];
      const int when = when_add[i];
      Rule r;
      r.name = name_add[i];
      r.predecessor = from[i];
      r.cond = [when](const Context& c) { return c.stage == when; };
      r.successor = [next](Context& c) {
        c.add(Face_selector::end(), next);
      };
      g.add_rule(r);
    }

    const char* glue[3] = { "EN", "ME", "EX" };
    const int when_glue[3] = { 3, 5, 7 };
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
    g.derive(map, step < 1 ? steps(p) : step);
    return map;
  }
}

#endif //GMAP_FRUITS_CHERRY_H