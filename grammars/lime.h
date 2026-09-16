#ifndef GMAP_FRUITS_LIME_H
#define GMAP_FRUITS_LIME_H

#include "../lmap.h"

#include <string>

namespace lime_grammar {
  struct Parameters {
    double core_height = 64;
    double fruit_radius = 46;
    double flatness = 2.3;
    double angletop = 46;
    double anglebase = 46;
    int precision = 15;

    unsigned int seed_no = 7;
    double core_size = 10;
    double coretop_height = 9;

    double seed_ratio = 0.17;
    double albedo_size = 3.0;
    double peel_size = 1.5;
  };

  inline const CGAL::IO::Color core_color{240, 238, 220};
  inline const CGAL::IO::Color seed_color{186, 176, 125};
  inline const CGAL::IO::Color pulp_color{197, 230, 99};
  inline const CGAL::IO::Color albedo_color{228, 226, 202};
  inline const CGAL::IO::Color flavedo_color{144, 186, 32};

  inline int steps(const Parameters& p) { return p.precision + 6; }

  inline lmap::LMap build_lime(const Parameters& p = {}, int step = 0) {
    using namespace lmap;

    const int S = p.precision;

    Grammar g;

    g.define("precision", S);
    g.define("sector_no", p.seed_no * 2);
    g.define("fruit_radius", p.fruit_radius);
    g.define("flatness", p.flatness);
    g.define("core_radius", p.core_size / 2);
    g.define("core_disc_height", p.core_height / (S - 1));
    g.define("angle_step", (p.angletop + p.anglebase) / S);
    g.define("angletop", p.angletop);

    g.define("seed_ratio", p.seed_ratio);
    g.define("segment_ratio", 1 / p.seed_ratio - 1);
    g.define("albedo_size", p.albedo_size);
    g.define("peel_size", p.peel_size);
    g.define("wall_size", p.albedo_size + p.peel_size);

    g.define("seed_glue_stage", S);
    g.define("segment_add_stage", S + 1);
    g.define("segment_glue_stage", S + 2);
    g.define("albedo_add_stage", S + 3);
    g.define("albedo_glue_stage", S + 4);
    g.define("peel_add_stage", S + 5);
    g.define("peel_glue_stage", S + 6);

    const Volume_attributes unit(4, 0, 0, 0, 1, 1, 1);
    g.define_volume("SE", unit, seed_color);
    g.define_volume("PU", unit, pulp_color);
    g.define_volume("SG", unit, pulp_color);
    g.define_volume("AL", unit, albedo_color);
    g.define_volume("FL", unit, flavedo_color);
    g.define_volume("CO",
                    Volume_attributes(p.seed_no * 2, 0, 0, 0, p.coretop_height, p.core_size, p.core_size),
                    core_color);

    Frame axiom_frame = Frame::from_ez_ex(Point(0, 0, 0), Vector(0, 0, -1), Vector(1, 0, 0));
    g.set_axiom("CO",
                Volume_attributes(p.seed_no * 2, 0, 0, 0, p.coretop_height, p.core_size, p.core_size),
                axiom_frame);

    Rule p01;
    p01.name = "p01";
    p01.predecessor = "CO";
    p01.cond = [](const Context& c) { return c.stage < static_cast<int>(c.get("precision")); };
    p01.successor = [](Context& c) {
      c.add(Face_selector::end(), "CO", Attribute_overload{}.with_H(c.get("core_disc_height")));
    };
    g.add_rule(p01);

    Rule p02;
    p02.name = "p02";
    p02.predecessor = "CO";
    p02.cond = [](const Context& c) { return c.stage < static_cast<int>(c.get("precision")); };
    p02.block2 = [](Context& c) {
      const double u = 2.0 * (c.stage - 0.5) / (c.get("precision") - 1) - 1.0;
      const double n = c.get("flatness");
      const double shape = std::pow(std::max(0.0, 1.0 - std::pow(std::abs(u), n)), 1.0 / n);
      c.set("a", c.get("angletop") - c.stage * c.get("angle_step"));
      c.set("l", std::max(2.0, c.get("fruit_radius") * shape - c.get("core_radius")));
    };
    p02.successor = [](Context& c) {
      const int order = static_cast<int>(c.get("sector_no"));
      const double l = c.get("l");
      const double inside = std::max(0.25 * l, l - c.get("wall_size"));
      const double t = inside * c.get("seed_ratio");
      for (int i = 1; i <= order; ++i) {
        const bool seed = ((c.stage + i) % 2) == 0;
        c.add(Face_selector::side(static_cast<std::size_t>(i)), seed ? "SE" : "PU",
              Attribute_overload{}.with_aty(c.get("a")).with_H(t));
      }
    };
    g.add_rule(p02);

    const char* seed_band[2] = { "SE", "PU" };
    for (int i = 0; i < 2; ++i) {
      Rule r;
      r.name = "p04";
      r.predecessor = seed_band[i];
      r.cond = [](const Context& c) { return c.stage == static_cast<int>(c.get("segment_add_stage")); };
      r.successor = [](Context& c) {
        const double t = c.map->volume(c.self).attr.H * c.get("segment_ratio");
        c.add(Face_selector::end(), "SG", Attribute_overload{}.with_H(t));
      };
      g.add_rule(r);
    }

    const char* wall_from[2] = { "SG", "AL" };
    const char* wall_to[2] = { "AL", "FL" };
    const char* wall_thick[2] = { "albedo_size", "peel_size" };
    const char* wall_when[2] = { "albedo_add_stage", "peel_add_stage" };
    const char* wall_name[2] = { "p06", "p08" };
    for (int i = 0; i < 2; ++i) {
      const std::string to = wall_to[i], thick = wall_thick[i], when = wall_when[i];
      Rule r;
      r.name = wall_name[i];
      r.predecessor = wall_from[i];
      r.cond = [when](const Context& c) { return c.stage == static_cast<int>(c.get(when)); };
      r.successor = [to, thick](Context& c) {
        c.add(Face_selector::end(), to, Attribute_overload{}.with_H(c.get(thick)));
      };
      g.add_rule(r);
    }

    struct Glue {
      const char* name;
      const char* label;
      const char* other;
      const char* when;
    };
    const Glue glue[] = {
      { "p03a", "SE", "SE", "seed_glue_stage" },
      { "p03b", "SE", "PU", "seed_glue_stage" },
      { "p03c", "PU", "PU", "seed_glue_stage" },
      { "p05", "SG", "SG", "segment_glue_stage" },
      { "p07", "AL", "AL", "albedo_glue_stage" },
      { "p09", "FL", "FL", "peel_glue_stage" }
    };

    for (const Glue& G : glue) {
      const std::string other = G.other, when = G.when;
      Rule r;
      r.name = G.name;
      r.predecessor = G.label;
      r.cond = [when](const Context& c) { return c.stage == static_cast<int>(c.get(when)); };
      r.successor = [other](Context& c) {
        c.adjacency(Face_selector::all_sides(), other, Face_selector::all_sides());
      };
      g.add_rule(r);
    }

    LMap map;
    g.derive(map, step < 1 ? steps(p) : step);
    return map;
  }
}

#endif //GMAP_FRUITS_LIME_H