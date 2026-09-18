#include "lmap.h"
#include "grammars/tomato.h"
#include "grammars/cherry.h"
#include "grammars/lime.h"

#include <CGAL/Graphics_scene_options.h>
#include <CGAL/draw_linear_cell_complex.h>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <QApplication>

using Gso = CGAL::Graphics_scene_options<
  lmap::LCC, lmap::LCC::Dart_const_descriptor, lmap::LCC::Dart_const_descriptor,
  lmap::LCC::Dart_const_descriptor, lmap::LCC::Dart_const_descriptor>;

class Viewer : public CGAL::Qt::Basic_viewer {
public:
    using Basic_viewer::Basic_viewer;
    void redraw() override { Basic_viewer::redraw(); pin(); }

protected:
    void init() override { Basic_viewer::init(); pin(); }

private:
  void pin() { size_edges(15.f); size_rays(2.f); size_lines(2.f); size_vertices(4.f); }
};

std::function<bool(const lmap::LCC&, lmap::CDart)> lower_half(const lmap::LMap& model) {
  const auto bb = model.bbox();
  const double at = (bb.first.z() + bb.second.z()) / 2;
  return [at](const lmap::LCC& lcc, lmap::CDart dh) {
    double z = 0;
    std::size_t n = 0;
    for (auto it = lcc.one_dart_per_incident_cell<0, 3>(dh).begin(),
              e = lcc.one_dart_per_incident_cell<0, 3>(dh).end(); it != e; ++it) {
      z += lcc.point(it).z();
      ++n;
    }
    return n == 0 || z / n <= at;
  };
}

int show(int argc, char** argv, const lmap::LMap& model, const std::string& title, bool cut) {
  Gso gso;
  gso.colored_volume = [](const lmap::LCC&, lmap::CDart) { return true; };
  gso.volume_color = [](const lmap::LCC& l, lmap::CDart dh) { return l.info<3>(dh).color; };
  if (cut) gso.draw_volume = lower_half(model);

  CGAL::Graphics_scene buffer;
  CGAL::add_to_graphics_scene(model.lcc, buffer, gso);

  CGAL::Qt::init_ogl_context(4, 3);
  QApplication app(argc, argv);
  Viewer viewer(app.activeWindow(), buffer, title.c_str());
  viewer.show();
  return app.exec();
}

int pick(int step, int max, const std::string& condition = "") {
  if (step > max)
    throw std::runtime_error("This fruit's grammar only has 1-" + std::to_string(max) + " steps" + condition + ".");
  return step < 1 ? max : step;
}

int main(int argc, char** argv) {
  const bool cut = argc > 1 && std::string(argv[argc - 1]) == "cut";
  const int n = cut ? argc - 1 : argc;

  const std::string fruit = (n > 1) ? argv[1] : "";
  const int step = (n > 2) ? std::atoi(argv[2]) : 0;

  if (fruit == "tomato") {
    grammar_1::Parameters p;
    const std::string condition = " at precision " + std::to_string(p.precision);
    const int s = pick(step, grammar_1::steps(p), condition);
    return show(argc, argv, grammar_1::build_round_tomato(p, s), "Grammar 1, Fig. 4 - step", cut);
  }

  if (fruit == "cherry") {
    grammar_2::Parameters p;
    const int s = pick(step, grammar_2::steps(p));
    return show(argc, argv, grammar_2::build_cherry(p, s), "Grammar 2, Fig. 11 - step", cut);
  }

  if (fruit == "lime") {
    lime_grammar::Parameters p;
    const std::string condition = " at precision " + std::to_string(p.precision);
    const int s = pick(step, lime_grammar::steps(p), condition);
    return show(argc, argv, lime_grammar::build_lime(p, s), "Fig. 13, lime - step", cut);
  }

  throw std::runtime_error("Pick one of the fruits to display: tomato, cherry, lime.");
}