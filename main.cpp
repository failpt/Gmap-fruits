#include "lmap.h"
#include "tomato/grammar_1.h"

#include <CGAL/Graphics_scene_options.h>
#include <CGAL/draw_linear_cell_complex.h>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <QApplication>

using Tomato_gso = CGAL::Graphics_scene_options<
  lmap::LCC, lmap::LCC::Dart_const_descriptor, lmap::LCC::Dart_const_descriptor,
  lmap::LCC::Dart_const_descriptor, lmap::LCC::Dart_const_descriptor>;

Tomato_gso make_gso() {
    Tomato_gso gso;
    gso.colored_volume = [](const lmap::LCC&, lmap::LCC::Dart_const_descriptor) { return true; };
    gso.volume_color = [](const lmap::LCC& l, lmap::LCC::Dart_const_descriptor dh) { return l.template info<3>(dh).color; };
    return gso;
}

class Tomato_viewer : public CGAL::Qt::Basic_viewer {
public:
    using Base = CGAL::Qt::Basic_viewer;
    using Base::Base;

    void redraw() override { Base::redraw(); pin_sizes(); }

protected:
    void init() override { Base::init(); pin_sizes(); }

private:
    void pin_sizes() {
        size_edges(4.f);
        size_rays(4.f);
        size_lines(4.f);
        size_vertices(6.f);
    }
};

int main(int argc, char** argv) {
    int step = 3;

    const tomato_g_1::Parameters params;
    const tomato_g_1::Step which = static_cast<tomato_g_1::Step>(step);

    lmap::LMap model = tomato_g_1::build_round_tomato(params, which);

    CGAL::Graphics_scene buffer;
    CGAL::add_to_graphics_scene(model.lcc, buffer, make_gso());

    CGAL::Qt::init_ogl_context(4, 3);
    QApplication app(argc, argv);
    Tomato_viewer viewer(app.activeWindow(), buffer, "Grammar 1, Fig. 4");
    viewer.show();
    return app.exec();
}