#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Cell_attribute_with_point.h>
#include <CGAL/IO/Color.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace lmap {
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Traits = CGAL::Linear_cell_complex_traits<3, Kernel>;

  struct Vertex_info {
    unsigned int weight = 1;
  };

  struct Barycentric_merge {
    template <class Cell_attribute>
    void operator()(Cell_attribute& a, Cell_attribute& b) const {
      using Pt = typename Cell_attribute::Point;
      const double wa = static_cast<double>(a.info().weight);
      const double wb = static_cast<double>(b.info().weight);
      const auto& pa = a.point();
      const auto& pb = b.point();
      a.point() = Pt((wa * pa.x() + wb * pb.x()) / (wa + wb),
                     (wa * pa.y() + wb * pb.y()) / (wa + wb),
                     (wa * pa.z() + wb * pb.z()) / (wa + wb));
      a.info().weight = static_cast<unsigned int>(wa + wb);
    }
  };

  struct Volume_info {
    std::string label;
    CGAL::IO::Color color{200, 200, 200};
    int stage = 0;
  };

  struct LMap_items {
    template <class Refs>
    struct Dart_wrapper {
      using Vertex_attribute =
        CGAL::Cell_attribute_with_point<Refs, Vertex_info, CGAL::Tag_true,
                                        Barycentric_merge>;
      using Volume_attribute =
        CGAL::Cell_attribute<Refs, Volume_info, CGAL::Tag_true>;
      using Attributes = std::tuple<Vertex_attribute, void, void, Volume_attribute>;
    };
  };

  using LCC = CGAL::Linear_cell_complex_for_generalized_map<3, 3, Traits, LMap_items>;
  using Dart = LCC::Dart_descriptor;
  using CDart = LCC::Dart_const_descriptor;
  using Point = LCC::Point;
  using Vector = LCC::Vector;

  namespace geom {
    inline double norm(const Vector& v) { return std::sqrt(CGAL::to_double(v.squared_length())); }

    inline Vector normalized(const Vector& v) {
      const double n = norm(v);
      if (n < 1e-14) return Vector(0, 0, 0);
      return v / n;
    }

    inline Vector newell_normal(const std::vector<Point>& p) {
      double nx = 0, ny = 0, nz = 0;
      const std::size_t m = p.size();
      for (std::size_t i = 0; i < m; ++i) {
        const Point& a = p[i];
        const Point& b = p[(i + 1) % m];
        nx += (a.y() - b.y()) * (a.z() + b.z());
        ny += (a.z() - b.z()) * (a.x() + b.x());
        nz += (a.x() - b.x()) * (a.y() + b.y());
      }
      return Vector(nx, ny, nz);
    }
  }

  struct Frame {
    Point origin{0, 0, 0};
    Vector ex{1, 0, 0};
    Vector ey{0, 1, 0};
    Vector ez{0, 0, 1};

    static Frame from_ez_ex(const Point& o, const Vector& ez_in, const Vector& ex_hint) {
      Frame f;
      f.origin = o;
      f.ez = geom::normalized(ez_in);
      Vector ex = ex_hint - (ex_hint * f.ez) * f.ez;
      if (geom::norm(ex) < 1e-9) {
        Vector alt(1, 0, 0);
        ex = alt - (alt * f.ez) * f.ez;
        if (geom::norm(ex) < 1e-9) { alt = Vector(0, 1, 0); ex = alt - (alt * f.ez) * f.ez; }
      }
      f.ex = geom::normalized(ex);
      f.ey = CGAL::cross_product(f.ez, f.ex);
      return f;
    }
  };

  struct Volume_attributes {
    unsigned int order = 4;
    double atx = 0, aty = 0, atz = 0;
    double H = 1, L = 1, W = 1;

    Volume_attributes() = default;
    Volume_attributes(unsigned int O_, double ax, double ay, double az,
                      double H_, double L_, double W_)
      : order(O_), atx(ax), aty(ay), atz(az), H(H_), L(L_), W(W_) {}
  };

  struct Attribute_overload {
    std::optional<unsigned int> order;
    std::optional<double> atx, aty, atz, H, L, W;

    Attribute_overload& with_order(unsigned int v) { order = v; return *this; }
    Attribute_overload& with_atx(double v) { atx = v; return *this; }
    Attribute_overload& with_aty(double v) { aty = v; return *this; }
    Attribute_overload& with_atz(double v) { atz = v; return *this; }
    Attribute_overload& with_H(double v) { H = v; return *this; }
    Attribute_overload& with_L(double v) { L = v; return *this; }
    Attribute_overload& with_W(double v) { W = v; return *this; }

    void apply_to(Volume_attributes& a) const {
      if (order) a.order = *order;
      if (atx) a.atx = *atx;
      if (aty) a.aty = *aty;
      if (atz) a.atz = *atz;
      if (H) a.H = *H;
      if (L) a.L = *L;
      if (W) a.W = *W;
    }
  };

  enum class Face_role { Origin, End, Side };

  using Volume_id = std::size_t;
  inline constexpr Volume_id null_volume = static_cast<Volume_id>(-1);

  struct Face_ref {
    Volume_id volume = null_volume;
    Face_role role = Face_role::Origin;
    std::size_t side_index = 0;

    bool operator==(const Face_ref& o) const {
      return volume == o.volume && role == o.role && side_index == o.side_index;
    }
  };

  class LMap;

  // 4.1.2: A_O, A_E, A_{C1}, A_{C*}, A_{C*-2}, A_*
  class Face_selector {
  public:
    static Face_selector origin() { Face_selector s; s.m_kind = Kind::Origin; return s; }
    static Face_selector end() { Face_selector s; s.m_kind = Kind::End; return s; }

    static Face_selector side(std::size_t i) {
      Face_selector s; s.m_kind = Kind::One_side; s.m_index = i; return s;
    }

    static Face_selector all_sides() {
      Face_selector s; s.m_kind = Kind::All_sides; return s;
    }

    static Face_selector all_sides_except(std::initializer_list<std::size_t> ex) {
      Face_selector s; s.m_kind = Kind::All_sides; s.m_excluded.assign(ex); return s;
    }

    static Face_selector all() { Face_selector s; s.m_kind = Kind::All; return s; }

    std::vector<Face_ref> resolve(const LMap& map, Volume_id v) const;

    std::string to_string() const {
      switch (m_kind) {
      case Kind::Origin: return "O";
      case Kind::End: return "E";
      case Kind::One_side: return "C" + std::to_string(m_index);
      case Kind::All: return "*";
      case Kind::All_sides: {
          std::string s = "C*";
          for (std::size_t i : m_excluded) s += "-" + std::to_string(i);
          return s;
      }
      }
      return "?";
    }

  private:
    enum class Kind { Origin, End, One_side, All_sides, All };
    Kind m_kind = Kind::All;
    std::size_t m_index = 1;
    std::vector<std::size_t> m_excluded;
  };

  struct Volume {
    Volume_id id = null_volume;
    std::string label;
    Volume_attributes attr;
    Frame frame;
    Dart d_origin{};
    Dart d_end{};
    std::vector<Dart> d_side;
    int stage = 0;
    CGAL::IO::Color color{200, 200, 200};
  };

  class LMap {
  public:
    LCC lcc;

    std::size_t nb_volumes() const { return m_volumes.size(); }
    const Volume& volume(Volume_id v) const { return m_volumes.at(v); }
    Volume& volume(Volume_id v) { return m_volumes.at(v); }
    const std::vector<Volume>& volumes() const { return m_volumes; }

    std::vector<Volume_id> volumes_with_label(const std::string& lab) const {
      std::vector<Volume_id> r;
      for (const Volume& v : m_volumes) if (v.label == lab) r.push_back(v.id);
      return r;
    }

    Dart dart_of(const Face_ref& f) const {
      const Volume& v = m_volumes.at(f.volume);
      switch (f.role) {
      case Face_role::Origin: return v.d_origin;
      case Face_role::End: return v.d_end;
      case Face_role::Side:
        if (f.side_index == 0 || f.side_index > v.d_side.size())
          throw std::out_of_range("lmap: side face index out of range");
        return v.d_side[f.side_index - 1];
      }
      throw std::logic_error("lmap: bad face role");
    }

    std::vector<Point> face_loop(CDart d) const {
      std::vector<Point> p;
      CDart c = d;
      do { p.push_back(lcc.point(c)); c = lcc.next(c); } while (c != d);
      return p;
    }

    Vector face_outward_normal(CDart d) const {
      return geom::normalized(geom::newell_normal(face_loop(d)));
    }

    bool is_glued(const Face_ref& f) const { return !lcc.is_free<3>(dart_of(f)); }

    bool are_adjacent(const Face_ref& f, const Face_ref& g) const {
      const Dart df = dart_of(f), dg = dart_of(g);
      if (lcc.belong_to_same_cell<2>(df, dg)) return false;
      return find_common_edge_darts(df, dg).has_value();
    }

    bool are_adjacent(Volume_id a, Volume_id b) const {
      if (a == b) return false;
      for (const Face_ref& f : Face_selector::all().resolve(*this, a))
        for (const Face_ref& g : Face_selector::all().resolve(*this, b))
          if (are_adjacent(f, g)) return true;
      return false;
    }

    Volume_id add_axiom(const std::string& label,
                        const Volume_attributes& attr_in,
                        const Frame& frame_in,
                        const CGAL::IO::Color& color = CGAL::IO::Color(200, 200, 200),
                        int stage = 0) {
      Volume_attributes attr = attr_in;
      const Frame f = rotate_frame(frame_in, attr);
      const std::vector<Point> base = regular_base(f, attr.order, attr.L, attr.W);
      std::vector<Point> top;
      top.reserve(base.size());
      for (const Point& p : base) top.push_back(p + attr.H * f.ez);
      return emplace_prism(base, top, label, attr, f, color, stage);
    }

    std::vector<Volume_id> add(Volume_id support,
                               const Face_selector& on_faces,
                               const std::string& new_label,
                               const Volume_attributes& defaults,
                               const Attribute_overload& overload = {},
                               const CGAL::IO::Color& color = CGAL::IO::Color(200, 200, 200),
                               int stage = 0) {
      std::vector<Volume_id> created;
      for (const Face_ref& f : on_faces.resolve(*this, support)) {
        const Dart ds = dart_of(f);

        if (!lcc.is_free<3>(ds)) continue;

        const std::vector<Point> P = face_loop(ds);
        const std::size_t m = P.size();

        const Vector ez = geom::normalized(geom::newell_normal(P));
        Frame base_frame = Frame::from_ez_ex(P[0], ez, Vector(0, 0, 1));

        Volume_attributes attr = defaults;
        overload.apply_to(attr);
        attr.order = static_cast<unsigned int>(m);

        const Frame f_rot = rotate_frame(base_frame, attr);

        if (std::abs(attr.H) < 1e-12) {
          std::cerr << "lmap: skipping degenerate addition of '" << new_label
                    << "' (H == 0) on " << on_faces.to_string() << "\n";
          continue;
        }

        const std::vector<Point>& base = P;
        std::vector<Point> top;
        top.reserve(m);
        for (const Point& p : base) top.push_back(p + attr.H * f_rot.ez);

        const Volume_id nid =
          emplace_prism(base, top, new_label, attr, f_rot, color, stage);

        const Dart d2 = lcc.alpha<1>(m_volumes[nid].d_origin);

        CGAL_assertion(lcc.point(ds) == lcc.point(d2));
        CGAL_assertion(lcc.point(lcc.alpha<0>(ds)) == lcc.point(lcc.alpha<0>(d2)));

        if (!lcc.is_sewable<3>(ds, d2)) {
          std::cerr << "lmap: addition of '" << new_label
                    << "' on volume " << support << ' ' << on_faces.to_string()
                    << " is not 3-sewable.\n";
          continue;
        }
        lcc.sew<3>(ds, d2);

        CGAL_assertion(equal_loops(face_loop(ds), P));

        created.push_back(nid);
      }
      return created;
    }

    struct Adjacency_candidate {
      Face_ref f, g;
      Dart df, dg;
    };

    std::vector<Adjacency_candidate>
    match_adjacencies(Volume_id a, const Face_selector& fa,
                      Volume_id b, const Face_selector& fb) const {
      std::vector<Adjacency_candidate> out;
      if (a == b) return out;
      if (!share_a_vertex(a, b)) return out;

      const std::vector<Face_ref> A = fa.resolve(*this, a);
      const std::vector<Face_ref> B = fb.resolve(*this, b);
      for (const Face_ref& f : A) {
        const Dart df = dart_of(f);
        if (!lcc.is_free<3>(df)) continue;
        for (const Face_ref& g : B) {
          const Dart dg = dart_of(g);
          if (!lcc.is_free<3>(dg)) continue;
          if (lcc.belong_to_same_cell<2>(df, dg)) continue;

          const auto pair = find_common_edge_darts(df, dg);
          if (!pair) continue;

          CGAL_assertion(face_outward_normal(df) * face_outward_normal(dg) < 0);

          out.push_back(Adjacency_candidate{f, g, pair->first, pair->second});
        }
      }
      return out;
    }

    std::size_t apply_adjacencies(const std::vector<Adjacency_candidate>& cands) {
      std::size_t glued = 0;
      for (const Adjacency_candidate& c : cands) {
        if (!lcc.is_free<3>(c.df)) continue;
        if (!lcc.is_free<3>(c.dg)) continue;

        if (lcc.is_sewable<3>(c.df, c.dg)) { lcc.sew<3>(c.df, c.dg); ++glued; } else { ++m_failed_gluings; }
      }
      return glued;
    }

    std::size_t create_adjacency(Volume_id a, const Face_selector& fa, Volume_id b, const Face_selector& fb) {
      return apply_adjacencies(match_adjacencies(a, fa, b, fb));
    }

    std::size_t create_adjacency_between_labels(const std::string& la, const Face_selector& sa, const std::string& lb, const Face_selector& sb) {
      const std::vector<Volume_id> A = volumes_with_label(la);
      const std::vector<Volume_id> B = volumes_with_label(lb);
      std::vector<Adjacency_candidate> all;
      for (Volume_id a : A)
        for (Volume_id b : B) {
          if (a == b) continue;
          if (la == lb && b <= a) continue;
          const auto c = match_adjacencies(a, sa, b, sb);
          all.insert(all.end(), c.begin(), c.end());
        }
      dedupe_candidates(all);
      return apply_adjacencies(all);
    }

    static void dedupe_candidates(std::vector<Adjacency_candidate>& v) {
      auto key = [](const Face_ref& r) {
        return std::make_tuple(r.volume, static_cast<int>(r.role), r.side_index);
      };
      std::set<std::pair<decltype(key(std::declval<Face_ref>())), decltype(key(std::declval<Face_ref>()))>> seen;
      std::vector<Adjacency_candidate> out;
      out.reserve(v.size());
      for (const Adjacency_candidate& c : v) {
        auto kf = key(c.f), kg = key(c.g);
        auto pr = (kg < kf) ? std::make_pair(kg, kf) : std::make_pair(kf, kg);
        if (seen.insert(pr).second) out.push_back(c);
      }
      v.swap(out);
    }

    std::size_t failed_gluings() const { return m_failed_gluings; }

    std::map<std::string, std::size_t> volume_count_per_label() const {
      std::map<std::string, std::size_t> r;
      for (const Volume& v : m_volumes) ++r[v.label];
      return r;
    }

    std::size_t nb_boundary_faces() const {
      std::size_t n = 0;
      for (auto it = lcc.one_dart_per_cell<2>().begin(), e = lcc.one_dart_per_cell<2>().end(); it != e; ++it)
        if (lcc.is_free<3>(it)) ++n;
      return n;
    }

    std::pair<Point, Point> bbox() const {
      double lo[3] = { std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max() };
      double hi[3] = { std::numeric_limits<double>::lowest(),
                       std::numeric_limits<double>::lowest(),
                       std::numeric_limits<double>::lowest() };
      for (auto it = lcc.one_dart_per_cell<0>().begin(), e = lcc.one_dart_per_cell<0>().end(); it != e; ++it) {
        const Point& p = lcc.point(it);
        const double c[3] = { p.x(), p.y(), p.z() };
        for (int i = 0; i < 3; ++i) { lo[i] = std::min(lo[i], c[i]); hi[i] = std::max(hi[i], c[i]); }
      }
      return { Point(lo[0], lo[1], lo[2]), Point(hi[0], hi[1], hi[2]) };
    }

  private:
    std::vector<Volume> m_volumes;
    std::size_t m_failed_gluings = 0;

    static std::vector<Point> regular_base(const Frame& f, unsigned int n, double L, double W) {
      if (n < 3) throw std::invalid_argument("lmap: prism order must be >= 3");
      const double pi = std::acos(-1.0);
      double bx = 0, by = 0;
      std::vector<double> cs(n), sn(n);
      for (unsigned int k = 0; k < n; ++k) {
        const double phi = 2.0 * pi * k / n + pi / n;
        cs[k] = std::cos(phi);
        sn[k] = std::sin(phi);
        bx = std::max(bx, 2.0 * std::abs(cs[k]));
        by = std::max(by, 2.0 * std::abs(sn[k]));
      }
      std::vector<Point> base;
      base.reserve(n);
      for (unsigned int k = 0; k < n; ++k)
        base.push_back(f.origin + (L / bx) * cs[k] * f.ex + (W / by) * sn[k] * f.ey);
      return base;
    }

    static Frame rotate_frame(const Frame& f, const Volume_attributes& a) {
      const double d2r = std::acos(-1.0) / 180.0;
      const double ca = std::cos(a.atx * d2r), sa = std::sin(a.atx * d2r);
      const double cb = std::cos(a.aty * d2r), sb = std::sin(a.aty * d2r);
      const double cc = std::cos(a.atz * d2r), sc = std::sin(a.atz * d2r);

      const double X[3] = { cc * cb, sc * cb, -sb };
      const double Y[3] = { cc * sb * sa - sc * ca, sc * sb * sa + cc * ca, cb * sa };
      const double Z[3] = { cc * sb * ca + sc * sa, sc * sb * ca - cc * sa, cb * ca };

      Frame r;
      r.origin = f.origin;
      r.ex = geom::normalized(X[0] * f.ex + X[1] * f.ey + X[2] * f.ez);
      r.ey = geom::normalized(Y[0] * f.ex + Y[1] * f.ey + Y[2] * f.ez);
      r.ez = geom::normalized(Z[0] * f.ex + Z[1] * f.ey + Z[2] * f.ez);
      return r;
    }

    Volume_id emplace_prism(const std::vector<Point>& base,
                            const std::vector<Point>& top,
                            const std::string& label,
                            const Volume_attributes& attr,
                            const Frame& frame,
                            const CGAL::IO::Color& color,
                            int stage) {
      const std::size_t n = base.size();
      if (n < 3 || top.size() != n)
        throw std::invalid_argument("lmap: bad prism base/top");

      CGAL::Linear_cell_complex_incremental_builder_3<LCC> builder(lcc);
      for (std::size_t i = 0; i < n; ++i) builder.add_vertex(base[i]);
      for (std::size_t i = 0; i < n; ++i) builder.add_vertex(top[i]);

      builder.begin_surface();

      Volume v;
      v.id = m_volumes.size();
      v.label = label;
      v.attr = attr;
      v.frame = frame;
      v.stage = stage;
      v.color = color;

      builder.begin_facet();
      builder.add_vertex_to_facet(0);
      for (std::size_t i = n - 1; i >= 1; --i) builder.add_vertex_to_facet(i);
      v.d_origin = builder.end_facet();

      builder.begin_facet();
      for (std::size_t i = 0; i < n; ++i) builder.add_vertex_to_facet(n + i);
      v.d_end = builder.end_facet();

      v.d_side.reserve(n);
      for (std::size_t k = 0; k < n; ++k) {
        builder.begin_facet();
        builder.add_vertex_to_facet(k);
        builder.add_vertex_to_facet((k + 1) % n);
        builder.add_vertex_to_facet(n + ((k + 1) % n));
        builder.add_vertex_to_facet(n + k);
        v.d_side.push_back(builder.end_facet());
      }

      builder.end_surface();

      lcc.set_attribute<3>(v.d_origin, lcc.create_attribute<3>(Volume_info{label, color, stage}));

      m_volumes.push_back(v);
      return v.id;
    }

    std::optional<std::pair<Dart, Dart>>
    find_common_edge_darts(Dart df, Dart dg) const {
      LCC& m = const_cast<LCC&>(lcc);
      auto rf = m.darts_of_cell<2, 3>(df);
      auto rg = m.darts_of_cell<2, 3>(dg);
      for (auto i1 = rf.begin(), e1 = rf.end(); i1 != e1; ++i1)
        for (auto i2 = rg.begin(), e2 = rg.end(); i2 != e2; ++i2)
          if (lcc.belong_to_same_cell<1>(i1, i2) && lcc.belong_to_same_cell<0>(i1, i2))
            return std::make_pair(Dart(i1), Dart(i2));
      return std::nullopt;
    }

    bool share_a_vertex(Volume_id a, Volume_id b) const {
      LCC& m = const_cast<LCC&>(lcc);
      auto ra = m.darts_of_cell<3>(m_volumes[a].d_origin);
      auto rb = m.darts_of_cell<3>(m_volumes[b].d_origin);
      std::unordered_set<const void*> va;
      for (auto it = ra.begin(), e = ra.end(); it != e; ++it)
        va.insert(static_cast<const void*>(&lcc.point(CDart(it))));
      for (auto it = rb.begin(), e = rb.end(); it != e; ++it)
        if (va.count(static_cast<const void*>(&lcc.point(CDart(it))))) return true;
      return false;
    }

    static bool equal_loops(const std::vector<Point>& a, const std::vector<Point>& b) {
      if (a.size() != b.size()) return false;
      for (std::size_t i = 0; i < a.size(); ++i)
        if ((a[i] - b[i]).squared_length() > 1e-18) return false;
      return true;
    }
  };

  inline std::vector<Face_ref> Face_selector::resolve(const LMap& map, Volume_id v) const {
    const Volume& vol = map.volume(v);
    const std::size_t n = vol.d_side.size();
    std::vector<Face_ref> r;

    auto push_sides = [&]() {
      for (std::size_t i = 1; i <= n; ++i) {
        if (std::find(m_excluded.begin(), m_excluded.end(), i) != m_excluded.end())
          continue;
        r.push_back(Face_ref{v, Face_role::Side, i});
      }
    };

    switch (m_kind) {
    case Kind::Origin: r.push_back(Face_ref{v, Face_role::Origin, 0}); break;
    case Kind::End: r.push_back(Face_ref{v, Face_role::End, 0}); break;
    case Kind::One_side:
      if (m_index >= 1 && m_index <= n)
        r.push_back(Face_ref{v, Face_role::Side, m_index});
      break;
    case Kind::All_sides: push_sides(); break;
    case Kind::All:
      r.push_back(Face_ref{v, Face_role::Origin, 0});
      r.push_back(Face_ref{v, Face_role::End, 0});
      push_sides();
      break;
    }
    return r;
  }

  using Variables = std::map<std::string, double>;

  class Grammar;

  struct Context {
    LMap* map = nullptr;
    Grammar* grammar = nullptr;
    Volume_id self = null_volume;
    int stage = 0;
    Variables* vars = nullptr;

    double get(const std::string& n) const {
      auto it = vars->find(n);
      if (it == vars->end())
        throw std::runtime_error("lmap: undefined variable '" + n + "'");
      return it->second;
    }
    void set(const std::string& n, double v) { (*vars)[n] = v; }

    // self -> self[label]_sel
    std::vector<Volume_id> add(const Face_selector& sel, const std::string& label, const Attribute_overload& overload = {}) const;

    // self_sa | other_sb
    std::size_t adjacency(const Face_selector& sa, const std::string& other_label, const Face_selector& sb) const;
  };

  struct Rule {
    std::string predecessor;
    std::function<void(Context&)> block1;
    std::function<bool(const Context&)> cond;
    std::function<void(Context&)> block2;
    std::function<void(Context&)> successor;
    std::string name;
  };

  struct Volume_prototype {
    Volume_attributes attr;
    CGAL::IO::Color color{200, 200, 200};
  };

  class Grammar {
  public:
    void define(const std::string& name, double value) { m_vars[name] = value; }
    double value(const std::string& name) const {
      auto it = m_vars.find(name);
      if (it == m_vars.end())
        throw std::runtime_error("lmap: undefined variable '" + name + "'");
      return it->second;
    }
    Variables& variables() { return m_vars; }

    void define_volume(const std::string& label, const Volume_attributes& a,
                       const CGAL::IO::Color& c = CGAL::IO::Color(200, 200, 200)) {
      m_protos[label] = Volume_prototype{a, c};
    }

    const Volume_prototype& prototype(const std::string& label) const {
      auto it = m_protos.find(label);
      if (it == m_protos.end())
        throw std::runtime_error("lmap: no #define for volume label '" + label + "'");
      return it->second;
    }

    void set_axiom(const std::string& label, const Volume_attributes& a, const Frame& f) {
      m_axiom_label = label; m_axiom_attr = a; m_axiom_frame = f; m_has_axiom = true;
    }

    void add_rule(Rule r) { m_rules.push_back(std::move(r)); }

    void derive(LMap& map, int nb_stages) {
      if (!m_has_axiom) throw std::runtime_error("lmap: grammar has no axiom");
      if (map.nb_volumes() == 0) {
        const Volume_prototype& p = prototype(m_axiom_label);
        map.add_axiom(m_axiom_label, m_axiom_attr, m_axiom_frame, p.color, 0);
      }

      for (int stage = 1; stage <= nb_stages; ++stage) {
        const std::size_t vol_before = map.nb_volumes();
        m_gluings_this_stage = 0;

        std::vector<Volume_id> snapshot;
        snapshot.reserve(map.nb_volumes());
        for (const Volume& v : map.volumes()) snapshot.push_back(v.id);

        m_pending.clear();

        for (const Rule& r : m_rules) {
          for (Volume_id v : snapshot) {
            if (map.volume(v).label != r.predecessor) continue;
            Context ctx{&map, this, v, stage, &m_vars};
            if (r.block1) r.block1(ctx);
            if (!r.cond || r.cond(ctx)) {
              if (r.block2) r.block2(ctx);
              if (r.successor) r.successor(ctx);
            }
          }
        }

        LMap::dedupe_candidates(m_pending);
        count_gluings(map.apply_adjacencies(m_pending));
        m_pending.clear();
      }
    }

    void count_gluings(std::size_t n) { m_gluings_this_stage += n; m_total_gluings += n; }
    std::size_t total_gluings() const { return m_total_gluings; }

    void queue_adjacencies(const std::vector<LMap::Adjacency_candidate>& c) {
      m_pending.insert(m_pending.end(), c.begin(), c.end());
    }

  private:
    Variables m_vars;
    std::map<std::string, Volume_prototype> m_protos;
    std::vector<Rule> m_rules;
    std::string m_axiom_label;
    Volume_attributes m_axiom_attr;
    Frame m_axiom_frame;
    bool m_has_axiom = false;
    std::size_t m_gluings_this_stage = 0;
    std::size_t m_total_gluings = 0;
    std::vector<LMap::Adjacency_candidate> m_pending;
  };

  inline std::vector<Volume_id>
  Context::add(const Face_selector& sel, const std::string& label, const Attribute_overload& overload) const {
    const Volume_prototype& p = grammar->prototype(label);
    return map->add(self, sel, label, p.attr, overload, p.color, stage);
  }

  inline std::size_t
  Context::adjacency(const Face_selector& sa, const std::string& other_label, const Face_selector& sb) const {
    std::size_t found = 0;
    for (Volume_id other : map->volumes_with_label(other_label)) {
      if (other == self) continue;
      const auto c = map->match_adjacencies(self, sa, other, sb);
      found += c.size();
      grammar->queue_adjacencies(c);
    }
    return found;
  }
}