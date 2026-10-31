// Drawings to scale of a setup, for slides, from the table of placed volumes that /tof/det/print writes into a log
// (macro/vis_czt.mac's log; scripts/geometry_czt.sh runs both). One drawing per file, PDF and PNG, in outDir:
//   geometry_side        a section through the whole setup at y = yCut, x to the right and z up, labelled
//   geometry_side_zoom   the same section around the target, with its dimensions
//   geometry_top         the setup seen from above, without the lead layer, the steel plate and the Al plates
// A volume is drawn as its extent in the table (axis-aligned: exact for the boxes of the setup); the corner beams,
// L profiles with 1 cm legs (DetectorConstruction.cc), are L shapes in the view from above and miss the section.
// The PNG files are the PDF ones rasterised by pdftoppm (smooth lines, 300 dpi), ROOT's own if pdftoppm is missing.
// Texts are plain Helvetica with the signs x and mu as the WinAnsi bytes \327 and \265: TLatex's #times and #mu, set
// in a symbol font, come out with uneven spaces in the PDF.
//   root -l -b -q 'analysis/geometry_schematic.C("output/geometry_czt/geometry.log", "output/geometry_czt")'
// yCut (cm): through a row of CZT pixels and the middle bar of the top and bottom TOF by default. Colours as in
// macro/vis_czt.mac: lead grey, steel light grey, Al plates light blue, plastic green, target orange, CZT purple
// (Si red and YSO cyan for the Si/YSO stack of the default setup).

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>

namespace {
// a placed volume of the table: centre and full size (cm)
struct Volume {
  std::string name, material;
  int copy = 0;
  double c[3] = {0., 0., 0.}, s[3] = {0., 0., 0.};
  double Lo(int a) const { return c[a] - 0.5 * s[a]; }
  double Hi(int a) const { return c[a] + 0.5 * s[a]; }
};

enum Part { kSkip, kLead, kWall, kSteel, kBeam, kPlate, kBar, kSi, kYSO, kTarget, kCZT, kOther };

Part PartOf(const Volume& v)
{
  const auto starts = [&](const char* prefix) { return v.name.rfind(prefix, 0) == 0; };
  if (starts("WorldPhys") || starts("StationPhys") || starts("CameraEnvelope")) return kSkip;   // air envelopes
  if (starts("AbsorberPhys")) return kLead;
  if (starts("LeadWall")) return kWall;
  if (starts("FrameworkFe")) return kSteel;
  if (starts("SupportBeam")) return kBeam;
  if (starts("FrameworkAl")) return kPlate;
  if (starts("Scintillator")) return kBar;
  if (starts("SiliconLayer")) return kSi;
  if (starts("YSOLayer")) return kYSO;
  if (v.name == "Target") return kTarget;
  if (starts("CZTPixel")) return kCZT;
  return kOther;
}

int ColourOf(Part p)
{
  switch (p) {
    case kLead: case kWall: return TColor::GetColor("#5a5f66");
    case kSteel: case kBeam: return TColor::GetColor("#b5bac2");
    case kPlate: return TColor::GetColor("#a8cce8");
    case kBar: return TColor::GetColor("#4db050");
    case kSi: return TColor::GetColor("#d9534f");
    case kYSO: return TColor::GetColor("#17a2b8");
    case kTarget: return TColor::GetColor("#f29c12");
    case kCZT: return TColor::GetColor("#7a40a1");
    default: return TColor::GetColor("#9a9a9a");
  }
}

std::string Short(const std::string& material)
{
  if (material == "G4_STAINLESS-STEEL") return "steel";
  if (material == "G4_PLASTIC_SC_VINYLTOLUENE") return "plastic";
  return material.rfind("G4_", 0) == 0 ? material.substr(3) : material;
}

// the last table of placed volumes in the log
std::vector<Volume> ReadTable(const char* path)
{
  std::vector<Volume> table;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("=== Placed volumes", 0) != 0) continue;
    table.clear();
    std::getline(in, line);   // column titles
    while (std::getline(in, line)) {
      std::istringstream is(line);
      std::string nameCopy;
      Volume v;
      if (!(is >> nameCopy >> v.material >> v.c[0] >> v.c[1] >> v.c[2] >> v.s[0] >> v.s[1] >> v.s[2])) break;
      const auto colon = nameCopy.rfind(':');
      if (colon == std::string::npos) break;
      v.name = nameCopy.substr(0, colon);
      v.copy = std::atoi(nameCopy.c_str() + colon + 1);
      table.push_back(v);
    }
  }
  return table;
}

// a canvas of `width` pixels showing x1-x2 by y1-y2 (cm), one scale for both
TCanvas* ScaledCanvas(const char* name, double x1, double y1, double x2, double y2, int width)
{
  const int height = std::lround(width * (y2 - y1) / (x2 - x1));
  auto c = new TCanvas(name, name, width, height);
  c->SetCanvasSize(width, height);
  c->SetMargin(0., 0., 0., 0.);
  c->SetFillColor(kWhite);
  c->Range(x1, y1, x2, y2);
  return c;
}

// a filled box with a darker outline; the part outside the view [vx1, vx2] x [vy1, vy2] is cut off, and so are the
// outline's edges there
void Box(double x1, double y1, double x2, double y2, int colour, double alpha = 1., int lineStyle = 1,
         double vx1 = -1e9, double vy1 = -1e9, double vx2 = 1e9, double vy2 = 1e9)
{
  const double a1 = std::max(x1, vx1), b1 = std::max(y1, vy1), a2 = std::min(x2, vx2), b2 = std::min(y2, vy2);
  if (a1 >= a2 || b1 >= b2) return;
  if (alpha > 0.) {
    auto fill = new TBox(a1, b1, a2, b2);
    if (alpha < 1.) fill->SetFillColorAlpha(colour, alpha);
    else fill->SetFillColor(colour);
    fill->SetLineWidth(0);
    fill->Draw();
  }
  const int dark = TColor::GetColorDark(colour);
  auto edge = [&](double xa, double ya, double xb, double yb) {
    auto l = new TLine(xa, ya, xb, yb);
    l->SetLineColor(dark);
    l->SetLineStyle(lineStyle);
    l->SetLineWidth(alpha > 0. ? 1 : 2);
    l->Draw();
  };
  if (y1 >= vy1) edge(a1, y1, a2, y1);
  if (y2 <= vy2) edge(a1, y2, a2, y2);
  if (x1 >= vx1) edge(x1, b1, x1, b2);
  if (x2 <= vx2) edge(x2, b1, x2, b2);
}

TLatex* Text(double x, double y, const char* s, int size, int align = 12, int colour = kBlack, int font = 43,
             double angle = 0.)
{
  auto t = new TLatex(x, y, s);
  t->SetTextFont(font);
  t->SetTextSize(size);
  t->SetTextAlign(align);
  t->SetTextColor(colour);
  t->SetTextAngle(angle);
  t->Draw();
  return t;
}

// a label in a column beside the drawing: a bold title over a grey detail line, a leader to its anchor (ax, ay)
struct Label {
  double ax, ay;
  std::string title, detail;
  double y = 0.;
};

// puts the labels' y near their anchors', at least `gap` apart and within ymin-ymax: neighbours that would overlap
// form a cluster, stacked around the mean of its anchors
void Spread(std::vector<Label>& labels, double gap, double ymin, double ymax)
{
  std::sort(labels.begin(), labels.end(), [](const Label& a, const Label& b) { return a.ay > b.ay; });
  std::vector<std::pair<int, int>> clusters;   // first label, number of labels
  for (int i = 0; i < int(labels.size()); ++i) clusters.push_back({i, 1});
  auto top = [&](const std::pair<int, int>& cl) {
    double mean = 0.;
    for (int k = 0; k < cl.second; ++k) mean += labels[cl.first + k].ay / cl.second;
    const double half = 0.5 * (cl.second - 1) * gap;
    return std::min(std::max(mean, ymin + half), ymax - half) + half;
  };
  for (bool merged = true; merged;) {
    merged = false;
    for (std::size_t k = 0; k + 1 < clusters.size() && !merged; ++k) {
      const double bottom = top(clusters[k]) - (clusters[k].second - 1) * gap;
      if (bottom - top(clusters[k + 1]) < gap) {
        clusters[k].second += clusters[k + 1].second;
        clusters.erase(clusters.begin() + k + 1);
        merged = true;
      }
    }
  }
  for (const auto& cl : clusters)
    for (int k = 0; k < cl.second; ++k) labels[cl.first + k].y = top(cl) - k * gap;
}

// the labels at x = xText; cm: centimetres per pixel of the canvas
void DrawLabels(const std::vector<Label>& labels, double xText, double cm)
{
  const int grey = TColor::GetColor("#52514e");
  for (const auto& l : labels) {
    const double yTitle = l.y + 15. * cm, yDetail = l.y - 15. * cm;   // the pair of lines centred on l.y
    auto leader = new TLine(l.ax, l.ay, xText - 30. * cm, yTitle);
    auto elbow = new TLine(xText - 30. * cm, yTitle, xText - 8. * cm, yTitle);
    for (auto line : {leader, elbow}) {
      line->SetLineColor(grey);
      line->SetLineWidth(1);
      line->Draw();
    }
    auto dot = new TMarker(l.ax, l.ay, 20);
    dot->SetMarkerColor(grey);
    dot->SetMarkerSize(0.7);
    dot->Draw();
    Text(xText, yTitle, l.title.c_str(), 30, 12, kBlack, 63);
    if (!l.detail.empty()) Text(xText, yDetail, l.detail.c_str(), 26, 12, grey);
  }
}

// a dimension: arrow between (x1, y1) and (x2, y2) with its text; along a line the text sits on a white patch in the
// middle, otherwise beside it (side: alignment of the text, 12 to the right of a vertical arrow, 32 to its left)
void Dimension(double x1, double y1, double x2, double y2, const char* text, double cm, int side = 12)
{
  auto arrow = new TArrow(x1, y1, x2, y2, 0.008, "<|>");
  arrow->SetLineWidth(1);
  arrow->SetAngle(40.);
  arrow->SetLineColor(kBlack);
  arrow->SetFillColor(kBlack);
  arrow->Draw();
  const double xm = 0.5 * (x1 + x2), ym = 0.5 * (y1 + y2);
  if (std::abs(y2 - y1) < 1e-9) {   // horizontal: text on the line
    const double half = 0.5 * (std::strlen(text) * 13. + 14.) * cm;
    auto patch = new TBox(xm - half, ym - 14. * cm, xm + half, ym + 14. * cm);
    patch->SetFillColor(kWhite);
    patch->SetLineWidth(0);
    patch->Draw();
    Text(xm, ym, text, 26, 22);
  } else {
    Text(xm + (side == 32 ? -10. : 10.) * cm, ym, text, 26, side);
  }
}

// thin line from (x1, y1) to (x2, y2): extension lines of the dimensions
void Thin(double x1, double y1, double x2, double y2)
{
  auto l = new TLine(x1, y1, x2, y2);
  l->SetLineColor(kGray + 2);
  l->SetLineWidth(1);
  l->Draw();
}

TGaxis* Axis(double x1, double y1, double x2, double y2, double w1, double w2, const char* title, const char* opt,
             double offset)
{
  auto axis = new TGaxis(x1, y1, x2, y2, w1, w2, 510, opt);
  axis->SetLabelFont(43);
  axis->SetLabelSize(24);
  axis->SetTitleFont(43);
  axis->SetTitleSize(26);
  axis->SetTitle(title);
  axis->SetTitleOffset(offset);
  axis->CenterTitle();
  axis->SetTickSize(0.015);
  axis->SetLineColor(kGray + 2);
  axis->Draw();
  return axis;
}

void Save(TCanvas* c, const std::string& base)
{
  c->SaveAs((base + ".pdf").c_str());
  char* tool = gSystem->Which(gSystem->Getenv("PATH"), "pdftoppm");
  const bool rastered =
      tool && gSystem->Exec(Form("pdftoppm -png -singlefile -cropbox -r 300 %s.pdf %s", base.c_str(), base.c_str())) == 0;
  delete[] tool;
  if (rastered) printf("%s.png made from the PDF\n", base.c_str());
  else c->SaveAs((base + ".png").c_str());
}

// the stations of plastic bars: their z (cm), from the top down
std::vector<double> StationZ(const std::vector<Volume>& table)
{
  std::vector<double> z;
  for (const auto& v : table)
    if (PartOf(v) == kBar && std::none_of(z.begin(), z.end(), [&](double s) { return std::abs(s - v.c[2]) < 0.01; }))
      z.push_back(v.c[2]);
  std::sort(z.rbegin(), z.rend());
  return z;
}

std::string StationName(const std::vector<double>& stations, double z)
{
  const char* names[3] = {"Top TOF", "Middle TOF", "Bottom TOF"};
  for (std::size_t k = 0; k < stations.size(); ++k)
    if (std::abs(stations[k] - z) < 0.01) return stations.size() == 3 ? names[k] : Form("TOF station %zu", k);
  return "TOF";
}

// bars of the station at z: number, and whether they lie along x
std::pair<int, bool> StationBars(const std::vector<Volume>& table, double z)
{
  int n = 0;
  bool alongX = true;
  for (const auto& v : table)
    if (PartOf(v) == kBar && std::abs(v.c[2] - z) < 0.01) {
      ++n;
      alongX = v.s[0] > v.s[1];
    }
  return {n, alongX};
}

// pixels per side and full size of the CZT plane
struct Plane { int n = 0; double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9}, pixel = 0.; };
Plane CZTPlane(const std::vector<Volume>& table)
{
  Plane p;
  int pixels = 0;
  for (const auto& v : table)
    if (PartOf(v) == kCZT) {
      ++pixels;
      p.pixel = v.s[0];
      for (int a = 0; a < 3; ++a) { p.lo[a] = std::min(p.lo[a], v.Lo(a)); p.hi[a] = std::max(p.hi[a], v.Hi(a)); }
    }
  p.n = std::lround(std::sqrt(pixels));
  return p;
}

const Volume* Find(const std::vector<Volume>& table, Part part)
{
  for (const auto& v : table)
    if (PartOf(v) == part) return &v;
  return nullptr;
}

// ---------------------------------------------------------------- section through the whole setup
void DrawSide(const std::vector<Volume>& table, double yCut, const std::string& base)
{
  std::vector<const Volume*> cut;   // the volumes the plane y = yCut goes through
  double x0 = 1e9, x1 = -1e9, z0 = 1e9, z1 = -1e9;
  for (const auto& v : table) {
    if (PartOf(v) == kSkip || v.Lo(1) >= yCut || v.Hi(1) <= yCut) continue;
    cut.push_back(&v);
    x0 = std::min(x0, v.Lo(0)); x1 = std::max(x1, v.Hi(0));
    z0 = std::min(z0, v.Lo(2)); z1 = std::max(z1, v.Hi(2));
  }
  const int width = 1800;
  const double left = x0 - 15., right = x1 + 78., bottom = z0 - 16., topEdge = z1 + 22.;
  const double cm = (right - left) / width, xText = x1 + 9.;
  auto c = ScaledCanvas("side", left, bottom, right, topEdge, width);
  const int grey = TColor::GetColor("#52514e");

  // the ground
  auto ground = new TBox(x0 - 4., z0 - 2.5, x1 + 4., z0);
  ground->SetFillColor(kGray + 1);
  ground->SetFillStyle(3354);
  ground->Draw();
  auto groundLine = new TLine(x0 - 4., z0, x1 + 4., z0);
  groundLine->SetLineColor(kGray + 2);
  groundLine->Draw();

  for (auto v : cut) Box(v->Lo(0), v->Lo(2), v->Hi(0), v->Hi(2), ColourOf(PartOf(*v)));

  // cosmic muons
  for (double dx : {-6., 0., 6.}) {
    auto mu = new TArrow(-28. + dx + 4., z1 + 18., -28. + dx, z1 + 3., 0.012, "|>");
    mu->SetLineWidth(2);
    mu->SetLineColor(grey);
    mu->SetFillColor(grey);
    mu->Draw();
  }
  Text(-18., z1 + 11., "Cosmic \265\261", 30, 12, grey, 63);

  // labels in the lead
  if (auto lead = Find(table, kLead))
    Text(lead->c[0], lead->c[2], Form("Pb  %g \327 %g \327 %g cm", lead->s[0], lead->s[1], lead->s[2]), 30, 22,
         kWhite, 63);
  const Volume* wall = nullptr;   // the wall on the left, which no leader crosses
  for (auto v : cut)
    if (PartOf(*v) == kWall && (!wall || v->c[0] < wall->c[0])) wall = v;
  if (wall)
    Text(wall->c[0], 0.5 * (wall->Lo(2) + std::min(wall->Hi(2), z1 - 12.)),
         Form("Pb wall, %g cm", std::min(wall->s[0], wall->s[1])), 30, 22, kWhite, 63, 90.);

  // labels beside the setup, each pointing at the right end of its part in the section
  std::vector<Label> labels;
  auto rightEnd = [&](auto select, double& ax, double& ay) {
    ax = -1e9; ay = 0.;
    int n = 0;
    for (auto v : cut)
      if (select(*v)) { ax = std::max(ax, v->Hi(0)); ay += v->c[2]; ++n; }
    if (n) ay /= n;
    return n > 0;
  };
  double ax, ay;
  if (auto steel = Find(table, kSteel); steel && rightEnd([](const Volume& v) { return PartOf(v) == kSteel; }, ax, ay))
    labels.push_back({ax, ay, "Steel plate", Form("%g \327 %g \327 %g cm", steel->s[0], steel->s[1], steel->s[2])});
  const auto stations = StationZ(table);
  for (double z : stations)
    if (rightEnd([&](const Volume& v) { return PartOf(v) == kBar && std::abs(v.c[2] - z) < 0.01; }, ax, ay)) {
      const auto bars = StationBars(table, z);
      labels.push_back({ax, ay, StationName(stations, z), Form("%d plastic bars along %s", bars.first,
                                                               bars.second ? "x" : "y")});
    }
  if (auto target = Find(table, kTarget); target && rightEnd([](const Volume& v) { return PartOf(v) == kTarget; }, ax, ay))
    labels.push_back({ax, ay, Short(target->material) + " target",
                      Form("%g \327 %g \327 %g cm", target->s[0], target->s[1], target->s[2])});
  const Plane czt = CZTPlane(table);
  if (czt.n && rightEnd([](const Volume& v) { return PartOf(v) == kCZT; }, ax, ay))
    labels.push_back({ax, ay, "CZT plane", Form("%d \327 %d pixels, %g \327 %g \327 %g cm", czt.n, czt.n,
                                                czt.hi[0] - czt.lo[0], czt.hi[1] - czt.lo[1], czt.hi[2] - czt.lo[2])});
  if (rightEnd([](const Volume& v) { return PartOf(v) == kSi || PartOf(v) == kYSO; }, ax, ay))
    labels.push_back({ax, ay, "Si/YSO stack", "Compton detector"});
  for (auto v : cut)
    if (PartOf(*v) == kOther) labels.push_back({v->Hi(0), v->c[2], v->name, Short(v->material)});
  // one label for the Al plates, at the plate farthest from the other labels
  const Volume* plate = nullptr;
  double clearance = -1.;
  for (auto v : cut) {
    if (PartOf(*v) != kPlate) continue;
    double d = 1e9;
    for (const auto& l : labels) d = std::min(d, std::abs(l.ay - v->c[2]));
    if (d > clearance) { clearance = d; plate = v; }
  }
  if (plate)
    labels.push_back({plate->Hi(0), plate->c[2], "Al plates", Form("%g \327 %g \327 %g cm", std::max(plate->s[0], plate->s[1]),
                                                                  std::min(plate->s[0], plate->s[1]), plate->s[2])});
  Spread(labels, 78. * cm, z0 + 2., z1 - 2.);
  DrawLabels(labels, xText, cm);

  // scales
  Axis(x0 - 5., z0, x0 - 5., 10. * std::floor(z1 / 10.), z0, 10. * std::floor(z1 / 10.), "z (cm)", "-", 1.6);
  const double xs = 10. * std::floor(std::min(-x0, x1) / 10.);
  Axis(-xs, z0 - 6.5, xs, z0 - 6.5, -xs, xs, "x (cm)", "", 1.2);
  std::string note = Form("Section at y = %g cm, to scale.", yCut);
  if (auto tof = Find(table, kBar))
    note += Form(" TOF bars %g \327 %g \327 %g cm.", std::max(tof->s[0], tof->s[1]), std::min(tof->s[0], tof->s[1]),
                 tof->s[2]);
  Text(xText - 30. * cm, z0 - 11., note.c_str(), 24, 12, grey);
  Save(c, base);
}

// ---------------------------------------------------------------- the same section around the target
void DrawZoom(const std::vector<Volume>& table, double yCut, const std::string& base)
{
  const Volume* target = Find(table, kTarget);
  if (!target) {
    printf("no target in the table: no zoomed section\n");
    return;
  }
  const double vx1 = target->Lo(0) - 5.5, vx2 = target->Hi(0) + 4.5, vz1 = target->Lo(2) - 4., vz2 = target->Hi(2) + 5.;
  const int width = 2000;
  const double cm = (vx2 - vx1) / width;
  auto c = ScaledCanvas("zoom", vx1, vz1, vx2, vz2, width);
  const int grey = TColor::GetColor("#52514e");

  // the section cut off 2.8 cm from the left edge, clear of the z scale
  std::vector<const Volume*> above;   // in the section, over the target's width and above it, from the bottom up
  for (const auto& v : table) {
    if (PartOf(v) == kSkip || v.Lo(1) >= yCut || v.Hi(1) <= yCut) continue;
    Box(v.Lo(0), v.Lo(2), v.Hi(0), v.Hi(2), ColourOf(PartOf(v)), 1., 1, vx1 + 2.8, vz1, vx2, vz2);
    if (v.Lo(2) >= target->Hi(2) - 1e-6 && v.Lo(2) < vz2 && v.Hi(0) > target->Lo(0) && v.Lo(0) < target->Hi(0))
      above.push_back(&v);
  }
  std::sort(above.begin(), above.end(), [](const Volume* a, const Volume* b) { return a->Lo(2) < b->Lo(2); });

  // the target: name inside, thickness and gap to the part above on its right, width below
  Text(target->c[0], target->c[2], Form("%s target  %g \327 %g \327 %g cm", Short(target->material).c_str(),
                                         target->s[0], target->s[1], target->s[2]), 32, 22, kBlack, 63);
  const double xr = target->Hi(0) + 1.1;
  Thin(target->Hi(0), target->Lo(2), xr + 0.4, target->Lo(2));
  Thin(target->Hi(0), target->Hi(2), xr + 0.4, target->Hi(2));
  Dimension(xr, target->Lo(2), xr, target->Hi(2), Form("%g cm", target->s[2]), cm);
  if (!above.empty()) {
    const Volume* next = above[0];
    Thin(target->Hi(0), next->Lo(2), xr + 0.4, next->Lo(2));
    Dimension(xr, target->Hi(2), xr, next->Lo(2), Form("%g cm", next->Lo(2) - target->Hi(2)), cm);
    if (PartOf(*next) == kPlate)
      Text(0.5 * (target->Lo(0) + target->Hi(0)), next->Hi(2) + 0.55, Form("Al plate, %g cm", next->s[2]), 28, 22, grey,
           63);
    // the gap between that part and the next one up, on the left
    for (std::size_t k = 1; k < above.size(); ++k)
      if (above[k]->Lo(2) > next->Hi(2) + 1e-6) {
        const double xl = target->Lo(0) - 1.1;
        Thin(xl - 0.4, next->Hi(2), target->Lo(0) + 2., next->Hi(2));
        Thin(xl - 0.4, above[k]->Lo(2), target->Lo(0) + 2., above[k]->Lo(2));
        Dimension(xl, next->Hi(2), xl, above[k]->Lo(2), Form("%g cm", above[k]->Lo(2) - next->Hi(2)), cm, 32);
        if (PartOf(*above[k]) == kBar) {
          const auto bars = StationBars(table, above[k]->c[2]);
          Text(0.5 * (target->Lo(0) + target->Hi(0)), above[k]->Hi(2) + 0.5,
               Form("%s: %d plastic bars along %s, %g \327 %g cm in this section",
                    StationName(StationZ(table), above[k]->c[2]).c_str(), bars.first, bars.second ? "x" : "y",
                    above[k]->s[0], above[k]->s[2]),
               28, 22, grey, 63);
        }
        break;
      }
  }

  // the CZT plane under it: name inside, thickness and pixels on the right, width below
  const Plane czt = CZTPlane(table);
  double zw = target->Lo(2) - 1.;   // height of the first width dimension
  if (czt.n) {
    const double czx = 0.5 * (czt.lo[0] + czt.hi[0]), czz = 0.5 * (czt.lo[2] + czt.hi[2]), dz = czt.hi[2] - czt.lo[2];
    Text(czx, czz, "CZT", 30, 22, kWhite, 63);
    const double xc = czt.hi[0] + 0.9;
    Thin(czt.hi[0], czt.lo[2], xc + 0.4, czt.lo[2]);
    Dimension(xc, czt.lo[2], xc, czt.hi[2], Form("%g cm", dz), cm);
    // its pixels right of the target's edge, clear of the extension lines of the widths
    Text(target->Hi(0) + 0.35, czz + 0.25, Form("CZT: %d \327 %d pixels", czt.n, czt.n), 26, 12, grey);
    Text(target->Hi(0) + 0.35, czz - 0.25, Form("of %g \327 %g \327 %g cm", czt.pixel, czt.pixel, dz), 26, 12, grey);
    zw = czt.lo[2] - 0.9;
    Thin(czt.lo[0], czt.lo[2], czt.lo[0], zw - 0.3);
    Thin(czt.hi[0], czt.lo[2], czt.hi[0], zw - 0.3);
    Dimension(czt.lo[0], zw, czt.hi[0], zw, Form("%g cm", czt.hi[0] - czt.lo[0]), cm);
    zw -= 1.1;
  }
  Thin(target->Lo(0), target->Lo(2), target->Lo(0), zw - 0.3);
  Thin(target->Hi(0), target->Lo(2), target->Hi(0), zw - 0.3);
  Dimension(target->Lo(0), zw, target->Hi(0), zw, Form("%g cm", target->s[0]), cm);

  // z scale on the left
  const double za = std::ceil(vz1 + 0.3), zb = std::floor(vz2 - 0.3);
  Axis(vx1 + 2., za, vx1 + 2., zb, za, zb, "z (cm)", "-", 1.3);
  Text(vx2 - 0.3, vz1 + 0.35, Form("Section at y = %g cm, to scale", yCut), 24, 32, grey);
  Save(c, base);
}

// ---------------------------------------------------------------- seen from above
void DrawTop(const std::vector<Volume>& table, const std::string& base)
{
  double x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
  for (const auto& v : table) {
    const Part p = PartOf(v);
    if (p == kSkip || p == kLead || p == kSteel || p == kPlate) continue;
    x0 = std::min(x0, v.Lo(0)); x1 = std::max(x1, v.Hi(0));
    y0 = std::min(y0, v.Lo(1)); y1 = std::max(y1, v.Hi(1));
  }
  const int width = 1900;
  const double left = x0 - 15., right = x1 + 80., bottom = y0 - 15., topEdge = y1 + 5.;
  const double cm = (right - left) / width, xText = x1 + 11.;
  auto c = ScaledCanvas("top", left, bottom, right, topEdge, width);
  const int grey = TColor::GetColor("#52514e");
  std::vector<Label> labels;

  const Volume* wall = nullptr;   // the wall on the right
  for (const auto& v : table)
    if (PartOf(v) == kWall) {
      Box(v.Lo(0), v.Lo(1), v.Hi(0), v.Hi(1), ColourOf(kWall));
      if (!wall || v.c[0] > wall->c[0]) wall = &v;
    }
  if (wall)
    labels.push_back({wall->Hi(0), wall->c[1] + 0.3 * wall->s[1], "Pb walls",
                      Form("%g cm thick, %g cm high", std::min(wall->s[0], wall->s[1]), wall->s[2])});

  // L beams: legs 1 cm thick along the outer edges of their extent
  const Volume* beam = nullptr;   // the one at +x, +y
  for (const auto& v : table) {
    if (PartOf(v) != kBeam) continue;
    const double cx = v.c[0] > 0. ? v.Hi(0) : v.Lo(0), cy = v.c[1] > 0. ? v.Hi(1) : v.Lo(1);
    const double sx = v.c[0] > 0. ? 1. : -1., sy = v.c[1] > 0. ? 1. : -1., leg = v.s[0], t = 1.;
    double px[7] = {cx, cx, cx - sx * t, cx - sx * t, cx - sx * leg, cx - sx * leg, cx};
    double py[7] = {cy, cy - sy * leg, cy - sy * leg, cy - sy * t, cy - sy * t, cy, cy};
    auto fill = new TPolyLine(7, px, py);
    fill->SetFillColor(ColourOf(kBeam));
    fill->Draw("f");
    auto edge = new TPolyLine(7, px, py);
    edge->SetLineColor(TColor::GetColorDark(ColourOf(kBeam)));
    edge->Draw();
    if (v.c[0] > 0. && v.c[1] > 0.) beam = &v;
  }
  if (beam) labels.push_back({beam->Hi(0), beam->Hi(1) - 0.5, "Steel L beams", Form("%g cm legs, 1 cm thick", beam->s[0])});

  // the target, the Si/YSO stack or the CZT plane under it (drawn over it)
  const Volume* target = Find(table, kTarget);
  if (target) {
    Box(target->Lo(0), target->Lo(1), target->Hi(0), target->Hi(1), ColourOf(kTarget));
    labels.push_back({target->Hi(0), target->Lo(1) + 0.3 * target->s[1], Short(target->material) + " target",
                      Form("%g \327 %g cm, %g cm thick", target->s[0], target->s[1], target->s[2])});
  }
  const Plane czt = CZTPlane(table);
  for (const auto& v : table)
    if (PartOf(v) == kCZT) Box(v.Lo(0), v.Lo(1), v.Hi(0), v.Hi(1), ColourOf(kCZT));
  if (czt.n)
    labels.push_back({czt.hi[0], czt.lo[1] + 0.6 * (czt.hi[1] - czt.lo[1]),
                      target ? "CZT plane, under the target" : "CZT plane",
                      Form("%g \327 %g cm, %d \327 %d pixels", czt.hi[0] - czt.lo[0], czt.hi[1] - czt.lo[1], czt.n,
                           czt.n)});
  if (auto si = Find(table, kSi); si && !czt.n)
    labels.push_back({si->Hi(0), si->c[1], target ? "Si/YSO stack, under the target" : "Si/YSO stack",
                      Form("%g \327 %g cm", si->s[0], si->s[1])});

  // the bars as outlines: along y solid, along x dashed (stations along x lie over each other)
  const auto stations = StationZ(table);
  std::string alongX, alongY;
  double axX = -1e9, ayX = 0., axY = -1e9, ayY = 0.;
  for (double z : stations) {
    const auto bars = StationBars(table, z);
    std::string& names = bars.second ? alongX : alongY;
    std::string name = StationName(stations, z);
    if (!names.empty()) { name[0] = std::tolower(name[0]); names += " and "; }
    names += name;
  }
  double maxY = 0.;   // centre of the outermost bar along x
  for (const auto& v : table)
    if (PartOf(v) == kBar && v.s[0] > v.s[1]) maxY = std::max(maxY, v.c[1]);
  for (const auto& v : table) {
    if (PartOf(v) != kBar) continue;
    const bool x = v.s[0] > v.s[1];
    Box(v.Lo(0), v.Lo(1), v.Hi(0), v.Hi(1), ColourOf(kBar), 0., x ? 2 : 1);
    // anchors: the right end of the bar along x nearest to 0.45 of maxY, the right edge of the outermost bar along y
    if (x && (axX < -1e8 || std::abs(v.c[1] - 0.45 * maxY) < std::abs(ayX - 0.45 * maxY))) { axX = v.Hi(0); ayX = v.c[1]; }
    if (!x && v.Hi(0) > axY) { axY = v.Hi(0); ayY = v.Lo(1) + 0.25 * v.s[1]; }
  }
  auto firstBar = Find(table, kBar);
  const std::string size =
      firstBar ? Form(", %g \327 %g cm", std::max(firstBar->s[0], firstBar->s[1]), std::min(firstBar->s[0], firstBar->s[1]))
               : "";
  if (!alongX.empty() && axX > -1e9) labels.push_back({axX, ayX, alongX, "bars along x (dashed)" + size});
  if (!alongY.empty() && axY > -1e9) labels.push_back({axY, ayY, alongY, "bars along y" + size});

  Spread(labels, 78. * cm, y0 + 9., y1 - 2.);   // above the note at the bottom
  DrawLabels(labels, xText, cm);
  const double xs = 10. * std::floor(std::min(-x0, x1) / 10.), ys = 10. * std::floor(std::min(-y0, y1) / 10.);
  Axis(-xs, y0 - 5., xs, y0 - 5., -xs, xs, "x (cm)", "", 1.2);
  Axis(x0 - 5., -ys, x0 - 5., ys, -ys, ys, "y (cm)", "-", 1.6);
  Text(xText - 30. * cm, y0 + 3.5, "Seen from above, to scale; the Pb layer,", 24, 12, grey);
  Text(xText - 30. * cm, y0 + 0.5, "the steel plate and the Al plates not drawn.", 24, 12, grey);
  Save(c, base);
}
}  // namespace

void geometry_schematic(const char* log = "output/geometry_czt/geometry.log", const char* outDir = "output/geometry_czt",
                        double yCut = 1.25)
{
  const auto table = ReadTable(log);
  if (table.empty()) {
    printf("no table of placed volumes (/tof/det/print) in %s\n", log);
    return;
  }
  printf("%zu placed volumes read from %s\n", table.size(), log);
  gROOT->SetBatch(true);
  const std::string out = outDir;
  DrawSide(table, yCut, out + "/geometry_side");
  DrawZoom(table, yCut, out + "/geometry_side_zoom");
  DrawTop(table, out + "/geometry_top");
}
