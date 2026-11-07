// A cut profile (2026-11-05): the analysis cuts of one detector configuration, kept in a text file so that each
// configuration has its own and an analysis chooses one by name: analysis/cuts/<name>.cuts (czt: the CZT plane;
// lyso: the LYSO panel of 20 x 20 pixels, 12.69 x 12.69 cm), or a path to a file. Each line "<variable> <lo> <hi>", "-" for an open side, "#"
// starts a comment; a cut keeps lo <= value <= hi, a variable not in the file is not cut on. The variables, all
// measured (the skims' columns of muonic_xray_stack_skim.C):
//   track   the track at the target's centre plane, max(|xTgt|, |yTgt|) (cm)
//   zenith  the track's zenith angle, atan(hypot(slopeX, slopeY)) (deg)
//   topE    the top-TOF energy, eTopM (MeV)          midE  the middle-TOF energy, eMidM (MeV)
//   tkin    the TOF kinetic energy (MeV)             entryE  the pixel entry's energy, e (keV)
//   otherE  the other pixels' energy within 100 ns, eOther (MeV)
//   npix    the pixels hit in the event within 100 ns (a whole number: "npix - 4" keeps at most 4)
// Used by analysis/muonic_xray_target_select.C, analysis/panel_intrinsic_check.C and the work folder's event_flow.C and
// background_composition.C; analysis/cut_profile_optimize.C writes new ones.

#ifndef CutProfile_h
#define CutProfile_h 1

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// the variables of one pixel entry, from the skim's columns (as muonic_xray_target_select.C computes them)
struct CutVars {
  double track = 0., zenith = 0., topE = 0., midE = 0., tkin = 0., entryE = 0., otherE = 0., npix = 0.;
  static CutVars FromSkim(double xTgt, double yTgt, double slopeX, double slopeY, double eTopM, double eMidM, double tkin,
                          double e, double eOther, int npix)
  {
    CutVars v;
    v.track = std::max(std::abs(xTgt), std::abs(yTgt));
    v.zenith = std::atan(std::hypot(slopeX, slopeY)) * 180. / M_PI;
    v.topE = eTopM; v.midE = eMidM; v.tkin = tkin; v.entryE = e; v.otherE = eOther / 1000.; v.npix = npix;
    return v;
  }
  double Get(const std::string& k) const
  {
    return k == "track" ? track : k == "zenith" ? zenith : k == "topE" ? topE : k == "midE" ? midE : k == "tkin" ? tkin
         : k == "entryE" ? entryE : k == "otherE" ? otherE : npix;
  }
};

struct CutProfile {
  std::string name, path;
  std::vector<std::string> notes;                           // the file's comment lines
  std::map<std::string, std::pair<double, double>> range;   // variable: lo, hi (+-1e9: open)

  static const std::vector<std::string>& Variables()
  {
    static const std::vector<std::string> v = {"track", "zenith", "topE", "midE", "tkin", "entryE", "otherE", "npix"};
    return v;
  }
  // a name (analysis/cuts/<name>.cuts from the working directory, the build or work folder) or a path
  bool Load(const std::string& nameOrPath)
  {
    path = (nameOrPath.find('/') != std::string::npos || nameOrPath.find(".cuts") != std::string::npos)
             ? nameOrPath : "analysis/cuts/" + nameOrPath + ".cuts";
    const auto slash = path.find_last_of('/');
    name = path.substr(slash == std::string::npos ? 0 : slash + 1);
    if (name.size() > 5 && name.compare(name.size() - 5, 5, ".cuts") == 0) name.resize(name.size() - 5);
    std::ifstream in(path);
    if (!in) { printf("CutProfile: cannot read %s\n", path.c_str()); return false; }
    range.clear(); notes.clear();
    std::string line;
    while (std::getline(in, line)) {
      const auto hash = line.find('#');
      if (hash != std::string::npos) {
        std::string note = line.substr(hash + 1);
        if (!note.empty() && note[0] == ' ') note.erase(0, 1);
        if (hash == 0) notes.push_back(note);
        line.resize(hash);
      }
      std::istringstream words(line);
      std::string key, lo, hi;
      if (!(words >> key)) continue;
      if (!(words >> lo >> hi)) { printf("CutProfile: %s: \"%s\" needs <lo> <hi>\n", path.c_str(), key.c_str()); return false; }
      bool known = false;
      for (const auto& v : Variables()) known = known || v == key;
      if (!known) { printf("CutProfile: %s: unknown variable %s\n", path.c_str(), key.c_str()); return false; }
      range[key] = {lo == "-" ? -1e9 : std::stod(lo), hi == "-" ? 1e9 : std::stod(hi)};
    }
    return true;
  }
  bool Has(const std::string& k) const { return range.count(k) > 0; }
  double Lo(const std::string& k) const { auto it = range.find(k); return it == range.end() ? -1e9 : it->second.first; }
  double Hi(const std::string& k) const { auto it = range.find(k); return it == range.end() ? 1e9 : it->second.second; }
  bool Keeps(const std::string& k, double v) const { return v >= Lo(k) && v <= Hi(k); }
  // every cut at once
  bool Pass(const CutVars& v) const
  {
    for (const auto& r : range)
      if (!(v.Get(r.first) >= r.second.first && v.Get(r.first) <= r.second.second)) return false;
    return true;
  }
  bool Pass(double xTgt, double yTgt, double slopeX, double slopeY, double eTopM, double eMidM, double tkin, double e,
            double eOther, int npix) const
  {
    return Pass(CutVars::FromSkim(xTgt, yTgt, slopeX, slopeY, eTopM, eMidM, tkin, e, eOther, npix));
  }
  // the variables cut on, in the order of Variables()
  std::vector<std::string> Keys() const
  {
    std::vector<std::string> k;
    for (const auto& v : Variables()) if (Has(v)) k.push_back(v);
    return k;
  }
  // a cut's range in words: "1.6-3.5 MeV", "> 1.6 MeV", "< 8 cm"
  std::string Range(const std::string& k, const char* unit) const
  {
    const double lo = Lo(k), hi = Hi(k);
    if (lo > -1e8 && hi < 1e8) return Text("%g-%g %s", lo, hi, unit);
    if (lo > -1e8) return Text("> %g %s", lo, unit);
    return Text("< %g %s", hi, unit);
  }
  // a cut in words, long ("track at the target within 8 cm", as event_flow.C's table) or short ("track < 8 cm", as the
  // plots' notes); det names the pixel entry ("CZT": "CZT entry < 600 keV")
  std::string Words(const std::string& k, bool longForm, const std::string& det = "") const
  {
    const std::string entry = det.empty() ? "entry" : det + " entry", tag = det.empty() ? "" : " (" + det + " entry)";
    const bool upper = Hi(k) < 1e8, lower = Lo(k) > -1e8;
    if (k == "track")
      return upper && !lower ? (longForm ? Text("track at the target within %g cm", Hi(k)) : Text("track < %g cm", Hi(k)))
                             : (longForm ? "track at the target " : "track ") + Range(k, "cm");
    if (k == "zenith")
      return upper && !lower ? (longForm ? Text("within %g deg of the vertical", Hi(k)) : Text("zenith < %g deg", Hi(k)))
                             : "zenith " + Range(k, "deg");
    if (k == "topE") return (longForm ? "top-TOF energy " : "top-TOF ") + Range(k, "MeV");
    if (k == "midE") return (longForm ? "middle-TOF energy " : "middle-TOF ") + Range(k, "MeV");
    if (k == "tkin") return (longForm ? "TOF kinetic energy " : "kinetic energy ") + Range(k, "MeV");
    if (k == "entryE") return longForm ? "entry energy " + Range(k, "keV") + tag : entry + " " + Range(k, "keV");
    if (k == "otherE") return "other pixels " + Range(k, "MeV") + (longForm ? tag : "");
    // npix
    const std::string n = upper && !lower ? Text("at most %g pixels hit", Hi(k))
                                          : lower && !upper ? Text("at least %g pixels hit", Lo(k))
                                                            : Text("%g-%g pixels hit", Lo(k), Hi(k));
    return n + (longForm ? tag : "");
  }
  // all the cuts in short words, e.g. "track < 8 cm, zenith < 35 deg, top-TOF 1.6-3.5 MeV, ..."
  std::string Describe(const std::string& det = "") const
  {
    std::string s;
    for (const auto& k : Keys()) s += (s.empty() ? "" : ", ") + Words(k, false, det);
    return s;
  }
  template <typename... A>
  static std::string Text(const char* format, A... a)   // printf into a string
  {
    char text[200];
    snprintf(text, sizeof(text), format, a...);
    return text;
  }
  // the profile as a file, with its comment lines first
  void Write(const std::string& file) const
  {
    std::ofstream out(file);
    for (const auto& n : notes) out << "# " << n << "\n";
    for (const auto& v : Variables()) {
      if (!Has(v)) continue;
      auto side = [](double x) { return std::abs(x) > 1e8 ? std::string("-") : std::to_string(x).substr(0, std::to_string(x).find_last_not_of('0') + 1); };
      std::string lo = side(Lo(v)), hi = side(Hi(v));
      if (!lo.empty() && lo.back() == '.') lo.pop_back();
      if (!hi.empty() && hi.back() == '.') hi.pop_back();
      char line[96];
      snprintf(line, sizeof(line), "%-8s %8s %8s\n", v.c_str(), lo.c_str(), hi.c_str());
      out << line;
    }
  }
};

#endif
