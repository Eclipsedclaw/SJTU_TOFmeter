// Merges skims of analysis/muonic_xray_stack_skim.C made in parts into one file, as hadd does (the trees "cand" and
// "info" chained in the order given; the TParameters keep the first part's values, see the skim), with their
// conditions merged too: the lines all parts share once, the others part by part (each part's processes, exposure and
// seed), under a line that names the parts. The parts are left in place.
//   root -l -b -q 'analysis/skim_merge.C("xray_skim.root", "part_0.root part_1.root part_2.root")'

#include <TFileMerger.h>

#include <sstream>

void skim_merge(const char* out, const char* parts)
{
  std::vector<std::string> files;
  std::istringstream list(parts);
  for (std::string f; list >> f;) files.push_back(f);
  if (files.empty()) {
    printf("no parts given\n");
    return;
  }
  std::vector<std::vector<std::string>> lines;   // every part's conditions, line by line
  for (const auto& f : files) {
    TFile file(f.c_str());
    if (file.IsZombie()) {
      printf("cannot open %s: nothing merged\n", f.c_str());
      return;
    }
    std::vector<std::string> l;
    if (auto text = dynamic_cast<TNamed*>(file.Get("conditions"))) {
      std::istringstream is(text->GetTitle());
      for (std::string s; std::getline(is, s);) l.push_back(s);
    } else {
      l.push_back("no conditions recorded");
    }
    lines.push_back(l);
  }
  TFileMerger merger(false);
  merger.OutputFile(out, "RECREATE");
  for (const auto& f : files) merger.AddFile(f.c_str());
  if (!merger.Merge()) {
    printf("merging into %s FAILED\n", out);
    return;
  }

  std::ostringstream c;
  c << "merged by analysis/skim_merge.C from " << files.size() << " parts (";
  for (std::size_t k = 0; k < files.size(); ++k) c << (k ? ", " : "") << "part " << k << " " << files[k];
  c << "); written " << TDatime().AsSQLString() << "\n";
  const std::size_t n = lines[0].size();
  for (std::size_t i = 0; i < n; ++i) {
    bool same = true;
    for (const auto& l : lines) same = same && l.size() == n && l[i] == lines[0][i];
    if (same) c << lines[0][i] << "\n";
    else
      for (std::size_t k = 0; k < lines.size(); ++k) c << "part " << k << ": " << (i < lines[k].size() ? lines[k][i] : "") << "\n";
  }
  TFile merged(out, "UPDATE");
  auto info = dynamic_cast<TTree*>(merged.Get("info"));
  double live = 0., days = 0.;
  long long muons = 0, totalMuons = 0;
  int processes = 0, totalProcesses = 0;
  if (info) {
    info->SetBranchAddress("live_time", &live); info->SetBranchAddress("muons", &muons); info->SetBranchAddress("files", &processes);
    for (Long64_t i = 0; i < info->GetEntries(); ++i) {
      info->GetEntry(i);
      days += live / 86400.; totalMuons += muons; totalProcesses += processes;
    }
  }
  c << "in all: " << totalProcesses << " processes, " << totalMuons << " muons, live time " << days << " days\n";
  TNamed("conditions", c.str().c_str()).Write("conditions", TObject::kOverwrite);
  printf("%s", c.str().c_str());
}
