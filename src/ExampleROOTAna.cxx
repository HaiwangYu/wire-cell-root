#include "WireCellRoot/ExampleROOTAna.h"
#include "WireCellIface/ITrace.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"

#include "WireCellIface/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"

#include <string>
#include <vector>

WIRECELL_FACTORY(ExampleROOTAna, WireCell::Root::ExampleROOTAna,
                 WireCell::IFrameSink, WireCell::IConfigurable)

using namespace WireCell;

Root::ExampleROOTAna::ExampleROOTAna() : log(Log::logger("ana")) {}

Root::ExampleROOTAna::~ExampleROOTAna() {}

void Root::ExampleROOTAna::configure(const WireCell::Configuration &cfg) {
  std::string fn;

  fn = cfg["output_filename"].asString();
  if (fn.empty()) {
    THROW(ValueError() << errmsg{
              "Must provide output filename to ExampleROOTAna"});
  }

  auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
  m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

  m_cfg = cfg;

  recreate_out_file();
}

WireCell::Configuration Root::ExampleROOTAna::default_configuration() const {
  Configuration cfg;

  cfg["anode"] = "AnodePlane";

  // Name of ROOT file to write.
  cfg["output_filename"] = "";

  cfg["nrebin"] = 1;

  return cfg;
}

namespace {
std::vector<WireCell::Binning> collate_byplane(const ITrace::vector &traces,
                                               const IAnodePlane::pointer anode,
                                               ITrace::vector byplane[]) {
  std::vector<int> uvwt[4];
  for (auto trace : traces) {
    const int chid = trace->channel();
    auto wpid = anode->resolve(chid);
    const int iplane = wpid.index();
    if (iplane < 0 || iplane >= 3) {
      THROW(RuntimeError() << errmsg{"Illegal wpid"});
    }
    uvwt[iplane].push_back(chid);
    byplane[iplane].push_back(trace);
    uvwt[3].push_back(trace->tbin());
    uvwt[3].push_back(trace->tbin() + trace->charge().size());
  }

  std::vector<Binning> binnings(4);
  for (int ind = 0; ind < 4; ++ind) {
    auto const &one = uvwt[ind];
    // std::cerr << "[wgu] get ind: " << ind << " size: " << one.size() <<
    // std::endl;
    if (one.empty()) {
      // THROW(ValueError() << errmsg{"ExampleROOTAna: bogus bounds"});
      std::cerr << "[wgu] plane: " << ind << " has not traces. " << std::endl;
    } else {
      auto mme = std::minmax_element(one.begin(), one.end());
      const int vmin = *mme.first;
      const int vmax = *mme.second;
      if (ind == 3) {
        const int n = vmax - vmin;
        // binnings.push_back(Binning(n, vmin, vmax););
        binnings.at(ind) = Binning(n, vmin, vmax);
      } else {
        // Channel-centered binning
        const double diff = vmax - vmin;
        // binnings.push_back(Binning(diff+1, vmin-0.5, vmax+0.5));
        binnings.at(ind) = Binning(diff + 1, vmin - 0.5, vmax + 0.5);
      }
    }
  }
  return binnings;
}
} // namespace

void Root::ExampleROOTAna::recreate_out_file() const {
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = "RECREATE";
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());
  output_tf->Close("R");
  delete output_tf;
  output_tf = nullptr;
}

void Root::ExampleROOTAna::peak_frame(const IFrame::pointer &frame) const {

  std::string frame_tags = "";
  for (auto tag : frame->frame_tags()) {
    frame_tags += " ";
    frame_tags += tag;
  }

  for (auto tag : frame->trace_tags()) {
    log->info(
        "ExampleROOTAna: trace tag: \"{}\" in frame: tags \"{}\" - ident: {}",
        tag, frame_tags, frame->ident());
  }
}

void Root::ExampleROOTAna::fill_hist(const IFrame::pointer &frame,
                                     TFile *output_tf) const {

  for (auto tag : frame->trace_tags()) {
    ITrace::vector traces = FrameTools::tagged_traces(frame, tag);
    if (traces.empty()) {
      log->warn("ExampleROOTAna: no tagged traces for \"{}\"", tag);
      continue;
    }

    log->debug("ExampleROOTAna: tag: \"{}\" with {} traces", tag,
               traces.size());

    ITrace::vector traces_byplane[3];
    auto binnings = collate_byplane(traces, m_anode, traces_byplane);
    Binning tbin = binnings[3];
    for (int iplane = 0; iplane < 3; ++iplane) {
      if (traces_byplane[iplane].empty())
        continue;
      const std::string name = Form("h%c_%s", 'u' + iplane, tag.c_str());
      Binning cbin = binnings[iplane];
      std::stringstream ss;
      ss << "ExampleROOTAna:"
         << " cbin:" << cbin.nbins() << "[" << cbin.min() << "," << cbin.max()
         << "]"
         << " tbin:" << tbin.nbins() << "[" << tbin.min() << "," << tbin.max()
         << "]";
      log->debug(ss.str());

      // consider to add nrebin ...
      int nbins = tbin.nbins();

      TH2F *hist =
          new TH2F(name.c_str(), name.c_str(), cbin.nbins(), cbin.min(),
                   cbin.max(), nbins, tbin.min(), tbin.max());

      hist->SetDirectory(output_tf);

      for (auto trace : traces_byplane[iplane]) {
        const int tbin1 = trace->tbin();
        const int ch = trace->channel();
        auto const &charges = trace->charge();
        for (size_t itick = 0; itick < charges.size(); ++itick) {

          int ibin = (tbin1 - tbin.min() + itick);

          hist->SetBinContent(
              cbin.bin(ch) + 1, ibin + 1,
              charges[itick] + hist->GetBinContent(cbin.bin(ch) + 1, ibin + 1));
        }
      }
    }
  }
}

bool Root::ExampleROOTAna::operator()(const IFrame::pointer &frame) {

  if (!frame) {
    // eos
    log->debug("ExampleROOTAna: EOS");
    return true;
  }
  if (frame->traces()->empty()) {
    log->debug("ExampleROOTAna: passing through empty frame ID {}",
               frame->ident());
    return true;
  }

  /// open output file
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = "UPDATE";
  log->debug("ExampleROOTAna: opening for output: {} with mode {}", ofname,
             mode);
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());

  peak_frame(frame);

  fill_hist(frame, output_tf);

  /// write and close output file
  auto count = output_tf->Write();
  log->debug("ExampleROOTAna: closing output file {}, wrote {} bytes", ofname,
             count);
  output_tf->Close();
  delete output_tf;
  output_tf = nullptr;
  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
