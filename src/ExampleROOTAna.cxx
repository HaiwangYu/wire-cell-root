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

Root::ExampleROOTAna::ExampleROOTAna()
    : m_nrebin(1), log(Log::logger("ana")) {}

Root::ExampleROOTAna::~ExampleROOTAna() {}

void Root::ExampleROOTAna::configure(const WireCell::Configuration &cfg) {
  std::string fn;

  fn = cfg["input_filename"].asString();
  if (fn.empty() and !cfg["shunt"].empty()) {
    log->error("ExampleROOTAna: asked to shunt but not given input file name");
    THROW(ValueError() << errmsg{"ExampleROOTAna: must provide input filename "
                                 "to shunt objects to output"});
  }

  fn = cfg["output_filename"].asString();
  if (fn.empty()) {
    THROW(ValueError() << errmsg{
              "Must provide output filename to ExampleROOTAna"});
  }

  auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
  m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

  m_cfg = cfg;

  m_nrebin = get<int>(cfg, "nrebin", m_nrebin);

  create_file();
}

WireCell::Configuration Root::ExampleROOTAna::default_configuration() const {
  Configuration cfg;

  cfg["anode"] = "AnodePlane";

  // fixme: this TOTALLY violates the design of wire cell DFP.
  cfg["input_filename"] = "";

  // List of TObjects to copy from input file to output file.
  cfg["shunt"] = Json::arrayValue;

  // Name of ROOT file to write.
  cfg["output_filename"] = "";

  // A list of trace tags defining which waveforms are saved to Magnify
  // histograms.
  cfg["frames"] = Json::arrayValue;

  // If no tags for traces, i.e. trace_has_tag=false in a frame,
  // set desired tag to ""
  // as FrameTool::tagged_traces(frame, "") calls untagged_traces(frame).
  cfg["trace_has_tag"] = true;

  // A list of pairs mapping a cmm key name to a ttree name.
  cfg["cmmtree"] = Json::arrayValue;

  // The ROOT file mode with which to open the file.  Use "RECREATE"
  // to overrite an existing file.  This might be useful for the
  // first ExampleROOTAna in a chain.  Use "UPDATE" for subsequent
  // sinks that add to the file.
  cfg["root_file_mode"] = "RECREATE";

  // If runinfo is given it should be a JSON object and its values
  // will be copied into the Trun tree.  If instead it is null AND
  // an input file is given AND it contains a Trun tree, it will be
  // copied to output.
  cfg["runinfo"] = Json::nullValue;

  cfg["nrebin"] = 1;

  // List tagged traces from which to save the "trace summary"
  // vector into a 1D histogram which will be named after the tag.
  // See "summary_operator".
  cfg["summaries"] = Json::arrayValue;

  // An object mapping tags to operators for aggregating trace
  // summary values on the same channel.  Operator may be "sum" to
  // add up all values on the same channel or "set" to assign values
  // to the channel bin (last one wins).  If a tag is not found, the
  // default operator is "sum".
  cfg["summary_operator"] = Json::objectValue;

  return cfg;
}

namespace {
typedef std::unordered_set<std::string> string_set_t;
string_set_t getset(const WireCell::Configuration &cfg) {
  string_set_t ret;
  for (auto jone : cfg) {
    ret.insert(jone.asString());
  }
  return ret;
}

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
}

void Root::ExampleROOTAna::create_file() {
  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = "RECREATE";
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());
  output_tf->Close("R");
  delete output_tf;
  output_tf = nullptr;
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

  const std::string ofname = m_cfg["output_filename"].asString();
  const std::string mode = m_cfg["root_file_mode"].asString();
  log->debug("ExampleROOTAna: opening for output: {} with mode {}", ofname,
             mode);
  TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());

  for (auto tag : getset(m_cfg["frames"])) {

    auto trace_tag = tag;
    auto trace_has_tag = m_cfg["trace_has_tag"].asBool();
    if (!trace_has_tag) {
      trace_tag = "";
      log->debug("ExampleROOTAna: set desired trace tag to \"\" as "
                 "cfg::trace_has_tag=false");
    }

    ITrace::vector traces_byplane[3],
        traces = FrameTools::tagged_traces(frame, trace_tag);
    if (traces.empty()) {
      log->warn("ExampleROOTAna: no tagged traces for \"{}\"", tag);
      continue;
    }

    log->debug("ExampleROOTAna: tag: \"{}\" with {} traces", tag,
               traces.size());

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
      int nbins = tbin.nbins() / m_nrebin;

      TH2F *hist =
          new TH2F(name.c_str(), name.c_str(), cbin.nbins(), cbin.min(),
                   cbin.max(), nbins, tbin.min(), tbin.max());

      hist->SetDirectory(output_tf);

      for (auto trace : traces_byplane[iplane]) {
        const int tbin1 = trace->tbin();
        const int ch = trace->channel();
        auto const &charges = trace->charge();
        for (size_t itick = 0; itick < charges.size(); ++itick) {

          int ibin = (tbin1 - tbin.min() + itick) / m_nrebin;

          hist->SetBinContent(
              cbin.bin(ch) + 1, ibin + 1,
              charges[itick] + hist->GetBinContent(cbin.bin(ch) + 1, ibin + 1));
        }
      }
    }
  }

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
