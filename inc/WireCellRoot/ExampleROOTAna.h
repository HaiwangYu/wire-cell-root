/** Sink data to a file format used by the "Magnify" GUI .
 *
 * This is technically a "filter" as it passes on its input.  This
 * allows an instance of the sink to sit in the middle of some longer
 * chain.
 *
 * FIXME: currently this class TOTALLY violates the encapsulation of
 * DFP by requiring the input file in order to transfer input data out
 * of band of the flow.
 */

#ifndef WIRECELLROOT_EXAMPLEROOTANA
#define WIRECELLROOT_EXAMPLEROOTANA

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellUtil/Logging.h"

class TFile;

namespace WireCell {
namespace Root {

class ExampleROOTAna : public IFrameFilter, public IConfigurable {
public:
  ExampleROOTAna();
  virtual ~ExampleROOTAna();

  /// working operation - interface from IFrameFilter
  /// executed when called by pgrapher
  virtual bool operator()(const IFrame::pointer &in, IFrame::pointer& out);

  /// interfaces from IConfigurable

  /// exeexecuted once at node creation
  virtual WireCell::Configuration default_configuration() const;

  /// executed once after node creation
  virtual void configure(const WireCell::Configuration &config);

private:
  Configuration m_cfg; /// copy of configuration
  IAnodePlane::pointer m_anode; /// pointer to some APA, needed to associate chnnel ID to planes

  /// recreate output file
  void recreate_out_file() const;

  /// print trace tags and frame tags
  void peak_frame(const IFrame::pointer &frame) const;

  /// fill time-channel 2D hist for each trace tag
  void fill_hist(const IFrame::pointer &frame, TFile *output_tf) const;

  /// fill time-channel 2D hist for each trace tag
  IFrame::pointer fft_frame(const IFrame::pointer &in) const;

  /// SPD logger
  Log::logptr_t log;
};
} // namespace Root
} // namespace WireCell

#endif
