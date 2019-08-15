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
#include "WireCellUtil/Logging.h"

class TFile;

namespace WireCell {
namespace Root {

class ExampleROOTAna : public IFrameSink, public IConfigurable {
public:
  ExampleROOTAna();
  virtual ~ExampleROOTAna();

  /// IFrameFilter
  virtual bool operator()(const IFrame::pointer &in);

  /// IConfigurable
  virtual WireCell::Configuration default_configuration() const;
  virtual void configure(const WireCell::Configuration &config);

private:
  Configuration m_cfg;
  IAnodePlane::pointer m_anode;

  int m_nrebin;
  void create_file();

  Log::logptr_t log;
};
} // namespace Root
} // namespace WireCell

#endif
