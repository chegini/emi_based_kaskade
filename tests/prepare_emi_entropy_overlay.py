#!/usr/bin/env python3
"""Create a diagnostic-only symlink overlay; original Kaskade tree is untouched."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
source=a.source.resolve();target=a.output.resolve();target.mkdir(parents=True,exist_ok=False)
for f in source.iterdir():
 if f.name not in ('.git','mg'): (target/f.name).symlink_to(f,target_is_directory=f.is_dir())
(target/'mg').mkdir()
for f in (source/'mg').iterdir():(target/'mg'/f.name).symlink_to(f,target_is_directory=f.is_dir())
modified={}
def replace(name,fn):
 old=(source/'mg'/name).read_text();new=fn(old)
 if new==old:raise RuntimeError('No overlay edit: '+name)
 out=target/'mg'/name;out.unlink();out.write_text(new)
 modified[name]=dict(original_sha256=hashlib.sha256(old.encode()).hexdigest(),overlay_sha256=hashlib.sha256(new.encode()).hexdigest())
def topology(s):
 include='#include "mg/bddcMessageEntropy.hh"\n';s=include+s
 old='        nodeTopology.configure();';assert s.count(old)==1
 return s.replace(old,old+'''\n        if (!MessageEntropy::prefix().empty()) {
          std::vector<int> entropyRoutes(mpiSize());
          for(int p=0;p<mpiSize();++p) entropyRoutes[p]=nodeTopology.route(mpiRank(),p);
          std::vector<int> entropyOwners(size());
          for(int id=0;id<size();++id) entropyOwners[id]=mpiOwner(id);
          MessageEntropy::topology(mpiRank(),mpiSize(),entropyRoutes,entropyOwners);
        }''')
replace('bddc.hh',topology)
def routed(s):
 old='        packetSentBytes+=destinations[target].size();';assert s.count(old)==2
 pieces=s.split(old)
 return pieces[0]+old+'\n        MessageEntropy::routed("restrict",nodeTopology.route(mpiRank(),target),destinations[target]);'+pieces[1]+old+'\n        MessageEntropy::routed("prolongate",nodeTopology.route(mpiRank(),target),destinations[target]);'+pieces[2]
replace('bddc.hpp',routed)
def coefficients(s):
 old='if (!symbolHistogramDirectory.empty())';assert s.count(old)==6
 s=s.replace(old,'if (!symbolHistogramDirectory.empty() || !MessageEntropy::prefix().empty())')
 for phase in ('restrict','prolongate'):
  old='SymbolHistogram::write(symbolHistogramDirectory,"'+phase+'",'
  assert s.count(old)==1
  s=s.replace(old,'SymbolHistogram::write("subdomain"+std::to_string(this->id),"'+phase+'",')
 return s
replace('bddcSpaceTransfer.hpp',coefficients)

f=target/'mg/bddcSymbolHistogram.hh';f.unlink();f.write_text('''// Diagnostic adapter: full-run bounded entropy counters, no per-packet CSVs.
#ifndef KASKADE_BDDC_SYMBOL_HISTOGRAM_HH
#define KASKADE_BDDC_SYMBOL_HISTOGRAM_HH
#include "mg/bddcMessageEntropy.hh"
namespace Kaskade::BDDC::SymbolHistogram {
inline bool physicalEnabled(std::uint64_t){return false;}
template<class... T> void writePhysical(T const&...){throw std::runtime_error("Physical histogram disabled in entropy build");}
template<class U> void write(std::string const& directory,std::string const& phase,
 std::uint64_t,int peer,int bits,std::size_t,std::vector<U> const& quantized,
 std::vector<U> const& mapped,std::vector<U> const& lengths,bool encoded,bool floating,bool integer){
 if(!directory.empty())MessageEntropy::coefficients(directory,phase,peer,bits,quantized,mapped,lengths,encoded,floating,integer);
}
}
#endif
''')
header=Path(__file__).with_name('emi_message_entropy.hh');(target/'mg/bddcMessageEntropy.hh').write_bytes(header.read_bytes())
(target/'entropy-overlay.json').write_text(json.dumps(dict(source=str(source),changes=modified,
 notes=['Observer-only overlay; production sources unchanged.','Observer reads uint32 coefficients directly, including 32-bit precision; legacy histogram CLI remains disabled.','Existing full-numerical/codec implementation retained.']),indent=2)+'\n')
print('Entropy overlay:',target)
