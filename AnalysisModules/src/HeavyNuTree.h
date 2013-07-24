#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"

#ifndef HeavyNuTree_h_included
#define HeavyNuTree_h_included

class HeavyNuTree
{
private:
   TDirectory& m_file;
   TTree* m_tree;
   TBranch* branch_;
   int nevt;

public:

   struct HNuSlopeFitInfo
   {
      Float_t mlljj, mll, weight, l1pt, l1eta, l1phi, l1jdR, l2pt, l2eta, l2phi, l2jdR, j1pt, j1eta, j1phi, j2pt, j2eta, j2phi;
      Short_t flavor, l1id, l2id, cutlevel, n_pileup, n_primaryVertex;
      Float_t  rhE1, rhE2, sE1, sE2, met;
      Short_t j1B, j2B;
      Int_t run, ls, event;
   } event_;

   HeavyNuTree(TDirectory& f, bool writable) : m_file(f)
   {
      nevt = 0;
      if(writable)
      {
         m_file.cd();
         m_tree = new TTree("HeavyNuTuple", "HeavyNuTuple");
         branch_ = m_tree->Branch("slopefit", &event_, "mlljj/F:mll:weight:l1pt:l1eta:l1phi:l1jdR:l2pt:l2eta:l2phi:l2jdR:j1pt:j1eta:j1phi:j2pt:j2eta:j2phi:flavor/S:l1id:l2id:cutlevel:npu:npv:rhE1/F:rhE2:sE1:sE2:met:j1B/S:j2B:run/I:ls:event");
      }
      else
      {
         m_tree = (TTree*) m_file.Get("HeavyNuTuple");
         m_tree->SetBranchAddress("slopefit", &event_);
      }
   }

   int Entries()
   {
      return m_tree->GetEntries();
   }

   bool GetNextEvent()
   {
      if(nevt == Entries())
      {
         return false;
      }
      m_tree->GetEntry(nevt);
      nevt++;
      return true;
   }

   void clear()
   {
      event_.mlljj = event_.mll = event_.weight = 0.0;
      event_.flavor = event_.cutlevel = event_.n_pileup = event_.n_primaryVertex = 0;
   }

   void fill()
   {
      m_tree->Fill();
   }
};

#endif
