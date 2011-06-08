#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "TVector3.h"

void HeavyNuEvent::regularize() {
  mu[0]=mu1;
  if (!mu2.isNull()) mu[1] = mu2 ; 
  if (!e1.isNull()) e[1] = e1 ; 
  // mu[1]=mu2;
  j[0]=j1;
  j[1]=j2;
  tjV[0]=tjV1;
  tjV[1]=tjV2;
}

static double planeCosAngle(const reco::Particle::Vector& plane1,
			    const reco::Particle::Vector& plane2,
			    const reco::Particle::Vector& vect) {
  reco::Particle::Vector cp=plane1.unit().Cross(plane2.unit()).unit();
  double rv=fabs(vect.unit().Dot(cp));
  //  std::cout << plane1 << " " << plane2 << " " << vect << " " << rv << std::endl;
  return rv;
}

void HeavyNuEvent::calculateMuMu(double muptfactor) {
  MESscale = muptfactor;
  vMuMu    = MESscale*(mu1->p4()+mu2->p4());
  mMuMu    = vMuMu.M();
  czeta_mumu   =mu1->momentum().unit().Dot(mu2->momentum().unit());
  ctheta_mumu  =planeCosAngle(mu1->momentum(),mu2->momentum(),reco::Particle::Vector(0,0,1));
}

void HeavyNuEvent::calculateMuE(double muptfactor,double elefactor) {
  MESscale = muptfactor;
  EEScale  = elefactor;
  vMuMu    = MESscale*mu1->p4() + EEScale*e1->p4();
  mMuMu    = vMuMu.M();
}

void HeavyNuEvent::calculate(int nMu) {

  reco::Particle::LorentzVector j1p4 = j1->p4();
  reco::Particle::LorentzVector j2p4 = j2->p4();

  reco::Particle::LorentzVector mu1p4 = mu1->p4();
  reco::Particle::LorentzVector mu2p4 = ( (nMu == 2) ? mu2->p4() : e1->p4() );
  // reco::Particle::LorentzVector mu2p4 = mu2->p4();

  // if doing JECU studies, apply scaling factor here
  //
  if( j1scale != 1.0 ) j1p4 *= j1scale;
  if( j2scale != 1.0 ) j2p4 *= j2scale;
 
  reco::Particle::Vector j1mom = j1p4.Vect();
  reco::Particle::Vector j2mom = j2p4.Vect();

  // if doing MES studies, apply scaling factor here
  //
  if ( nMu == 2 ) {
    if ( MESscale != 1.0 ) { mu1p4 *= MESscale; mu2p4 *= MESscale; }
  } else { 
    if ( MESscale != 1.0 ) mu1p4 *= MESscale;  
    if ( EEScale != 1.0 )  mu2p4 *= EEScale; 
  }

  reco::Particle::Vector mu1mom = mu1p4.Vect();
  reco::Particle::Vector mu2mom = mu2p4.Vect();

  vMuMu  = mu1p4+mu2p4;
  vJJ    = j1p4+j2p4;
  lv_evt = vMuMu+vJJ;

  czeta_mumu   =fabs(mu1mom.unit().Dot(mu2mom.unit()));

  ctheta_mumu  =planeCosAngle(mu1mom,mu2mom,reco::Particle::Vector(0,0,1));
  ctheta_jj    =planeCosAngle(j1mom,j2mom,reco::Particle::Vector(0,0,1));
  ctheta_mu1_jj=planeCosAngle(j1mom,j2mom,mu1mom);
  ctheta_mu2_jj=planeCosAngle(j1mom,j2mom,mu2mom);

  // LorentzVector of just the Z deboost.
  reco::Particle::LorentzVector deboostz(0,0,-lv_evt.pz(),lv_evt.pz());
  
  reco::Particle::LorentzVector mu1z=mu1p4+deboostz;
  reco::Particle::LorentzVector mu2z=mu2p4+deboostz;
  reco::Particle::LorentzVector j1z=j1p4+deboostz;
  reco::Particle::LorentzVector j2z=j2p4+deboostz;
  
  cthetaz_mumu   = planeCosAngle(mu1z.Vect(),mu2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_jj     = planeCosAngle(j1z.Vect(),j2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_mu1_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),mu1z.Vect());
  cthetaz_mu2_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),mu2z.Vect());

  float dRmu1jet1 = deltaR( mu1->eta(), mu1->phi(), j1->eta(), j1->phi() ) ; 
  float dRmu1jet2 = deltaR( mu1->eta(), mu1->phi(), j2->eta(), j2->phi() ) ; 
  float dRmu2jet1 = (( e1.isNull() ) ? 
 		     (deltaR( mu2->eta(), mu2->phi(), j1->eta(), j1->phi() )) : 
 		     (deltaR( e1->eta(), e1->phi(), j1->eta(), j1->phi() ))) ; 
  float dRmu2jet2 = (( e1.isNull() ) ? 
 		     (deltaR( mu2->eta(), mu2->phi(), j2->eta(), j2->phi() )) : 
 		     (deltaR( e1->eta(), e1->phi(), j2->eta(), j2->phi() ))) ; 
  
  // find the closest jets
  dRminMu1jet = std::min( dRmu1jet1,dRmu1jet2 );
  dRminMu2jet = std::min( dRmu2jet1,dRmu2jet2 );

  // what are the muon transverse momenta relative to the closest jets?
  reco::Particle::Vector jmom4mu1 = (dRminMu1jet == dRmu1jet1) ? j1mom : j2mom;
  reco::Particle::Vector jmom4mu2 = (dRminMu2jet == dRmu2jet1) ? j1mom : j2mom;

  TVector3 mu1vec( mu1mom.X(), mu1mom.Y(), mu1mom.Z() );
  TVector3 mu2vec( mu2mom.X(), mu2mom.Y(), mu2mom.Z() );

  TVector3 jt4mu1vec( jmom4mu1.X(), jmom4mu1.Y(), jmom4mu1.Z() );
  TVector3 jt4mu2vec( jmom4mu2.X(), jmom4mu2.Y(), jmom4mu2.Z() );

  ptrelMu1 = mu1vec.Perp( jt4mu1vec );
  ptrelMu2 = mu2vec.Perp( jt4mu2vec );

  // Composite objects
  mJJ   = vJJ.M();
  mMuMu = vMuMu.M();

  mWR   = lv_evt.M();

  mNuR1 = (vJJ + mu1p4).M();
  mNuR2 = (vJJ + mu2p4).M();
}

void HeavyNuEvent::decayID(const HepMC::GenEvent& genE) {

  HepMC::GenEvent::vertex_const_iterator vtex;
  HepMC::GenVertex::particles_out_const_iterator Pout;

  mc_class=0;

  for (vtex=genE.vertices_begin();vtex!=genE.vertices_end();vtex++){
    if(((*vtex)->particles_in_size())==1){ // We want a decay, not collision
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==9900024){ // Is a WR
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 1;
            break;
	    
          case 13:	// mu
            mc_class = 2;
            break;
	    
          case 15:  // tau
            mc_class = 3;
            break;
	    
          default:	// else
	    break;
	  }
	  if (mc_class!=0) break;
	
	}
	break; //end of id loop
      }
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==23){ // Is a Z
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 11;
            break;
	    
          case 13:	// mu
            mc_class = 12;
            break;
	    
          case 15:  // tau
            mc_class = 13;
            break;
	    
          default:	// else
            break;
	  }
	  if (mc_class!=0) break;	
	  
	}
	break; //end of id loop
      }
    }

  }
}
