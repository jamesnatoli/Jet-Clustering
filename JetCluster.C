#define JetCLuster_cxx
#include "JetCLuster.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <string>
#include <strstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>

void JetCLuster::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L JetCLuster.C
//      Root > JetCLuster t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry < nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}  

//void JetCLuster::Algorithm(){

//Create Canvas for Histogram
TCanvas *c1 = new TCanvas("c1", "demo", 200, 10, 700, 500);
c1 -> SetFillColor(42);

TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);

Long64_t nentries = fChain->GetEntriesFast();
//Declare a struct (P) to hold arrays which will hold all of the values. The object particle of type P can access all of these values
     struct P
     {
       float pt[nentries];
       float phi[nentries];
       float eta[nentries];
       float mass[nentries];
     } particles, jets;

     //solve the phi wrap problem?

     //Fill the array "particles" with the values from the vectors
for (int wh = 0; wh < (particles.pt.size() - 1); wh++)
       {
	 particles.pt[wh] = pt->GetEntry(wh);
	 particles.phi[wh] = phi->GetEntry(wh);
	 particles.eta[wh] = eta->GetEntry(wh);
	 particles.mass[wh] = mass->GetEntry(wh);
       }
     
     //These variables are used in the calculations
     float Dib, Dij, Rij, deltaphi, deltaeta, Dijmin, Dibmin, minindex_i, minindex_j, minindex_jB;

     //This value can be changed easily and it denotes the size of the cones
     float Rrr = .4;
     
     //This is to make sure that the first value is smaller than this
     Dijmin = 100000;
     Dibmin = 100000;

     //Loop over every pair of particle
for (int ispot = 0; ispot < (particles.pt.size() - 1); ispot++)
  {
    for (int jspot = ispot + 1; jspot < (particles.pt.size() - 1); jspot++)
      {
	Dib = pow(particles.pt[jspot], 2);
        deltaphi = particles.phi[jspot] - particles.pho[ispot];
        deltaeta = particles.eta[jspot] - particles.eta[ispot];
        Rij = hypot( deltaphi, deltaeta);
        Dij = min( pow( particles.pt[jspot], 2), pow( particles.pt[ispot]), 2);
        Dij = Dij * pow( ( Rij / R), 2);
        
        //Determine if this is the smallest so far
        if (Dij < Dijmin)
          {
   	 Dijmin = Dij;
   	 //Save the particle index so they can be removed late
       	 minindex_i = ispot;
       	 minindex_j = jspot;
          }
            
        //Determine if this is the smallest so far
        if (Dib < Dibmin)
          {
	    Dibmin = Dib;
	    //Save the particle index so it can be removed later
	    minindex_jB = jspot; 
	  }
      }
  }
    //If the smallest is a beam, add to beam list and remove from particle list
    if (Dibmin < Dijmin)
      {
	jets.pt[0] = particles.pt[minindex_jB];
	jets.eta[0] = particles.eta[minindex_jB];
	jets.phi[0] = particles.phi[minindex_jB];
	jets.mass[0] = particles.mass[minindex_jB];
        for( int kspot = minindex_jB; kspot < (particles.pt.size() - 2); kspot++)
	  {
	    particles.pt[kspot] = particles.pt[kspot + 1];
	    particles.eta[kspot] = particles.eta[kspot + 1];
	    particles.phi[kspot] = particles.phi[kspot + 1];
	    particles.mass[kspot] = particles.mass[kspot + 1];	    
	  }
      }
    //If the smallest is not a beam, add momenta???, add to list, and remove other two particles
    else
      {
	//add momenta
	for( int kspot = minindex_j; kspot < (particles.pt.size() - 2); kspot++)
          {
            particles.pt[kspot] = particles.pt[kspot + 1];
            particles.eta[kspot] = particles.eta[kspot + 1];
            particles.phi[kspot] = particles.phi[kspot + 1];
            particles.mass[kspot] = particles.mass[kspot + 1];
	  }
	particles.pt[minindex_i] = ;//added momenta?
	particles.eta[minindex_i] = ;//??
	particles.phi[minindex_i] = ;//??
	particles.mass[minindex_i] = ;//??
      }	     

     //float pti, ptj, ptBeam, phii, phij, massi, massj, etai, etaj;

     TBranch *beampt = T->Branch("pt", &pt, "pt/F");
     T-> SetBranchAddress("pt", &pt);
     T-> SetBranchAddress("eta", &eta);
     T-> SetBranchAddress("phi", &phi);
     T-> SetBranchAddress("mass", &mass);

     TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);
     TH1F* histo2 = new TH1F("histo2", "Test for pT", 100, 0, 200);
     float Dib, Dij, Rij;
     float min = 1000000;
     float Rrr = .4;
     for (Long64_t ispot = 0; ispot < nentries ; ispot++)
       {
    	 for (Long64_t jspot = ispot + 1; jspot < nentries; jspot++)
	   {
	     //Anti -kT algorithm                                                            
	     Dib = pow( particle.pt[jspot], -2);
	     Rij = hypot( (particle.phi[jspot] - phij), (etai - etaj) );
	     Dij = min( pow( ptj, -2), pow( pti, -2)) * ( pow( (Rij/R), 2));
	     
	 
//Attempt to solve phi wrap problem
	     //if ( phij > 180)
	     //{
	     //	 float temp = phij - 360;
	     //	 phi.assign (jspot, temp);
	     //}
	     
	     //Assign these temporary variables to the values from the vector
	     phii = phi-> GetEntry(i);
	     phij = phi-> GetEntry(j);
	     etai = eta-> GetEntry(i);
	     etaj = eta-> GetEntry(j);
	     pti = pt-> GetEntry(i);
	     ptj = pt-> GetEntry(j);

	     //Anti -kT algorithm
	     Dib = pow( ptj], -2);
	     Rij = hypot( (phii - phij), (etai - etaj) );
	     Dij = min( pow( ptj, -2), pow( pti, -2)) * ( pow( (Rij/R), 2));
	     if ( min( Dij, Dib) == Dib)
	       {
		 //remove jspot from particle list 
		 //add to beam list?
	       }
	     else if ( min (Dif, Dib) == Dij)
	       {
		 //add 4 momenta? How do I access px, py, pz
		 //remove jspot and ispot and add this new particle, but to where? fchain? or individually from every tree?
	       }
       }
	   }
