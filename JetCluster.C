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
#include <list>
#include <map>

using namespace std;

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
   //TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry < nentries;jentry++) 
     {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     }
  

//void JetCLuster::Algorithm(){

//These variables are used in the calculations                                                       
   float Dib, Dij, Rij, deltaphi, deltaeta, Dijmin, Dibmin;
   //Int minindex_i, minindex_j, minindex_jB;
   float px, py, pz, newphi, neweta, newmass, newtheta, totalpm, newpT;

//This value can be changed easily and it denotes the size of the cones                             
   float Rrr = .4;

//Set Speed of light for Energy Calc
   const float c = 299792458;

//Function to find distance between particles
   float Distance_P (float aphi, float bphi, float ceta, float deta, float ept, float fpt)
   {
     deltaphi = aphi - bphi;
     deltaeta = ceta - deta;
     Rij = hypot( deltaphi, deltaeta);
     Dij = min( pow( ept, -2), pow( fpt, -2);
     Dij = Dij * pow( ( Rij / Rrr), 2);
     return Dij;
   }

//Function for Beam Distance
   float Distance_J (float gpt)
   {
     Dib = pow( gpt, 2);
     return Dib;
   }

//Function for Phi Wrap Soln

//Function for Adding Momentum
   

//Create Canvas for Histogram
 TCanvas *c1 = new TCanvas("c1", "demo", 200, 10, 700, 500);
 c1 -> SetFillColor(42);

 TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);
 int numentry = pt->size();

 //Declare a struct (P) to hold arrays which will hold all of the values. The object particle of type P can access all of these values
 struct P 
 {
   float pt;
   float phi;
   float eta;
   float mass;
 } item, minindex_i, minindex_j, minindex_jB;
 
 list<P> Particles;
 list<P> Jets;
 
 //Fill the array "particles" with the values from the vectors
 //Check pointer to array syntax for accesing element
 for (int wh = 0; wh < (particles.pt.size() - 1); wh++)
   {
     item.pt = pt[wh];
     item.phi = phi[wh];
     item.eta = eta[wh];
     item.mass = mass[wh];
     Particles.push_back( item );
   }
 
 //This is to make sure that the first value checked is assigned to the smallest
 Dijmin = 100000;
 Dibmin = 100000;
 
 /* struct Dist_ij
 {
   float dist;
   P iparticle;
   P jparticle;
 } dist_ij;
 
 list<Dist_ij> ParticleDist;
 
 struct Dist_jb
 {
   float dist;
   P jbeamparticle;
 } dist_jb;
 
 list<Dist_jb> BeamDist;
 list<float> AllParticleDists;
 list<float> AllBeamDists;
 */

 //Loop over every pair of particle 
 list <P>::iterator it;//Create iterators to loop through lists
 list <P>::jterator jt;
 while ( !Particles.empty() ) //Stop condition
   {
     for( it = Particles.begin(); it != Particles.end(); it++)
       {
	 for (jt = it + 1; jt != Particles.end(); jt++)
	   {
	     Dij = Distance_P( (*it).phi, (*jt).phi, (*it).eta, (*jt).eta, (*it).pt, (*jt).pt );
	     Dib = Distance_J( (*jt).pt );
	     if (Dij < Dijmin) //Determine if this is the smallest so far     
	       { 
		 Dijmin = Dij;
		 minindex_i = *it;                 	 
		 minindex_j = *jt;                                                                                      
	       }
	     if (jt = (it + 1))
	       {
		 if (Dib < Dibmin)
		   {
		     Dibmin = Dib;
		     minindex_jB = *jt; //Save the particle so it can be removed later                            
		   }
	       }
	   }
       }
     
     //If the smallest is a beam, add to beam list and remove from particle list
     if (Dibmim < Dijmin)
       {
	 Jets.push_back( minindex_jB );
	 Particles.remove
     for( int kspot = minindex_jB; kspot < (particles.pt.size() - 2); kspot++)
       {
	 particles.pt[kspot] = particles.pt[kspot + 1];
	 particles.eta[kspot] = particles.eta[kspot + 1];
	 particles.phi[kspot] = particles.phi[kspot + 1];
	 particles.mass[kspot] = particles.mass[kspot + 1];	    
	 pleasestop--;
       }
     
     //If the smallest is not a beam, add momenta, add to list, and remove other two particles
     
     else
       {
	 //ADDING MOMENTUM
	 
	 //Adds the components of momentum in the x, y, and z directions
	 px = (particles.pt[minindex_i] * cos(particles.phi[minindex_i])) + (particles.pt[minindex_j] * cos(particles.phi[minindex_j]));
	 py = (particles.pt[minindex_i] * sin(particles.phi[minindex_i])) + (particles.pt[minindex_j] * sin(particles.phi[minindex_j]));
	 pz = (particles.pt[minindex_i] * sinh(particles.eta[minindex_i])) + (particles.pt[minindex_i] * sinh(particles.eta[minindex_j]));
	 
	 //Calculate new values for new paricle
	 newpT = hypot(px, py);
	 totalpm = sqrt((px * px) + (py * py) + (pz * pz));
	 newtheta = asin( newpT / totalpm);
	 neweta = -log( tan( newtheta/ 2));
	 newphi = acos( px / (totalpm * sin(newtheta)));
	 //***Add energies not mass***
	 newmass = particles.mass[minindex_i] + particles.mass[minindex_j];
	 
	 //Remove the j particle and shift everything down
	 for( int kspot = minindex_j; kspot < pleasestop; kspot++)
	   {
	     particles.pt[kspot] = particles.pt[kspot + 1];
	     particles.eta[kspot] = particles.eta[kspot + 1];
	     particles.phi[kspot] = particles.phi[kspot + 1];
	     particles.mass[kspot] = particles.mass[kspot + 1];
	   }
	 //Assign new particle values to the i particle;
	 particles.pt[minindex_i] = newpT;
	 particles.eta[minindex_i] = neweta;
	 particles.phi[minindex_i] = newphi;
	 particles.mass[minindex_i] = newmass;
       }
   }
   }
   //Fill and Draw Histograms
   histo1 -> SetMarkerStyle(21);
   histo1 -> Fill(jet.pt.size());
   histo1 -> Draw("");
   c1 -> SaveAs("c1.gif");
}
