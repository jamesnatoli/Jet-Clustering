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
   Int minindex_i, minindex_j, minindex_jB;
   float px, py, pz, newphi, neweta, newmass, newtheta, totalpm, newpT;

//This value can be changed easily and it denotes the size of the cones                             
   float Rrr = .4;

//Set Speed of light for Energy Calc
   const float c = 299792458;

//Function to find distance between particles
   float Distance_P (int a, int b)
   {
     deltaphi = particles.phi[a] - particles.phi[b];
     deltaeta = particles.eta[a] - particles.eta[b];
     Rij = hypot( deltaphi, deltaeta);
     Dij = min( pow( particles.pt[a], -2), pow( particles.pt[b]), -2);
     Dij = Dij * pow( ( Rij / Rrr), 2);
     return Dij;
   }

//Function for Beam Distance
   float Distance_J (int a)
   {
     Dib = pow(particles.pt[jspot], 2);
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
 } item;

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

//Loop over every pair of particle
 while ( !Particles.empty() )
   {
     for
//Stop condition while list is not empty
 while ( !Particles.empty() )
  {
    for (int ispot = 0; ispot < pleasestop; ispot++)
      {
	for (int jspot = ispot + 1; jspot < pleasestop; jspot++)
	  {
	    Dij = Distance_P(ispot, jspot);
	    Dib = Distance_J(jspot);
	  } 
	    //Determine if this is the smallest so far
	    if (Dij < Dijmin)
	      {
		Dijmin = Dij;
		minindex_i = ispot; //Save the particle index so they can be removed later
		minindex_j = jspot;
	      }    
	    //Determine if this is the smallest so far
	    if (Dib < Dibmin)
	      {
		Dibmin = Dib;
		minindex_jB = jspot; //Save the particle index so it can be removed later
	      }
	  }
      }
    //If the smallest is a beam, add to beam list and remove from particle list
    if (Dibmin < Dijmin)
      {
	Jets.push_back(

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
	    pleasestop--;
	  }
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
//Fill and Draw Histograms
histo1 -> SetMarkerStyle(21);
histo1 -> Fill(jet.pt.size());
histo1 -> Draw("");
c1 -> SaveAs("c1.gif");
}
