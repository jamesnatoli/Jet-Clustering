#define JetCluster_cxx
#include "JetCluster.h"
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
#include <vector>
#include <iterator>

using namespace std;

//These variables are used in the calculations                                        
float Dib, Dij, Rij, deltaphi, deltaeta, Dijmin, Dibmin;

//This value can be changed easily and it denotes the size of the cones                
float Rrr = .4;

//Set Speed of light for Energy Calc
const float csq = pow(299792458, 2);

//Set Pi for Phi Wrap
const double PI = 4 * atan(1);

//Function to find distance between particles
float Distance_P (float aphi, float bphi, float ceta, float deta, float ept, float fpt)
{
  deltaphi = aphi - bphi;
  deltaeta = ceta - deta;
  Rij = hypot( deltaphi, deltaeta);
  Dij = min( pow( ept, -2), pow( fpt, -2));
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
float PhiWrap( float val )
{
  float output;
  if( val > PI)
    output = ( val - (2 * PI));
  else
    output = val;
  return output;
}

//Function for Adding Momentum
// float Addmomentum ( list<P>iterator:: it1, list<P>iterator:: it2)
//{

//}

//Declare a struct (P) to hold arrays which will hold all of the values. The object particle of type P can acess all of these values
struct JVN_particle
{
  float pt;
  float phi;
  float eta;
  float mass;
  float energy;
};

JVN_particle item;

list <JVN_particle> Particles;
list <JVN_particle> Jets;

//Function to get next iterator       
vector <JVN_particle>::iterator Next ( vector <JVN_particle>::iterator nxt )
{
  vector<JVN_particle>::iterator junk;
  junk++;
  return junk;
}

//Alternate Declaration for next                      
vector <JVN_particle>::iterator Next ( list<JVN_particle>::iterator nxt, int n )
{
  vector<JVN_particle>::iterator junk;
  junk = nxt;
  for (int bip = 0; bip < n; bip++)
    junk++;
  return junk;
}

  vector <JVN_particles>::iterator it, jt, minindex_i, minindex_j, minindex_jB;//Create iterators to loop through list 

//Create Instance Variables                                    
float px, py, pz, newphi, neweta, newmass, newtheta, totalpm, newpT;

void JetCluster::Loop()
{

//   In a ROOT session, you can do:
//      Root > .L JetCluster.C
//      Root > JetCluster t
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

//Create Canvas for Histogram                                                                        
  TCanvas *c1 = new TCanvas("c1", "demo", 200, 10, 700, 500);
  c1 -> SetFillColor(42);
  
  TH1F* histo1 = new TH1F("histo1", "Number of Jets", 100, 0, 200);
  histo1 -> SetMarkerStyle(21);

  if (fChain == 0) 
    return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  int numentry = pt->size();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      {
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
	// if (Cut(ientry) < 0) continue;
      }
      
      //Fill the array "particles" with the values from the vectors
      //Check pointer to array syntax for accesing element
      for (int wh = 0; wh < numentry - 1; wh++)
	{
	  item.pt = (*pt)[wh];
	  item.phi = PhiWrap( (*phi)[wh] );
	  item.eta = (*eta)[wh];
	  item.mass = (*mass)[wh];
	  item.energy = (*mass)[wh] * csq;
	  Particles.push_back( item );
	}
      
      //This is to make sure that the first value checked is assigned to the smallest
      Dijmin = 100000;
      Dibmin = 100000;
      
      while ( !Particles.empty() ) //Stop condition
	 {
	   for( it = Particles.begin(); it != Particles.end(); it++)
	     {
	       for (jt = Next(it); jt != Particles.end(); jt++)
		 {
		   Dij = Distance_P( (*it).phi, (*jt).phi, (*it).eta, (*jt).eta, (*it).pt, (*jt).pt );
		   Dib = Distance_J( (*jt).pt );
		   if (Dij < Dijmin) //Determine if this is the smallest so far     
		     { 
		       Dijmin = Dij;
		       minindex_i = it;  
		       minindex_j = jt; 
		     }
		   if (distance( it, jt) >= 0)
		     {
		       if (Dib < Dibmin)
			 {
			   Dibmin = Dib;
			   minindex_jB = jt; //Save the particle so it can be removed later
			 }
		     }
		 }
	     }
	   
	   //If the smallest is a beam, add to beam list and remove from particle list
	   if (Dibmin < Dijmin)
	     {
	       Jets.push_back( (*minindex_jB) );
	       Particles.erase( minindex_j );
	     }
	   
	   //If the smallest is not a beam, add momenta, add to list, and remove other two particles     
	   else
	     {
	       //ADDING MOMENTUM
	       
	       //Adds the components of momentum in the x, y, and z directions
	       px = ( ( (*minindex_i).pt * cos((*minindex_i).phi) ) + ((*minindex_j).pt * cos((*minindex_j).phi) ) );
	       py = ( ( (*minindex_i).pt * sin((*minindex_i).phi) ) + ((*minindex_j).pt * sin((*minindex_j).phi) ) );
	       pz = ( ( (*minindex_i).pt * sinh((*minindex_i).eta) ) + ((*minindex_i).pt * sinh((*minindex_j).eta) ) );
	             
	       //Calculate new values for new paricle
	       newpT = hypot(px, py);
	       totalpm = sqrt((px * px) + (py * py) + (pz * pz));
	       newtheta = asin( newpT / totalpm);
	       neweta = -log( tan( newtheta/ 2));
	       newphi = acos( px / (totalpm * sin(newtheta)));
	       //Add energies not mass
	       newmass = ( ( (*minindex_i).energy + (*minindex_j).energy ) / csq);
	             
	       //Build New Struct with all new values
	       item.pt = newpT;
	       item.eta = neweta;
	       item.phi = newphi;
	       item.mass = newmass;
	       item.energy = (newmass * csq);
	       
	       //Remove the Particles and add the new one
	       Particles.erase( minindex_j );
	       Particles.push_back( item );
	     }
	 }
      histo1 -> Fill( Jets.size() );
    }
  histo1 -> Draw("");
  c1 -> SaveAs("c1.gif");
}
