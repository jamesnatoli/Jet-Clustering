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

int main()
{
  struct fruit
  {
    int name;
    int color;
  }produce<int>, vegetables<int>;

  produce.insert( 0, 7);
  produce.insert( 1, 8);
  vegetables.insert( 0, 4);
}

int printstruct (struct bum)
{
  for (int i = 0; i < bum.size(); i ++)
    {
      bum[i].name
