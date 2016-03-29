#include <vector>
#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <cmath>
#include <iterator>
using namespace std;

//This program is used to test the workings of Structs                                                    
struct P
{
  int height;
  int weight;
  bool gender;
} item, part;

list <P>::iterator Next ( list<P>::iterator nxt )
{
  list<P>::iterator junk;
  junk = nxt;
  junk++;
  return junk;
}

list <P>::iterator Next ( list<P>::iterator nxt, int n )  
{
  list<P>::iterator junk;
  junk = nxt;
  for (int bip = 0; bip < n; bip++)
    junk++;
  return junk;
}


int main()
{

  list <P> Person;

  item.height = 60;
  item.weight = 160;
  item.gender = true;
  Person.push_back(item);

  item.height = 54;
  item.weight = 145;
  item.gender = false;
  Person.push_back(item);

  item.height = 48;
  item.weight = 105;
  item.gender = false;
  Person.push_back(item);

  list <P>::iterator it, bit;
  it = Person.begin();
  Person.erase(it);
  for ( it = Person.begin(); it != Person.end(); it++)
    {
      bit = Next( it ); 
      cout << "First Height = " << it -> height << endl;
      part = *it;
      cout << "Copied Height = " << part.height << endl;
    }

  /*Vector Practice
  vector<float> *myvec = (10, 0);
  int sz = (*myvec).size();
  cout << sz << endl;

  float map;
  for (int i = 0; i < sz; i++)
    {
      map = (float)i;
      (*myvec)[i] = pow(map, 2);
    }
  
  float bob;
  cout << "The vector contains: ";
  for (int j = 0; j < sz; j++)
    {
      item.height =  (*myvec)[j];
      cout << " " << item.height;
    }
  */
  return 0;
}

