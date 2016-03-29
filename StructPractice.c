#include <vector>
#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

//This program is used to test the workings of Structs                                                    
struct P
{
  int height;
  int weight;
  bool gender;
} item, part;

  
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
      cout << "First Height = " << (*it).height << endl;
      part = *it;
      cout << "Copied Height = " << part.height << endl;
    }

  //Vector Practice
  vector<float> *myvec;
  int sz = (*myvec).size();
  cout << sz << endl;

  float map;
  for (int i = 0; i < sz; i++)
    {
      map = (float)i;
      (*myvec[i]) = pow(map, 2);
    }
  
  float bob;
  cout << "The vector contains: ";
  for (int j = 0; j < sz; j++)
    {
      item.height =  (*myvec)[j];
      cout << " " << item.height;
    }

  return 0;
}

