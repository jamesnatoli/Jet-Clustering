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
} item;

int main()
{

  list <P> Person;

  item.height = 60;
  item.weight = 160;
  item.gender = true;

  Person.push_front(item);

  item.height = 54;
  item.weight = 145;
  item.gender = false;

  Person.push_front(item);

  list <P>::iterator it;
  for ( it = Person.begin(); it != Person.end(); it++)
    cout << (*it).height << endl ;
  
  return 0;
}

