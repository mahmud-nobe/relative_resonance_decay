// Implementation of our vec2d class for ROOT
// now we have to include our header
// file !!!

#include "vec2d.h"

#include <iostream>  

ClassImp(vec2d)

using namespace std;

vec2d vec2d::operator+(vec2d  &v)
{
  return vec2d(v.GetX()+GetX(),v.GetY()+GetY());
}

vec2d vec2d::operator+=(vec2d &v)
{
  //DEBUG:
  //cout<<x<<" "<<v.GetX()<<endl;
  //cout<<y<<" "<<v.GetY()<<endl;
  
  x += v.GetX(); y +=v.GetY();

  //DEBUG:
  //cout<<"After Add : x = "<<x<<" "<<y<<endl;

  return *this;
}

vec2d::vec2d()
{
  //default constructor
  //DEBUG:
  //cout<<" vec2d: default constructor called"<<endl;

  x=0;y=0;
}

vec2d::vec2d(double mx, double my)
{
  //user defined constructor
  //cout<<" vec2d: user defined constructor called"<<endl;

  x=mx; y=my;
}

vec2d::~vec2d() 
{
  //cout<<" vec2d: default destructor called"<<endl;
}; 

void vec2d::Print()
{
  cout<<"Vec2d Clas:"<<endl;
  cout<<"x = "<<x<<endl;
  cout<<"y = "<<y<<endl;
}

