#include <fstream>
#include <random>
#include "Particles.h"

namespace HACCabana
{

Particles::Particles() 
{
  ;
}

Particles::~Particles() 
{
  ;
}

void Particles::convert_phys2grid(int ng, float rL, float a)
{
  auto velocity = Cabana::slice<Velocity>(aosoa_host, "velocity");

  const float phys2grid_pos = ng/rL;
  const float phys2grid_vel = phys2grid_pos/100.0;
  const float scaling = phys2grid_vel*a*a;
  for (int i=0; i<num_p; ++i)
  {
    for (int j=0; j<3; ++j) {
      velocity(i,j) *= scaling;
    }
  }
}

// 
void Particles::generateData(const int np, const float rl, const float mean_vel)
{
  num_p = np*np*np;
  aosoa_host = aosoa_host_type("aosoa_host", num_p);

  auto id = Cabana::slice<ParticleID>(aosoa_host, "id");
  auto position = Cabana::slice<Position>(aosoa_host, "position");

  const float delta = rl/np;

  // generate data from a grid and offset from a normal random distribution
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<float> d1{0.0, 0.05}; // mu=0.0 sigma=0.05

  for (int i=0; i<num_p; ++i)
  {
    id(i) = i;
    position(i,0) = (float)(i % np) * delta + d1(gen);
    position(i,1) = (float)(i / np % np) * delta + d1(gen);
    position(i,2) = (float)(i/(np*np)) * delta + d1(gen);
  }

  const float vel_1d = mean_vel/sqrt(3.0);

  std::normal_distribution<float> d2{0.0, 1.0};
  auto velocity = Cabana::slice<Velocity>(aosoa_host, "velocity");

  for (int i=0; i<num_p; ++i)
  {
    for (int j=0; j<3; ++j) {
      velocity(i,j) = vel_1d * d2(gen);
    }
  }

  this->begin = 0;
  this->end = num_p;
  reorder();
}

void Particles::readRawData(std::string file_name) 
{
  std::ifstream infile(file_name, std::ifstream::binary);

  // the first int has the number of particles
  infile.read((char*)&num_p, sizeof(int));

  aosoa_host = aosoa_host_type("aosoa_host", num_p);

  auto id = Cabana::slice<ParticleID>(aosoa_host, "id");
  auto position = Cabana::slice<Position>(aosoa_host, "position");
  auto velocity = Cabana::slice<Velocity>(aosoa_host, "velocity");

  for (int i=0; i<num_p; ++i)
    infile.read((char*)&id(i),sizeof(int64_t));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&position(i,0),sizeof(float));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&position(i,1),sizeof(float));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&position(i,2),sizeof(float));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&velocity(i,0),sizeof(float));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&velocity(i,1),sizeof(float));
  for (int i=0; i<num_p; ++i)
    infile.read((char*)&velocity(i,2),sizeof(float));

  infile.close();

  this->begin = 0;
  this->end = num_p;
  reorder();
}

void Particles::reorder()
{
  auto id = Cabana::slice<ParticleID>(aosoa_host, "id");
  auto position = Cabana::slice<Position>(aosoa_host, "position");
  auto velocity = Cabana::slice<Velocity>(aosoa_host, "velocity");

  // Relocate any particle outside of the boundary to the end of the 
  // aosoa -- outside particles start at this->end until the end of the aosoa.
  for (int i=this->begin; i<this->end; ++i) {
    if (position(i,0) < MIN_POS || position(i,1) < MIN_POS || position(i,2) < MIN_POS ||
        position(i,0) >= MAX_POS || position(i,1) >= MAX_POS || position(i,2) >= MAX_POS)
    {
      for (int j=0; j<3; ++j)
      {
        float tmp;
        tmp = position(i,j);
        position(i,j) = position(this->end-1,j);
        position(this->end-1,j) = tmp;
        tmp = velocity(i,j);
        velocity(i,j) = velocity(this->end-1,j);
        velocity(this->end-1,j) = tmp;
      }
      int64_t tmp2 = id(i);
      id(i) = id(this->end-1);
      id(this->end-1) = tmp2;

      --this->end;
      --i;
    }
  }
}

}
