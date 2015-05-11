#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "array3_utils.h"
#include "faithful_surfacing3.h"
#include "SDFGen.h"
using namespace std;

//int grid_resolution = 128;
int gridX = 192, gridY = 192, gridZ = 96;
float timestep = 0.01f;
int frame = 0;

float grid_width = 1;

FluidSim sim;
SDFGenerator g_sdfSolid;
SDFGenerator g_water;

using namespace std;
float box_phi(const Vec3f& position, const Vec3f & centre, Vec3f & b)
{
	//vec3 d = abs(p) - b;
	//return min(max(d.x,max(d.y,d.z)),0.0) +
	//	length(max(d,0.0));


	Vec3f p = position - centre;
	Vec3f d = fabs(p) - b;
	return min(max(d[0], max(d[1],d[2])),0.0f) 
		+ dist(Vec3f(max(d[0],0.0f),max(d[1],0.0f),max(d[2],0.0f)),Vec3f(0,0,0));
}

float sphere_phi(const Vec3f& position, const Vec3f& centre, float radius) {
   return (dist(position,centre) - radius);
}

Vec3f ball_pos(0.2,0.175,0.5);
Vec3f ball_vel(0.8,0.0,0.0);
Vec3f c0(2.0f,0.5f,0.5f);
float rad0 = 0.35f;
float ball_radius = 0.05;
float moving_phi(const Vec3f& position)
{
	return box_phi(position,ball_pos,Vec3f(0.02,0.07,0.06));
}
float moving_phi2(const Vec3f& position)
{
	return box_phi(position,ball_pos,Vec3f(0.00,0.0,0.0));
}
float boundary_phi(const Vec3f& position) {
    //return -sphere_phi(position, c0, rad0);
	float d1 = box_phi(position, Vec3f(0.5,0.5,0.25),Vec3f(0.45,0.45,0.225));
	//float d2 = g_sdfSolid.getPhi(position);
	//float d = max(-d2,d1);
	return -d1;
	
}

float liquid_phi(const Vec3f& position) {
    //return sphere_phi(position, Vec3f(0.55f, 0.55f, 0.4f), 0.23f);
	//float d1 = box_phi(position,ball_pos,Vec3f(0.02,0.07,0.06));
	float d1 = box_phi(position,Vec3f(0.2,0.3,0.25),Vec3f(0.16,0.26,0.225));
	//float d2 = box_phi(position,Vec3f(130,0.2,0.25),Vec3f(120,0.16,0.23));
	//return max(-d1,d2);
	return d1;
	//return 0;//g_water.getPhi(position);
}

void export_particles(string path, int frame, const std::vector<FLIP_particle>& particles, float radius);

//Main testing code
int advection_type=0;
//-------------
int main(int argc, char **argv)
{
   //g_sdfSolid.dx = 1.25*grid_width/(float)gridX;
   //g_water.dx = g_sdfSolid.dx;
   //g_sdfSolid.objToSDF("C:/Users/xinxin/Desktop/testScene2.obj",Vec3f(0,0,0), Vec3f(2,0.5,2));
   //g_water.objToSDF("C:/Users/xinxin/Desktop/water2.obj",Vec3f(0,0,0), Vec3f(2,1,1));
   if(argc!=3){
      cerr << "The first parameter should be the folder to write the output liquid meshes into. (eg. c:\\output\\)"<<"advection_method" << endl;
      return 1;
   }

   string outpath(argv[1]);
   sscanf(argv[2],"%d", &advection_type);
   
   printf("Initializing data\n");
   sim.initialize(grid_width, gridX, gridY, gridZ);
   
   printf("Initializing boundary\n");
   sim.set_boundary(boundary_phi);
   
   printf("Initializing liquid\n");
   sim.set_liquid(liquid_phi);
      
   printf("Exporting initial data\n");
   //export_particles(outpath, 0, sim.particles, sim.particle_radius);
   //sim.setEmitter(Vec3f(0.12,0,0), Vec3f(0.14,0.27,0));
   for(frame = 1; frame <=300; ++frame) {
      printf("--------------------\nFrame %d\n", frame);
      
      //Simulate
      printf("Simulating liquid\n");
	  //if(frame<21)
	  //{

	  ////ball_vel = Vec3f(0,0.0,0);
		 // for(int i=0;i<10;i++){

			//  sim.set_boundary(boundary_phi);
			//  sim.set_forceStr(0.0);
			//  sim.set_moving_boundary(moving_phi, ball_vel,true);
			//  sim.advance(0.002, advection_type);

			//  

			//  ball_pos = ball_pos + 0.002f * ball_vel;
			//	  // // + 0.5f * (float)(timestep * timestep)*Vec3f(0,-9.81,0);
			//	  // //ball_vel = ball_vel + Vec3f(0, -9.81f * timestep ,0 ) ;

			//  
		 //}
	  //}
	  //else if (frame<61)
	  //{
		 // ball_vel = Vec3f(0,0.2,0);
		 // sim.set_boundary(boundary_phi);
		 // sim.set_moving_boundary(moving_phi, ball_vel,false);
		 // sim.set_forceStr(0);
		 // sim.advance(0.02, advection_type);
		 // ball_pos = ball_pos + 0.02f * ball_vel;
	  //}
	  //else
	  
	  {
		  //ball_vel = Vec3f(0,0.0,0);
		  sim.set_boundary(boundary_phi);
		  //sim.set_moving_boundary(moving_phi, ball_vel,false);
		  sim.advance(0.02, advection_type);
	  }
	  
      
      printf("Exporting particle data\n");
      export_particles(outpath, frame-1, sim.particles, sim.particle_radius);
   }

   return 0;
}


void export_particles(string path, int frame, const std::vector<FLIP_particle>& particles_in, float radius) {
   //Write the output
   
   //std::stringstream strout;
   //strout << path << "particles_" << frame << ".txt";
   //string filepath = strout.str();
   //
   //ofstream outfile(filepath.c_str());
   ////write vertex count and particle radius
   //outfile << particles.size() << " " << radius << std::endl;
   ////write vertices
   //for(unsigned int i = 0; i < particles.size(); ++i)
   //   outfile << particles[i].pos[0] << " " << particles[i].pos[1] << " " << particles[i].pos[2] << std::endl;
   //outfile.close();

	FaithfulSurfacing3 g_mesher;
	vector<Vec3f> particles;
	particles.resize(particles_in.size());
	tbb::parallel_for(0,(int)particles.size(),1,[&](int p)
	{
		particles[p] = particles_in[p].pos;
	});
	

   g_mesher.grid_dx = radius/1.001*2.0/sqrt(3.0);
   g_mesher.particle_x = particles;
   g_mesher.inner_radius = 0.95*radius;
   g_mesher.outer_radius = 1.05*radius;

   g_mesher.shrink_steps = 0;
   g_mesher.smooth_steps = 0;
   g_mesher.free_smooth_steps = 4;
   g_mesher.run_surfacing();


   ostringstream strout;
   strout << path << "liquidmesh_" << frame << ".obj";

   string filepath = strout.str();

   ofstream outfile(filepath.c_str());

   //write vertices
   for(unsigned int i = 0; i < g_mesher.x.size(); ++i)
	   outfile<<"v"<<" "<< g_mesher.x[i].v[0] << " " << g_mesher.x[i].v[1] << " " << g_mesher.x[i].v[2] << std::endl;
   for(unsigned int i = 0; i < g_mesher.tri.size(); ++i)
	   outfile<<"f"<<" "<< g_mesher.tri[i].v[0]+ 1<< " " << g_mesher.tri[i].v[1] + 1 << " " << g_mesher.tri[i].v[2]+ 1<< std::endl;
   outfile.close();


   particles.resize(0);
   particles.shrink_to_fit();

   



}


