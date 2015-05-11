#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array3.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"  
#include <vector>
typedef unsigned int uint;
using namespace std;
struct FLIP_particle
{
	Vec3f pos;
	Vec3f vel;
	FLIP_particle()
	{
		for (int i=0; i<3;i++)
		{
			pos[i] = 0;
			vel[i] = 0;
		}
		
	}
	~FLIP_particle()
	{
	}
	FLIP_particle(const FLIP_particle &p)
	{
		pos = p.pos;
		vel = p.vel;
	}
	FLIP_particle(const Vec3f & p,
	              const Vec3f & v)
	{
		pos = p;
		vel = v;
	}
};
struct emitter
{
	Vec3f vel;
	Vec3f plan;
	emitter()
	{
	}
	~emitter(){}
	emitter(const emitter &e)
	{
		vel = e.vel;
		plan = e.plan;
	}
	emitter(const Vec3f &_vel,
		    const Vec3f &_plan)
	{
		vel = _vel;
		plan = _plan;
	}
	
};
class FluidSim {

public:
   void initialize(float width, int ni_, int nj_, int nk_);
   void setEmitter(Vec3f &vel, Vec3f &plan){ simple_emitters.push_back(emitter(vel,plan)); }
   void set_boundary(float (*phi)(const Vec3f&));
   void set_moving_boundary(float (*phi)(const Vec3f&), Vec3f & vel_solid, bool b_moving);
   void set_liquid(float (*phi)(const Vec3f&));
   void add_particle(const Vec3f& pos);
   void set_forceStr(float str){_external_str = str;}

   void advance(float dt, int adv_type);

   //Grid dimensions
   int ni,nj,nk;
   float dx;
   
   //Fluid velocity
   vector<emitter> simple_emitters;
   Array3f u, v, w;
   Array3f u_solid, v_solid, w_solid;
   Array3f u_save, v_save, w_save;
   Array3f u_coef, v_coef, w_coef;
   Array3f temp_u, temp_v, temp_w;
   
   //Static geometry representation
   Array3f nodal_solid_phi;
   Array3f u_weights, v_weights, w_weights;
   Array3c u_valid, v_valid, w_valid;//, p_valid;
   Array3c marker_cell;

   std::vector<FLIP_particle> particles;
   std::vector<Vec3f> vort;
   float particle_radius;

   Array3f liquid_phi;

   //Array3f potential, temp_potential;

   //Data arrays for extrapolation
   Array3c valid, old_valid;

   //Solver data
   PCGSolver<double> solver;
   SparseMatrixd matrix;
   std::vector<double> rhs;
   std::vector<double> pressure;
   
   Vec3f get_velocity(const Vec3f& position);

   float _external_str;

private:
	void emitFluids(float dt);
	Vec3f back_trace_pos(float dt, Vec3f & pos, float _pcfl);

	 

	float H(float r)
	{
		float res = 0;
		if(r>=-1 && r<0) res = 1+r;
		if(r>=0 && r<1) res = 1-r;
		return res;
	}
	float compute_weight(Vec3f gp, Vec3f pp)
	{
		//k(x,y,z) = H(dx/hx)H(dy/hx)H(dz/hx)
		//H(r) = 1-r 0<=r<1  1+r -1<=r<0 0 else;
		Vec3f dd = gp - pp;
		return H(dd[0]/dx)*H(dd[1]/dx)*H(dd[2]/dx);
	}
	void divide_weight(Array3f &u, Array3f & u_coef);
	Vec3f get_dvelocity(const Vec3f & position);
	void compute_delta(Array3f & u, Array3f &u_old, Array3f &u_temp);
	void particle_interpolate(float alpha);
	void FLIP_advection(float dt); 

   Vec3f trace_rk2(const Vec3f& position, float dt);

   float cfl();

   void advect_particles(float dt);
   void advect(float dt);
   void add_force(float dt);
   void add_external_force(float dt, float str);
   void project(float dt);
   void constrain_velocity();

   //helpers for pressure projection
   void compute_weights();
   void solve_pressure(float dt);
   void compute_phi(vector<vector<int>> & particle_hash, vector<Vec3i> cell_list);
   //void computeGradPhi(Array3f & u_temp, int dir);
   //void advectPotential(float dt);



};


#endif