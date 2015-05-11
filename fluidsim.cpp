#include "fluidsim.h"

#include "array3_utils.h"
#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "tbb/tbb.h"
#include "AlgebraicMultigrid.h"
void extrapolate(Array3f& grid, Array3c& valid);

void FluidSim::initialize(float width, int ni_, int nj_, int nk_) {
   ni = ni_;
   nj = nj_;
   nk = nk_;
   dx = width / (float)ni;
   u.resize(ni+1,nj,nk); temp_u.resize(ni+1,nj,nk); u_weights.resize(ni+1,nj,nk); u_valid.resize(ni+1,nj,nk);
   u_save.resize(ni+1,nj,nk); u_coef.resize(ni+1,nj,nk);
   v.resize(ni,nj+1,nk); temp_v.resize(ni,nj+1,nk); v_weights.resize(ni,nj+1,nk); v_valid.resize(ni,nj+1,nk);
   v_save.resize(ni,nj+1,nk); v_coef.resize(ni,nj+1,nk);
   w.resize(ni,nj,nk+1); temp_w.resize(ni,nj,nk+1); w_weights.resize(ni,nj,nk+1); w_valid.resize(ni,nj,nk+1);
   w_save.resize(ni,nj,nk+1); w_coef.resize(ni,nj,nk+1);

   particle_radius = (float)(dx * 1.01*sqrt(3.0)/2.0); 
   //make the particles large enough so they always appear on the grid

   u.set_zero();
   v.set_zero();
   w.set_zero();
   nodal_solid_phi.resize(ni+1,nj+1,nk+1);
   valid.resize(ni+1, nj+1, nk+1);
   old_valid.resize(ni+1, nj+1, nk+1);
   liquid_phi.resize(ni,nj,nk);
   marker_cell.resize(ni,ni,nk);


 


}



//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(const Vec3f&)) {

   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec3f pos(i*dx,j*dx,k*dx);
      nodal_solid_phi(i,j,k) = phi(pos);
   }

}

void FluidSim::set_liquid(float (*phi)(const Vec3f&)) {
   //surface.reset_phi(phi, dx, Vec3f(0.5f*dx,0.5f*dx,0.5f*dx), ni, nj, nk);
   
   //initialize particles
   int seed = 0;
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
	   for(int sub_p=0; sub_p<6; sub_p++)
	   {

		   Vec3f pos(i*dx,j*dx,k*dx);
		   float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
		   pos += dx * Vec3f(a,b,c);

		   if(phi(pos) <= -particle_radius) {
			   float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
			   if(solid_phi >= 0)
				   particles.push_back(FLIP_particle(pos,Vec3f(0,0,0)));
		   }
       }
   }
}
void FluidSim::emitFluids(float dt)
{
	int seed = 0;
	for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
		for(int sub_p=0; sub_p<6; sub_p++)
		{

			Vec3f pos(i*dx,j*dx,k*dx);
			float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
			pos += dx * Vec3f(a,b,c);

			for (int e=0;e<simple_emitters.size();e++)
			{
				if(interpolate_value(pos/dx,liquid_phi)>=0)
				{
					Vec3f back_pos = pos - dt*simple_emitters[e].vel;
					if (back_pos[0]<simple_emitters[e].plan[0] 
					    && back_pos[0]>0.11 && back_pos[1]<simple_emitters[e].plan[1] 
						&&back_pos[1]>0.13)
					{
						float solid_phi = interpolate_value(pos/dx, nodal_solid_phi);
						if(solid_phi >= 0)
							particles.push_back(FLIP_particle(pos,simple_emitters[e].vel));
					}
					
				}
			}
			
		}
	}
	printf("num %d particles\n",particles.size());
}
//The main fluid simulation step
void FluidSim::advance(float dt, int adv_type) {
   float t = 0;

   while(t < dt) {
      float substep = 10*cfl();   
      if(t + substep > dt)
         substep = dt - t;
      printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t+substep)/dt);
      
      //printf(" Surface (particle) advection\n");
      //advect_particles(substep);

      //printf(" Velocity advection\n");
      //Advance the velocity
      //advect(substep);
	  printf("emission\n");
	  emitFluids(substep);

	   
		  printf("FLIP advection\n");
		  FLIP_advection(substep);
	   


      add_force(substep);

      printf(" Pressure projection\n");
      project(substep); 
       
      //Pressure projection only produces valid velocities in faces with non-zero associated face area.
      //Because the advection step may interpolate from these invalid faces, 
      //we must extrapolate velocities from the fluid domain into these invalid faces.
      printf(" Extrapolation\n");
      extrapolate(u, u_valid);
      extrapolate(v, v_valid);
      extrapolate(w, w_valid);
    
      //For extrapolated velocities, replace the normal component with
      //that of the object.
      printf(" Constrain boundary velocities\n");
      constrain_velocity();

      t+=substep;
   }
}


float FluidSim::cfl() {

	float maxvel = 0;
	for(unsigned int i = 0; i < u.a.size(); ++i)
		maxvel = max(maxvel, fabs(u.a[i]));
	for(unsigned int i = 0; i < v.a.size(); ++i)
		maxvel = max(maxvel, fabs(v.a[i]));
	for(unsigned int i = 0; i < w.a.size(); ++i)
		maxvel = max(maxvel, fabs(w.a[i]));
	float dt = dx / maxvel;

	if (dt<0.0005)
	{
		dt = 0.0005;
	}

	return dt;
}

void FluidSim::add_particle(const Vec3f& pos) {
   particles.push_back(FLIP_particle(pos,Vec3f(0,0,0)));
}

void FluidSim::add_force(float dt) { //needs tbb

   //gravity
	
   //for(int k = 0;k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
   //   v(i,j,k) -= 9.81f * dt;
   //}

	//int compute_num = nk*(nj)*ni;
	//int slice = ni*(nj);
	//tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	//{
	//	int k = thread_idx/slice;
	//	int j = (thread_idx%slice)/ni;
	//	int i = thread_idx%ni;

	//	for (int e=0;e<simple_emitters.size();e++)
	//	{
	//		Vec3f pos = Vec3f((float)i+0.5, (float)j+0.5,(float)k+0.5)*dx;
	//		if (pos[0]>=simple_emitters[e].plan[0] 
	//		&& pos[0]<simple_emitters[e].plan[0]+dt*simple_emitters[e].vel[0]
	//		  &&pos[1]>0.13  && pos[1]<simple_emitters[e].plan[1])
	//		{
	//			if(liquid_phi(i,j,k)<0){

	//				Vec3f emit_vel = simple_emitters[e].vel;
	//				u(i,j,k) = emit_vel[0];
	//				u(i+1,j,k) = emit_vel[0];
	//				v(i,j,k) = emit_vel[1];
	//				v(i,j+1,k) = emit_vel[1];
	//				w(i,j,k) = emit_vel[2];
	//				w(i,j,k+1) = emit_vel[2];

	//			}		
	//		  }
	//		
	//	}
	//});


   int compute_num = nk*(nj+1)*ni;
   int slice = ni*(nj+1);
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/ni;
	   int i = thread_idx%ni;

	   v(i,j,k) -= 9.81f * dt;
   });

}



//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::constrain_velocity() {
   temp_u = u;
   temp_v = v;
   temp_w = w;

   //(At lower grid resolutions, the normal estimate from the signed
   //distance function can be poor, so it doesn't work quite as well.
   //An exact normal would do better if we had it for the geometry.)

   
   //constrain u
   //needs tbb
   //for(int k = 0; k < u.nk;++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
   //   if(u_weights(i,j,k) == 0) {
   //      //apply constraint
   //      Vec3f pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
   //      Vec3f vel = get_velocity(pos);
   //      Vec3f normal(0,0,0);
   //      interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
   //      normalize(normal);
   //      float perp_component = dot(vel, normal);
   //      vel -= perp_component*normal;
   //      temp_u(i,j,k) = vel[0];
   //   }
   //}
   int compute_num; int slice;
   compute_num = u.nk*u.nj*u.ni;
   slice = u.ni*u.nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/u.ni;
	   int i = thread_idx%u.ni;
	   if(u_weights(i,j,k) == 0) {
		   //apply constraint
		   Vec3f pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
		   Vec3f vel = get_velocity(pos);
		   Vec3f normal(0,0,0);
		   interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
		   normalize(normal);
		   float perp_component = dot(vel, normal);
		   vel -= perp_component*normal;
		   temp_u(i,j,k) = vel[0];
	   }
   });

   //constrain v
   //needs tbb
   //for(int k = 0; k < v.nk;++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
   //   if(v_weights(i,j,k) == 0) {
   //      //apply constraint
   //      Vec3f pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
   //      Vec3f vel = get_velocity(pos);
   //      Vec3f normal(0,0,0);
   //      interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
   //      normalize(normal);
   //      float perp_component = dot(vel, normal);
   //      vel -= perp_component*normal;
   //      temp_v(i,j,k) = vel[1];
   //   }
   //}
   compute_num = v.nk*v.nj*v.ni;
   slice = v.ni*v.nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/v.ni;
	   int i = thread_idx%v.ni;

	   if(v_weights(i,j,k) == 0) {
		   //apply constraint
		   Vec3f pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
		   Vec3f vel = get_velocity(pos);
		   Vec3f normal(0,0,0);
		   interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
		   normalize(normal);
		   float perp_component = dot(vel, normal);
		   vel -= perp_component*normal;
		   temp_v(i,j,k) = vel[1];
	   }
   });

   //constrain w
   //needs tbb
   //for(int k = 0; k < w.nk;++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
   //   if(w_weights(i,j,k) == 0) {
   //      //apply constraint
   //      Vec3f pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
   //      Vec3f vel = get_velocity(pos);
   //      Vec3f normal(0,0,0);
   //      interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
   //      normalize(normal);
   //      float perp_component = dot(vel, normal);
   //      vel -= perp_component*normal;
   //      temp_w(i,j,k) = vel[2];
   //   }
   //}
   compute_num = w.nk*w.nj*w.ni;
   slice = w.ni*w.nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/w.ni;
	   int i = thread_idx%w.ni;
	   if(w_weights(i,j,k) == 0) {
		   //apply constraint
		   Vec3f pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
		   Vec3f vel = get_velocity(pos);
		   Vec3f normal(0,0,0);
		   interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
		   normalize(normal);
		   float perp_component = dot(vel, normal);
		   vel -= perp_component*normal;
		   temp_w(i,j,k) = vel[2];
	   }
   });

   //update
   u = temp_u;
   v = temp_v;
   w = temp_w;

}

void FluidSim::advect_particles(float dt) { 
	//needs tbb
   //for(unsigned int p = 0; p < particles.size(); ++p) {
   //   particles[p] = trace_rk2(particles[p], dt);
   //
   //   //check boundaries and project exterior particles back in
   //   float phi_val = interpolate_value(particles[p]/dx, nodal_solid_phi); 
   //   if(phi_val < 0) {
   //      Vec3f grad;
   //      interpolate_gradient(grad, particles[p]/dx, nodal_solid_phi);
   //      if(mag(grad) > 0)
   //         normalize(grad);
   //      particles[p] -= phi_val * grad;
   //   }
   //}

   int num = particles.size();
   tbb::parallel_for(0,num, 1, [&](int p){
	   particles[p].pos = trace_rk2(particles[p].pos, dt);

	   //check boundaries and project exterior particles back in
	   float phi_val = interpolate_value(particles[p].pos/dx, nodal_solid_phi); 
	   if(phi_val < 0) {
		   Vec3f grad;
		   interpolate_gradient(grad, particles[p].pos/dx, nodal_solid_phi);
		   if(mag(grad) > 0)
			   normalize(grad);
		   particles[p].pos -= phi_val * grad;
	   }

   });
   

}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::advect(float dt) {

   temp_u.assign(0);
   temp_v.assign(0);
   temp_w.assign(0);

   //semi-Lagrangian advection on u-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec3f pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
      pos = trace_rk2(pos, -dt);
      temp_u(i,j,k) = get_velocity(pos)[0];  
   }

   //semi-Lagrangian advection on v-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      Vec3f pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
      pos = trace_rk2(pos, -dt);
      temp_v(i,j,k) = get_velocity(pos)[1];
   }

   //semi-Lagrangian advection on w-component of velocity
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vec3f pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
      pos = trace_rk2(pos, -dt);
      temp_w(i,j,k) = get_velocity(pos)[2];
   }

   //move update velocities into u/v vectors
   u = temp_u;
   v = temp_v;
   w = temp_w;
}

void FluidSim::compute_phi(vector<vector<int>> & particle_hash, vector<Vec3i> cell_list) {
   
   //grab from particles 
	//keep for now, parallel later
   liquid_phi.assign(3*dx);
   marker_cell.assign(0);

   int compute_num = cell_list.size();
   int slice = ni*nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = cell_list[thread_idx].v[2];
	   int j = cell_list[thread_idx].v[1];
	   int i = cell_list[thread_idx].v[0];
	   if( i<ni && j<nj && k<nk)
	   {
		   
			Vec3f sample_pos((i+0.5f)*dx, (j+0.5f)*dx,(k+0.5f)*dx);
			for(int kk = max(0,k-1); kk<=min(k+1, nk-1);kk++)
			for(int jj = max(0,j-1); jj<=min(j+1, nj-1);jj++)
			for(int ii = max(0,i-1); ii<=min(i+1, ni-1);ii++)
			{
				int index = kk*ni*nj + jj*ni + ii;
				for (int p=0; p<particle_hash[index].size();p++)
				{
					Vec3f pos = particles[particle_hash[index][p]].pos;
					float test_val = dist(sample_pos, pos) - particle_radius;
					if(test_val < liquid_phi(i,j,k))
						liquid_phi(i,j,k) = test_val;
				}

			}
		   
	   }
   });


   //for(unsigned int p = 0; p < particles.size(); ++p) {
   //   Vec3i cell_ind(particles[p].pos / dx);
	  //marker_cell(cell_ind[0], cell_ind[1],cell_ind[2]) = 1;
   //   for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,nk-1); ++k) {
   //      for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,nj-1); ++j) {
   //         for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,ni-1); ++i) {
   //            Vec3f sample_pos((i+0.5f)*dx, (j+0.5f)*dx,(k+0.5f)*dx);
   //            float test_val = dist(sample_pos, particles[p].pos) - particle_radius;
   //            if(test_val < liquid_phi(i,j,k))
   //               liquid_phi(i,j,k) = test_val;
   //         }
   //      }
   //   }
   //}
   Array3c marker_temp = marker_cell;
   //eat small bubble;
   compute_num = ni*nj*nk;
   slice = ni*nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/ni;
	   int i = thread_idx%ni;
	   if(i>0&&i<ni-1&&j>0&&j<nj-1&&k>0&&k<nk-1)
	   {
		   if(marker_cell(i,j,k)==0)
		   {
			   if(marker_cell(i-1,j,k)&&marker_cell(i+1,j,k)
				   &&marker_cell(i,j-1,k)&&marker_cell(i,j+1,k)
				   &&marker_cell(i,j,k-1)&&marker_cell(i,j,k+1)
				   )
			   { marker_temp(i,j,k) = 1; }
		   }
	   }
   });
   marker_cell = marker_temp;


   //extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
   //needs parallel
   Array3f phi_temp = liquid_phi;
   //for(int k = 0; k < nk; ++k) {
   //   for(int j = 0; j < nj; ++j) {
   //      for(int i = 0; i < ni; ++i) {
   //         if(liquid_phi(i,j,k) < 0.5*dx) {
   //            float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
   //               + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
   //            if(solid_phi_val < 0)
   //               phi_temp(i,j,k) = -0.5f*dx;
   //         }
   //      }
   //   }
   //}
   compute_num = ni*nj*nk;
   slice = ni*nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/ni;
	   int i = thread_idx%ni;
	   if(liquid_phi(i,j,k) < 0.5*dx) {
		   float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
			   + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
		   if(solid_phi_val < 0)
			   phi_temp(i,j,k) = -0.5f*dx;
	   }

   });
   liquid_phi = phi_temp;
   

   //phi_temp = liquid_phi;
   //compute_num = ni*nj*nk;
   //slice = ni*nj;
   //tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   //{
	  // int k = thread_idx/slice;
	  // int j = (thread_idx%slice)/ni;
	  // int i = thread_idx%ni;
	  // if(liquid_phi(i,j,k) < 0.5*dx) {
		 //  float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
			//   + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
		 //  if(solid_phi_val > 0){

			//   if(marker_cell(i,j,k) == 0)
			//   {
			//	   liquid_phi(i,j,k) += 0.5*dx;
			//   }
		 //  }
	  // }

   //});


   //compute_num = ni*nj*nk;
   //slice = ni*nj;
   //tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   //{
	  // int k = thread_idx/slice;
	  // int j = (thread_idx%slice)/ni;
	  // int i = thread_idx%ni;
	  // float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
		 //  + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
	  // if(liquid_phi(i,j,k) < 0 || solid_phi_val <0 ) {
			//	marker_cell(i,j,k) = 1;
	  // }

   //});
   

}
 
void FluidSim::project(float dt) {

   
   
   //Compute finite-volume type face area weight for each velocity sample.
   compute_weights();

   //Set up and solve the variational pressure solve.
   solve_pressure(dt);
   
}


//Apply RK2 to advect a point in the domain.
Vec3f FluidSim::trace_rk2(const Vec3f& position, float dt) {
	float c1 = 2.0/9.0*dt, c2 = 3.0/9.0 * dt, c3 = 4.0/9.0 * dt;
   Vec3f input = position;
   Vec3f velocity1 = get_velocity(input);
   Vec3f midp1 = input + ((float)(0.5*dt))*velocity1;
   Vec3f velocity2 = get_velocity(midp1);
   Vec3f midp2 = input + ((float)(0.75*dt))*velocity2;
   Vec3f velocity3 = get_velocity(midp2);
   //velocity = get_velocity(input + 0.5f*dt*velocity);
   //input += dt*velocity;
   input = input + c1*velocity1 + c2*velocity2 + c3*velocity3;
   return input;
}

//Interpolate velocity from the MAC grid.
Vec3f FluidSim::get_velocity(const Vec3f& position) {

   //Interpolate the velocity from the u and v grids
   float u_value = interpolate_value(position / dx - Vec3f(0, 0.5f, 0.5f), u);
   float v_value = interpolate_value(position / dx - Vec3f(0.5f, 0, 0.5f), v);
   float w_value = interpolate_value(position / dx - Vec3f(0.5f, 0.5f, 0), w);

   return Vec3f(u_value, v_value, w_value);
}

Vec3f FluidSim::get_dvelocity(const Vec3f& position) {

	//Interpolate the velocity from the u and v grids
	float u_value = interpolate_value(position / dx - Vec3f(0, 0.5f, 0.5f), temp_u);
	float v_value = interpolate_value(position / dx - Vec3f(0.5f, 0, 0.5f), temp_v);
	float w_value = interpolate_value(position / dx - Vec3f(0.5f, 0.5f, 0), temp_w);

	return Vec3f(u_value, v_value, w_value);
}

//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::compute_weights() {

   //Compute face area fractions (using marching squares cases).
	int compute_num; int slice;
	//parallel
   //for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
   //   u_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,j,  k),
   //                                          nodal_solid_phi(i,j+1,k),
   //                                          nodal_solid_phi(i,j,  k+1),
   //                                          nodal_solid_phi(i,j+1,k+1));
   //   u_weights(i,j,k) = clamp(u_weights(i,j,k),0.0f,1.0f);
   //}
	compute_num = u.ni*u.nj*u.nk;
	slice = u.ni*u.nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/u.ni;
		int i = thread_idx%u.ni;
		u_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,j,  k),
			nodal_solid_phi(i,j+1,k),
			nodal_solid_phi(i,j,  k+1),
			nodal_solid_phi(i,j+1,k+1));
		u_weights(i,j,k) = clamp(u_weights(i,j,k),0.0f,1.0f);
	});
   //parallel
   //for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
   //   v_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,k),
   //                                          nodal_solid_phi(i,  j,k+1),
   //                                          nodal_solid_phi(i+1,j,k),
   //                                          nodal_solid_phi(i+1,j,k+1));
   //   v_weights(i,j,k) = clamp(v_weights(i,j,k),0.0f,1.0f);
   //}
	compute_num = v.ni*v.nj*v.nk;
	slice = v.ni*v.nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/v.ni;
		int i = thread_idx%v.ni;
		v_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,k),
			nodal_solid_phi(i,  j,k+1),
			nodal_solid_phi(i+1,j,k),
			nodal_solid_phi(i+1,j,k+1));
		v_weights(i,j,k) = clamp(v_weights(i,j,k),0.0f,1.0f);
	});
   //parallel
   //for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
   //   w_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,  k),
   //                                          nodal_solid_phi(i,  j+1,k),
   //                                          nodal_solid_phi(i+1,j,  k),
   //                                          nodal_solid_phi(i+1,j+1,k));
   //   w_weights(i,j,k) = clamp(w_weights(i,j,k),0.0f,1.0f);
   //}
	compute_num = w.ni*w.nj*w.nk;
	slice = w.ni*w.nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/w.ni;
		int i = thread_idx%w.ni;
		w_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,  k),
			nodal_solid_phi(i,  j+1,k),
			nodal_solid_phi(i+1,j,  k),
			nodal_solid_phi(i+1,j+1,k));
		w_weights(i,j,k) = clamp(w_weights(i,j,k),0.0f,1.0f);
	});


}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim::solve_pressure(float dt) {


   int ni = v.ni;
   int nj = u.nj;
   int nk = u.nk;
   int compute_num = ni*nj*nk;
   int slice = ni*nj;
   //vector<int> index_table;
   //vector<char> mask;
   //mask.resize(ni*nj*nk);
   //mask.assign(mask.size(),0);
   //index_table.resize(ni*nj*nk);
   //index_table.assign(index_table.size(),0);
   //tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   //{
	  // int k = thread_idx/slice;
	  // int j = (thread_idx%slice)/ni;
	  // int i = thread_idx%ni;
	  // if(i>=1 && i<ni-1 && j>=1 && j<nj-1 && k>=1 && k<nk-1)
	  // {
		 //  int index = i + ni*j + ni*nj*k;

		 //  float centre_phi = liquid_phi(i,j,k);
		 //  if(centre_phi < 0) 
		 //  {
			//   mask[index]=1;
		 //  }
	  // }
   //});
   //index_table[0] = mask[0] - 1;
   //for (int i = 1;i<index_table.size();i++)
   //{
	  // index_table[i]=index_table[i-1]+mask[i];
   //}
   //

   //int system_size = index_table[index_table.size()-1] + 1;
   pressure.resize(ni*nj*nk);
   pressure.assign(ni*nj*nk,0);
   vector<double> x;
   vector<Vec3i> dof_ijk;
   x.resize(0);
   dof_ijk.resize(0);


   Array3ui dof_index;
   dof_index.resize(ni,nj,nk);
   dof_index.assign(0);

   for (int k=0;k<nk;k++)for(int j=0;j<nj;j++)for(int i=0;i<ni;i++)
   {
	   if (liquid_phi(i,j,k)<0)
	   {
		   dof_index(i,j,k) = x.size();
		   dof_ijk.push_back(Vec3i(i,j,k));
		   x.push_back(0);

	   }
	   
   }
   int system_size = x.size();
   


   //if(rhs.size() != system_size) {
      rhs.resize(system_size);
      //x.resize(system_size);
      matrix.resize(system_size);
   //}
   
   matrix.zero();
   rhs.assign(rhs.size(), 0);
   x.assign(x.size(), 0);

   
   compute_num = ni*nj*nk;
   slice = ni*nj;
   tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/ni;
	   int i = thread_idx%ni;
	   if(i>=1 && i<ni-1 && j>=1 && j<nj-1 && k>=1 && k<nk-1)
	   {
		   int index = i + ni*j + ni*nj*k;

		   //rhs[index_table[index]] = 0;
		   //pressure[index_table[index]] = 0;
		   float centre_phi = liquid_phi(i,j,k);
		   uint dof_idx = dof_index(i,j,k);
		   if(centre_phi < 0) 
		   {
			   rhs[dof_idx] = 0;
			   //right neighbour
			   float term = u_weights(i+1,j,k) * dt / sqr(dx);
			   float right_phi = liquid_phi(i+1,j,k);
			   if(right_phi < 0) {
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i+1,j,k), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, right_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term/theta);
			   }
			   rhs[dof_idx] -= u_weights(i+1,j,k)*u(i+1,j,k) / dx;

			   //left neighbour
			   term = u_weights(i,j,k) * dt / sqr(dx);
			   float left_phi = liquid_phi(i-1,j,k);
			   if(left_phi < 0) {
				   matrix.add_to_element(dof_idx,
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i-1,j,k), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, left_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term/theta);
			   }
			   rhs[dof_idx] += u_weights(i,j,k)*u(i,j,k) / dx;

			   //top neighbour
			   term = v_weights(i,j+1,k) * dt / sqr(dx);
			   float top_phi = liquid_phi(i,j+1,k);
			   if(top_phi < 0) {
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i,j+1,k), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, top_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term/theta);
			   }
			   rhs[dof_idx] -= v_weights(i,j+1,k)*v(i,j+1,k) / dx;

			   //bottom neighbour
			   term = v_weights(i,j,k) * dt / sqr(dx);
			   float bot_phi = liquid_phi(i,j-1,k);
			   if(bot_phi < 0) {
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i,j-1,k), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, bot_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					  dof_idx, term/theta);
			   }
			   rhs[dof_idx] += v_weights(i,j,k)*v(i,j,k) / dx;


			   //far neighbour
			   term = w_weights(i,j,k+1) * dt / sqr(dx);
			   float far_phi = liquid_phi(i,j,k+1);
			   if(far_phi < 0) {
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i,j,k+1), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, far_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term/theta);
			   }
			   rhs[dof_idx] -= w_weights(i,j,k+1)*w(i,j,k+1) / dx;

			   //near neighbour
			   term = w_weights(i,j,k) * dt / sqr(dx);
			   float near_phi = liquid_phi(i,j,k-1);
			   if(near_phi < 0) {
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term);
				   matrix.add_to_element(dof_idx, 
					   dof_index(i,j,k-1), -term);
			   }
			   else {
				   float theta = fraction_inside(centre_phi, near_phi);
				   if(theta < 0.01f) theta = 0.01f;
				   matrix.add_to_element(dof_idx, 
					   dof_idx, term/theta);
			   }
			   rhs[dof_idx] += w_weights(i,j,k)*w(i,j,k) / dx;
		   }
	   }
   });

   //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

    double tolerance;
    int iterations;
    //solver.set_solver_parameters(1e-18, 1000);
    //bool success = solver.solve(matrix, rhs, x, tolerance, iterations);
    //printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
    //bool success = AMGPCGSolveCompressed(matrix,rhs,x,mask,index_table, 1e-15,50,tolerance,iterations,ni,nj,nk);
	bool success = AMGPCGSolveSparse(matrix,rhs,x,dof_ijk,1e-6,50,tolerance,iterations,ni,nj,nk);
    printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
    if(!success) {
        printf("WARNING: Pressure solve failed!************************************************\n");
    }
	tbb::parallel_for(0,ni*nj*nk,1,[&](int idx)
	{
		if(liquid_phi.a[idx]<0)
			pressure[idx] = x[dof_index.a[idx]];
	});
   //Apply the velocity update
   u_valid.assign(0);
   //needs parallel
   //for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i) {
   //   int index = i + j*ni + k*ni*nj;
   //   if(u_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i-1,j,k) < 0)) {
   //      float theta = 1;
   //      if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
   //         theta = fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
   //      if(theta < 0.01f) theta = 0.01f;
   //      u(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-1]) / dx / theta; 
   //      u_valid(i,j,k) = 1;
   //   }
   //}
   compute_num = u.ni*u.nj*u.nk;
   slice = u.ni*u.nj;
   tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/u.ni;
	   int i = thread_idx%u.ni;
	   if(k<u.nk && j<u.nj && i<u.ni-1 && i>0)
	   {
		   int index = i + j*ni + k*ni*nj;
		   if(u_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i-1,j,k) < 0)) {
			   float theta = 1;
			   if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
				   theta = fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
			   if(theta < 0.01f) theta = 0.01f;
			   u(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-1]) / dx / theta; 
			   u_valid(i,j,k) = 1;
		   }

	   }
   });
   
   v_valid.assign(0);
   //needs parallel
   //for(int k = 0; k < v.nk; ++k) for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i) {
   //   int index = i + j*ni + k*ni*nj;
   //   if(v_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j-1,k) < 0)) {
   //      float theta = 1;
   //      if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
   //         theta = fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
   //      if(theta < 0.01f) theta = 0.01f;
   //      v(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni]) / dx / theta; 
   //      v_valid(i,j,k) = 1;
   //   }
   //}
   compute_num = v.ni*v.nj*v.nk;
   slice = v.ni*v.nj;
   tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/v.ni;
	   int i = thread_idx%v.ni;
	   if(k<v.nk && j>0 && j<v.nj-1 && i<v.ni)
	   {
		   int index = i + j*ni + k*ni*nj;
		   if(v_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j-1,k) < 0)) {
			   float theta = 1;
			   if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
				   theta = fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
			   if(theta < 0.01f) theta = 0.01f;
			   v(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni]) / dx / theta; 
			   v_valid(i,j,k) = 1;
		   }

	   }
   });

   w_valid.assign(0);
   //needs parallel
   //for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 1; i < w.ni-1; ++i) {
   //   int index = i + j*ni + k*ni*nj;
   //   if(w_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j,k-1) < 0)) {
   //      float theta = 1;
   //      if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
   //         theta = fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
   //      if(theta < 0.01f) theta = 0.01f;
   //      w(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni*nj]) / dx / theta; 
   //      w_valid(i,j,k) = 1;
   //   }
   //}
   compute_num = w.ni*w.nj*w.nk;
   slice = w.ni*w.nj;
   tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
   {
	   int k = thread_idx/slice;
	   int j = (thread_idx%slice)/w.ni;
	   int i = thread_idx%w.ni;
	   if(k>0 && k<w.nk-1 && j<w.nj && i<w.ni)
	   {
		   int index = i + j*ni + k*ni*nj;
		   if(w_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j,k-1) < 0)) {
			   float theta = 1;
			   if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
				   theta = fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
			   if(theta < 0.01f) theta = 0.01f;
			   w(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni*nj]) / dx / theta; 
			   w_valid(i,j,k) = 1;
		   }

	   }
   });
 
   //needs parallel
   int num;
   //for(unsigned int i = 0; i < u_valid.a.size(); ++i)
   //   if(u_valid.a[i] == 0)
   //      u.a[i] = 0;
   num = u_valid.a.size();
   tbb::parallel_for(0,num,1,[&](int i)
   {
	   if(u_valid.a[i] == 0)
		   u.a[i] = 0;
   });
   //for(unsigned int i = 0; i < v_valid.a.size(); ++i)
   //   if(v_valid.a[i] == 0)
   //      v.a[i] = 0;

   num = v_valid.a.size();
   tbb::parallel_for(0,num,1,[&](int i)
   {
	   if(v_valid.a[i] == 0)
		   v.a[i] = 0;
   });

   //for(unsigned int i = 0; i < w_valid.a.size(); ++i)
   //   if(w_valid.a[i] == 0)
   //      w.a[i] = 0;

   num = w_valid.a.size();
   tbb::parallel_for(0,num,1,[&](int i)
   {
	   if(w_valid.a[i] == 0)
		   w.a[i] = 0;
   });
}

void FluidSim::particle_interpolate(float alpha)
{
	int num = particles.size();
	tbb::parallel_for(0,num, 1, [&](int p){
		
		particles[p].vel = alpha * get_velocity(particles[p].pos) + (1-alpha)*(particles[p].vel + get_dvelocity(particles[p].pos));

	});
}

void FluidSim::compute_delta( Array3f & u, Array3f &u_old, Array3f &u_temp )
{
	int compute_num = u_temp.ni*u_temp.nj*u_temp.nk;
	int slice = u_temp.ni*u_temp.nj;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/u_temp.ni;
		int i = thread_idx%u_temp.ni;
		if( k<u_temp.nk && j<u_temp.nj && i<u_temp.ni)
		{
			u_temp(i,j,k) = u(i,j,k) - u_old(i,j,k);
		}
	});
}
void FluidSim::divide_weight(Array3f &u, Array3f & u_coef)
{
	int compute_num = u.ni*u.nj*u.nk;
	int slice = u.ni*u.nj;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/u.ni;
		int i = thread_idx%u.ni;
		if( k<u.nk && j<u.nj && i<u.ni)
		{
			u(i,j,k) = u(i,j,k)/(u_coef(i,j,k)+1e-5);
		}
	});
}
void FluidSim::FLIP_advection(float dt)
{
	//
	temp_u.set_zero();
	temp_v.set_zero();
	temp_w.set_zero();
	u_coef.set_zero();
	v_coef.set_zero();
	w_coef.set_zero();

	//compute delta, temp_U = U - U_save
	compute_delta(u, u_save,temp_u);
	compute_delta(v, v_save,temp_v);
	compute_delta(w, w_save,temp_w);

	//for each particle, p.vel = alpha*Interp(U) + (1-alpha)*(U_p + Interp(dU))
	particle_interpolate(0.03);
	//move particle
	float t = 0;
	float substep = cfl();
	while(t < dt) {

		if(t + substep > dt)
			substep = dt - t;
		advect_particles(substep);
		t+=substep;
	}

	//particle to grid
	u.set_zero();
	v.set_zero();
	w.set_zero();

	u_coef.assign(1e-7);
	v_coef.assign(1e-7);
	w_coef.assign(1e-7);
	marker_cell.assign(0);
	tbb::parallel_for(0,(int)(particles.size()),1,[&](int tidx)
	{
		Vec3i cell_ind(particles[tidx].pos / dx);
		int i = cell_ind[0];
		int j = cell_ind[1];
		int k = cell_ind[2];
		for(int kk = max(0,k-2); kk<=min(k+2, nk-1);kk++)
		for(int jj = max(0,j-2); jj<=min(j+2, nj-1);jj++)
		for(int ii = max(0,i-2); ii<=min(i+2, ni-1);ii++)
		{
			marker_cell(ii,jj,kk) = 1;
		}
	});
	vector<Vec3i> cell_list;
	cell_list.resize(0);
	for (int k=0; k<nk;k++)for(int j=0;j<nj;j++)for(int i=0;i<ni;i++)
	{
		if(marker_cell(i,j,k)==1)
		{
			cell_list.push_back(Vec3i(i,j,k));
		}
	}
	
	//sort particles into list
	vector<vector<int>> particle_hash;
	particle_hash.resize(ni*nj*nk);
	for (int p=0;p<particles.size();p++)
	{
		Vec3i cell_ind(particles[p].pos / dx);
		int index = ni*nj*cell_ind[2] + ni*cell_ind[1] + cell_ind[0];
		particle_hash[index].push_back(p);
	}
	//int compute_num = u.ni*u.nj*u.nk;
	//int slice = u.ni*u.nj;
	int compute_num = cell_list.size();
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = cell_list[thread_idx].v[2];
		int j = cell_list[thread_idx].v[1];
		int i = cell_list[thread_idx].v[0];
		if( k<nk && j<nj && i<ni)
		{
			Vec3f sample_pos0((float)i*dx, ((float)j+0.5f)*dx,((float)k+0.5f)*dx);
			Vec3f sample_pos1(((float)i+0.5f)*dx, ((float)j)*dx,((float)k+0.5f)*dx);
			Vec3f sample_pos2(((float)i+0.5f)*dx, ((float)j+0.5)*dx,((float)k)*dx);
			for(int kk = max(0,k-1); kk<=min(k+1, nk-1);kk++)
				for(int jj = max(0,j-1); jj<=min(j+1, nj-1);jj++)
					for(int ii = max(0,i-1); ii<=min(i+1, ni-1);ii++)
					{
						int index = kk*ni*nj + jj*ni + ii;
						for (int p=0; p<particle_hash[index].size();p++)
						{
							Vec3f pos = particles[particle_hash[index][p]].pos;
							Vec3f vel = particles[particle_hash[index][p]].vel;
							float weight0 = compute_weight(sample_pos0,pos);
							float weight1 = compute_weight(sample_pos1,pos);
							float weight2 = compute_weight(sample_pos2,pos);
							u(i,j,k) += weight0 * vel[0];
							u_coef(i,j,k)   += weight0;
							v(i,j,k) += weight1 * vel[1];
							v_coef(i,j,k)   += weight1;
							w(i,j,k) += weight2 * vel[2];
							w_coef(i,j,k)   += weight2;
						}

					}
		}
	});

	//compute_num = v.ni*v.nj*v.nk;
	//slice = v.ni*v.nj;
	//tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	//{
	//	int k = thread_idx/slice;
	//	int j = (thread_idx%slice)/v.ni;
	//	int i = thread_idx%v.ni;
	//	if( k<v.nk && j<v.nj && i<v.ni)
	//	{
	//		Vec3f sample_pos(((float)i+0.5f)*dx, ((float)j)*dx,((float)k+0.5f)*dx);
	//		for(int kk = max(0,k-1); kk<=min(k+1, nk-1);kk++)
	//			for(int jj = max(0,j-1); jj<=min(j+1, nj-1);jj++)
	//				for(int ii = max(0,i-1); ii<=min(i+1, ni-1);ii++)
	//				{
	//					int index = kk*ni*nj + jj*ni + ii;
	//					for (int p=0; p<particle_hash[index].size();p++)
	//					{
	//						Vec3f pos = particles[particle_hash[index][p]].pos;
	//						Vec3f vel = particles[particle_hash[index][p]].vel;
	//						float weight = compute_weight(sample_pos,pos);
	//						v(i,j,k) += weight * vel[1];
	//						v_coef(i,j,k)   += weight;
	//					}

	//				}
	//	}
	//});


	//compute_num = w.ni*w.nj*w.nk;
	//slice = w.ni*w.nj;
	//tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	//{
	//	int k = thread_idx/slice;
	//	int j = (thread_idx%slice)/w.ni;
	//	int i = thread_idx%w.ni;
	//	if( k<w.nk && j<w.nj && i<w.ni)
	//	{
	//		Vec3f sample_pos(((float)i+0.5f)*dx, ((float)j+0.5)*dx,((float)k)*dx);
	//		for(int kk = max(0,k-1); kk<=min(k+1, nk-1);kk++)
	//			for(int jj = max(0,j-1); jj<=min(j+1, nj-1);jj++)
	//				for(int ii = max(0,i-1); ii<=min(i+1, ni-1);ii++)
	//				{
	//					int index = kk*ni*nj + jj*ni + ii;
	//					for (int p=0; p<particle_hash[index].size();p++)
	//					{
	//						Vec3f pos = particles[particle_hash[index][p]].pos;
	//						Vec3f vel = particles[particle_hash[index][p]].vel;
	//						float weight = compute_weight(sample_pos,pos);
	//						w(i,j,k) += weight * vel[2];
	//						w_coef(i,j,k)   += weight;
	//					}

	//				}
	//	}
	//});
	//Estimate the liquid signed distance
	compute_phi(particle_hash,cell_list);
	particle_hash.resize(0);
	particle_hash.shrink_to_fit();
	cell_list.resize(0);
	cell_list.shrink_to_fit();

	//for (int p=0; p<particles.size();p++)
	//{
	//	Vec3i cell_ind(particles[p].pos / dx);
	//	//u
	//	for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,nk-1); ++k) {
	//		for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,nj-1); ++j) {
	//			for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,ni); ++i) {
	//				Vec3f sample_pos((float)i*dx, ((float)j+0.5f)*dx,((float)k+0.5f)*dx);
	//				float weight = compute_weight(sample_pos, particles[p].pos);
	//				u(i,j,k) += weight * (particles[p].vel[0]);
	//				u_coef(i,j,k) += weight;
	//			}
	//		}
	//	}
	//	//v
	//	for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,nk-1); ++k) {
	//		for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,nj); ++j) {
	//			for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,ni-1); ++i) {
	//				Vec3f sample_pos(((float)i+0.5f)*dx, ((float)j)*dx,((float)k+0.5f)*dx);
	//				float weight = compute_weight(sample_pos, particles[p].pos);
	//				v(i,j,k) += weight * (particles[p].vel[1]);
	//				v_coef(i,j,k) += weight;
	//			}
	//		}
	//	}
	//	//w
	//	for(int k = max(0,cell_ind[2] - 1); k <= min(cell_ind[2]+1,nk); ++k) {
	//		for(int j = max(0,cell_ind[1] - 1); j <= min(cell_ind[1]+1,nj-1); ++j) {
	//			for(int i = max(0,cell_ind[0] - 1); i <= min(cell_ind[0]+1,ni-1); ++i) {
	//				Vec3f sample_pos(((float)i+0.5f)*dx, ((float)j+0.5)*dx,((float)k)*dx);
	//				float weight = compute_weight(sample_pos, particles[p].pos);
	//				w(i,j,k) += weight * (particles[p].vel[2]);
	//				w_coef(i,j,k) += weight;
	//			}
	//		}
	//	}
	//}

	divide_weight(u,u_coef);
	divide_weight(v,v_coef);
	divide_weight(w,w_coef);

	u_save = u;
	v_save = v;
	w_save = w;

	//modifyFreeSurfaceVelocity();
}



//Apply several iterations of a very simple propagation of valid velocity data in all directions
void extrapolate(Array3f& grid, Array3c& valid) {

   Array3f temp_grid = grid;
   Array3c old_valid(valid.ni,valid.nj,valid.nk);
   for(int layers = 0; layers < 10; ++layers) {
      old_valid = valid;

	  int num = grid.ni*grid.nj*grid.nk;
	  int slice =  grid.ni*grid.nj;
	  tbb::parallel_for(0,num, 1, [&](int thread_idx)
	  {
		  int k = thread_idx/slice;
		  int j = (thread_idx%slice)/grid.ni;
		  int i = thread_idx%grid.ni;
      //for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) {
		  if(i>0&&i<grid.ni-1 &&j>0&&j<grid.nj-1 &&k>0&&k<grid.nk-1)
		  {
			  float sum = 0;
			  int count = 0;

			  if(!old_valid(i,j,k)) 
			  {

				  if(old_valid(i+1,j,k))
				  {
					  sum += grid(i+1,j,k);
					  ++count;
				  }
				  if(old_valid(i-1,j,k)) 
				  {
					  sum += grid(i-1,j,k);
					  ++count;
				  }
				  if(old_valid(i,j+1,k)) 
				  {
					  sum += grid(i,j+1,k);
					  ++count;
				  }
				  if(old_valid(i,j-1,k)) 
				  {
					  sum += grid(i,j-1,k);
					  ++count;
				  }
				  if(old_valid(i,j,k+1))
				  {
					  sum += grid(i,j,k+1);
					  ++count;
				  }
				  if(old_valid(i,j,k-1)) 
				  {
					  sum += grid(i,j,k-1);
					  ++count;
				  }

				  //If any of neighbour cells were valid, 
				  //assign the cell their average value and tag it as valid
				  if(count > 0) 
				  {
					  temp_grid(i,j,k) = sum /(float)count;
					  valid(i,j,k) = 1;
				  }
			  }

         }
      });
      grid = temp_grid;

   }

}
