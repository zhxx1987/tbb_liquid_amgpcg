#include <iostream>
#include "faithful_surfacing3.h"

FaithfulSurfacing3::
	FaithfulSurfacing3(void)
	: inner_radius(0.5f),
	outer_radius(1),
	grid_dx(0.3f),
	gs_damping(0.5f),
	shrink_steps(5),
	smooth_steps(25),
	free_smooth_steps(3),
	tri(march.tri),
	x(march.x)
{}

void FaithfulSurfacing3::
	run_surfacing(void)
{
	//std::cout<<"INIT"<<std::endl;
	init();
	//std::cout<<"BUILDING PHI AND BUCKETS"<<std::endl;
	build_phi_and_particle_buckets();
	//std::cout<<"CONTOURING"<<std::endl;
	march.contour_grid();
	//std::cout<<"BUILDING_A"<<std::endl;
	build_A();
	for(int s=0; s<shrink_steps; ++s){
		//std::cout<<"SHRINK "<<s<<std::endl;
		shrink_mesh(true);
	}
#if 0
	//std::cout<<"BUILDING_B"<<std::endl;
	build_B();
	for(int s=0; s<smooth_steps; ++s){
		//std::cout<<"SMOOTH "<<s<<std::endl;
		smooth_mesh(true);
	}
	for(int s=0; s<free_smooth_steps; ++s){
		//std::cout<<"FREE SMOOTH "<<s<<std::endl;
		smooth_mesh(false);
	}
#else
	for(int s=0; s<free_smooth_steps; ++s){
		//std::cout<<"FREE SHRINK "<<s<<std::endl;
		shrink_mesh(false);
	}
#endif
}

void FaithfulSurfacing3::
	init(void)
{
	// find bounding box
	Vec3f xmin(particle_x[0]), xmax(particle_x[0]);
	for(unsigned int p=1; p<particle_x.size(); ++p) update_minmax(particle_x[p], xmin, xmax);
	// set up grid sizes
	origin=xmin-Vec3f(outer_radius+2*grid_dx);
	nx=8+4*(int)std::floor((xmax[0]+outer_radius+2*grid_dx-origin[0])/(4*grid_dx));
	ny=8+4*(int)std::floor((xmax[1]+outer_radius+2*grid_dx-origin[1])/(4*grid_dx));
	nz=8+4*(int)std::floor((xmax[2]+outer_radius+2*grid_dx-origin[2])/(4*grid_dx));
	//std::cout<<"   grid size "<<nx<<" "<<ny<<" "<<nz<<std::endl;
	march.origin=origin;
	march.dx=grid_dx;
}

void FaithfulSurfacing3::
	build_phi_and_particle_buckets(void)
{
	march.phi.assign(nx, ny, nz, 2*sqr(outer_radius));
	for(unsigned int p=0; p<particle_x.size(); ++p){
		int i0=(int)((particle_x[p][0]-origin[0]-(outer_radius+1.5*grid_dx))/grid_dx); if(i0<0) i0=0;
		int i1=(int)((particle_x[p][0]-origin[0]+(outer_radius+1.5*grid_dx))/grid_dx); if(i1>=nx) i1=nx-1;
		int j0=(int)((particle_x[p][1]-origin[1]-(outer_radius+1.5*grid_dx))/grid_dx); if(j0<0) j0=0;
		int j1=(int)((particle_x[p][1]-origin[1]+(outer_radius+1.5*grid_dx))/grid_dx); if(j1>=ny) j1=ny-1;
		int k0=(int)((particle_x[p][2]-origin[2]-(outer_radius+1.5*grid_dx))/grid_dx); if(k0<0) k0=0;
		int k1=(int)((particle_x[p][2]-origin[2]+(outer_radius+1.5*grid_dx))/grid_dx); if(k1>=nz) k1=nz-1;
		for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j) for(int i=i0; i<=i1; ++i){
			float d2=sqr(particle_x[p][0] - (origin[0]+grid_dx*i))
				+sqr(particle_x[p][1] - (origin[1]+grid_dx*j))
				+sqr(particle_x[p][2] - (origin[2]+grid_dx*k));
			if(d2<march.phi(i,j,k)) march.phi(i,j,k)=d2;
		}
		//if(p%10000==0) //std::cout<<" particle "<<p<<std::endl;
	}
	for(int k=0; k<nz; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx; ++i){
		march.phi(i,j,k)=std::sqrt(march.phi(i,j,k))-outer_radius;
	}
	particle_bucket.resize(0, 0, 0);
	particle_bucket.resize(nx, ny, nz);
	for(unsigned int p=0; p<particle_x.size(); ++p){
		int i=(int)((particle_x[p][0]-origin[0])/grid_dx); assert(i>=0 && i<nx);
		int j=(int)((particle_x[p][1]-origin[1])/grid_dx); assert(j>=0 && j<ny);
		int k=(int)((particle_x[p][2]-origin[2])/grid_dx); assert(k>=0 && k<nz);
		particle_bucket(i,j,k).push_back(p);
	}
}

void FaithfulSurfacing3::
	build_A(void)
{
	A.zero();
	A.resize((unsigned int)x.size());
	for(unsigned int t=0; t<tri.size(); ++t){
		int i, j, k; assign(tri[t], i, j, k);
		A.add_to_element(i, j, -0.5);
		A.add_to_element(i, k, -0.5);
		A.add_to_element(j, i, -0.5);
		A.add_to_element(j, k, -0.5);
		A.add_to_element(k, i, -0.5);
		A.add_to_element(k, j, -0.5);
	}
	// set diagonal to get zero row sums
	for(unsigned int i=0; i<x.size(); ++i){
		float s=0;
		for(unsigned int a=0; a<A.index[i].size(); ++a) s+=A.value[i][a];
		A.add_to_element(i, i, -s);
	}
}

void FaithfulSurfacing3::
	shrink_mesh(bool enforce_fidelity)
{
	for(unsigned int i=0; i<x.size(); ++i){
		Vec3f xnew(0.f, 0.f, 0.f);
		double d;
		for(unsigned int a=0; a<A.index[i].size(); ++a){
			if(A.index[i][a]==i) d=A.value[i][a];
			else xnew-=A.value[i][a]*x[A.index[i][a]];
		}
		xnew/=(float)d;
		x[i]+=gs_damping*(xnew-x[i]);
		if(enforce_fidelity) make_point_faithful(x[i]);
	}
}

void FaithfulSurfacing3::
	build_B(void)
{
	SparseMatrixf W((unsigned int)x.size(),9);
	std::vector<float> D(x.size(), 0);
	for(unsigned int t=0; t<tri.size(); ++t){
		int i, j, k; assign(tri[t], i, j, k);
		float area=(float)(mag(cross(x[j]-x[i], x[k]-x[i]))+1e-30);
		float coti=dot(x[j]-x[i], x[k]-x[i])/area; if(!(coti>0.01f)) coti=0.01f;
		float cotj=dot(x[i]-x[j], x[k]-x[j])/area; if(!(cotj>0.01f)) cotj=0.01f;
		float cotk=dot(x[i]-x[k], x[j]-x[k])/area; if(!(cotk>0.01f)) cotk=0.01f;
		W.add_to_element(i, i, cotj+cotk);
		W.add_to_element(i, j, -cotk);
		W.add_to_element(i, k, -cotj);
		W.add_to_element(j, i, -cotk);
		W.add_to_element(j, j, coti+cotk);
		W.add_to_element(j, k, -coti);
		W.add_to_element(k, i, -cotj);
		W.add_to_element(k, j, -coti);
		W.add_to_element(k, k, coti+cotj);
		float minarea=0.05f*mag(x[j]-x[i])*mag(x[k]-x[i]);
		if(area<minarea) area=minarea;
		D[i]+=area;
		D[j]+=area;
		D[k]+=area;
	}
	for(unsigned int i=0; i<D.size(); ++i){
		D[i]=(float)(1./D[i]);
	}
	multiply_sparse_matrices_with_diagonal_weighting(W, D, W, B);
}

void FaithfulSurfacing3::
	smooth_mesh(bool enforce_fidelity)
{
	for(unsigned int i=0; i<x.size(); ++i){
		Vec3f xnew(0.f, 0.f, 0.f);
		float d=0;
		for(unsigned int a=0; a<B.index[i].size(); ++a){
			if(B.index[i][a]==i) d=B.value[i][a];
			else xnew-=B.value[i][a]*x[B.index[i][a]];
		}
		assert(d>0 && d<1e15);
		assert(xnew[0]==xnew[0]);
		assert(xnew[1]==xnew[1]);
		assert(xnew[2]==xnew[2]);
		xnew/=d;
		x[i]=xnew;
		if(enforce_fidelity) make_point_faithful(x[i]);
	}
}

void FaithfulSurfacing3::
	make_point_faithful(Vec3f& xp)
{
	// first find closest particle
	int nearest_q=-1;
	float dist2_q=1e30f; // should be an upper bound!
	float search_band=outer_radius;
find_nearest_particle:
	int i0=(int)((xp[0]-search_band-origin[0])/grid_dx); if(i0<0) i0=0; else if(i0>nx-1) i0=nx-1;
	int i1=(int)((xp[0]+search_band-origin[0])/grid_dx); if(i1<0) i1=0; else if(i1>nx-1) i1=nx-1;
	int j0=(int)((xp[1]-search_band-origin[1])/grid_dx); if(j0<0) j0=0; else if(j0>ny-1) j0=ny-1;
	int j1=(int)((xp[1]+search_band-origin[1])/grid_dx); if(j1<0) j1=0; else if(j1>ny-1) j1=ny-1;
	int k0=(int)((xp[2]-search_band-origin[2])/grid_dx); if(k0<0) k0=0; else if(k0>nz-1) k0=nz-1;
	int k1=(int)((xp[2]+search_band-origin[2])/grid_dx); if(k1<0) k1=0; else if(k1>nz-1) k1=nz-1;
	for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j) for(int i=i0; i<=i1; ++i){
		const std::vector<unsigned int>& bucket=particle_bucket(i,j,k);
		for(unsigned int a=0; a<bucket.size(); ++a){
			unsigned int q=bucket[a];
			float d2=dist2(xp, particle_x[q]);
			if(d2<dist2_q){
				dist2_q=d2;
				nearest_q=q;
			}
		}
	}
	if(nearest_q==-1){
		search_band+=outer_radius;
		goto find_nearest_particle;
	}
	// and now make sure xp is within limits
	float d=std::sqrt(dist2_q);
	if(d<inner_radius){
		xp+=((inner_radius-d)/d)*(xp-particle_x[nearest_q]);
	}else if(d>outer_radius){
		xp-=((d-outer_radius)/d)*(xp-particle_x[nearest_q]);
	}
}