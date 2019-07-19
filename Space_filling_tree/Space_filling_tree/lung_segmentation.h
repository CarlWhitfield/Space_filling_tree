#ifndef LUNG_SEGMENTATION_H
#define LUNG_SEGMENTATION_H

#include<stdio.h>
#include"stl_reader.h"
#include<network_3D.h>
#include<list_template.h>

class PointCloud
{
protected:
	std::vector<network::Position*> cloud;
	network::Position com;
	inline void PointCloud::calc_com()
	{
		using namespace network;
		this->com = Position(0,0,0);
		for(size_t npt = 0; npt < this->cloud.size(); npt++)
		{
			this->com = this->com + *(this->cloud[npt]);
		}
		this->com = this->com/((double) this->cloud.size());
	}
public:
	PointCloud(const std::vector<network::Position*> & cld)
	{ 
		this->cloud = cld; 
		this->calc_com();
	}
	std::vector<PointCloud*> partition_point_cloud(const std::vector<network::Position> & pos_to_be_assigned) const;
	std::vector<PointCloud*> split_point_cloud(const network::Position & plane_pt, const network::Position & plane_norm) const;

	inline network::Position get_com() const { return (this->com); }
	inline size_t point_count() const { return this->cloud.size(); } 

	void remove_closest_point( const network::Position & pos);
};

class LobeMesh: public stl_reader::StlMesh<double, unsigned int>
{
protected:
	network::Position minpos, maxpos;
	double bin_spacing[2], volume;
	std::vector<network::Position> min_xyz, max_xyz;
	std::vector<std::vector<std::vector<size_t>>> xy_binned_face_indices;
	std::string id;
public:
	LobeMesh(const std::string & fname, const size_t & Nbins);
	PointCloud *pt_cloud;

	bool is_in(const network::Position & pos);

	inline network::Position get_tri_corner_coords(const size_t & ntri, const size_t & ncorner) const 
	{
		const double *c = this->tri_corner_coords(ntri, ncorner);
		return network::Position(c[0], c[1], c[2]);
	}
	inline network::Position get_tri_normal(const size_t & ntri) const 
	{
		const double *c = this->tri_normal(ntri);
		network::Position p(c[0], c[1], c[2]);
		p.normalise();
		return p;
	}
	inline std::vector<size_t> get_bin_coord(const network::Position & pos) const 
	{
		std::vector<size_t> coord;
		coord.resize(2);
		for(size_t dim = 0; dim < 2; dim++) coord[dim] = ((size_t) ((pos.x[dim] - this->minpos.x[dim])/bin_spacing[dim]));
		return coord;  
	}
	inline void initialise_point_cloud(const std::vector<network::Position*> & pv)
	{
		this->pt_cloud = new PointCloud(pv);
	}
	inline double get_volume() const { return this->volume; }
	inline network::Position get_minpos() const { return this->minpos; }
	inline network::Position get_maxpos() const { return this->maxpos; }
	inline std::string get_id() const { return this->id; }
};

class LungSegNetwork: public network::Network<network::Node,network::Edge<network::Node>> 
{
public:
	LungSegNetwork(std::map<long int, network::Node*> &node, std::map<long int, network::Edge<network::Node>*> &edge)
		:Network<network::Node,network::Edge<network::Node>>(node, edge){};

	size_t grow_lobe(LobeMesh*, const double & cloud_h, const std::vector<network::Node*> &, const bool & prune,
		           const bool & reassign, const double & l_factor, const double & max_angle,
				   const double & l_cutoff, const size_t & prev_iter);
	void print_step_by_step(const std::vector<size_t> & end_nodes, const std::vector<PointCloud*> & pt_clouds, 
		                          const size_t & num, const std::string lobe_id, const double & cloud_h) const;
	void fill_in_radii();
};

class LungSegOptions: public inlist::OptionList<bool,char>
{
protected:

public:
	LungSegOptions();
};

class LungSegParams: public inlist::ParameterList<int,double>
{
protected:

public:
	LungSegParams();
};

class LungSegmentation
{
protected:
	std::vector<LobeMesh*> lobe_meshes;
	double volume;
	network::Position min, max;

	inline void vol_calc()	
	{
		this->volume = 0;
		for(size_t n = 0; n < this->lobe_meshes.size(); n++) this->volume += this->lobe_meshes[n]->get_volume();
	}

	inline void min_calc()
	{
		this->min = this->lobe_meshes[0]->get_minpos();
		for(size_t n = 1; n < this->lobe_meshes.size(); n++) 
		{
			for(size_t dim = 0; dim < 3; dim++)
			{
				if(this->lobe_meshes[n]->get_minpos().x[dim] < this->min.x[dim])
				{
					this->min.x[dim] = this->lobe_meshes[n]->get_minpos().x[dim];
				}
			}
		}
	}

	inline void max_calc()
	{
		this->max = this->lobe_meshes[0]->get_maxpos();
		for(size_t n = 1; n < this->lobe_meshes.size(); n++) 
		{
			for(size_t dim = 0; dim < 3; dim++)
			{
				if(this->lobe_meshes[n]->get_maxpos().x[dim] > this->max.x[dim]) 
				{
					this->max.x[dim] = this->lobe_meshes[n]->get_maxpos().x[dim];
				}
			}
		}
	}

	inline void initialise()
	{
		this->vol_calc();
		this->min_calc();
		this->max_calc();
	}
public:
	LungSegNetwork *tree;
	LungSegOptions *options;
	LungSegParams *params;

	LungSegmentation(){};
	LungSegmentation(int argc, char *argv[]);

	LungSegNetwork* read_PTK_node_and_element_files(const std::string & node_fname, const std::string & el_fname);
	void build_tree();

	inline size_t count_lobes() const { return (this->lobe_meshes.size()); }
	inline LobeMesh* get_lobe(const size_t & n) const { return (this->lobe_meshes[n]); }
	inline double get_lung_volume() const { return this->volume; }

	inline network::Position get_minpos() const { return this->min; }
	inline network::Position get_maxpos() const { return this->max; }
};

network::Position new_angle_calc(const network::Position & old_vec, const network::Position & parent_dir, const double & max_angle);

double calc_radius(const size_t & pgen, const double & prad, const size_t & tgen);

#endif