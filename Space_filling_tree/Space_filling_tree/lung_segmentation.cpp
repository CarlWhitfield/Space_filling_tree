#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
	#ifndef _XOPEN_SOURCE
		#define _XOPEN_SOURCE
	#endif
#endif

#include"lung_segmentation.h"
#include"file_manip.h"
#include"globals.h"
#include<math.h>
#include"define_seg_codes.h"
#include<string>
#include<stdlib.h>
#include<time.h>

//need to comment everything

void PointCloud::remove_closest_point( const network::Position & pos)
{
	//removes point closest to pos from the point cloud
	size_t i_min = 0;
	double dist_min = this->cloud[0]->distance(pos);
	for(size_t ipt = 1; ipt < this->cloud.size(); ipt++)
	{
		double disth = this->cloud[ipt]->distance(pos);
		if(disth < dist_min)
		{
			dist_min = disth;
			i_min = ipt;
		}
	}
	this->cloud.erase(this->cloud.begin() + i_min);
}

LobeMesh::LobeMesh(const std::string & fname, const size_t & Nbins):StlMesh<double, unsigned int>(fname)
{
	//initialise lobe mesh
	this->min_xyz.clear();
	this->min_xyz.resize(this->num_tris());
	this->max_xyz.clear();
	this->max_xyz.resize(this->num_tris());

	std::vector<std::string> fsplit = string_split(fname, ".");
	this->id = string_split(fsplit[fsplit.size()-2],"_").back();

	network::Position p = this->get_tri_corner_coords(0, 0);
	//initialise
	this->minpos = p;
	this->maxpos = p;
	this->volume = 0;
	for(size_t ntri = 0; ntri < this->num_tris(); ntri++)   //loop over stl file triangles
	{
		double x_centroid = 0;
		//record the smallest x y and z values
		for(size_t corner = 0; corner < 3; corner++)
		{
			network::Position p = this->get_tri_corner_coords(ntri, corner);
			if(corner==0) 
			{
				this->min_xyz[ntri] = p;
				this->max_xyz[ntri] = p;
			}
			else
			{
				for(size_t d = 0; d < 3; d++)
				{
					if(p.x[d] < this->min_xyz[ntri].x[d]) this->min_xyz[ntri].x[d] = p.x[d];
					if(p.x[d] > this->max_xyz[ntri].x[d]) this->max_xyz[ntri].x[d] = p.x[d];
				}
			}
			x_centroid += p.x[0] / 3.0;
		}
		for(size_t d = 0; d < 3; d++)
		{
			if(this->min_xyz[ntri].x[d] < this->minpos.x[d]) this->minpos.x[d] = this->min_xyz[ntri].x[d];
			if(this->max_xyz[ntri].x[d] > this->maxpos.x[d]) this->maxpos.x[d] = this->max_xyz[ntri].x[d];
		}

		network::Position norm = this->get_tri_normal(ntri);
		//get side lengths
		double l[3];
		l[0] = (this->get_tri_corner_coords(ntri,1) - this->get_tri_corner_coords(ntri,0)).magnitude();
		l[1] = (this->get_tri_corner_coords(ntri,2) - this->get_tri_corner_coords(ntri,0)).magnitude();
		l[2] = (this->get_tri_corner_coords(ntri,2) - this->get_tri_corner_coords(ntri,1)).magnitude();
		double s = 0.5*(l[0] + l[1] + l[2]);
		//get area
		double area = sqrt(s*(s - l[0])*(s - l[1])*(s - l[2]));
		//add to volume sum
		volume += norm.x[0] * x_centroid * area;
	}

	//bin points into a grid in x and y
	this->bin_spacing[0] = (1.001*this->maxpos.x[0] - this->minpos.x[0])/double(Nbins);
	this->bin_spacing[1] = (1.001*this->maxpos.x[1] - this->minpos.x[1])/double(Nbins);
	this->xy_binned_face_indices.clear();
	this->xy_binned_face_indices.resize(Nbins);
	for(size_t i = 0; i < Nbins; i++)
	{
		xy_binned_face_indices[i].resize(Nbins);
	}

	for(size_t ntri = 0; ntri < this->num_tris(); ntri++)   //bin stl file triangles
	{
		std::vector<size_t> min_ind = this->get_bin_coord(this->min_xyz[ntri]);
		std::vector<size_t> max_ind = this->get_bin_coord(this->max_xyz[ntri]);

		for(size_t i = min_ind[0]; i <= max_ind[0]; i++)
		{
			for(size_t j = min_ind[1]; j <= max_ind[1]; j++)
			{
				this->xy_binned_face_indices[i][j].push_back(ntri);
			}
		}
	}
}

bool LobeMesh::is_in(const network::Position & pos)
{
	//find if a point is inside a particular lobe
	bool check = true;
	for(size_t dim = 0; dim < 3; dim++)   //check if it is inside the binned grid
	{
		if(pos.x[dim] < this->minpos.x[dim] || pos.x[dim] > this->maxpos.x[dim]) check = false;
		//std::cout << this->minpos.x[dim] << ' ' << this->maxpos.x[dim] << '\n';
	}
	size_t faces_crossed = 0;
	if(check)
	{
		std::vector<size_t> p_ind = this->get_bin_coord(pos);
		std::vector<size_t> reduced_tri_list;
		reduced_tri_list.reserve(this->xy_binned_face_indices[p_ind[0]][p_ind[1]].size());
		//get list of triangles that it might cross
		for(size_t il = 0; il < this->xy_binned_face_indices[p_ind[0]][p_ind[1]].size(); il++)
		{
			size_t ntri = this->xy_binned_face_indices[p_ind[0]][p_ind[1]][il];
			bool tri_check = true;
			size_t dim = 0;
			while(tri_check == true && dim < 3)
			{
				//check position is not less that min triangle coord (except in z)
				if(dim < 2 && pos.x[dim] < this->min_xyz[ntri].x[dim]) tri_check = false;
				//check position is not more that max triangle coord 
				if(pos.x[dim] > this->max_xyz[ntri].x[dim]) tri_check = false;
				dim++;
			}
			if(tri_check) reduced_tri_list.push_back(ntri);
		}
		//loop over reduced list and cast point onto tangent surface of triangles
		std::vector<double> final_zpts_on_corner;
		for(size_t il = 0; il < reduced_tri_list.size(); il++)
		{
			size_t ntri = reduced_tri_list[il];
			network::Position norm = this->get_tri_normal(ntri);
			if(norm.x[2]*norm.x[2] > 0)    //if z component of normal is non-zero
			{
				network::Position c[3];
				for(size_t nc = 0; nc < 3; nc++)
				{
					c[nc] = this->get_tri_corner_coords(ntri, nc);
				}
				network::Position pos_cast(pos.x[0], pos.x[1], c[0].x[2] + (norm.x[0]*(c[0].x[0] - pos.x[0]) 
					                                          + norm.x[1]*(c[0].x[1] - pos.x[1]))/norm.x[2]);
				if(pos_cast.x[2] > pos.x[2])
				{
					//check if crossing point is inside triangle
					network::Position u[3];
					for(size_t nc = 0; nc < 3; nc++) u[nc] = pos_cast - c[nc];
					if(u[0].magnitude() == 0 || u[1].magnitude() == 0 || u[2].magnitude() == 0)
					{
						bool duplicate_check = false;
						for(size_t ifl = 0; ifl < final_zpts_on_corner.size(); ifl++)
						{
							if(pos_cast.x[2] == final_zpts_on_corner[ifl]) duplicate_check = true;
						}

						if(!duplicate_check)
						{
							faces_crossed++;
							final_zpts_on_corner.push_back(pos_cast.x[2]);
						}
					}
					else
					{
						double tot_angle = 0;
						for(size_t nc = 0; nc < 3; nc++) tot_angle += u[nc].angle(u[(nc+1)%3]);
						//std::cout << tot_angle << '\n';
						if(tot_angle > 2*M_PI - 1E-03 && tot_angle < 2*M_PI + 1E-03)
						{
							faces_crossed++;
							/*for(size_t nc = 0; nc < 3; nc++) std::cout << c[nc].x[0] << ' ' << c[nc].x[1] << ' ' << c[nc].x[2] << '\n';
							std::cout << "\n";*/
						}
					}
				}
			}
		}
	}
	//std::cout << faces_crossed << '\n';
	if(faces_crossed % 2) return true;
	else return false;
}

LungSegOptions::LungSegOptions():OptionList<bool,char>()
{
	//initialise options for running the code
	//bool
	std::shared_ptr<inlist::Option<bool>> pt = std::make_shared<inlist::Option<bool>>(false, std::string(PRUNE_TERM_BRANCHES_KEY));
	std::shared_ptr<inlist::Option<bool>> ro = std::make_shared<inlist::Option<bool>>(true, std::string(REASSIGN_KEY));
	this->add(PRUNE_TERM_BRANCHES_KEY, pt);
	this->add(REASSIGN_KEY, ro);
}

LungSegParams::LungSegParams():ParameterList<int,double>()
{
	//initialise code parameters

	//int
	std::shared_ptr<inlist::Parameter<int>> np = std::make_shared<inlist::Parameter<int>>(100, std::string(NUMBER_MESH_BINS_KEY));
	std::shared_ptr<inlist::Parameter<int>> ac = std::make_shared<inlist::Parameter<int>>(30000, std::string(ACINI_COUNT_KEY));
	this->add(NUMBER_MESH_BINS_KEY, np);
	this->add(ACINI_COUNT_KEY, ac);

	//double
	std::shared_ptr<inlist::Parameter<double>> dcom = std::make_shared<inlist::Parameter<double>>(0.4, std::string(DISTANCE_TO_COM_FACTOR_KEY));
	std::shared_ptr<inlist::Parameter<double>> mba = std::make_shared<inlist::Parameter<double>>(M_PI/3.0, std::string(MAX_BRANCHING_ANGLE_KEY));
	std::shared_ptr<inlist::Parameter<double>> lc = std::make_shared<inlist::Parameter<double>>(1.0, std::string(LENGTH_CUTOFF_MM_KEY));
	this->add(DISTANCE_TO_COM_FACTOR_KEY, dcom);
	this->add(MAX_BRANCHING_ANGLE_KEY, mba);
	this->add(LENGTH_CUTOFF_MM_KEY, lc);
}

LungSegmentation::LungSegmentation(int argc, char *argv[])
{
	//initialise segmentation algorithm
	using namespace stl_reader;

	this->options = std::make_shared<LungSegOptions>();
	this->params = std::make_shared<LungSegParams>();

	bool node_file_done = false, element_file_done = false;
	size_t node_file_index, el_file_index;
	std::vector<std::string> ext;
	ext.resize(4);
	ext[0] = "txt";
	ext[1] = "stl";
	ext[2] = "options";
	ext[3] = "params";

	this->options->get_filenames_from_args(ext, argc, argv);
	//extract all filenames
	//find node and element files from PTK and process
	for(size_t n = 0; n < this->options->count_files_with_ext("txt"); n++)
	{
		std::vector<std::string> fname_split = string_split(this->options->get_filename("txt", n), ".");
		std::string fhead = fname_split[fname_split.size()-2];
		std::vector<std::string> fhead_split = string_split(fhead, "_");
		if (fhead_split.back() == "node" && !node_file_done) //process node file
		{
			node_file_index = n;
			node_file_done = true;
		}
		if (fhead_split.back() == "element" && !element_file_done)  //process element file
		{
			el_file_index = n;
			element_file_done = true;
		}
	}
	//find and process stl files
	unsigned stl_files_done = 0;
	for(size_t n = 0; n < this->options->count_files_with_ext("stl"); n++)
	{
		this->lobe_meshes.push_back(std::make_shared<LobeMesh>(this->options->get_filename("stl", n), 
									((size_t) this->params->get_param<int>(NUMBER_MESH_BINS_KEY)->get_value())));
		stl_files_done++;
	}

	if(node_file_done)
	{
		if(element_file_done)
		{
			if(stl_files_done > 0)
			{
				//build tree
				this->tree = this->read_PTK_node_and_element_files(this->options->get_filename("txt", node_file_index),
					                                                this->options->get_filename("txt", el_file_index));
			}
			else
			{
				std::cout << "Error, no stl files found.\n";
				abort_on_failure();
			}
		}
		else
		{
			std::cout << "Error, no element file found.\n";
			abort_on_failure();
		}
	}
	else
	{
		std::cout << "Error, no node file found.\n";
		abort_on_failure();
	}
	//read options and param files
	if(this->options->filename_exists("options")) this->options->read_file(this->options->get_filename("options"));
	if(this->options->filename_exists("params")) this->params->read_file(this->options->get_filename("params"));
	//set up segmentation
	this->initialise();
}

std::shared_ptr<LungSegNetwork> LungSegmentation::read_PTK_node_and_element_files(const std::string & node_fname, const std::string & el_fname)
{
	//get info from PTK output files
	using namespace std;
	
	vector<vector<double>> node_file_lines;
	if(!check_infile(node_fname))
	{ 
		node_file_lines = parse_csv_file<double>(node_fname);
	}
	else
	{
		cout << "Problem reading node file, aborting.\n";
		abort_on_failure();
	}

	//first line is number of nodes - 1
	size_t Nnodes = ((size_t) node_file_lines[0][0]) + 1;
	if(node_file_lines.size() != Nnodes+1)
	{
		cout << "Problem: node file has wrong number of lines.\n";
	}
	map<long int, std::shared_ptr<network::Node>> node_map;
	map<long int, double> node_radii;
	for(size_t l = 1; l <= Nnodes; l++)
	{
		node_map[((long int) node_file_lines[l][0])] = std::make_shared<network::Node>(node_file_lines[l][1], node_file_lines[l][2],
			                                                             node_file_lines[l][3]);
		node_radii[((long int) node_file_lines[l][0])] = node_file_lines[l][4];
	}

	vector<vector<long int>> edge_file_lines;
	if(!check_infile(el_fname))
	{ 
		edge_file_lines = parse_csv_file<long int>(el_fname);
	}
	else
	{
		cout << "Problem reading edge file, aborting.\n";
		abort_on_failure();
	}
	size_t Nedges = ((size_t) edge_file_lines[0][0]);
	if(Nedges != Nnodes-1)
	{
		cout << "Problem, node and edge files do not match.\n";
		abort_on_failure();
	}
	map<long int, std::shared_ptr<network::Edge<network::Node>>> edge_map;
	for(long int l = 1; l < ((long int) edge_file_lines.size()); l++)
	{
		if(edge_file_lines[l][0] != edge_file_lines[l][1]) //ignore loop edges 
		{
			edge_map[l] = std::make_shared<network::Edge<network::Node>>(node_map[edge_file_lines[l][0]],
				                                           node_map[edge_file_lines[l][1]], 1.0,
														   node_radii[edge_file_lines[l][1]]);
		}
	}

	return (std::make_shared<LungSegNetwork>(node_map,edge_map));
}

void LungSegmentation::build_tree()
{
	//code to build whole airway tree using PTK outputs
	using namespace std;

	//assign terminal edges to segmented regions
	vector<vector<network::Node*>> lobe_termnodes;
	lobe_termnodes.resize(this->count_lobes());
	for(size_t k = this->tree->get_first_term_index(); k < this->tree->count_nodes(); k++)
	{
		bool assigned = false;
		int iterations = 0;
		while(!assigned && iterations < 10)
		{
			for(size_t n = 0; n < this->count_lobes(); n++)
			{
				if(this->get_lobe(n)->is_in(this->tree->get_node(k)->get_pos()))
				{
					lobe_termnodes[n].push_back(tree->get_node(k));
					assigned = true;
				}
			}
			if(assigned == false)
			{
				cout << "Warning, node " << k << " not assigned to lobe, extending branch.\n";
				size_t j = this->tree->get_edge_in_index(k,0);   //find edge in
				size_t k_in =  this->tree->get_node_in_index(j);  //find node in
				network::Position sep = this->tree->get_node(k)->get_pos() - this->tree->get_node(k_in)->get_pos();
				this->tree->get_node(k)->set_pos(this->tree->get_node(k_in)->get_pos() + sep*1.01);  //move along by 10%
				this->tree->get_edge(j)->update_length_from_nodes(); //update branch length
				cout << "New pos of node " << k << " is (" << this->tree->get_node(k)->get_pos().x[0] << ',' 
					 << this->tree->get_node(k)->get_pos().x[1] << ',' << this->tree->get_node(k)->get_pos().x[2] << ").\n";
					 
			}
			iterations++;
		}
		if(!assigned)
		{
			cout << "Error, node " << k << " not assigned to lobe.\n";
		}
	}

	vector<vector<std::shared_ptr<network::Position>>> point_clouds;
	point_clouds.resize(this->count_lobes());
	for(size_t nl = 0; nl < this->count_lobes(); nl++)
	{
		point_clouds[nl].reserve(((size_t) ((this->get_lobe(nl)->get_volume() / this->get_lung_volume()) 
			                               *((double) this->params->get_param<int>(ACINI_COUNT_KEY)->get_value()))));
	}

	double h = pow(this->get_lung_volume() / ((double) this->params->get_param<int>(ACINI_COUNT_KEY)->get_value()), 1.0/3.0);
	size_t Nx = ((size_t)   ((this->get_maxpos().x[0] - this->get_minpos().x[0])/h)) + 1;
	size_t Ny = ((size_t)   ((this->get_maxpos().x[1] - this->get_minpos().x[1])/h)) + 1;
	size_t Nz = ((size_t)   ((this->get_maxpos().x[2] - this->get_minpos().x[2])/h)) + 1;
	srand (time(NULL));
	for(size_t i = 0; i < Nx; i++)
	{
		for(size_t j = 0; j < Ny; j++)
		{
			for(size_t k = 0; k < Nz; k++)
			{
				double rand0 = 2*(rand() - 0.5) / ((double) RAND_MAX);
				double rand1 = 2*(rand() - 0.5) / ((double) RAND_MAX);
				double rand2 = 2*(rand() - 0.5) / ((double) RAND_MAX);
				std::shared_ptr<network::Position> pt = std::make_shared<network::Position>(this->get_minpos() + 
					                   network::Position((i + 0.01*rand0)*h, (j + 0.01*rand1)*h, (k + 0.01*rand2)*h));
				size_t nl = 0;
				bool not_found = true;
				while(not_found && nl < this->count_lobes())
				{
					if(this->get_lobe(nl)->is_in(*pt))
					{
						point_clouds[nl].push_back(pt);
						not_found = false;
					}
					nl++;
				}
			}
		}
	}
	for(size_t nl = 0; nl < this->count_lobes(); nl++)
	{
		this->lobe_meshes[nl]->initialise_point_cloud(point_clouds[nl]);
	}

	size_t pi = 0;
	for(size_t n = 0; n < this->count_lobes(); n++)
	{
		cout << "Growing lobe " << this->get_lobe(n)->get_id() << endl;
		pi += this->tree->grow_lobe(this->lobe_meshes[n].get(), h, lobe_termnodes[n], this->options->get_option<bool>(PRUNE_TERM_BRANCHES_KEY)->get_value(),
			                  this->options->get_option<bool>(REASSIGN_KEY)->get_value(), this->params->get_param<double>(DISTANCE_TO_COM_FACTOR_KEY)->get_value(),
							  this->params->get_param<double>(MAX_BRANCHING_ANGLE_KEY)->get_value(), this->params->get_param<double>(LENGTH_CUTOFF_MM_KEY)->get_value(), pi);
	}
	this->tree->fill_in_radii();
}

size_t LungSegNetwork::grow_lobe(LobeMesh* mesh, const double & cloud_h, const std::vector<network::Node*> & tnodes, 
							   const bool & prune, const bool & reassign, const double & l_factor,
							   const double & max_angle, const double & l_cutoff, const size_t & prev_iter)
{
	//function for growing airways in a signle lobe
	using namespace network;

	//prune terminal edges as their lengths are unreliable
	std::vector<size_t> start_node;
	std::vector<Position> sn_position;
	std::cout << tnodes.size() << std::endl;
	start_node.resize(tnodes.size());
	sn_position.resize(tnodes.size());
	for(size_t kt = 0; kt < tnodes.size(); kt++)
	{
		size_t k = this->get_node_index(tnodes[kt]);
		start_node[kt] = k;
		
		if(prune)
		{
			size_t j = this->get_edge_in_index(k,0);
			size_t k_in = this->get_node_in_index(j);
			size_t j_in = this->get_edge_in_index(j,0);
			//change term node position -- cap to length of parent branch		
			Position node_in_pos = this->get_node(k_in)->get_pos();
			Position dir = (this->get_node(k)->get_pos() - node_in_pos);
			dir.normalise();
			if(this->get_edge(j)->get_geom()->get_length() > this->get_edge(j_in)->get_geom()->get_length()) //longer than parent -- cutoff
			{
				double new_length = this->get_edge(j_in)->get_geom()->get_length();
				this->get_node(k)->set_pos(node_in_pos + dir * new_length);
				this->get_edge(j)->update_length();
			}
		}
		sn_position[kt] = this->get_node(k)->get_pos();
	}

	PointCloud total_pt_cloud = *(mesh->pt_cloud);
	std::vector<std::shared_ptr<PointCloud>> sn_point_cloud;
	total_pt_cloud.partition_point_cloud(sn_position, sn_point_cloud);
	size_t iter = 0;
	while(start_node.size() > 0)
	{
		std::vector<size_t> next_start_node;
		std::vector<Position> next_sn_position;
		next_start_node.reserve(2*start_node.size());
		std::vector<std::shared_ptr<PointCloud>> next_point_clouds;
		if(!reassign) next_point_clouds.reserve(2*start_node.size());

		std::vector<Position> sn_parent_dir, sn_grandparent_dir;
		sn_parent_dir.resize(start_node.size());
		sn_grandparent_dir.resize(start_node.size());
		for(size_t isn = 0; isn < start_node.size(); isn++)
		{
			sn_position[isn] = this->get_node(start_node[isn])->get_pos();
			size_t k_parent = this->get_node_in_index(this->get_edge_in_index(start_node[isn],0));
			sn_parent_dir[isn] = sn_position[isn] - this->get_node(k_parent)->get_pos();
			sn_parent_dir[isn].normalise();
			size_t k_grandparent = this->get_node_in_index(this->get_edge_in_index(k_parent,0));
			sn_grandparent_dir[isn] = this->get_node(k_parent)->get_pos() - this->get_node(k_grandparent)->get_pos();
			sn_grandparent_dir[isn].normalise();
		}

		for(size_t isn = 0; isn < start_node.size(); isn++)
		{
			Position com_dir = sn_point_cloud[isn]->get_com() - sn_position[isn];
			com_dir.normalise();
			Position plane_norm;
			if(reassign && com_dir.distance(sn_parent_dir[isn]) > 1E-09)
			{
				plane_norm = com_dir.cross(sn_parent_dir[isn]);
			}
			else plane_norm = com_dir.cross(sn_grandparent_dir[isn]);
			plane_norm.normalise();

			std::vector<std::shared_ptr<PointCloud>> split_cloud;
			sn_point_cloud[isn]->split_point_cloud(sn_point_cloud[isn]->get_com(), plane_norm, split_cloud);
			for(size_t ic = 0; ic < split_cloud.size(); ic++)
			{

				Position new_com_vect = split_cloud[ic]->get_com() - sn_position[isn];
				//add new node in com direction
				std::shared_ptr<Node> no = std::make_shared<Node>(sn_position[isn] +  new_com_vect * l_factor);
				this->NodeVec.push_back(no);
				//add connecting edge
				std::shared_ptr<Edge<Node>> ed = std::make_shared<Edge<Node>>(this->NodeVec[start_node[isn]], this->NodeVec.back(), 1, -1.0);
				this->EdgeVec.push_back(ed);   //set rad at -1.0 as flag
				//check branching angle
				double br = new_com_vect.angle(sn_parent_dir[isn]);
				if(this->NodeVec.back()->get_pos().x[0] != this->NodeVec.back()->get_pos().x[0])
				{
					std::cout << "Pause\n";
				}
				if(br*br > max_angle*max_angle)
				{
					Position old_vec = this->NodeVec.back()->get_pos() - sn_position[isn];
					Position new_vec = new_angle_calc(old_vec, sn_parent_dir[isn], max_angle);
					this->NodeVec.back()->set_pos(sn_position[isn] + new_vec);
				}

				//check it is inside
				if(!mesh->is_in(this->NodeVec.back()->get_pos()))
				{
					Position old_vec = this->NodeVec.back()->get_pos() - sn_position[isn];
					double dl = 0.01;
					while(!mesh->is_in(this->NodeVec.back()->get_pos()) && dl < 0.995)  //if not shorted it
					{
						Position new_vec = old_vec*(1 - dl);
						this->NodeVec.back()->set_pos(sn_position[isn] + new_vec);
						dl += 0.01;
					}
					this->EdgeVec.back()->update_length();
				}

				//determine if terminal
				if(this->EdgeVec.back()->get_geom()->get_length() >= l_cutoff && (reassign || split_cloud[ic]->point_count() > 1))
				{
					next_start_node.push_back(this->NodeVec.size()-1);
					next_sn_position.push_back(this->NodeVec.back()->get_pos());
					if(!reassign) next_point_clouds.push_back(split_cloud[ic]);
				}
				else total_pt_cloud.remove_closest_point(sn_position[isn]);
			}
		}
		start_node = next_start_node;
		sn_position = next_sn_position;
		//if reassigning -> next pt_clouds come from partition algorithm
		//otherwise -> next point clouds are those from splitting algorithm (already filled)
		if(reassign) 
		{
			total_pt_cloud.partition_point_cloud(sn_position, sn_point_cloud);
			//if using reassign option it is possible that a cloud can have < 2 pts -> make these terminal 
			bool any_removed = true;
			while(any_removed)
			{
				any_removed = false;
				for(long int ipc = ((long int) sn_point_cloud.size())-1; ipc >= 0; ipc--)
				{
					if(sn_point_cloud[ipc]->point_count() < 2)
					{
						//any nodes with not enough points are term nodes -- remove points from cloud
						total_pt_cloud.remove_closest_point(sn_position[ipc]);
						start_node.erase(start_node.begin() + ipc);
						sn_position.erase(sn_position.begin() + ipc);
						if(!any_removed) any_removed = true;
					}
				}
				if(any_removed) total_pt_cloud.partition_point_cloud(sn_position, sn_point_cloud);
			}
		}
		else sn_point_cloud = next_point_clouds;

		this->update_node_edge_maps();

		//print step-by-step output
		//this->print_step_by_step(start_node, sn_point_cloud, prev_iter + iter, mesh->get_id(), cloud_h);
		iter++;
	}
	
	this->reorder_network();

	return iter;
}

void LungSegNetwork::print_step_by_step(const size_t & num) const
{
	std::stringstream ss;
	ss << "tree_iter" << num;
	this->print_vtk(ss.str().c_str(), 1.0);
}

void LungSegNetwork::fill_in_radii()
{
	//fill in all branch radii for airways descended from CT airways
	//loop over weibel orders, track order and radius of last ancestor that had one
	std::vector<size_t> parent_horsfield_gen;
	parent_horsfield_gen.resize(this->count_edges());
	std::fill(parent_horsfield_gen.begin(), parent_horsfield_gen.end(), 0);
	std::vector<double> parent_rad;
	parent_rad.resize(this->count_edges());
	std::fill(parent_rad.begin(), parent_rad.end(), 0);

	for(size_t wo = 0; wo < this->count_weibel_orders(); wo++)
	{
		for(size_t iedge = 0; iedge < this->count_edges_in_weibel_order(wo); iedge++)
		{
			size_t j = this->get_edge_index_from_weibel_order(wo,iedge);   //j is edge index
			if(this->get_edge(j)->get_geom()->get_inner_radius() > 100)
			{
					std::cout << j << ' ' << this->get_edge(j)->get_geom()->get_inner_radius() << ' ' 
				                  << this->get_weibel_order(j) << ' '
								  << this->get_horsfield_order(j) << '\n';
			}

			if(this->get_edge(j)->get_geom()->get_inner_radius() > 0)   //radius already defined
			{			
				for(size_t ieo = 0; ieo < this->count_edges_out(this->get_node_out_index(j)); ieo++)
				{
					//branch out index
					size_t jo = this->get_edge_out_index(this->get_node_out_index(j), ieo);
					//assign parent branches as this one
					parent_horsfield_gen[jo] = this->get_horsfield_order(j);
					parent_rad[jo] = this->get_edge(j)->get_geom()->get_inner_radius();
				}
			}
			else
			{
				this->get_edge(j)->update_geometry(calc_radius(parent_horsfield_gen[j], parent_rad[j], 
					                                           this->get_horsfield_order(j)),
											       this->get_edge(j)->get_geom()->get_length());
				for(size_t ieo = 0; ieo < this->count_edges_out(this->get_node_out_index(j)); ieo++)
				{
					//branch out index
					size_t jo = this->get_edge_out_index(this->get_node_out_index(j), ieo);
					//assign parent branches same as this one (i.e. last branch from CT)
					parent_horsfield_gen[jo] = parent_horsfield_gen[j];
					parent_rad[jo] = parent_rad[j];
				}

			}

		}
	}

}

void PointCloud::partition_point_cloud(const std::vector<network::Position> & pos_to_be_assigned, 
									   std::vector<std::shared_ptr<PointCloud>> & partitioned_cloud) const
{
	//assign points in cloud to nearest point in the list "pos_to_be_assigned"
	using namespace std;

	vector<vector<std::shared_ptr<network::Position>>> ptd_point_cloud;
	ptd_point_cloud.resize(pos_to_be_assigned.size());
	
	for(size_t icloud = 0; icloud < this->cloud.size(); icloud++)
	{
		double mindist;
		size_t ipt_closest;
		for(size_t ipt = 0; ipt < pos_to_be_assigned.size(); ipt++)
		{
			double dist = this->cloud[icloud]->distance(pos_to_be_assigned[ipt]);
			if(ipt == 0 || dist < mindist)
			{
				mindist = dist;
				ipt_closest = ipt;
			}
		}
		ptd_point_cloud[ipt_closest].push_back(this->cloud[icloud]);
	}

	partitioned_cloud.resize(pos_to_be_assigned.size());
	for(size_t ipt = 0; ipt < pos_to_be_assigned.size(); ipt++)
	{
		partitioned_cloud[ipt] = std::make_shared<PointCloud>(ptd_point_cloud[ipt]);
	}
}

void PointCloud::split_point_cloud(const network::Position & plane_pt, const network::Position & plane_norm,
								   std::vector<std::shared_ptr<PointCloud>> & partitioned_cloud) const
{
	//Split the point cloud in two by the plane given
	using namespace std;

	vector<vector<std::shared_ptr<network::Position>>> split_pt_cloud;
	split_pt_cloud.resize(2);
	if(this->cloud.size() == 2)
	{
		split_pt_cloud[0].push_back(this->cloud[0]);
		split_pt_cloud[1].push_back(this->cloud[1]);
	}
	else
	{
		for(size_t icloud = 0; icloud < this->cloud.size(); icloud++)
		{
			network::Position dir = *(this->cloud[icloud]) - plane_pt;
			if(dir.dot(plane_norm) >= 0) split_pt_cloud[0].push_back(this->cloud[icloud]);
			else split_pt_cloud[1].push_back(this->cloud[icloud]);
		}
	}

	partitioned_cloud.resize(2);
	for(size_t ipt = 0; ipt < 2; ipt++)
	{
		partitioned_cloud[ipt] = std::make_shared<PointCloud>(split_pt_cloud[ipt]);
	}

	if(partitioned_cloud[0]->point_count() == 0 || partitioned_cloud[1]->point_count() == 0)
	{
		std::cout << "Pause\n";
	}
}

network::Position new_angle_calc(const network::Position & prev_vec, const network::Position & parent_dir, const double & max_angle)
{
	//rotate vector so that its angle with parent_dir is equal to the max_angle
	using namespace network;

	double theta = max_angle - prev_vec.angle(parent_dir);   //rotation angle (negative)
	Position u = parent_dir.cross(prev_vec);   //rotation axis
	u.normalise();
	double ct = cos(theta);
	double st = sin(theta);
	double pxnew = (ct + u.x[0]*u.x[0]*(1 - ct))*prev_vec.x[0] + (u.x[0]*u.x[1]*(1 - ct) - u.x[2]*st)*prev_vec.x[1] + (u.x[0]*u.x[2]*(1 - ct) + u.x[1]*st)*prev_vec.x[2];
	double pynew = (u.x[1]*u.x[0]*(1 - ct) + u.x[2]*st)*prev_vec.x[0] + (ct + u.x[1]*u.x[1]*(1 - ct))*prev_vec.x[1] + (u.x[1]*u.x[2]*(1 - ct) - u.x[0]*st)*prev_vec.x[2];
	double pznew = (u.x[2]*u.x[0]*(1 - ct) - u.x[1]*st)*prev_vec.x[0] + (u.x[2]*u.x[1]*(1 - ct) + u.x[0]*st)*prev_vec.x[1] + (ct + u.x[2]*u.x[2]*(1 - ct))*prev_vec.x[2];

	Position pnew = Position(pxnew, pynew, pznew);
	double check = pnew.angle(parent_dir);

	return pnew;
}

double calc_radius(const size_t & pgen, const double & prad, const size_t & tgen)
{
	//radius scaling rule: function of ancestor generation and radius (from CT) and gen number here (tgen)
	return (prad*pow(1.15, double(tgen) - double(pgen)));
}