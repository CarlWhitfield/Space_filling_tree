#ifndef DEFINE_SEG_CODES_H
#define DEFINE_SEG_CODES_H

#define PRUNE_TERM_BRANCHES_KEY "Prune"   //whether to prune terminal branches of original CT image and change length
#define REASSIGN_KEY "Reassign"  //after each iteration: reassign points to terminal nodes of tree

#define ACINI_COUNT_KEY "Total_acini"     
#define DISTANCE_TO_COM_FACTOR_KEY "Distance_factor"    //each new branch goes x*distance to c-o-m of assigned point cloud where x is distance factor 0<x<1
#define MAX_BRANCHING_ANGLE_KEY "Max_branching_angle"
#define LENGTH_CUTOFF_MM_KEY "Length_cutoff"    //terminal branch size cutoff -> terminates when it reaches this value
#define NUMBER_MESH_BINS_KEY "Bins"     //affects speed of computation - number of bins to split point clouds into

#endif