//main function for code to generate tree from segmented CT data
#include"lung_segmentation.h"

int main(int argc, char *argv[])
{
	LungSegmentation seg(argc, argv);   //construct simulation environment from command line arguments
	
	seg.build_tree();
	
	if(seg.tree->print_vtk("test", 1.0)) return 1;

	return 0;
}

