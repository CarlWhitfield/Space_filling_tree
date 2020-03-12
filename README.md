# Space_filling_tree

Command-line based tool to generate an airway tree given lobe surface meshes (.stl format) and CT-based centreline model (output from PulmonaryToolkit (PTK) by Tom Doel). To run the executable from cmd terminal in Windows the following command-line arguments need to be provided in any order
- Paths to the .txt files containing the PTK output
- Paths to the .stl files containing the lobe meshes
- (Optional) Paths to .options and .params files containing lists of options and parameters for the algorithm (see examples). If not given, default options and parameters are used. 
