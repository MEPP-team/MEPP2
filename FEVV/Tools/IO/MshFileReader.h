// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

/*
 * gmsh ref:
 * - http://geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
 */

#include <iostream>

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"

//Comment this one if you want a strict usage of msh file with only one Nodes and one Elements block.
#define USE_MULTI_BLOC

namespace FEVV {
namespace IO {


using namespace StrUtils;
using namespace FileUtils;	
	//Return the dimension of an element type
  template< typename IndexType >
	inline IndexType dimensions(IndexType type)
	{
		//If it is a 2D element
		if (type==2 || type==3 || type==9)
			return 2;
		//If it is a 3D element
		else if (type==4 || type==5)
			return 3;
		//If we don't find a known element we return 0
		return 0;
	}

	//Get every node of an element
  template< typename IndexType >
	std::vector<IndexType> get_index(char* buffer, IndexType d)
	{
		//We split the line at every whitespace
		char* token = strtok(buffer, " ");
		token = strtok(NULL, " ");
		//We get the type of the element
		unsigned long type = static_cast<unsigned long>(atol(token));
		unsigned long temp;
		//If the dimensions of the element is unknown, meaning that the element is not implemented yet
		if (dimensions(type)==0)
			throw "MSH Reader: Error there is an unknown element or a not implemented yet element";
		//If the dimensions of the element is different from the other elements
		if (dimensions(type)!=d)
			throw "MSH Reader: Error there is 2D elements with 3D Elements";

		//We get the number of nodes composing the element 
		switch(type)
		{
			case 2:
				temp = 3;
				break;
			case 3:
			case 4:
				temp = 4;
				break;
			case 9:
				temp = 6;
				break;
			case 5:
				temp = 8;
				break;
			default:
				temp=0;
		}
		//We get the next element of the buffer
		token = strtok(NULL, " ");
		std::vector<IndexType> indexes;
		//We convert it to a long and we got the number of tags of the element
		unsigned long nb = atol(token);
		//If there is a wrong number of tag
		if (nb > 4 || nb < 2)
			throw "MSH Reader: Wrong number of tags, it should be between 2 and 4 included";
		//We loop until the end of the line
		for(unsigned long i(0); i<nb+temp; i++)
		{
			token = strtok(NULL, " ");
			//If we finished the tags we push back in the vector the value of the node
			if (i>=nb)
				indexes.push_back(static_cast<IndexType>(atol(token)-1));
		}
		return indexes;
	}

	//Initialize vectors and the dimension of the file
	//It will read the file and get the number of nodes/elements to
	//resize the vectors and get the dimension of the file
template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >   	
	bool init_vectors(	const std::string &file_path,
						std::vector< std::vector<CoordType> >& points_coords,
						std::vector< std::vector<CoordNType> >& normals_coords,
						std::vector< std::vector<CoordCType> >& vertex_color_coords,
						std::vector< std::vector<IndexType> >& line_indices,
						std::vector< std::vector<CoordCType> >& lines_color_coords,
						std::vector< std::vector<IndexType> >& face_indices,
						std::vector< std::vector<CoordCType> >& face_color_coords,						
            IndexType& dim)
	{
		FILE* file;

		file = fopen(file_path.c_str(), "r");
		if (file==NULL)
		{
			std::cerr << "MSH Reader : Unable to open " << file_path << std::endl;
			return false;
		}

		long nb_points(0);
    long nb_p_color(0);
    long nb_p_normal(0);
    long nb_e_color(0);
    long nb_elements(0);
		int nb_param_reads;
		long number;

#ifndef USE_MULTI_BLOC
		bool found_multi_nodes = false;
		bool found_multi_elements = false;
#endif
		const unsigned int T_MAX = 256;
		char buffer[T_MAX];

		while(fgets(buffer, T_MAX, file) != NULL)
		{
			//If we found a Node block
			if (strcmp(buffer, "$Nodes\n")==0)
			{
#ifndef USE_MULTI_BLOC
				if (!found_multi_nodes)
					found_multi_nodes = true;
				else
				{
					std::cerr << "You can't have multiple Nodes blocks" << std::endl;
					return false;
				}
#endif
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file) == NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndNodes missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				//We get the number of nodes declared in this block
				nb_param_reads = sscanf(buffer, "%ld", &number);
				//If there is a wrong declaration
				if (nb_param_reads < 1)
				{
					std::cerr << "MSH Reader: Error in the file structure ($Nodes tag with no number associated)" << std::endl;
					fclose(file);
					return false;
				}
				//We add the number of points found in this block to the total amount
				nb_points += number;
			}

			//If we found a Element block
			else if (strcmp(buffer, "$Elements\n")==0)
			{
#ifndef USE_MULTI_BLOC
				if (!found_multi_elements)
					found_multi_elements = true;
				else
				{
					std::cerr << "You can't have multiple Nodes blocks" << std::endl;
					return false;
				}
#endif
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file) == NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndElements missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				//We get the number of elements declared in this block
				nb_param_reads = sscanf(buffer, "%ld", &number);
				//If there is a wrong declaration
				if (nb_param_reads < 1)
				{
					std::cerr << "MSH Reader: Error in the file structure ($Nodes tag with no number associated)" << std::endl;
					fclose(file);
					return false;
				}
				//We add the nomber of elements declared in this block to the total amount of elements
				nb_elements += number;
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file) == NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndElements missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				//If the dimension hasn't been changed until the beginning
				if (dim==0)
				{
					//We read the type of the element
					//The first time we use number variables but we don't need the first value
					nb_param_reads = sscanf(buffer, "%ld %ld", &number, &number);
					if (nb_param_reads != 2)
					{
						std::cerr << "MSH Reader: Error in the element declaration" << std::endl;
						fclose(file);
						return false;
					}
					//We get the dimension of the element
					dim = dimensions(number);
					//If we get a 0 then there is an unexpected element and we explain the error
					if (dim==0)
					{
						std::cerr << "MSH Reader: Error the element is not recognized or not supported yet" << std::endl;
						fclose(file);
						return false;
					}
				}
			}

			//If we found a data block
			else if (!(strcmp(buffer, "$NodeData\n")!=0 && strcmp(buffer, "$ElementData\n")!=0))
			{
				const char* types[] = {"NodeData", "ElementData"};
				unsigned short type;
				//If it is a NodeData block
				if (strcmp(buffer, "$NodeData\n")==0)
				{
					type = 0;
				}
				//Otherwise it is an ElementData block
				else
				{
					type = 1;
				}
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file) == NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($End" << types[type] << " missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				//We get the number of string tags refering to this block
				nb_param_reads = sscanf(buffer, "%ld", &number);
				//If the file is bad structured
				if (nb_param_reads < 1)
				{
					std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " tag with no string tag associated)" << std::endl;
					fclose(file);
					return false;
				}
				bool found_color = false;
				bool found_normal = false;
				//We loop through every string tags searching for a color or a normal tag
				for(int i=0; i<number; i++)
				{
					//If the file is bad structured
					if (fgets(buffer, T_MAX, file) == NULL)
					{
						std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " with wrong number of string datas)" << std::endl;
						fclose(file);
						return false;
					}
					//If we found a color tag
					if (strcmp(buffer, "\"Color\"\n")==0 ||
						strcmp(buffer, "\"color\"\n")==0 ||
						strcmp(buffer, "\"COLOR\"\n")==0 ||
						strcmp(buffer, "\"Colour\"\n")==0 ||
						strcmp(buffer, "\"colour\"\n")==0 ||
						strcmp(buffer, "\"COLOUR\"\n")==0 ||
						strcmp(buffer, "\"Colors\"\n")==0 ||
						strcmp(buffer, "\"colors\"\n")==0 ||
						strcmp(buffer, "\"COLORS\"\n")==0 )
							found_color = true;
					//If we found a normal tag
					if (strcmp(buffer, "\"Normal\"\n")==0 ||
						strcmp(buffer, "\"normal\"\n")==0 ||
						strcmp(buffer, "\"NORMAL\"\n")==0 ||
						strcmp(buffer, "\"Normals\"\n")==0 ||
						strcmp(buffer, "\"normals\"\n")==0 ||
						strcmp(buffer, "\"NORMALS\"\n")==0)
							found_normal = true;
				}
				//We can't have both color and normal under the same tag
				if (found_color && found_normal)
				{
					std::cerr << "MSH Reader: Color and normal values declared in the same data block" << std::endl;
					fclose(file);
					return false;
				}
				//Then if one of those two are present
				if (found_color || found_normal)
				{						
					//If the file is bad structured
					if (fgets(buffer, T_MAX, file) == NULL)
					{
						std::cerr << "MSH Reader: Error in the file structure ($End" << types[type] << " missing ?)" << std::endl;
						fclose(file);
						return false;
					}
					//We get the number of real parameters
					nb_param_reads = sscanf(buffer, "%ld", &number);
					//If the file is bad structured
					if (nb_param_reads != 1)
					{
						std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " tag with no real tag associated)" << std::endl;
						fclose(file);
						return false;
					}
					//We loop through all real tags
					for(int i=0; i<number; i++)
					{
						//If the file is bad structured
						if (fgets(buffer, T_MAX, file) == NULL)
						{
							std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " with wrong number of real datas)" << std::endl;
							fclose(file);
							return false;
						}
					}
					//If the file is bad structured
					if (fgets(buffer, T_MAX, file) == NULL)
					{
						std::cerr << "MSH Reader: Error in the file structure ($End" << types[type] << " missing ?)" << std::endl;
						fclose(file);
						return false;
					}
					//We get the number of integer tags
					nb_param_reads = sscanf(buffer, "%ld", &number);
					//If the file is bad structured
					if (nb_param_reads != 1)
					{
						std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " tag with no integer tag associated)" << std::endl;
						fclose(file);
						return false;
					}
					//We loop throught the integer tags
					for(int i=0; i<number; i++)
					{
						//If the file is bad structured
						if (fgets(buffer, T_MAX, file) == NULL)
						{
							std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " with wrong number of real datas)" << std::endl;
							fclose(file);
							return false;
						}
					}
					//We get the number of datas declared in this block
					nb_param_reads = sscanf(buffer, "%ld", &number);
					//If the file is bad structured
					if (nb_param_reads != 1)
					{
						std::cerr << "MSH Reader: Error in the file structure ($" << types[type] << " tag with no integer tag associated)" << std::endl;
						fclose(file);
						return false;
					}
					//If we are in a NodeData block
					if (type==0)
					{
						//If we found a color string in this block
						if (found_color)
							nb_p_color += number;
						//Otherwise it is a normal data
						else
							nb_p_normal += number;
					}
					//If we are in an ElementData block
					else
						nb_e_color += number;
				}
			}
		}
		fclose(file);

		//If only a part of points or elements are colored or only a part of point have a normal
		if ((nb_points!=nb_p_color && nb_p_color!=0) || (nb_elements!=nb_e_color && nb_e_color!=0) || (nb_points!=nb_p_normal && nb_p_normal!=0))
		{
			std::cerr << "MSH Reader: Error in the file, all nodes or elements should have or not color" << std::endl;
			return false;
		}
		//We resize the vectors
		points_coords.resize(nb_points);
		normals_coords.resize(nb_p_normal);		
		vertex_color_coords.resize(nb_p_color);
		//If the dimension is 2 we resize only the vectors for the 2 dimensions elements
		if (dim==2)
		{
			face_indices.resize(nb_elements);
			face_color_coords.resize(nb_e_color);
		}
		//Otherwise we resize the vectors for the 3 dimensions elements
		else
		{
			line_indices.resize(nb_elements);
			lines_color_coords.resize(nb_e_color);
		}
		return true;
	}
   
   /**
   * Read gmsh (.msh) file and store mesh data in an intermediate vector representation.
   */
template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >   
	bool read_gmsh_file(const std::string &file_path,
						std::vector< std::vector<CoordType> >& points_coords,
						std::vector< std::vector<CoordNType> >& normals_coords,
						std::vector< std::vector<CoordCType> >& vertex_color_coords,
						std::vector< std::vector<IndexType> >& line_indices,
						std::vector< std::vector<CoordCType> >& lines_color_coords,
						std::vector< std::vector<IndexType> >& face_indices,
						std::vector< std::vector<CoordCType> >& face_color_coords,
						std::vector< std::vector < std::vector<double> > >& field_attributes,
						std::vector< std::string >& field_names)
	{
		//We check if the file has the good extension
		if (!FileUtils::has_extension(file_path, ".msh"))
		{
			std::cerr << "MSH Reader: File extension is inappropriate" << std::endl;
		}

		points_coords.clear();
		normals_coords.clear();
		vertex_color_coords.clear();		
		face_indices.clear();
		face_color_coords.clear();
		line_indices.clear();
		lines_color_coords.clear();
		field_names.clear();
		field_attributes.clear();

		//Dimension of the elements
    IndexType dim(0);

		if (!init_vectors(file_path, points_coords, normals_coords, vertex_color_coords, line_indices, lines_color_coords, face_indices, face_color_coords, dim))
		{
			std::cerr << "Error during initialization" << std::endl;
			return false;
		}

		const long T_MAX = 256;
		char buffer[T_MAX];
		long number_param;
		
		FILE* file = fopen(file_path.c_str(), "r");
		
		if (file==NULL)
		{
			std::cerr << "MSH Reader: Unable to open file : " << file_path << std::endl;
			return false;
		}

		while( fgets(buffer, T_MAX, file) != NULL)
		{
			//------------------$MeshFormat-------------------
			if (strcmp(buffer, "$MeshFormat\n")==0)
			{
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file) == NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndMeshFormat missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				float version;
				unsigned int type;
				unsigned int size;

				number_param = sscanf(buffer, "%f %u %u", &version, &type, &size);
				//Currently only the version 2.2 is supported by Gmsh
				if (version != 2.2f)
				{
					std::cerr << "MSH Reader: The " << version << " version is not supported yet" << std::endl;
					fclose(file);
					return false;
				}
				//We check if the type of the file is not ASCI (i.e : a type of 0 is an ASCI file)
				if (type!=0)
				{
					std::cerr << "MSH Reader: Only the ASCII file type is supported" << std::endl;
					fclose(file);
					return false;
				}
				//We check if the size of the elements is good or not
				//Actually Gmsh only supports sizeof(double)
				if (size != sizeof(double))
				{
					std::cerr << "MSH Reader: The size of the elements should be " << sizeof(double) << " and not " << size << std::endl;
					fclose(file);
					return false;
				}
				//We check that we got the $EndMeshFormat tag
				if (fgets(buffer, T_MAX, file) == NULL || strcmp(buffer, "$EndMeshFormat\n")!=0)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndMeshFormat missing)" << std::endl;
					fclose(file);
					return false;
				}
			}
			//--------------------$EndMeshFormat---------------------

			//--------------------$Nodes-----------------------------
			else if (strcmp(buffer, "$Nodes\n")==0)
			{
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file)==NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndNodes missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				long temp;
				std::vector<double> temp_pts;
				
				//We get the number of nodes declared in this block
				number_param = sscanf(buffer, "%ld", &temp);
				if (number_param != 1)
				{
					std::cerr << "MSH Reader: Error in the number of nodes" << std::endl;
					fclose(file);
					return false;
				}
				long index;
				//Then we add all the nodes on the block
				for(long i(0); i<temp; ++i)
				{
					if (fgets(buffer, T_MAX, file) == NULL)
					{
						std::cerr << "MSH Reader: Error in the file, no node associated to a $Node tag" << std::endl;
						fclose(file);
						return false;
					}
					temp_pts.clear();
					temp_pts.resize(3);
					//We get the three coordinates
					number_param = sscanf(buffer, "%ld %lf %lf %lf", &index, &temp_pts[0], &temp_pts[1], &temp_pts[2]);
					if (number_param != 4)
					{
						std::cerr << "MSH Reader: Error in the file, wrong node declaration" << std::endl;
						fclose(file);
						return false;
					}
					//Then we assign it to the vector points_coords if the index is correct
					//Don't forget the index-1 because msh format starts array at 1 and not 0
					if (index-1 < static_cast<long>(points_coords.size()))
						points_coords[index-1] = std::vector<CoordType>(temp_pts.begin(), temp_pts.end());
					else
					{
						std::cerr << "MSH Reader: Error in the node index, " << points_coords.size() << " declared and there is the " << index << " declared." << std::endl;
						fclose(file);
						return false;
					}
				}
				//We check that there is the end block
				if (fgets(buffer, T_MAX, file)==NULL || (strcmp(buffer, "$EndNodes")!=0 && strcmp(buffer, "$EndNodes\n")!=0))
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndNodes missing)" << std::endl;
					fclose(file);
					return false;
				}
			}
			//--------------------$EndNodes--------------------------

			//--------------------$Elements--------------------------
			else if (strcmp(buffer, "$Elements\n")==0)
			{
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file)==NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndElements missing ?)" << std::endl;
					fclose(file);
					return false;
				}
        long temp;
				//We get the number of elements declared in this block
				number_param = sscanf(buffer, "%ld", &temp);
				if (number_param != 1)
				{
					std::cerr << "MSH Reader: Error in the number of elements" << std::endl;
					fclose(file);
					return false;
				}
        long index;
				//We get every elements
				for(long i(0); i<temp; ++i)
				{
					if (fgets(buffer, T_MAX, file) == NULL)
					{
						std::cerr << "MSH Reader: Error in the file, no element associated to a $Elements tag" << std::endl;
						fclose(file);
						return false;
					}
					//We get the number of the element
					number_param = sscanf(buffer, "%ld", &index);
					//If the file is bad structured
					if (number_param!=1)
					{
						std::cerr << "MSH Reader: Error in the file, wrong declaration of an element" << std::endl;
						fclose(file);
						return false;
					}
					try
					{
						//We get the number of the nodes composing the element
						std::vector<IndexType> temp = get_index<IndexType>(buffer, dim);
						//If it is a 2D file we add it to the 2D elements vector
						if (dim==2)
						{
							//If the index of the element isn't too big for the size declared
							//The file should NOT skip an index
							if (index-1 < static_cast<long>(face_indices.size()))
								face_indices[index-1] = temp;
							//Otherwise we throw an error
							else
							{
								std::cerr << "MSH Reader: Error in the 2D elements index, " << face_indices.size() << " declared and there is the " << index << " declared." << std::endl;
								fclose(file);
								return false;
							}
						}
						//Otherwise we add it to the 3D elements vector
						else
						{
							//If the index of the element isn't too big for the size declared
							//The file should NOT skip an index
							if (index-1 < static_cast<long>(line_indices.size()))
								line_indices[index-1] = temp;
							//Otherwise we throw an error
							else
							{
								std::cerr << "MSH Reader: Error in the 2D elements index, " << line_indices.size() << " declared and there is the " << index << " declared." << std::endl;
								fclose(file);
								return false;
							}
						}
					}
					catch (char* e)
					{
						std::cerr << e << std::endl;
						fclose(file);
						return false;
					}
					
				}
				buffer[strlen(buffer)-1]='\0';
				//If the file is bad structured
				if (fgets(buffer, T_MAX, file)==NULL || (strcmp(buffer, "$EndElements")!=0 && strcmp(buffer, "$EndElements\n")!=0))
				{
					std::cerr << "MSH Reader: Error in the file structure ($EndElements missing)" << std::endl;
					fclose(file);
					return false;
				}
			}
			//--------------------$EndElements-----------------------

			//--------------------$Data--------------------------
			else if (!(strcmp(buffer, "$NodeData\n")!=0 && strcmp(buffer, "$ElementData\n")!=0))
			{
				const char* types[] = {"NodeData", "ElementData"};
				unsigned short type, data = 2;
				//If we are in a NodeData block
				if (strcmp(buffer, "$NodeData\n")==0)
					type = 0;
				//If we are in a ElementData block 
				else
					type = 1;

				//String parameters
				if (fgets(buffer, T_MAX, file)==NULL)
				{
					std::cerr << "MSH Reader: Error in the file structure ($End" << types[type] << " missing ?)" << std::endl;
					fclose(file);
					return false;
				}
				long params;
				//We get the number of string parameters
				number_param = sscanf(buffer, "%ld", &params);
				if (number_param != 1)
				{
					std::cerr << "MSH Reader: Error in the file, wrong" << types[type] << " declaration in the number of string parameters" << std::endl;
					fclose(file);
					return false;
				}
				//We loop through all the string paramaters
				for(long i(0); i<=params; ++i)
				{
					if (fgets(buffer, T_MAX, file)==NULL)
					{
						std::cerr << "MSH Reader: Error in the file, wrong number of string parameters in " << types[type] << std::endl;
						fclose(file);
						return false;
					}
					if (i==0)
					{
						if (strcmp(buffer, "\"Color\"\n")==0 ||
						strcmp(buffer, "\"color\"\n")==0 ||
						strcmp(buffer, "\"COLOR\"\n")==0 ||
						strcmp(buffer, "\"Colour\"\n")==0 ||
						strcmp(buffer, "\"colour\"\n")==0 ||
						strcmp(buffer, "\"COLOUR\"\n")==0 ||
						strcmp(buffer, "\"Colors\"\n")==0 ||
						strcmp(buffer, "\"colors\"\n")==0 ||
						strcmp(buffer, "\"COLORS\"\n")==0 )
						{
							data = 0;
						}
						else if (strcmp(buffer, "\"Normal\"\n")==0 ||
						strcmp(buffer, "\"normal\"\n")==0 ||
						strcmp(buffer, "\"NORMAL\"\n")==0 ||
						strcmp(buffer, "\"Normals\"\n")==0 ||
						strcmp(buffer, "\"normals\"\n")==0 ||
						strcmp(buffer, "\"NORMALS\"\n")==0)
						{
							data = 1;
						}
						else
						{
							data = 2;
                            buffer[strlen(buffer)-1] = '\0';
                            std::string a = (type==0)?"POINT_DATA_":"ELEMENT_DATA_";
                            a += std::string(buffer).substr(1, std::string(buffer).size()-2);
                            field_names.push_back(a);
						}
					}
				}

				//Real paramaters
				//We get the number of real parameters
				number_param = sscanf(buffer, "%ld", &params);
				if (number_param!=1)
				{
					std::cerr << "MSH Reader: Error in the file, wrong " << types[type] << " declaration in the number of real parameters" << std::endl;
					fclose(file);
					return false;
				}
				//We loop through all the real paramaters
				for(long i(0); i<=params; ++i)
				{
					if (fgets(buffer, T_MAX, file)==NULL)
					{
						std::cerr << "MSH Reader: Error in the file, wrong number of real parameters in " << types[type] << std::endl;
						fclose(file);
						return false;
					}
				}

				//Integer parameters
				//We get the number of integer parameters
				number_param = sscanf(buffer, "%ld", &params);
				if (number_param!=1)
				{
					std::cerr << "MSH Reader: Error in the file, wrong " << types[type] << " declaration in the number of integer parameters" << std::endl;
					fclose(file);
					return false;
				}
				//We loop until the last integer parameter (the number of datas in the block)
				for(long i(0); i<params; ++i)
				{
					if (fgets(buffer, T_MAX, file)==NULL)
					{
						std::cerr << "MSH Reader: Error in the file, wrong number of integer parameters in " << types[type] << std::endl;
						fclose(file);
						return false;
					}
				}
				//We get the number of datas declared in this block
				number_param = sscanf(buffer, "%ld", &params);
				if (number_param!=1)
				{
					std::cerr << "MSH Reader: Error in the file, wrong number of integer parameters in " << types[type] << std::endl;
					fclose(file);
					return false;
				}
				
				std::vector<CoordCType> temp_vector_c;
				std::vector<CoordNType> temp_vector_n;
				std::vector< std::vector<double> > temp_vector_f;
				temp_vector_f.clear();
				double temp1, temp2, temp3;
				//We loop through every declared datas
				for(long i(0); i<params; i++)
				{
					if (fgets(buffer, T_MAX, file)==NULL)
					{
						std::cerr << "MSH Reader: Error in the file, wrong number of integer parameters in " << types[type] << std::endl;
						fclose(file);
						return false;
					}
					switch(data)
					{
					case 0:
						temp_vector_c.clear();
						temp_vector_c.reserve(3);
						break;
					case 1:
						temp_vector_n.clear();
						temp_vector_n.reserve(3);
						break;						
					}
					long index;
					//We get all the values concerning the data (the index and the data itself)
					//If it is a color we get three values
					//Duplicate this if block, create temp4 et check for 5 number_param for color if you want to use alpha
					number_param = sscanf(buffer, "%ld %lf %lf %lf", &index, &temp1, &temp2, &temp3);
					if (number_param<1 || number_param>4)
					{
						std::cerr << "MSH Reader: Error in the file, wrong number of arguments for integer parameters in " << type << std::endl;
						fclose(file);
						return false;
					}
					--number_param;
					switch(data)
					{
					case 0:
						temp_vector_c.push_back(static_cast<CoordCType>(temp1));
						if (number_param > 1)
							temp_vector_c.push_back(static_cast<CoordCType>(temp2));
						if (number_param > 2)
							temp_vector_c.push_back(static_cast<CoordCType>(temp3));
						break;
					case 1:
						temp_vector_n.push_back(static_cast<CoordNType>(temp1));
						if (number_param > 1)
							temp_vector_n.push_back(static_cast<CoordNType>(temp2));
						if (number_param > 2)
							temp_vector_n.push_back(static_cast<CoordNType>(temp3));
						break;
					case 2:
						if (temp_vector_f.size() < index)
							temp_vector_f.resize(index);
						temp_vector_f[index-1].push_back(temp1);
						if (number_param > 1)
							temp_vector_f[index-1].push_back(temp2);
						if (number_param > 2)
							temp_vector_f[index-1].push_back(temp3);
					}

					//Si ce n'est pas des points mais des éléments
					if (type!=0)
					{
						//Si ils sont de dimension 2 on les met dans le vector face_indices
						if (dim==2)
						{
							if (data == 0)
							{
								//On vérifie qu'on ne sort pas de la taille du tableau
								if (index - 1 < static_cast<long>(face_color_coords.size()))
									face_color_coords[index-1] = temp_vector_c;
								else
								{
									std::cerr << "MSH Reader: Wrong index usage in 2D element color, " << index << " used when the size is " << face_color_coords.size() << std::endl << buffer << std::endl;
									fclose(file);
									return false;
								}
							}
						}
						//Si ils sont de dimension 3 on les met dans le vecteur line_indices
						else
						{
							if (data==0)
							{
								//On vérifie qu'on ne sort pas de la taille du tableau
								if (index - 1 < static_cast<long>(lines_color_coords.size()))
									lines_color_coords[index-1] = temp_vector_c;
								else
								{
									std::cerr << "MSH Reader: Wrong index usage in 3D element color, " << index << " used when the size is " << lines_color_coords.size() << std::endl;
									fclose(file);
									return false;
								}
							}
						}
					}

					//Si ce sont des points
					else
					{
						if (data==0)
						{
							//On vérifie qu'on ne sort pas de la taille du tableau
							if (index - 1 < static_cast<long>(vertex_color_coords.size()))
								vertex_color_coords[index-1] = temp_vector_c;
							else
							{
								std::cerr << "MSH Reader: Wrong index usage in node color, " << index << " used when the size is " << vertex_color_coords.size() << std::endl;
								fclose(file);
								return false;
							}
						}
						else if (data==1)
						{
							//On vérifie qu'on ne sort pas de la taille du tableau
							if (index - 1 < static_cast<long>(normals_coords.size()))
								normals_coords[index-1] = temp_vector_n;
							else
							{
								std::cerr << "MSH Reader: Wrong index usage in node normal, " << index << " used when the size is " << normals_coords.size() << std::endl;
								fclose(file);
								return false;
							}
						}
					}
				}
				if (!temp_vector_f.empty())
					field_attributes.push_back(temp_vector_f);
				//We create a string with the end tag
				char s[20] = "$End";
				strcat(s, types[type]);
				//Then if the file is bad structured OR
				//If the line found is not an ending tag of the Data block
				//(i.e : $EndThingData or $EndThingData\n)
				if (fgets(buffer, T_MAX, file)==NULL || (strcmp(buffer, s)!=0 && strcmp(buffer, strcat(buffer, "\n"))!=0))
				{
					std::cerr << "MSH Reader: Error in the file, missing $End" << types[type] << std::endl;
					fclose(file);
					return false;
				}
			}
			//--------------------$EndData-----------------------
		}
    return true;
	}

}
}
