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

#include <iostream>
#include <cstdio>
#include <cstring>

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "FEVV/Tools/IO/StringUtilities.hpp"

namespace FEVV {
namespace IO {


using namespace StrUtils;
using namespace FileUtils;	

template< typename IndexType >
	inline IndexType get_type(const IndexType size, const IndexType dim)
	{
		switch(size)
		{
			case 3:
				return 2;
			case 4:
				return (dim==2)?3:4;
			case 6:
				return 9;
			case 8:
				return 5;
		}
		return 0;
	}

    /**
    * Write mesh data to gmsh (.msh) file.
    */
template< typename CoordType,
          typename CoordNType,
          typename CoordCType,
          typename IndexType >   	
	bool write_gmsh_file(	const std::string &file_path,
							const std::vector< std::vector<CoordType> >& points_coords,							
							const std::vector< std::vector<CoordNType> >& normals_coords,
							const std::vector< std::vector<CoordCType> >& vertex_color_coords,							
							const std::vector< std::vector<IndexType> >& line_indices,
							const std::vector< std::vector<CoordCType> >& lines_color_coords,
							const std::vector< std::vector<IndexType> >& face_indices,
							const std::vector< std::vector<CoordCType> >& face_color_coords,
							const std::vector< std::vector< std::vector<double> > >& field_attributes,
							const std::vector< std::string >& field_names)
	{
		//We check that the file path has the correct extension
		if(!(FileUtils::has_extension(file_path, ".msh")))
		{
			std::cerr << "MSH Writer: file extension is inappropriate" << std::endl;
			return false;
		}
		
		//We check that no node is colored or that every node are colored
		if (vertex_color_coords.size()!=0 && vertex_color_coords.size()!=points_coords.size())
		{
			std::cerr << "MSH Writer: Not every node are colored" << std::endl;
			return false;
		}
			
		//We check that no node has normal or that every node has normal
		if (normals_coords.size()>0 && normals_coords.size()!=points_coords.size()) 
		{
			std::cerr << "MSH Writer: Not every node are normalized" << std::endl;
			return false;
		}

		//We check that no element is colored or that every element is colored
		if (face_color_coords.size()!=0 && face_color_coords.size()!=face_indices.size())
		{
			std::cerr << "MSH Writer: Not every 2D element are colored" << std::endl;
			return false;
		}

		//We check that there is not 3D elements combined with 2D elements
		if (line_indices.size()!=0 && face_indices.size()!=0)
		{
			std::cerr << "MSH Writer: It is impossible to have 2D elements with 3D elements" << std::endl;
			return false;
		}
		std::vector<std::string> names;
		if (field_names.size() != field_attributes.size())
		{
			names.resize(field_attributes.size());
			// When data field names need to be set, we assume that nD (n>=2)
			// fields associated with positions are displacement fields while others
			// are understood as physical laws.
			for (unsigned long i(0); i<field_attributes.size(); ++i){
				if( (points_coords.size() == field_attributes[i].size()) && 
					((points_coords.size() != face_indices.size()) || (field_attributes[i][0].size()>1)))
					names[i] = std::string("Shifting_") + convert(i);  
				else
					names[i] = std::string("cell_law_")+ convert(i);
			}
		}
		else
			names = field_names;
    ///////////////////////////////////////////////////////////////////////////
		//If there is elements in the vector face_indices so it is the 2D vectors that we should use
		//Otherwise we should use the 3D vectors
    IndexType dim = (face_indices.size()) ? 2 : 3;
		const std::vector< std::vector<IndexType> >& elem = (face_indices.size()) ? face_indices : line_indices;
		const std::vector< std::vector<CoordCType> >& color = (face_indices.size()) ? face_color_coords : lines_color_coords;

		//We open the file
		FILE* file = fopen(file_path.c_str(), "w");
		if (file==NULL)
		{
			std::cerr << "MSH Reader: Unable to open file : " << file_path << std::endl;
			return false;
		}
		
		//We write the beginning of the file
		fprintf(file, "$MeshFormat\n2.2 0 %zd\n$EndMeshFormat\n$Nodes\n%zd\n", sizeof(double), points_coords.size());

		//If there is some points to write
		if (points_coords.size()!=0)
		{
			//We write every point while not forgetting the fact that index starts at 1 in msh file
			for(unsigned long i(0); i<points_coords.size(); i++)
			{
				fprintf(file, "%ld %.10g %.10g %.10g\n", i+1, static_cast<double>(points_coords[i][0]), static_cast<double>(points_coords[i][1]), static_cast<double>(points_coords[i][2]));
			}
			fprintf(file, "$EndNodes\n");

			//Then we write the element (at least points, msh file must at least have one Elements block)
			fprintf(file, "$Elements\n");
			//If there is elements different from points
			if (elem.size()!=0)
			{
				//Then we write the number of elements
				fprintf(file, "%zd\n", elem.size());
				//And we write all their attributes
        IndexType type;
				for(unsigned long i(0); i<elem.size(); i++)
				{
					type = get_type(static_cast<IndexType>(elem[i].size()), dim);
					//If we found an unknown element
					if (type==0)
					{
						std::cerr << "MSH Writer: Type of element is unknown" << std::endl;
						fclose(file);
						return false;
					}
					fprintf(file, "%ld %ld 2 99 2", i+1, static_cast<long>(type));
					for(unsigned long j(0); j<elem[i].size(); j++)
					{
						fprintf(file, " %s", convert(elem[i][j]+1).c_str());
					}
					fprintf(file, "\n");
				}
				fprintf(file, "$EndElements");
			}
			//If there is no element we write the points as elements
			else
			{
				fprintf(file, "%zd\n", points_coords.size());
				for(unsigned long i(0); i<points_coords.size(); i++)
				{
					fprintf(file, "%lu 15 2 99 2 %lu\n", i+1, i+1);
				}
				fprintf(file, "$EndElements");
			}

			//If there is point colors we write it
			if (vertex_color_coords.size()!=0)
			{
				fprintf(file, "\n$NodeData\n1\n\"Color\"\n0\n1\n%zd\n", vertex_color_coords.size());
				for(unsigned long i(0); i<vertex_color_coords.size(); i++)
				{
					if (vertex_color_coords[i].size()!=3)
					{
						std::cerr << "MSH Writer: Wrong number of values for declaration of colors" << std::endl;
						fclose(file);
						return false;
					}
					fprintf(file, "%ld %.10g %.10g %.10g\n", i+1, static_cast<double>(vertex_color_coords[i][0]), static_cast<double>(vertex_color_coords[i][1]), static_cast<double>(vertex_color_coords[i][2]));
				}
				fprintf(file, "$EndNodeData");
			}

			//If there is point normal we write it
			if (normals_coords.size()!=0)
			{
				fprintf(file, "\n$NodeData\n1\n\"Normal\"\n0\n1\n%zd\n", normals_coords.size());
				for(unsigned long i(0); i<normals_coords.size(); i++)
				{
					if (normals_coords[i].size()!=3)
					{
						std::cerr << "MSH Writer: Wrong number of values for declaration of normals" << std::endl;
						fclose(file);
						return false;
					}
					fprintf(file, "%ld %.10g %.10g %.10g\n", i+1, static_cast<double>(normals_coords[i][0]), static_cast<double>(normals_coords[i][1]), static_cast<double>(normals_coords[i][2]));
				}
				fprintf(file, "$EndNodeData");
			}

			//If there is element color we write it
			if (color.size()!=0)
			{
				fprintf(file, "\n$ElementData\n1\n\"Color\"\n0\n1\n%zd\n", color.size());
				for(unsigned long i(0); i<color.size(); i++)
				{
					if (color[i].size()!=3)
					{
						std::cerr << "MSH Writer: Wrong number of values for declaration of colors" << std::endl;
						fclose(file);
						return false;
					}
					fprintf(file, "%ld %.10g %.10g %.10g\n", i+1, static_cast<double>(color[i][0]), static_cast<double>(color[i][1]), static_cast<double>(color[i][2]));
				}
				fprintf(file, "$EndElementData");
			}

			if (names.size()>0)
			{
				bool isnode = true;
				for(unsigned long i(0); i<names.size(); ++i)
				{
					fprintf(file, "\n");
					if ( names[i].find("POINT_DATA_") != std::string::npos )
					{
						isnode = true;
						fprintf(file, "$NodeData\n1\n\"%s\"\n0\n1\n%zd\n", names[i].substr(11).c_str(), field_attributes[i].size()); 
					}
					else if ( (points_coords.size()==field_attributes[i].size()) && ( (points_coords.size()!=face_indices.size()) || (field_attributes[i][0].size()>1)) ) // the case of points_coords.size()=face_indices.size() is not well tackled, but for the time being we assume that in that case the attribute is a face attribute if only one value... [temporary solution; try to use the name of the field to decide?]
					{
						isnode = true;
						fprintf(file, "$NodeData\n1\n\"%s\"\n0\n1\n%zd\n", names[i].c_str(), field_attributes[i].size()); 
					}
					else if (names[i].find("ELEMENT_DATA_") != std::string::npos)
					{
						isnode = false;
						fprintf(file, "$ElementData\n1\n\"%s\"\n0\n1\n%zd\n", names[i].substr(13).c_str(), field_attributes[i].size()); 
					}
					else if (names[i].find("CELL_DATA_") != std::string::npos)
					{
						isnode = false;
						fprintf(file, "$ElementData\n1\n\"%s\"\n0\n1\n%zd\n", names[i].substr(10).c_str(), field_attributes[i].size()); 
					}					
					else
					{
						isnode = false;
						fprintf(file, "$ElementData\n1\n\"%s\"\n0\n1\n%zd\n", names[i].c_str(), field_attributes[i].size()); 
					}
					for(unsigned long j(0); j<field_attributes[i].size(); ++j)
					{
						fprintf(file, "%lu", j+1);
						for(unsigned long k(0); k<field_attributes[i][j].size(); ++k)
						{
							fprintf(file, " %.10g", field_attributes[i][j][k]);
						}
						fprintf(file, "\n");
					}
					if (isnode)
						fprintf(file, "$EndNodeData");
					else
						fprintf(file, "$EndElementData");
				}
			}
		}
		//If there is no point we write that there is nothing
		else
		{
			fprintf(file, "$EndNodes\n$Elements\n0\n$EndElements");
		}

		fclose(file);
		return true;
	}

}
}
