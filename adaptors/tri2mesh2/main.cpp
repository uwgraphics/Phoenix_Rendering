//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <fstream>

using namespace PhysBAM;


template<class T> void
Write_Mesh2_File(const TRIANGULATED_SURFACE<T>& tri_surface,const std::string& filename,const std::string& vectors_filename,const std::string& indices_filename)
{
    std::ofstream output(filename.c_str());
    if(!output.is_open())
        PHYSBAM_FATAL_ERROR("Could not open file "+filename+" for writing");
   
    bool using_uvs = false;

    std::ifstream vectors_in(vectors_filename.c_str());
    if(!vectors_in.is_open())
        {
            //PHYSBAM_WARNING("Could not open vector file '"+vectors_filename+"'. Assuming no uv vectors.");
        }
    std::ifstream indices_in(indices_filename.c_str());
    if(!indices_in.is_open())
        {
            //PHYSBAM_WARNING("Could not open indices file '"+indices_filename+"'. Assuming no uv indices");  
        }
    static const int d1=3;
    static const int d2=3;

	int index_min = 1;
	int index_max = tri_surface.particles.array_collection->Size();

		//Write Mesh2 Data

	output << "mesh2 { " <<  std::endl; 

		// Write Vertex Data
	output << "vertex_vectors { "<< std::endl;
	output << index_max << ", " << std::endl;

    for(int p=1;p<=tri_surface.particles.array_collection->Size();p++)
		{
			output<<"<";
			for(int i=1;i<=d1;i++)
				{
					output<<tri_surface.particles.X(p)(i);
					if(i!=d1) output << ",";
				}
				
			if( p != tri_surface.particles.array_collection->Size())
				output<<">, "<<std::endl;
			else
				output<<"> "<<std::endl;
				
		}

	output << "}" << std::endl;

		// Done with Verts

    // Check for UV vectors
    if(vectors_in.is_open())
        {
            using_uvs = true;
            std::string line;
            while( !vectors_in.eof() )
                {
                    getline(vectors_in, line, '\n');
                    output << line << std::endl;
                }
        }
    // Done with UV vectors

		// Write Face Data
	
	output << "face_indices {" << std::endl;
	output << tri_surface.mesh.elements.m <<", "<< std::endl;


    for(int e=1;e<=tri_surface.mesh.elements.m;e++)
		{
			output<<"<";
			for(int i=1;i<=d2;i++)
				{
					output<<tri_surface.mesh.elements(e)(i)-1;
					if(i!=d2) output << ",";
				}
			
			if( e != tri_surface.mesh.elements.m)
				output<<">, "<<std::endl;
			else
				output<<"> "<<std::endl;
		}

	output << "}" << std::endl;

		// Done with Faces


    // Check for UV indices
    if(indices_in.is_open())
        {
            using_uvs = true;
            std::string line;
            while( !indices_in.eof() )
                {
                    getline(indices_in, line, '\n');
                    output << line << std::endl;
                }
        }
    // Done with UV indices

    if( using_uvs )
        output << "uv_mapping" << std::endl;

	output << "}" << std::endl;
	
		// Done with Mesh2

    output.close();
    
}


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    static const int d=3;
    typedef VECTOR<T,d> TV;

		//LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","path","Input .tet filename");
    parse_args.Add_String_Argument("-output","","path","Output .tri filename");

    parse_args.Add_String_Argument("-uv_vectors","","path","Additional data for UVs");
    parse_args.Add_String_Argument("-uv_indices","","path","Additional data for UVs");

    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);

    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");

    std::string vectors_filename=parse_args.Get_String_Value("-uv_vectors");
    std::string indices_filename=parse_args.Get_String_Value("-uv_indices");


    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
	FILE_UTILITIES::Read_From_File(stream_type,input_filename,triangulated_surface);

	Write_Mesh2_File(triangulated_surface,output_filename,vectors_filename,indices_filename);

	delete &triangulated_surface;
 
		//LOG::Finish_Logging();

    return 0;
}
//#####################################################################
