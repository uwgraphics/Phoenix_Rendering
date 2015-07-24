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


template<class TV> void
Write_Obj_File(const SEGMENTED_CURVE<TV>& tri_surface,const std::string& filename)
{
    std::ofstream output(filename.c_str());
    if(!output.is_open())
        PHYSBAM_FATAL_ERROR("Could not open file "+filename+" for writing");

    static const int d1=3;
    static const int d2=2;

	int index_min = 1;
	int index_max = tri_surface.particles.array_collection->Size();

		//LOG::cout << "Min Index: " << index_min << std::endl;
		//LOG::cout << "Max Index: " << index_max << std::endl;
	

    for(int p=1;p<=tri_surface.particles.array_collection->Size();p++)
		{
			output<<"v";
			for(int i=1;i<=d1;i++)
				output<<" "<<tri_surface.particles.X(p)(i);
			if(d1==2) output<<" 0.0";
			output<<std::endl;
		}

    for(int e=1;e<=tri_surface.mesh.elements.m;e++)
		{
			output<<"f";
			for(int i=1;i<=d2;i++)
				{
					
						//PHYSBAM_ASSERT( tri_surface.mesh.elements(e)(i) >= index_min );
						//PHYSBAM_ASSERT( tri_surface.mesh.elements(e)(i) <= index_max );
					output<<" "<<tri_surface.mesh.elements(e)(i);
				}
			output<<std::endl;
		}
	
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
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);

    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");

    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create();
	FILE_UTILITIES::Read_From_File(stream_type,input_filename,segmented_curve);

	Write_Obj_File(segmented_curve,output_filename);

	delete &segmented_curve;
 
		//LOG::Finish_Logging();

    return 0;
}
//#####################################################################
