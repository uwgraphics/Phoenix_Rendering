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
    
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    FILE_UTILITIES::Read_From_File(stream_type,input_filename,tetrahedralized_volume);
    tetrahedralized_volume.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
    triangulated_surface.particles.array_collection->Append(*tetrahedralized_volume.particles.array_collection);
    triangulated_surface.mesh.Initialize_Mesh(tetrahedralized_volume.triangulated_surface->mesh);
    triangulated_surface.Discard_Valence_Zero_Particles_And_Renumber();
    FILE_UTILITIES::Write_To_File(stream_type,output_filename,triangulated_surface);
    delete &triangulated_surface;
    delete &tetrahedralized_volume;
 
		//LOG::Finish_Logging();

    return 0;
}
//#####################################################################
