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
#include <PhysBAM_Tools/Matrices/FRAME.h>

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
    parse_args.Add_String_Argument("-inputA","","path","Input .tet filename A");
	parse_args.Add_String_Argument("-inputB","","path","Input .tet filename B");
    parse_args.Add_String_Argument("-output","","path","Output .tri filename");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-inputA") ||
	   !parse_args.Is_Value_Set("-inputB") ||
	   !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);

    std::string inputA_filename=parse_args.Get_String_Value("-inputA");
    std::string inputB_filename=parse_args.Get_String_Value("-inputB");
    std::string output_filename=parse_args.Get_String_Value("-output");

    TRIANGULATED_SURFACE<T>& triangulated_surfaceA=*TRIANGULATED_SURFACE<T>::Create();
	FILE_UTILITIES::Read_From_File(stream_type,inputA_filename,triangulated_surfaceA);

    TRIANGULATED_SURFACE<T>& triangulated_surfaceB=*TRIANGULATED_SURFACE<T>::Create();
	FILE_UTILITIES::Read_From_File(stream_type,inputB_filename,triangulated_surfaceB);

	
	ARRAY<TRIANGULATED_SURFACE<T> *> surf_array;
	surf_array.Resize(2);
	surf_array(1) = &triangulated_surfaceA;
	surf_array(2) = &triangulated_surfaceB;

	ARRAY<FRAME<TV> > frames;
	frames.Resize( 2 );

    TRIANGULATED_SURFACE<T>& triangulated_surfaceC=*TRIANGULATED_SURFACE<T>::Union_Mesh_Objects_Relatively(surf_array,frames);
	

	FILE_UTILITIES::Write_To_File(stream_type,output_filename,triangulated_surfaceC);

	delete &triangulated_surfaceA;
	delete &triangulated_surfaceB;
	delete &triangulated_surfaceC;
 
		//LOG::Finish_Logging();

    return 0;
}
//#####################################################################
