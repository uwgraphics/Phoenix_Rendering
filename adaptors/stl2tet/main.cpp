//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>

#include <fstream>
#include <string>
#include <sstream>

#include "tetgen.h"

using namespace PhysBAM;

template<class T> void Convert_To_Tetrahedralized_Volume( tetgenio& out, 
                                                          TETRAHEDRALIZED_VOLUME<T>& volume )
{
    ARRAY<VECTOR<T,3> > X;
    for( int i=0; i < out.numberofpoints; i++){
        VECTOR<T,3> p;
        p.x = out.pointlist[ i*out.mesh_dim + 0 ];
        p.y = out.pointlist[ i*out.mesh_dim + 1 ];
        p.z = out.pointlist[ i*out.mesh_dim + 2 ];
        X.Append(p);
        //LOG::cout << "V" << i+1 << " : " << p << std::endl;
    }

    ARRAY<VECTOR<int,4> > elements;
    for( int i=0; i < out.numberoftetrahedra; i++){
        VECTOR<int,4> e;
        e(1) = out.tetrahedronlist[ i*out.numberofcorners + 0 ];
        e(2) = out.tetrahedronlist[ i*out.numberofcorners + 1 ];
        e(3) = out.tetrahedronlist[ i*out.numberofcorners + 2 ];
        e(4) = out.tetrahedronlist[ i*out.numberofcorners + 3 ];
        elements.Append(e);
        //LOG::cout << "T" << i+1 << " : " << e << std::endl;
    }

    int m = X.m;
    int n = elements.m;

    volume.mesh.Initialize_Mesh(m,elements);
    volume.particles.array_collection->Resize(m);
    volume.particles.X=X;
}


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    tetgenbehavior b;
    tetgenio in, out, addin, bgmin;

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","file","Input .tri filename");
    parse_args.Add_String_Argument("-output","","file","Input .tri_raw filename");
    parse_args.Add_String_Argument("-tetgen_flags","VpqmR","flags","Optional flags for TetGen");

    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output")){
        try{ b.usage(); } catch(int e){};
        parse_args.Print_Usage(true);
    }
    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");
    std::string flags=parse_args.Get_String_Value("-tetgen_flags");

   
    if(!b.parse_commandline( const_cast<char *> (flags.c_str())) )
        PHYSBAM_FATAL_ERROR( "Could not parse tetgen command switches." );

    if (!in.load_plc(const_cast<char *> (input_filename.c_str()), 4 /*STL TYPE*/)) {
        terminatetetgen(NULL, 10);
    }

    tetrahedralize(&b, &in, &out, &addin, &bgmin);

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    Convert_To_Tetrahedralized_Volume( out, tetrahedralized_volume );
    FILE_UTILITIES::Write_To_File(stream_type,output_filename,tetrahedralized_volume);
    delete &tetrahedralized_volume;

    LOG::Finish_Logging();

    return 0;
}
//#####################################################################
