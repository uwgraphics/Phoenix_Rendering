//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <fstream>

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","file","Input .tri filename");
    parse_args.Add_String_Argument("-output","","file","Input .tri_raw filename");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);
    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");

    ARRAY<VECTOR<int,3> > elements;
    ARRAY<VECTOR<T,3> > X;
    int n,m;

    std::istream* input=new std::ifstream(input_filename.c_str(),std::ios::in|std::ios::binary);
    Read_Binary<float>(*input,n,m);
    elements.Resize(n);
    X.Resize(m);
    Read_Binary_Array<float>(*input,&elements(1)(1),3*n);
    Read_Binary_Array<float>(*input,&X(1)(1),3*m);
    delete input;

    LOG::cout<<"Number of triangles = "<<n<<std::endl;
    LOG::cout<<"Number of particles = "<<m<<std::endl;

    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
    triangulated_surface.mesh.Initialize_Mesh(m,elements);
    triangulated_surface.particles.array_collection->Resize(m);
    triangulated_surface.particles.X=X;

    FILE_UTILITIES::Write_To_File(stream_type,output_filename,triangulated_surface);


    delete &triangulated_surface;

    LOG::Finish_Logging();

    return 0;
}
//#####################################################################
