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
#include <string>

using namespace PhysBAM;

template<class T> void ReadObj( std::istream* input, ARRAY<VECTOR<int,3> >& elements,  ARRAY<VECTOR<T,3> >& X )
{

    std::istream& in = (*input);

    std::string ind_type;
    
    T max_x=0.00000001, min_x=10000000;
    T max_y=0.00000001, min_y=10000000;
    T max_z=0.00000001, min_z=10000000;

    while( !in.eof() )
        {
            in >> ind_type;

            if( ind_type == std::string("f") )
                {
                    int v1, v2, v3;
                    in >> v1 >> v2 >> v3;
                    elements.Append( VECTOR<int,3>(v1, v2, v3) );
                }
            else if( ind_type == std::string("v") )
                {
                    T p1, p2, p3;
                    in >> p1 >> p2 >> p3;

                    if( max_x < p1 )
                        max_x = p1;

                    if( min_x > p1 )
                        min_x = p1;

                    if( max_y < p2 )
                        max_y = p2;

                    if( min_y > p2 )
                        min_y = p2;

                    if( max_z < p3 )
                        max_z = p3;

                    if( min_z > p3 )
                        min_z = p3;


                    X.Append( VECTOR<T,3>(p1, p2, p3) );
                }
            else if( ind_type == std::string("vn") )
                {
                    T p1, p2, p3;
                    in >> p1 >> p2 >> p3;
                }
            else
                LOG::cout << "Unrecognized field character " << ind_type << std::endl;
            
        }

    LOG::cout << "Max X : " << max_x << "    Min X : " << min_x << std::endl;
    LOG::cout << "Max Y : " << max_y << "    Min Y : " << min_y << std::endl;
    LOG::cout << "Max Z : " << max_z << "    Min Z : " << min_z << std::endl;


}


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

    std::istream* input=new std::ifstream(input_filename.c_str(),std::ios::in);

    ReadObj(input, elements, X );

    m = X.m;
    n = elements.m;

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
