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

using namespace PhysBAM;

template<class T> void ReadObj( std::istream* input, ARRAY<VECTOR<int,4> >& elements,  ARRAY<VECTOR<T,3> >& X )
{
    std::istream& in = (*input);

    std::string key;
    int value;

    // Read in header bits...
    in >> key;
    PHYSBAM_ASSERT( key == "MeshVersionFormatted" );
    in >> value;
    PHYSBAM_ASSERT( value == 1 );

    in >> key;
    PHYSBAM_ASSERT( key == "Dimension" );
    in >> value;
    PHYSBAM_ASSERT( value == 3 );

    while( !in.eof() )
    {
            in >> key;
            in >> value;

            if( key == "Vertices" ){
                for( int i = 0; i < value; i++ ){
                    T x, y, z, r;
                    in >> x >> y >> z >> r;
                    X.Append( VECTOR<T,3>( x, y, z ) );
                }
            }
            else if( key == "Tetrahedra" ){
                for( int i = 0; i < value; i++ ){
                    int a, b, c, d, r;
                    in >> a >> b >> c >> d >> r;
                    elements.Append( VECTOR<int,4>( a, b, c, d ) );
                }
            }
            else if( key == "Edges") {
                for( int i = 0; i < value; i++ ){
                    T a, b, r;
                    in >> a >> b >> r;
                }
            } 
            else if( key == "Triangles") {
                for( int i = 0; i < value; i++ ){
                    T a, b, c, r;
                    in >> a >> b >> c >> r;
                }
            }  
            else if( key == "Quadrilaterals") {
                for( int i = 0; i < value; i++ ){
                    T a, b, c, d, r;
                    in >> a >> b >> c >> d >> r;
                }
            }            
            else if( key == "Hexaedra") {
                for( int i = 0; i < value; i++ ){
                    T a, b, c, d, e, f, g, h, r;
                    in >> a >> b >> c >> d >> e >> f >> g >> h >> r;
                }
            }  
            else if( key == "End" )
                break;
            else
                PHYSBAM_FATAL_ERROR( "Unknown Field.");
    }

}


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","file","Input .mesh filename");
    parse_args.Add_String_Argument("-output","","file","Input .tet filename");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);
    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");

    ARRAY<VECTOR<int,4> > elements;
    ARRAY<VECTOR<T,3> > X;
    int n,m;

    std::istream* input=new std::ifstream(input_filename.c_str(),std::ios::in);

    ReadObj(input, elements, X );

    m = X.m;
    n = elements.m;

    delete input;

    LOG::cout<<"Number of tetrahedrons = "<<n<<std::endl;
    LOG::cout<<"Number of particles = "<<m<<std::endl;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    tetrahedralized_volume.mesh.Initialize_Mesh(m,elements);
    tetrahedralized_volume.particles.array_collection->Resize(m);
    tetrahedralized_volume.particles.X=X;

    FILE_UTILITIES::Write_To_File(stream_type,output_filename,tetrahedralized_volume);


    delete &tetrahedralized_volume;

    LOG::Finish_Logging();

    return 0;
}
//#####################################################################
