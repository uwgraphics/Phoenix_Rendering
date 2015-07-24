//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>

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

    in >> ind_type;

    while( !in.eof() )
        {
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
            else{
                //LOG::cout << "Unrecognized field character " << ind_type << std::endl;
            }

            in >> ind_type;            
        }

    //LOG::cout << "Max X : " << max_x << "    Min X : " << min_x << std::endl;
    //LOG::cout << "Max Y : " << max_y << "    Min Y : " << min_y << std::endl;
    //LOG::cout << "Max Z : " << max_z << "    Min Z : " << min_z << std::endl;


}

template<class T>
struct AnimationPointCloud{
    int num_frames;
    int num_points;
    ARRAY< VECTOR< T,3 > > rest_pose;
    ARRAY< ARRAY< VECTOR< T, 3> > > frames;
};


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","directory","Animation Directory");
    parse_args.Add_String_Argument("-output","","file","Output Pointcloud Path");
    parse_args.Add_Integer_Argument("-frames",1,"frames","Number of Frames to Process.");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);
    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");

    AnimationPointCloud<T> animation;

    ARRAY<VECTOR<int,3> > elements;
    ARRAY<VECTOR<T,3> > X;
    int n,m;
    std::istream* input;

    if(!FILE_UTILITIES::Directory_Exists(input_filename))
        PHYSBAM_FATAL_ERROR( "Animation Directory Doesn't Exist." );



    // Read in Rest Pose
    std::string rest_pose_filename = input_filename + "/" + "RestPose.obj";
    input=FILE_UTILITIES::Safe_Open_Input(rest_pose_filename);
    ReadObj(input, elements, X );
    m = X.m;
    n = elements.m;

    delete input;

    animation.num_frames = parse_args.Get_Integer_Value("-frames");
    animation.num_points = m;
    animation.rest_pose = X;

    LOG::cout << "Rest Pose: " << animation.num_points << " Points." << std::endl;

    // Read in animation frames
    animation.frames.Clean_Memory();
    for( int frame = 1; frame <= animation.num_frames; frame ++ ){
        LOG::cout << "Processing frame " << frame << std::endl;
        elements.Clean_Memory();
        X.Clean_Memory();

        std::string frame_filename = input_filename + "/" + STRING_UTILITIES::Value_To_String(frame) + ".obj";
        input=FILE_UTILITIES::Safe_Open_Input(frame_filename);
        ReadObj(input, elements, X );
        m = X.m;
        n = elements.m;
        
        PHYSBAM_ASSERT( m == animation.num_points );
        
        animation.frames.Append( X );

        delete input;
    }

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(output_filename);
    Write_Binary<T>(*output,animation.num_frames);
    Write_Binary<T>(*output,animation.num_points);
    Write_Binary<T>(*output,animation.rest_pose);
    Write_Binary<T>(*output,animation.frames);
    delete output;

    LOG::Finish_Logging();

    return 0;
}
//#####################################################################
