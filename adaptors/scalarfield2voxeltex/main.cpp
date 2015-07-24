//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>

#include <fstream>
#include <stdint.h> 

#include "minilzo/minilzo.h"

using namespace PhysBAM;

typedef union{
        char wide_block[4];        
        int  int_val;
        float float_val;
} CONVERTER;

void Write_Compressed( std::ofstream& outfile, char* data_in, int len_in, int mode ){
    CONVERTER cv;

    cv.wide_block[0] = mode;
    outfile.write( cv.wide_block, 1 );

    if( mode == 0 ){        
        cv.int_val = len_in;
        outfile.write( cv.wide_block, 4 );
        outfile.write( data_in, len_in );
    }

    if( mode == 1 ){
        char* compressed_data = new char[len_in];
        char* work_mem = new char[LZO1X_MEM_COMPRESS];
        unsigned long len_out;
        int status = lzo1x_1_compress( (unsigned char*)data_in,         len_in,
                                       (unsigned char*)compressed_data, &len_out,
                                       (unsigned char*)work_mem );
        PHYSBAM_ASSERT( status == 0 );
        cv.int_val = len_out;
        outfile.write( cv.wide_block, 4 );
        outfile.write( compressed_data, len_out );
        delete [] compressed_data;
        delete [] work_mem;
    }

}

void WriteBlenderCache( std::ofstream& outfile, int resolution, float dt, ARRAY<float> density, ARRAY<float> heat, bool use_heat ){

    CONVERTER cv;

    outfile << "BPHYSICS";
    cv.int_val = 3;
    outfile.write( cv.wide_block, 4 );
    cv.int_val = resolution*resolution*resolution;
    outfile.write( cv.wide_block, 4 );
    cv.int_val = 2;
    outfile.write( cv.wide_block, 4 );
    outfile << "1.04";
    
    if( use_heat ){
        cv.int_val = 1;
        outfile.write( cv.wide_block, 4 );
        cv.int_val = 9;
        outfile.write( cv.wide_block, 4 );
    }
    else{
        cv.int_val = 0;
        outfile.write( cv.wide_block, 4 );
        cv.int_val = 0;
        outfile.write( cv.wide_block, 4 );
    }

    cv.int_val = resolution;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 1.0f/resolution;
    outfile.write( cv.wide_block, 4 );  
    
    int compression_mode = 1;

    char* data_buffer= new char[ resolution*resolution*resolution * sizeof(float) ];
    memset(data_buffer, 0, resolution*resolution*resolution * sizeof(float) );

    // Shadow
    Write_Compressed( outfile,  data_buffer, resolution*resolution*resolution * sizeof(float), compression_mode);
    
    // Density
    PHYSBAM_ASSERT( density.m == resolution*resolution*resolution );
    Write_Compressed(outfile, (char*)density.base_pointer,resolution*resolution*resolution*sizeof(float),compression_mode);
    
    if( use_heat ){
        PHYSBAM_ASSERT( heat.m == resolution*resolution*resolution );

        // Heat
        Write_Compressed(outfile, (char*)heat.base_pointer,resolution*resolution*resolution*sizeof(float),compression_mode);
        // HeatOld
        Write_Compressed(outfile, (char*)heat.base_pointer,resolution*resolution*resolution*sizeof(float),compression_mode);
    }

    // vx
    Write_Compressed( outfile,  data_buffer, resolution*resolution*resolution * sizeof(float), compression_mode);
    // vy
    Write_Compressed( outfile,  data_buffer, resolution*resolution*resolution * sizeof(float), compression_mode);
    // vz
    Write_Compressed( outfile,  data_buffer, resolution*resolution*resolution * sizeof(float), compression_mode);

    // Obstacles
    Write_Compressed( outfile,  data_buffer, resolution*resolution*resolution * sizeof(char), compression_mode);

    // dt
    cv.float_val = dt;
    outfile.write( cv.wide_block, 4 );

    // dx
    cv.float_val = 1.0/resolution;
    outfile.write( cv.wide_block, 4 );

    // p0
    cv.float_val = 0.0f;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    
    // p1
    cv.float_val = 1.0;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    
    // dp0
    cv.float_val = 0.0;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // shift
    cv.int_val = 0;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // obj_shift
    cv.float_val = 0.0f;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // obj_mat
    cv.float_val = 1.0f;
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 0.0f;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    outfile.write( cv.wide_block, 4 );
    cv.float_val = 1.0f;
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 0.0f;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 1.0f;
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 0.0f;
    outfile.write( cv.wide_block, 4 );

    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    cv.float_val = 1.0f;
    outfile.write( cv.wide_block, 4 );

    // base res
    cv.int_val = resolution;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // min res
    cv.int_val = 0;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // max res
    cv.int_val = resolution;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );

    // color
    cv.float_val = 1.0;
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
    outfile.write( cv.wide_block, 4 );
}


int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    //LOG::Initialize_Logging();

    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-input","","path","Input .tet filename");
    parse_args.Add_String_Argument("-output","","path","Output .tri filename");
    parse_args.Add_Integer_Argument("-bitdepth",8,"depth","Bitdepth for voxel conversion [8 16 32 64]");    

    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-input") || !parse_args.Is_Value_Set("-output"))
        parse_args.Print_Usage(true);

    std::string input_filename=parse_args.Get_String_Value("-input");
    std::string output_filename=parse_args.Get_String_Value("-output");
    int bitdepth = parse_args.Get_Integer_Value("-bitdepth");
    if( bitdepth != 8 && bitdepth != 16 && bitdepth != 32 && bitdepth != 64 )
        parse_args.Print_Usage(true);

    
    ARRAY<T, T_INDEX> scalar_field;
    FILE_UTILITIES::Read_From_File(stream_type,input_filename,scalar_field);
    T_INDEX counts = scalar_field.counts;
    PHYSBAM_ASSERT( counts(1) == counts(2) && counts(1) == counts(3) );
    LOG::cout << "Scalar field ( " << counts << " ) read as input..." << std::endl;
    LOG::cout << "Max value: " << scalar_field.Maxabs() << std::endl;

    // Permute counts here?
    T_INDEX P( 1, 2, 3);

    ARRAY<float> linear_scalar_data;

    for( int i=1; i <= counts(P(1)); i++ )
        for( int j=1; j <=  counts(P(2)); j++ )
            for( int k=1; k <=  counts(P(3)); k++ ){
                T_INDEX index;
                index(P(1))=i; index(P(2))=j; index(P(3))=k;
                T value;
                value = scalar_field(index);
                linear_scalar_data.Append( value );
            }

    PHYSBAM_ASSERT( linear_scalar_data.m == counts(1)*counts(2)*counts(3) );

    std::ofstream outfile (output_filename.c_str(),std::ofstream::binary);
    WriteBlenderCache( outfile, counts(1), 1.0f, linear_scalar_data, linear_scalar_data, false );
    outfile.close();
 
    //LOG::Finish_Logging();

    return 0;
}
//#####################################################################
