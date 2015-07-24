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
#include <vector>

using namespace PhysBAM;

struct PC2
{
    typedef std::vector< float > SAMPLE;
    typedef std::vector< SAMPLE > FRAMES;
    
    typedef union {
        char BUFFER[32];
        struct{
            char    cacheSignature[12];   // Will be 'POINTCACHE2' followed by a trailing null character.
            int     fileVersion;          // Currently 1
            int     numPoints;            // Number of points per sample
            float   startFrame;           // Corresponds to the UI value of the same name.
            float   sampleRate;           // Corresponds to the UI value of the same name.
            int     numSamples;           // Defines how many samples are stored in the file.
        } FIELDS;
    } PC2_HEADER;
           
    FRAMES samples;

    int numPoints;
    float startFrame;
    float sampleRate;

    

    template<class T_STRUCTURE>
    void AddSample(T_STRUCTURE& structure){
        T_STRUCTURE* new_structure=T_STRUCTURE::Create();
        new_structure->particles.array_collection->Append(*structure.particles.array_collection);
        new_structure->mesh.Initialize_Mesh(structure.mesh);
        new_structure->Discard_Valence_Zero_Particles_And_Renumber();
        
        if(numPoints == -1)
            numPoints = new_structure->particles.X.m;
        else
            PHYSBAM_ASSERT(new_structure->particles.X.m == numPoints);
        
        SAMPLE temp;
        temp.reserve( new_structure->particles.X.m * 3 );
        for( int i=1; i <= new_structure->particles.X.m; i++)
            for( int w=1; w<=3; w++)
                temp.push_back(new_structure->particles.X(i)(w));

        samples.push_back( temp );
    }
    
    void Write_PC2(std::string filename){
        std::ofstream ofs( filename.c_str(), std::ofstream::binary );
        PC2_HEADER header;
        strcpy( header.FIELDS.cacheSignature, "POINTCACHE2" );
        header.FIELDS.fileVersion = 1;
        header.FIELDS.numPoints = numPoints;
        header.FIELDS.startFrame = startFrame;
        header.FIELDS.sampleRate = sampleRate;
        header.FIELDS.numSamples = samples.size();
        ofs.write( header.BUFFER, sizeof(char)*32 );

        LOG::cout << "Writing PC2 File." << std::endl;
        LOG::cout << "Number of Verts: " << numPoints << std::endl;
        LOG::cout << "Number of Samples: " << samples.size() << std::endl;

        char* _sample_buffer = new char[numPoints * sizeof(float) * 3 ];
        float* sample_buffer = reinterpret_cast<float*>( _sample_buffer); 
        for( FRAMES::iterator iter = samples.begin(); iter != samples.end(); iter++){
            for( int s=0; s< numPoints*3; s++)
                sample_buffer[s] = (*iter)[s];
            ofs.write( _sample_buffer, numPoints * sizeof(float) * 3 );            
        }
        delete _sample_buffer;
        
        ofs.close();
    }
    
};


template<class T_STRUCTURE>
void Write_Structure(const STREAM_TYPE stream_type,T_STRUCTURE& structure,const std::string& filename,const bool compacted)
{
    if(!compacted)
        FILE_UTILITIES::Write_To_File(stream_type,filename,structure);
    else{
        T_STRUCTURE* new_structure=T_STRUCTURE::Create();
        new_structure->particles.array_collection->Append(*structure.particles.array_collection);
        new_structure->mesh.Initialize_Mesh(structure.mesh);
        new_structure->Discard_Valence_Zero_Particles_And_Renumber();
        FILE_UTILITIES::Write_To_File(stream_type,filename,*new_structure);
        delete new_structure;
    }
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
    parse_args.Add_String_Argument("-sim_directory","","path","Simulation output directory");
    parse_args.Add_String_Argument("-name","Object","name","Name of the object");
    parse_args.Add_Integer_Argument("-start_frame",0,"n","Start Frame to export");
    parse_args.Add_Integer_Argument("-num_frames",1,"n","Number of frames to export");
    parse_args.Add_Integer_Argument("-object",0,"n","Object Number to generate PointCache.");
    parse_args.Add_Option_Argument("-convert_tets", "Convert Tet object to Surfaces.");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-sim_directory") || !parse_args.Is_Value_Set("-start_frame"))
        parse_args.Print_Usage(true);
    std::string sim_directory=parse_args.Get_String_Value("-sim_directory");
    std::string prefix=parse_args.Get_String_Value("-name");
    int start_frame=parse_args.Get_Integer_Value("-start_frame");
    int num_frames=parse_args.Get_Integer_Value("-num_frames");
    int export_object=parse_args.Get_Integer_Value("-object");
    bool convert_tets=parse_args.Is_Value_Set("-convert_tets");
    bool first_time = true;

    GEOMETRY_PARTICLES<TV> particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV> deformable_geometry(particles);

    //LOG::cout<<"Simulation directory : "<<sim_directory<<std::endl;
    //LOG::cout<<"Frame to export      : "<<frame<<std::endl;
    
    PC2 pc2;

    pc2.numPoints = -1;
    pc2.startFrame = 1.0;
    pc2.sampleRate = 1.0;

    for( int frame = start_frame; frame < num_frames+start_frame; frame++){
        
        {//gLOG::SCOPE scope("Reading frame");
            std::string frame_string=STRING_UTILITIES::string_sprintf("%s/%d/",sim_directory.c_str(),frame);
            int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:0;
            deformable_geometry.Read(STREAM_TYPE(RW()),sim_directory,sim_directory,frame,static_frame,true);
        }
        
        int number_of_structures=deformable_geometry.structures.m;
		if( frame == start_frame ) LOG::cout<<"Number of structures : "<<number_of_structures<<std::endl;
        
        if(number_of_structures < export_object)
            PHYSBAM_FATAL_ERROR("Requested object is not found in archive.");
               
        STRUCTURE<TV>* structure=deformable_geometry.structures(export_object);

        if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(structure)){
            if( frame == start_frame ) LOG::cout<<"Structure #"<<export_object<<" is a segmented curve"<<std::endl;
            pc2.AddSample(*segmented_curve);
        }
        else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(structure)){
            if( frame == start_frame ) LOG::cout<<"Structure #"<<export_object<<" is a point cloud"<<std::endl;
            PHYSBAM_FATAL_ERROR("This is not a supported object type for PointCache export.");
            //pc2.AddSample(*free_particles);
        }
        else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
            if( frame == start_frame ) LOG::cout<<"Structure #"<<export_object<<" is a triangulated surface"<<std::endl;
            pc2.AddSample(*triangulated_surface);
        }
        else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
            if( frame == start_frame ) LOG::cout<<"Structure #"<<export_object<<" is a tetrahedralized volume"<<std::endl;
            if(convert_tets){
                tetrahedralized_volume->Initialize_Triangulated_Surface();
                TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
                triangulated_surface.particles.array_collection->Append(*tetrahedralized_volume->particles.array_collection);
                triangulated_surface.mesh.Initialize_Mesh(tetrahedralized_volume->triangulated_surface->mesh);
                triangulated_surface.Discard_Valence_Zero_Particles_And_Renumber();
                pc2.AddSample(triangulated_surface);
            }           
            else
                pc2.AddSample(*tetrahedralized_volume);
        }
        else
            PHYSBAM_FATAL_ERROR("Unexpected structure type");
    


    }

    pc2.Write_PC2(prefix+".pc2");

    return 0;
}
//#####################################################################
