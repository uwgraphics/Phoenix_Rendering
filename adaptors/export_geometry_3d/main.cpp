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
    parse_args.Add_Integer_Argument("-frame",0,"n","Frame to export");
    parse_args.Add_String_Argument("-prefix","geometry","string","Prefix for output geometry files");
    parse_args.Add_Option_Argument("-compact_particles","Compact all particles prior to exporting geometry");
    parse_args.Parse(argc,argv);
    if(!parse_args.Is_Value_Set("-sim_directory") || !parse_args.Is_Value_Set("-frame"))
        parse_args.Print_Usage(true);
    std::string sim_directory=parse_args.Get_String_Value("-sim_directory");
    std::string prefix=parse_args.Get_String_Value("-prefix");
    int frame=parse_args.Get_Integer_Value("-frame");
    bool compact_particles=parse_args.Is_Value_Set("-compact_particles");

    GEOMETRY_PARTICLES<TV> particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV> deformable_geometry(particles);

		//LOG::cout<<"Simulation directory : "<<sim_directory<<std::endl;
		//LOG::cout<<"Frame to export      : "<<frame<<std::endl;

    {//LOG::SCOPE scope("Reading frame");
        // static bool first_time=true;
        std::string frame_string=STRING_UTILITIES::string_sprintf("%s/%d/",sim_directory.c_str(),frame);
        // std::string static_frame_string=frame_string;
        int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:0;
        // bool read_static_variables=static_frame!=-1 || first_time || !(deformable_geometry && deformable_geometry->structures.m);
        deformable_geometry.Read(STREAM_TYPE(RW()),sim_directory,sim_directory,frame,static_frame,true);
    }

    int number_of_structures=deformable_geometry.structures.m;
		//LOG::cout<<"Number of structures : "<<number_of_structures<<std::endl;

    for(int s=1;s<=number_of_structures;s++){

        STRUCTURE<TV>* structure=deformable_geometry.structures(s);
        std::string indexed_prefix=STRING_UTILITIES::string_sprintf("%s.%d",prefix.c_str(),s);

        if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(structure)){
				//LOG::cout<<"Structure #"<<s<<" is a segmented curve"<<std::endl;
            Write_Structure(stream_type,*segmented_curve,indexed_prefix+".curve",compact_particles);
        }
        else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(structure)){
				//LOG::cout<<"Structure #"<<s<<" is a point cloud"<<std::endl;
        }
        else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
				//LOG::cout<<"Structure #"<<s<<" is a triangulated surface"<<std::endl;
            Write_Structure(stream_type,*triangulated_surface,indexed_prefix+".tri",compact_particles);
        }
        else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
				//LOG::cout<<"Structure #"<<s<<" is a tetrahedralized volume"<<std::endl;
            Write_Structure(stream_type,*tetrahedralized_volume,indexed_prefix+".tet",compact_particles);
        }
        else
            PHYSBAM_FATAL_ERROR("Unexpected structure type");
    }

    // ARRAY<VECTOR<int,3> > elements;
    // ARRAY<VECTOR<T,3> > X;
    // int n,m;

    // std::istream* input=new std::ifstream(input_filename.c_str(),std::ios::in|std::ios::binary);
    // Read_Binary<float>(*input,n,m);
    // elements.Resize(n);
    // X.Resize(m);
    // Read_Binary_Array<float>(*input,&elements(1)(1),3*n);
    // Read_Binary_Array<float>(*input,&X(1)(1),3*m);
    // delete input;

    // LOG::cout<<"Number of triangles = "<<n<<std::endl;
    // LOG::cout<<"Number of particles = "<<m<<std::endl;

    // TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
    // triangulated_surface.mesh.Initialize_Mesh(m,elements);
    // triangulated_surface.particles.array_collection->Resize(m);
    // triangulated_surface.particles.X=X;

    // FILE_UTILITIES::Write_To_File(stream_type,output_filename,triangulated_surface);


    // delete &triangulated_surface;

    //LOG::Finish_Logging();

    return 0;
}
//#####################################################################
