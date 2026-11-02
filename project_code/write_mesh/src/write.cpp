#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "write.hpp"

void write_mesh()
{
    std::string mesh_type;

    std::cout<<"Input mesh type:"<<std::endl;
    std::cin>>mesh_type;

    if(mesh_type == "bw_step" || mesh_type == "bw_step\r")
    {
        write_bw_step_mesh();
    }
    else if(mesh_type == "cartesian" || mesh_type == "cartesian\r")
    {
        write_cartesian_mesh();
    }
    else
    {
        std::cerr<<"ERROR: UNKNOWN MESH TYPE"<<std::endl;
        std::abort();
    }
}

void write_cartesian_mesh()
{
    int n_x, n_y;   // mesh dimensions
    double x_min, x_max,
            y_min, y_max;

    int n_patches;  // number of patches

    std::string init_p, init_v, init_k, init_eps;
    std::string filepath = "../input/";
    std::string filename;

    std::cout<<"Input mesh dimensions:"<<std::endl;
    std::cin>>n_x>>n_y;

    std::cout<<"Input x boundaries:"<<std::endl;
    std::cin>>x_min>>x_max;
    std::cout<<"Input y boundaries:"<<std::endl;
    std::cin>>y_min>>y_max;

    write_points_cartesian(n_x, n_y, x_min, x_max, y_min, y_max);
    write_faces_cartesian(n_x, n_y);
    write_cells_cartesian(n_x, n_y);

    std::cout<<"Input number of patches:"<<std::endl;
    std::cin>>n_patches;

    write_patches_bcs_cartesian(n_x, n_y, n_patches);

    std::cout<<"Input initial conditions scenario for p:"<<std::endl;
    std::cin>>init_p;

    filename = "pressure.txt";
    write_scalar_field_cartesian(filepath, filename, n_x, n_y, init_p);

    std::cout<<"Input initial conditions scenario for v:"<<std::endl;
    std::cin>>init_v;

    filename = "velocity.txt";
    write_vector_field_cartesian(filepath, filename, n_x, n_y, x_min, x_max, y_min, y_max, init_v);
}



void write_points_cartesian(int n_x, int n_y, double x_min = 0, double x_max = 1, double y_min = 0, double y_max = 1)
{
    // initialisations
    std::string filename = "points.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate number of points
    int n_points_x = n_x+1;
    int n_points_y = n_y+1;
    int n_points = n_points_x * n_points_y;

    // calculate dimensions of each cell
    double dx = (x_max-x_min) / n_x;    // n_x as the number of segments is 1 less than the number of points
    double dy = (y_max-y_min) / n_y;

    // write total number of points to file
    file<<n_points<<std::endl;

    // open points
    file<<"("<<std::endl;

    // write all points
    for(int j=0; j<n_points_y; j++)
    {
        for(int i=0; i<n_points_x; i++)
        {
            file<<"\t("<<i*dx<<" "<<j*dy<<")"<<std::endl;
        }
    }

    // close points
    file<<")"<<std::endl;

    file.close();
}

void write_faces_cartesian(int n_x, int n_y)
{
    // initialisations
    std::string filename = "faces.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate number of points in x dimension for calculating point indices
    int n_points_x = n_x+1;

    // calculate total number of faces
    int n_bdr_faces = 2*(n_x + n_y);
    int n_internal_faces = (n_y-1)*n_x + (n_x-1)*n_y;
    int n_faces = n_bdr_faces + n_internal_faces;

    // write total number of faces to file
    file<<n_faces<<std::endl;

    // open faces
    file<<"("<<std::endl;

    // loop through internal faces
    int loop_size_y = n_y + n_y-1;
    int loop_size_x = n_x-1;

    for(int j=0; j<loop_size_y; j++)
    {
        if(j%2 == 0)    // filling faces parallel to y axis
        {
            for(int i=0; i<loop_size_x; i++)
            {
                int p0 = (j/2)*n_points_x + i + 1;
                int p1 = (j/2+1)*n_points_x + i + 1;
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
        }
        else            // filling faces parallel to x axis
        {
            for(int i=0; i<loop_size_x + 1; i++)
            {
                int p0 = (j/2+1)*n_points_x + i;
                int p1 = (j/2+1)*n_points_x + i + 1;
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
        }
    }

    // loop through boundary faces

    // j=n_y
    for(int i=0; i<n_x; i++)
    {
        int p0 = n_y * n_points_x + i;
        int p1 = n_y * n_points_x + i + 1;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // i=0
    for(int j=0; j<n_y; j++)
    {
        int p0 = j*n_points_x;
        int p1 = (j+1)*n_points_x;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // i=n_x
    for(int j=0; j<n_y; j++)
    {
        int p0 = j*n_points_x + n_points_x - 1;
        int p1 = (j+1)*n_points_x + n_points_x - 1;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // j=0
    for(int i=0; i<n_x; i++)
    {
        int p0 = i;
        int p1 = i + 1;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // close faces
    file<<")"<<std::endl;

    file.close();
}

void write_cells_cartesian(int n_x, int n_y)
{
    // initialisations
    std::string filename = "cells.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate the number of internal faces
    int n_internal_faces = (n_y-1)*n_x + (n_x-1)*n_y;

    // calculate total number of cells
    int n_cells = n_x * n_y;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write all cells
    for(int j=0; j<n_y; j++)
    {
        for(int i=0; i<n_x; i++)
        {
            std::array<int, 4> faces; // faces[0] - x_max, faces[1] - y_max, faces[2] - x_min, faces[3] - y_min

            if(i == 0)  // mesh x_min
            {
                faces[2] = n_internal_faces + n_x + j;
            }
            else
            {
                faces[2] = j*(n_x-1 + n_x) + i-1;
            }

            if(i == n_x-1) // mesh x_max
            {
                faces[0] = n_internal_faces + n_x + n_y + j;
            }
            else
            {
                faces[0] = j*(n_x-1 + n_x) + i;
            }

            if(j == 0)  // mesh y_min
            {
                faces[3] = n_internal_faces + n_x + 2*n_y + i;
            }
            else
            {
                faces[3] = (j-1)*(n_x-1 + n_x) + n_x-1 + i;
            }

            if(j == n_y-1)  // mesh y_max
            {
                faces[1] = n_internal_faces + i;
            }
            else
            {
                faces[1] = j*(n_x-1 + n_x) + n_x-1 + i;
            }

            // if mesh x_max, flip order to 1-0-3-2
            if(i == n_x-1)
            {
                int temp;   // helper variable

                // swap faces[0] and faces[1]
                temp = faces[0];
                faces[0] = faces[1];
                faces[1] = temp;

                // swap faces[2] and faces[3]
                temp = faces[2];
                faces[2] = faces[3];
                faces[3] = temp;
            }

            // write faces to file, in order
            file<<"\t4("<<faces[0]<<" "
                        <<faces[1]<<" "
                        <<faces[2]<<" "
                        <<faces[3]<<")"<<std::endl;
        }
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}

void write_patches_bcs_cartesian(int n_x, int n_y, int n_patches)
{
    // initialisations
    std::string patch_filename = "patches.txt";
    std::string bc_filename = "boundary_conditions.txt";

    std::string filepath = "../input/";
    std::ofstream patch_file(filepath + patch_filename);
    std::ofstream bc_file(filepath + bc_filename);

    int n_bdr_faces = 2*(n_x + n_y);
    int n_internal_faces = (n_y-1)*n_x + (n_x-1)*n_y;
    int n_faces = n_bdr_faces + n_internal_faces;

    int bdr_faces_used = 0;                 // variables used to make sure all boundary
    int min_face = n_faces, max_face = 0;   // faces were used in the patches file

    // write total number of patches to patch_file
    patch_file<<n_patches<<std::endl;

    // open patches
    patch_file<<"("<<std::endl;

    // write patch_file and bc_file in parallel
    for(int n=0; n<n_patches; n++)
    {
        // declarations
        std::string patch_name, patch_bc;
        int face0, face1;

        // ---------------------- INPUT ----------------------
        std::cout<<"Write patch name:"<<std::endl;
        std::cin>>patch_name;

        std::cout<<"Write the first and last face indices belonging to the patch:"<<std::endl;
        std::cin>>face0>>face1;

        std::cout<<"Write patch boundary condition:"<<std::endl;
        patch_bc = read_patch_bc();

        // swap if in wrong order
        if(face0 > face1)
        {
            int temp = face0;
            face0 = face1;
            face1 = temp;
        }

        // ---------------------- patch_file ----------------------

        // write patch name to patch_file
        patch_file<<"\t"<<patch_name<<std::endl;

        // write number of faces to patch_file
        patch_file<<"\t"<<face1-face0+1<<std::endl;

        // open patch
        patch_file<<"\t("<<std::endl;

        // write face indices to patch_file
        patch_file<<"\t\t"<<face0;
        for(int i=face0+1; i<=face1; i++)
            patch_file<<" "<<i;
        patch_file<<std::endl;

        // close patch
        patch_file<<"\t)"<<std::endl;

        // update bdr_faces_used, min_face and max_face
        bdr_faces_used += face1-face0+1;

        if(face0 < min_face)
            min_face = face0;

        if(face1 > max_face)
            max_face = face1;

        // ---------------------- bc_file ----------------------
        // write boundary condition to bc_file
        bc_file<<patch_name<<" "<<patch_bc<<std::endl;
    }

    // close patches
    patch_file<<")"<<std::endl;

    patch_file.close();

    // perform patch checks
    if(min_face > n_internal_faces)
        std::cerr<<"FACES FROM "<<n_internal_faces<<" TO "<<min_face - 1<<" WERE NOT USED!"<<std::endl;
    else if(min_face < n_internal_faces)
        std::cerr<<"MINIMUM FACE USED "<<min_face<<" IS NOT A BOUNDARY FACE!"<<std::endl;

    if(max_face < n_faces - 1)
        std::cerr<<"FACES FROM "<<max_face + 1<<" TO "<<n_faces - 1<<" WERE NOT USED!"<<std::endl;
    else if(max_face > n_faces - 1)
        std::cerr<<"MAXIMUM FACE USED "<<max_face<<" IS NOT A BOUNDARY FACE!"<<std::endl;

    if(bdr_faces_used != n_bdr_faces)
        std::cerr<<"NUMBER OF BOUNDARY FACES USED "<<bdr_faces_used
                <<" DOES NOT MATCH THE TOTAL NUMBER OF BOUNDARY FACES "<<n_bdr_faces<<"!"<<std::endl;
}



void write_scalar_field_cartesian(std::string& filepath, std::string& filename,
                        int n_x, int n_y,
                        std::string& init_field)
{
    // initialisations
    std::ofstream file(filepath + filename);

    // calculate total number of cells
    int n_cells = n_x * n_y;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write scalar to all cells
    if(init_field == "zero")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(0)"<<std::endl;
    }
    else if(init_field == "init1")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(1)"<<std::endl;
    }
    
    else
    {
        std::cerr<<"UNRECOGNISED INITIAL CONDITION CASE FOR "<<filename<<" : "<<init_field<<"!"<<std::endl;
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}

void write_vector_field_cartesian(std::string& filepath, std::string& filename,
                        int n_x, int n_y, double x_min, double x_max, double y_min, double y_max,
                        std::string& init_field)
{
    // initialisations
    std::ofstream file(filepath + filename);

    // calculate total number of cells
    int n_cells = n_x * n_y;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write vector to all cells
    if(init_field == "zero")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(0 0)"<<std::endl;
    }
    else if(init_field == "init1")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(1 0)"<<std::endl;
    }
    else if(init_field == "stepfunc")
    {
        for(int i=0; i<n_cells; i++)
        {
            if(i%n_x==1 && i>n_x && i<n_cells-n_x)
                file<<"\t(1 0)"<<std::endl;
            else
                file<<"\t(0 0)"<<std::endl;
        }
    }
    else if(init_field == "diffusion_x")    // set up diffusion test case to solve only velocity_systems[X]
    {
        // calculate constants
        double PI = 4 * atan(1);
        double dx = (x_max-x_min) / n_x;
        double dy = (y_max-y_min) / n_y;

        for(int i=0; i<n_x; i++)
        {
            double centroid_x = x_min + dx/2 + i*dx;

            for(int j=0; j<n_x; j++)
            {
                double centroid_y = y_min + dy/2 + j*dy;

                double diffusion_init_val = -2*PI*PI*sin(10*PI*centroid_x)*sin(10*PI*centroid_y);
                file<<"\t("<<diffusion_init_val<<" 0)"<<std::endl;
            }
        }
    }
    else if(init_field == "shearflow")
    {
        for(int i=0; i<n_cells; i++)
        {
            if(i<n_x*(n_y/2))
                file<<"\t(1 0)"<<std::endl;
            else
                file<<"\t(1.1 0)"<<std::endl;
        }
    }
    else if(init_field == "test_grad_vel_t")
    {
        double dx = (x_max-x_min)/n_x;
        double dy = (y_max-y_min)/n_y;

        for(int j=0; j<n_y; j++)
        {
            for(int i=0; i<n_x; i++)
            {
                file<<"\t("<<1+i*dx<<" 0)"<<std::endl;
            }
        }
    }
    
    else
    {
        std::cerr<<"UNRECOGNISED INITIAL CONDITION CASE FOR "<<filename<<" : "<<init_field<<"!"<<std::endl;
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}



void write_bw_step_mesh()
{
    int n_x_inlet, n_y_inlet,
        n_x_outlet, n_y_outlet;   // mesh dimensions
    double x_min_inlet, x_max_inlet,
            x_min_outlet, x_max_outlet,
            y_min_inlet, y_max_inlet,
            y_min_outlet, y_max_outlet;

    int n_patches;  // number of patches

    std::string init_p, init_v, init_k, init_eps;
    std::string filepath = "../input/";
    std::string filename;

    // inlet data
    std::cout<<"Input inlet dimensions:"<<std::endl;
    std::cin>>n_x_inlet>>n_y_inlet;

    std::cout<<"Input inlet x boundaries:"<<std::endl;
    std::cin>>x_min_inlet>>x_max_inlet;
    std::cout<<"Input inlet y boundaries:"<<std::endl;
    std::cin>>y_min_inlet>>y_max_inlet;

    // outlet data
    std::cout<<"Input outlet dimensions:"<<std::endl;
    std::cin>>n_x_outlet>>n_y_outlet;

    std::cout<<"Input outlet x boundaries:"<<std::endl;
    std::cin>>x_min_outlet>>x_max_outlet;
    std::cout<<"Input outlet y boundaries:"<<std::endl;
    std::cin>>y_min_outlet>>y_max_outlet;

    // prerequisite: x_max_inlet == x_min_outlet
    assert(x_max_inlet == x_min_outlet);

    write_points_bw_step(n_x_inlet, n_y_inlet, x_min_inlet, x_max_inlet, y_min_inlet, y_max_inlet,
                        n_x_outlet, n_y_outlet, x_min_outlet, x_max_outlet, y_min_outlet, y_max_outlet);
    write_faces_bw_step(n_x_inlet, n_y_inlet, n_x_outlet, n_y_outlet);
    write_cells_bw_step(n_x_inlet, n_y_inlet, n_x_outlet, n_y_outlet);

    std::cout<<"Input number of patches:"<<std::endl;
    std::cin>>n_patches;

    write_patches_bcs_bw_step(n_x_inlet, n_y_inlet, n_x_outlet, n_y_outlet, n_patches);

    std::cout<<"Input initial conditions scenario for p:"<<std::endl;
    std::cin>>init_p;

    filename = "pressure.txt";
    write_scalar_field_bw_step(filepath, filename, n_x_inlet, n_y_inlet, n_x_outlet, n_y_outlet, init_p);

    std::cout<<"Input initial conditions scenario for v:"<<std::endl;
    std::cin>>init_v;

    filename = "velocity.txt";
    write_vector_field_bw_step(filepath, filename, n_x_inlet, n_y_inlet, n_x_outlet, n_y_outlet, init_v);
}



void write_points_bw_step(int n_x_inlet, int n_y_inlet,
                            double x_min_inlet, double x_max_inlet, double y_min_inlet, double y_max_inlet,
                            int n_x_outlet, int n_y_outlet,
                            double x_min_outlet, double x_max_outlet, double y_min_outlet, double y_max_outlet)
{
    // initialisations
    std::string filename = "points.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate number of points,
    // assuming proper initial data
    int n_points = n_x_inlet * (n_y_inlet+1) + (n_x_outlet+1) * (n_y_outlet+1);

    // calculate dimensions of each cell
    double dx_inlet = (x_max_inlet-x_min_inlet) / n_x_inlet;    // n_x as the number of segments is 1 less than the number of points
    double dy_inlet = (y_max_inlet-y_min_inlet) / n_y_inlet;
    double dx_outlet = (x_max_outlet-x_min_outlet) / n_x_outlet;
    double dy_outlet = (y_max_outlet-y_min_outlet) / n_y_outlet;

    // assert(fabs(dx_inlet - dx_outlet) <= EPSILON);
    assert(fabs(dy_inlet - dy_outlet) <= EPSILON);

    // write total number of points to file
    file<<n_points<<std::endl;

    // open points
    file<<"("<<std::endl;

    // write inlet points
    for(int j=0; j<n_y_inlet+1; j++)
    {
        for(int i=0; i<n_x_inlet; i++)
        {
            file<<"\t("<<x_min_inlet + i*dx_inlet<<" "<<y_min_inlet + j*dy_inlet<<")"<<std::endl;
        }
    }

    // write outlet points
    for(int j=0; j<n_y_outlet+1; j++)
    {
        for(int i=0; i<n_x_outlet+1; i++)
        {
            file<<"\t("<<x_min_outlet + i*dx_outlet<<" "<<y_min_outlet + j*dy_outlet<<")"<<std::endl;
        }
    }

    // close points
    file<<")"<<std::endl;

    file.close();
}

void write_faces_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet)
{
    // initialisations
    std::string filename = "faces.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate total number of faces
    int n_internal_inlet = (n_y_inlet-1)*n_x_inlet + (n_x_inlet-1)*n_y_inlet;
    int n_internal_outlet = (n_y_outlet-1)*n_x_outlet + (n_x_outlet-1)*n_y_outlet;
    int n_inner_bdr = n_y_inlet;
    int n_bdr_faces = 2*(n_x_inlet + n_x_outlet + n_y_outlet);
    int n_faces = n_bdr_faces + n_internal_inlet + n_internal_outlet + n_inner_bdr;

    // calculate other parameters
    int n_points_inlet = (n_x_inlet)*(n_y_inlet+1);

    // write total number of faces to file
    file<<n_faces<<std::endl;

    // open faces
    file<<"("<<std::endl;

    // loop through internal faces
    // inlet
    {
        int loop_size_y = 2*n_y_inlet-1;
        int loop_size_x = n_x_inlet-1;

        for(int j=0; j<loop_size_y; j++)
        {
            if(j%2 == 0)    // filling faces parallel to y axis
            {
                for(int i=0; i<loop_size_x; i++)
                {
                    int p0 = (j/2)*n_x_inlet + i + 1;
                    int p1 = (j/2+1)*n_x_inlet + i + 1;
                    file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
                }
            }
            else            // filling faces parallel to x axis
            {
                for(int i=0; i<loop_size_x; i++)
                {
                    int p0 = (j/2+1)*n_x_inlet + i;
                    int p1 = (j/2+1)*n_x_inlet + i + 1;
                    file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
                }
                int p0 = (j/2+1)*n_x_inlet + loop_size_x;
                int p1 = n_points_inlet+(n_y_outlet-n_y_inlet+j/2+1)*(n_x_outlet+1);
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
        }
    }

    // inlet-outlet boundary
    for(int i=0; i<n_inner_bdr; i++)
    {
        int p0 = n_points_inlet+(n_y_outlet-n_y_inlet+i)*(n_x_outlet+1);
        int p1 = n_points_inlet+(n_y_outlet-n_y_inlet+i+1)*(n_x_outlet+1);
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // outlet
    {
        int loop_size_y = 2*n_y_outlet-1;
        int loop_size_x = n_x_outlet-1;

        for(int j=0; j<loop_size_y; j++)
        {
            if(j%2 == 0)    // filling faces parallel to y axis
            {
                for(int i=0; i<loop_size_x; i++)
                {
                    int p0 = n_points_inlet+(j/2)*(n_x_outlet+1) + i + 1;
                    int p1 = n_points_inlet+(j/2+1)*(n_x_outlet+1) + i + 1;
                    file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
                }
            }
            else            // filling faces parallel to x axis
            {
                for(int i=0; i<loop_size_x + 1; i++)
                {
                    int p0 = n_points_inlet+(j/2+1)*(n_x_outlet+1) + i;
                    int p1 = n_points_inlet+(j/2+1)*(n_x_outlet+1) + i + 1;
                    file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
                }
            }
        }
    }

    // loop through boundary faces

    // upper wall
    {
        // inlet
        {
            for(int i=0; i<n_x_inlet-1; i++)
            {
                int p0 = n_y_inlet*n_x_inlet + i;
                int p1 = n_y_inlet*n_x_inlet + i + 1;
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
            int p0 = n_points_inlet-1;
            int p1 = n_points_inlet + n_y_outlet*(n_x_outlet+1);
            file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
        }

        // outlet
        {
            for(int i=0; i<n_x_outlet; i++)
            {
                int p0 = n_points_inlet + n_y_outlet*(n_x_outlet+1) + i;
                int p1 = n_points_inlet + n_y_outlet*(n_x_outlet+1) + i + 1;
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
        }
    }


    // inlet
    {
        for(int j=0; j<n_y_inlet; j++)
        {
            int p0 = j*n_x_inlet;
            int p1 = (j+1)*n_x_inlet;
            file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
        }
    }

    // bw step
    {
        // parallel to x axis
        {
            for(int i=0; i<n_x_inlet-1; i++)
            {
                int p0 = i;
                int p1 = i + 1;
                file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
            }
            int p0 = n_x_inlet-1;
            int p1 = n_points_inlet + (n_y_outlet-n_y_inlet)*(n_x_outlet+1);
            file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
        }

        // parallel to y axis
        for(int j=0; j < n_y_outlet - n_y_inlet; j++)
        {
            int p0 = n_points_inlet + j*(n_x_outlet+1);
            int p1 = n_points_inlet + (j+1)*(n_x_outlet+1);
            file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
        }
    }

    // outlet
    for(int j=0; j<n_y_outlet; j++)
    {
        int p0 = n_points_inlet + j*(n_x_outlet+1) + n_x_outlet;
        int p1 = n_points_inlet + (j+1)*(n_x_outlet+1) + n_x_outlet;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // lower wall
    for(int i=0; i<n_x_outlet; i++)
    {
        int p0 = n_points_inlet + i;
        int p1 = n_points_inlet + i + 1;
        file<<"\t2("<<p0<<" "<<p1<<")"<<std::endl;
    }

    // close faces
    file<<")"<<std::endl;

    file.close();
}

void write_cells_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet)
{
    // initialisations
    std::string filename = "cells.txt";
    std::string filepath = "../input/";
    std::ofstream file(filepath + filename);

    // calculate parameters
    int n_internal_inlet = (n_y_inlet-1)*n_x_inlet + (n_x_inlet-1)*n_y_inlet;
    int n_internal_outlet = (n_y_outlet-1)*n_x_outlet + (n_x_outlet-1)*n_y_outlet;
    int n_inner_bdr = n_y_inlet;
    int n_internal_faces = n_internal_inlet + n_internal_outlet + n_inner_bdr;

    int n_cells = n_x_inlet * n_y_inlet + n_x_outlet * n_y_outlet;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write all cells
    // inlet
    for(int j=0; j<n_y_inlet; j++)
    {
        for(int i=0; i<n_x_inlet; i++)
        {
            std::array<int, 4> faces; // faces[0] - x_max, faces[1] - y_max, faces[2] - x_min, faces[3] - y_min

            if(i == 0)  // mesh x_min
            {
                faces[2] = n_internal_faces + (n_x_inlet + n_x_outlet) + j;
            }
            else
            {
                faces[2] = j*(n_x_inlet-1 + n_x_inlet) + i-1;
            }

            if(i == n_x_inlet-1) // mesh x_max
            {
                faces[0] = n_internal_inlet + j;
            }
            else
            {
                faces[0] = j*(n_x_inlet-1 + n_x_inlet) + i;
            }

            if(j == 0)  // mesh y_min
            {
                faces[3] = n_internal_faces + (n_x_inlet + n_x_outlet) + n_y_inlet + i;
            }
            else
            {
                faces[3] = (j-1)*(n_x_inlet-1 + n_x_inlet) + n_x_inlet-1 + i;
            }

            if(j == n_y_inlet-1)  // mesh y_max
            {
                faces[1] = n_internal_faces + i;
            }
            else
            {
                faces[1] = j*(n_x_inlet-1 + n_x_inlet) + n_x_inlet-1 + i;
            }

            // // if mesh x_max, flip order to 1-0-3-2
            // if(i == n_x_inlet-1)
            // {
            //     int temp;   // helper variable

            //     // swap faces[0] and faces[1]
            //     temp = faces[0];
            //     faces[0] = faces[1];
            //     faces[1] = temp;

            //     // swap faces[2] and faces[3]
            //     temp = faces[2];
            //     faces[2] = faces[3];
            //     faces[3] = temp;
            // }

            // write faces to file, in order
            file<<"\t4("<<faces[0]<<" "
                        <<faces[1]<<" "
                        <<faces[2]<<" "
                        <<faces[3]<<")"<<std::endl;
        }
    }

    // outlet
    for(int j=0; j<n_y_outlet; j++)
    {
        for(int i=0; i<n_x_outlet; i++)
        {
            std::array<int, 4> faces; // faces[0] - x_max, faces[1] - y_max, faces[2] - x_min, faces[3] - y_min

            if(i == 0)  // mesh x_min
            {
                if(j < n_y_outlet - n_y_inlet)
                {
                    faces[2] = n_internal_faces + (n_x_inlet+n_x_outlet) + n_y_inlet + n_x_inlet + j;
                }
                else
                {
                    faces[2] = n_internal_inlet + j - (n_y_outlet - n_y_inlet);
                }
            }
            else
            {
                faces[2] = n_internal_inlet + n_inner_bdr + j*(n_x_outlet-1 + n_x_outlet) + i-1;
            }

            if(i == n_x_outlet-1) // mesh x_max
            {
                faces[0] = n_internal_faces + (n_x_inlet + n_x_outlet) + n_y_outlet + n_x_inlet + j;
            }
            else
            {
                faces[0] = n_internal_inlet + n_inner_bdr + j*(n_x_outlet-1 + n_x_outlet) + i;
            }

            if(j == 0)  // mesh y_min
            {
                faces[3] = n_internal_faces + 2*n_x_inlet + n_x_outlet + 2*n_y_outlet + i;
            }
            else
            {
                faces[3] = n_internal_inlet + n_inner_bdr + (j-1)*(n_x_outlet-1 + n_x_outlet) + n_x_outlet-1 + i;
            }

            if(j == n_y_outlet-1)  // mesh y_max
            {
                faces[1] = n_internal_faces + n_x_inlet + i;
            }
            else
            {
                faces[1] = n_internal_inlet + n_inner_bdr + j*(n_x_outlet-1 + n_x_outlet) + n_x_outlet-1 + i;
            }

            // if mesh x_max, flip order to 1-0-3-2
            if(i == n_x_outlet-1)
            {
                int temp;   // helper variable

                // swap faces[0] and faces[1]
                temp = faces[0];
                faces[0] = faces[1];
                faces[1] = temp;

                // swap faces[2] and faces[3]
                temp = faces[2];
                faces[2] = faces[3];
                faces[3] = temp;
            }

            // write faces to file, in order
            file<<"\t4("<<faces[0]<<" "
                        <<faces[1]<<" "
                        <<faces[2]<<" "
                        <<faces[3]<<")"<<std::endl;
        }
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}

void write_patches_bcs_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet, int n_patches)
{
    // initialisations
    std::string patch_filename = "patches.txt";
    std::string bc_filename = "boundary_conditions.txt";

    std::string filepath = "../input/";
    std::ofstream patch_file(filepath + patch_filename);
    std::ofstream bc_file(filepath + bc_filename);

    int n_internal_inlet = (n_y_inlet-1)*n_x_inlet + (n_x_inlet-1)*n_y_inlet;
    int n_internal_outlet = (n_y_outlet-1)*n_x_outlet + (n_x_outlet-1)*n_y_outlet;
    int n_inner_bdr = n_y_inlet;
    int n_internal_faces = n_internal_inlet + n_inner_bdr+ n_internal_outlet;
    int n_bdr_faces = 2*(n_x_inlet + n_x_outlet + n_y_outlet);
    int n_faces = n_internal_faces + n_bdr_faces;

    int bdr_faces_used = 0;                 // variables used to make sure all boundary
    int min_face = n_faces, max_face = 0;   // faces were used in the patches file

    // write total number of patches to patch_file
    patch_file<<n_patches<<std::endl;

    // open patches
    patch_file<<"("<<std::endl;

    // write patch_file and bc_file in parallel
    for(int n=0; n<n_patches; n++)
    {
        // declarations
        std::string patch_name, patch_bc;
        int face0, face1;

        // ---------------------- INPUT ----------------------
        std::cout<<"Write patch name:"<<std::endl;
        std::cin>>patch_name;

        std::cout<<"Write the first and last face indices belonging to the patch:"<<std::endl;
        std::cin>>face0>>face1;

        std::cout<<"Write patch boundary condition:"<<std::endl;
        patch_bc = read_patch_bc();

        // swap if in wrong order
        if(face0 > face1)
        {
            int temp = face0;
            face0 = face1;
            face1 = temp;
        }

        // ---------------------- patch_file ----------------------

        // write patch name to patch_file
        patch_file<<"\t"<<patch_name<<std::endl;

        // write number of faces to patch_file
        patch_file<<"\t"<<face1-face0+1<<std::endl;

        // open patch
        patch_file<<"\t("<<std::endl;

        // write face indices to patch_file
        patch_file<<"\t\t"<<face0;
        for(int i=face0+1; i<=face1; i++)
            patch_file<<" "<<i;
        patch_file<<std::endl;

        // close patch
        patch_file<<"\t)"<<std::endl;

        // update bdr_faces_used, min_face and max_face
        bdr_faces_used += face1-face0+1;

        if(face0 < min_face)
            min_face = face0;

        if(face1 > max_face)
            max_face = face1;

        // ---------------------- bc_file ----------------------
        // write boundary condition to bc_file
        bc_file<<patch_name<<" "<<patch_bc<<std::endl;
    }

    // close patches
    patch_file<<")"<<std::endl;

    patch_file.close();

    // perform patch checks
    if(min_face > n_internal_faces)
        std::cerr<<"FACES FROM "<<n_internal_faces<<" TO "<<min_face - 1<<" WERE NOT USED!"<<std::endl;
    else if(min_face < n_internal_faces)
        std::cerr<<"MINIMUM FACE USED "<<min_face<<" IS NOT A BOUNDARY FACE!"<<std::endl;

    if(max_face < n_faces - 1)
        std::cerr<<"FACES FROM "<<max_face + 1<<" TO "<<n_faces - 1<<" WERE NOT USED!"<<std::endl;
    else if(max_face > n_faces - 1)
        std::cerr<<"MAXIMUM FACE USED "<<max_face<<" IS NOT A BOUNDARY FACE!"<<std::endl;

    if(bdr_faces_used != n_bdr_faces)
        std::cerr<<"NUMBER OF BOUNDARY FACES USED "<<bdr_faces_used
                <<" DOES NOT MATCH THE TOTAL NUMBER OF BOUNDARY FACES "<<n_bdr_faces<<"!"<<std::endl;
}



void write_scalar_field_bw_step(std::string& filepath, std::string& filename,
                        int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet,
                        std::string& init_field)
{
    // initialisations
    std::ofstream file(filepath + filename);

    // calculate total number of cells
    int n_cells = n_x_inlet * n_y_inlet + n_x_outlet * n_y_outlet;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write scalar to all cells
    if(init_field == "zero")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(0)"<<std::endl;
    }
    else if(init_field == "init1")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(1)"<<std::endl;
    }
    
    else
    {
        std::cerr<<"UNRECOGNISED INITIAL CONDITION CASE FOR "<<filename<<" : "<<init_field<<"!"<<std::endl;
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}

void write_vector_field_bw_step(std::string& filepath, std::string& filename,
                        int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet,
                        std::string& init_field)
{
    // initialisations
    std::ofstream file(filepath + filename);

    // calculate total number of cells
    int n_cells = n_x_inlet * n_y_inlet + n_x_outlet * n_y_outlet;

    // write total number of cells to file
    file<<n_cells<<std::endl;

    // open cells
    file<<"("<<std::endl;

    // write vector to all cells
    if(init_field == "zero")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(0 0)"<<std::endl;
    }
    else if(init_field == "init1")
    {
        for(int i=0; i<n_cells; i++)
            file<<"\t(1 0)"<<std::endl;
    }
    else if(init_field == "bw_step1")
    {
        int i=0;
        int start_lower_outlet_id = n_x_inlet * n_y_inlet;
        int start_higher_outlet_id = n_cells - n_x_outlet * (n_y_outlet - n_y_inlet);

        for(; i < start_lower_outlet_id; i++)
            file<<"\t(1 0)"<<std::endl;

        for(; i < start_higher_outlet_id; i++)
            file<<"\t(0 0)"<<std::endl;

        for(; i < n_cells; i++)
            file<<"\t(1 0)"<<std::endl;
    }
    // else if(init_field == "stepfunc")
    // {
    //     for(int i=0; i<n_cells; i++)
    //     {
    //         if(i%n_x==1 && i>n_x && i<n_cells-n_x)
    //             file<<"\t(1 0)"<<std::endl;
    //         else
    //             file<<"\t(0 0)"<<std::endl;
    //     }
    // }
    // else if(init_field == "diffusion_x")    // set up diffusion test case to solve only velocity_systems[X]
    // {
    //     // calculate constants
    //     double PI = 4 * atan(1);
    //     double dx = (x_max-x_min) / n_x;
    //     double dy = (y_max-y_min) / n_y;

    //     for(int i=0; i<n_x; i++)
    //     {
    //         double centroid_x = x_min + dx/2 + i*dx;

    //         for(int j=0; j<n_x; j++)
    //         {
    //             double centroid_y = y_min + dy/2 + j*dy;

    //             double diffusion_init_val = -2*PI*PI*sin(10*PI*centroid_x)*sin(10*PI*centroid_y);
    //             file<<"\t("<<diffusion_init_val<<" 0)"<<std::endl;
    //         }
    //     }
    // }
    // else if(init_field == "shearflow")
    // {
    //     for(int i=0; i<n_cells; i++)
    //     {
    //         if(i<n_x*(n_y/2))
    //             file<<"\t(1 0)"<<std::endl;
    //         else
    //             file<<"\t(1.1 0)"<<std::endl;
    //     }
    // }
    
    else
    {
        std::cerr<<"UNRECOGNISED INITIAL CONDITION CASE FOR "<<filename<<" : "<<init_field<<"!"<<std::endl;
    }

    // close cells
    file<<")"<<std::endl;

    file.close();
}


std::string read_patch_bc()
{
    std::string patch_bc;

    while(std::getline(std::cin, patch_bc))
    {
        if(patch_bc == "" || patch_bc == "\r")
            continue;
        else
        {
            if(patch_bc.back() == '\r')
                patch_bc.pop_back();

            break;
        }
    }

    return patch_bc;
}
