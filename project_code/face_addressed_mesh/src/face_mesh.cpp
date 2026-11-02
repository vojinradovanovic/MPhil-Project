#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "face_mesh.hpp"

// ----------------------------- FaceMesh methods -----------------------------

// public

// empty constructor
FaceMesh::FaceMesh()
{}

void FaceMesh::initialise()
{
    // read initial vectors from input, perform mesh checks
    read_from_input();
    check_mesh();

    calc_parameters();
    check_parameters();
}



int FaceMesh::get_num_points()
{
    return points.size();
}

int FaceMesh::get_num_faces()
{
    return faces.size();
}

int FaceMesh::get_num_cells()
{
    return cells.size();
}

int FaceMesh::get_num_patches()
{
    return patches.size();
}



std::vector<double> FaceMesh::get_cell_volumes()
{
    std::vector<double> cell_volumes;
    int num_cells = cells.size();
    cell_volumes.resize(num_cells);

    for(int i=0; i<num_cells; i++)
        cell_volumes[i] = cells[i].get_volume();

    return cell_volumes;
}

std::vector<Point> FaceMesh::get_face_area_vectors()
{
    std::vector<Point> face_area_vectors;
    int num_faces = faces.size();
    face_area_vectors.resize(num_faces);

    for(int i=0; i<num_faces; i++)
        face_area_vectors[i] = faces[i].get_area_vector();

    return face_area_vectors;
}

std::vector<Point> FaceMesh::get_cell_centroids()
{
    std::vector<Point> cell_centroids;
    int num_cells = cells.size();
    cell_centroids.resize(num_cells);

    for(int i=0; i<num_cells; i++)
        cell_centroids[i] = cells[i].get_centroid();

    return cell_centroids;
}

std::vector<Point> FaceMesh::get_face_centroids()
{
    std::vector<Point> face_centroids;
    int num_faces = faces.size();
    face_centroids.resize(num_faces);

    for(int i=0; i<num_faces; i++)
        face_centroids[i] = faces[i].get_centroid();

    return face_centroids;
}



std::vector<std::string> FaceMesh::get_patch_names()
{
    std::vector<std::string> patch_names(num_patches);

    for(int i=0; i<num_patches; i++)
        patch_names[i] = patches[i].get_name();

    return patch_names;
}

std::vector<int> FaceMesh::get_patch_face_ids(int patch_id)
{
    return patches[patch_id].get_face_ids();
}

std::vector<int> FaceMesh::get_patch_cell_ids(int patch_id)
{
    // initialisations
    std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
    int num_patch_faces = patch_face_ids.size();
    std::vector<int> patch_cell_ids(num_patch_faces);

    for(int i=0; i<num_patch_faces; i++)
        patch_cell_ids[i] = faces[patch_face_ids[i]].get_owner_id();

    return patch_cell_ids;
}

std::vector<int> FaceMesh::get_neighbouring_cell_ids(int cell_id)
{
    // initialisations
    std::vector<int> cell_face_ids = cells[cell_id].get_face_ids();
    int num_cell_faces = cell_face_ids.size();
    std::vector<int> neighbouring_cell_ids;
    neighbouring_cell_ids.reserve(num_cell_faces);

    for(int i=0; i<num_cell_faces; i++)
    {
        int face_owner_id = faces[cell_face_ids[i]].get_owner_id();
        int face_neighbour_id = faces[cell_face_ids[i]].get_neighbour_id();

        if(face_owner_id != cell_id)
        {
            if(face_owner_id >= 0)
                neighbouring_cell_ids.push_back(face_owner_id);
        }

        else if(face_neighbour_id != cell_id)
        {
            if(face_neighbour_id >= 0)
                neighbouring_cell_ids.push_back(face_neighbour_id);
        }
    }

    return neighbouring_cell_ids;
}



std::vector<int> FaceMesh::get_owner_ids()
{
    std::vector<int> owner_ids(num_faces);

    for(int i=0; i<num_faces; i++)
        owner_ids[i] = faces[i].get_owner_id();

    return owner_ids;
}

std::vector<int> FaceMesh::get_neighbour_ids()
{
    std::vector<int> neighbour_ids(num_faces);

    for(int i=0; i<num_faces; i++)
        neighbour_ids[i] = faces[i].get_neighbour_id();

    return neighbour_ids;
}



// protected

void FaceMesh::check_mesh()
{
    // initialisations
    bool abort_flag = false;

    // check for invalid point ids in faces: negative or over num_points
    for(int i=0; i<num_faces; i++)
    {
        for(int point_id: faces[i].get_point_ids())
            if(point_id < 0 || point_id >= num_points)
            {
                std::cerr<<"INVALID POINT ID:\t"<<point_id<<"\t(FACE "<<i<<")!"<<std::endl;
                abort_flag = true;
            }
    }

    // check for invalid face ids in cells: negative or over num_faces
    for(int i=0; i<num_cells; i++)
    {
        for(int face_id: cells[i].get_face_ids())
            if(face_id < 0 || face_id >= num_faces)
            {
                std::cerr<<"INVALID FACE ID:\t"<<face_id<<"\t(CELL "<<i<<")!"<<std::endl;
                abort_flag = true;
            }
    }

    // check for invalid face ids in patches: negative or over num_faces
    for(int i=0; i<num_patches; i++)
    {
        for(int face_id: patches[i].get_face_ids())
            if(face_id < 0 || face_id >= num_faces)
            {
                std::cerr<<"INVALID FACE ID:\t"<<face_id<<"\t(PATCH "<<i<<")!"<<std::endl;
                abort_flag = true;
            }
    }

    // check no two patches overlap
    std::vector<int> face_counts(num_faces, 0);

    for(int i=0; i<num_patches; i++)
    {
        for(int face_id: patches[i].get_face_ids())
            face_counts[face_id]++;
    }

    for(int i=0; i<num_faces; i++)
    {
        if(face_counts[i] > 1)
        {
            std::cerr<<"FACE "<<i<<" MENTIONED "<<face_counts[i]<<" TIMES IN PATCHES FILE!"<<std::endl;
            abort_flag = true;
        }
    }

    // topological conditions for cells and mesh boundary
    // are implemented independently for FaceMesh2D and FaceMesh3D

    // abort if flagged
    if(abort_flag)
        std::abort();
}

void FaceMesh::check_parameters()
{
    // initialisations
    bool abort_flag = false;

    // collect all boundary faces
    std::vector<int> boundary_face_ids;
    boundary_face_ids.reserve(num_faces);

    for(int i=0; i<num_patches; i++)
    {
        std::vector<int> patch_face_ids = patches[i].get_face_ids();
        boundary_face_ids.insert(boundary_face_ids.end(), patch_face_ids.begin(), patch_face_ids.end());
    }

    // ------------------------ CHECKS ------------------------

    // check all faces have assigned owners
    std::vector<int> owner_ids = get_owner_ids();

    for(int i=0; i<num_faces; i++)
    {
        if(owner_ids[i]<0)
        {
            std::cerr<<"INVALID OWNER ID:\t"<<owner_ids[i]<<"\t(FACE "<<i<<")!"<<std::endl;
            abort_flag = true;
        }
    }

    // check all internal faces have assigned neighbours
    std::vector<int> boundary_face_ids_copy(boundary_face_ids);

    for(int i=0; i<num_faces; i++)
    {
        // memorise iterator where i is found
        std::vector<int>::iterator find_i = std::find(boundary_face_ids_copy.begin(), boundary_face_ids_copy.end(), i);

        if(find_i == boundary_face_ids_copy.end())
        {
            if(faces[i].get_neighbour_id() == -1)
            {
                std::cerr<<"INTERNAL FACE "<<i<<" HAS -1 NEIGHBOUR!"<<std::endl;
                abort_flag = true;
            }
        }
        else
            boundary_face_ids_copy.erase(find_i);  // erase boundary face for easier lookup in later iterations
    }
    
    // no need to check that no face has more than 2 neighbours;
    // already performed when calling calc_owners_neighbours()

    // check for unconnected points;
    // since the topological checks were performed in check_mesh(), there can be no dangling points
    std::vector<int> point_counts(num_points, 0);

    for(int i=0; i<num_faces; i++)
    {
        for(int point_id: faces[i].get_point_ids())
            point_counts[point_id]++;
    }

    for(int i=0; i<num_points; i++)
    {
        if(point_counts[i] == 0)
        {
            std::cerr<<"UNCONNECTED POINT! ID:\t"<<i<<std::endl;
            abort_flag = true;
        }
    }

    // check all cells have total face area vector 0
    for(int i=0; i<num_cells; i++)
    {
        Point total_area_vector;

        for(int face_id: cells[i].get_face_ids())
        {
            Point face_area_vector = faces[face_id].get_area_vector();

            // negate if neighbour
            if(faces[face_id].get_neighbour_id() == i)
                face_area_vector.negate();

            total_area_vector.move_by(face_area_vector);
        }

        std::array<double, N_DIMS> total_area_vector_coords = total_area_vector.get_coordinates();

        if(total_area_vector.calc_modulus() > EPSILON)
        {
            std::cerr<<"CELL "<<i<<" HAS NON-ZERO TOTAL FACE AREA VECTOR: { ";
            for(Dimensions dim: ALL_DIMS)
                std::cerr<<total_area_vector_coords[dim]<<" ";
            std::cerr<<" }!"<<std::endl;

            std::cerr<<"TOLERANCE: "<<EPSILON<<std::endl;

            abort_flag = true;
        }
    }

    // check mesh boundary has total face area vector 0
    Point total_boundary_area_vector;

    for(int face_id: boundary_face_ids)
    {
        Point face_area_vector = faces[face_id].get_area_vector();

        // no need to negate as all face area vectors at boundary point away from mesh

        total_boundary_area_vector.move_by(face_area_vector);
    }

    std::array<double, N_DIMS> total_boundary_area_vector_coords = total_boundary_area_vector.get_coordinates();

    if(total_boundary_area_vector.calc_modulus() > EPSILON)
    {
        std::cerr<<"MESH BOUNDARY HAS NON-ZERO TOTAL FACE AREA VECTOR: { ";
        for(Dimensions dim: ALL_DIMS)
            std::cerr<<total_boundary_area_vector_coords[dim]<<" ";
        std::cerr<<" }!"<<std::endl;

        std::cerr<<"TOLERANCE: "<<EPSILON<<std::endl;

        abort_flag = true;
    }

    // abort if flagged
    if(abort_flag)
        std::abort();
}



std::string FaceMesh::trim(std::string& str, std::vector<char>& lchars, std::vector<char>& rchars)
{
    str = ltrim(str, lchars);
    str = rtrim(str, rchars);

    return str;
}

std::string FaceMesh::ltrim(std::string& str, std::vector<char>& chars)
{
    int begin, str_len = str.length();

    for(begin=0; begin<str_len; begin++)
        if(std::find(chars.begin(), chars.end(), str[begin]) == chars.end())
            break;

    return str.substr(begin, str_len-begin);
}

std::string FaceMesh::rtrim(std::string& str, std::vector<char>& chars)
{
    int end, str_len = str.length();

    for(end=str_len-1; end>=0; end--)
        if(std::find(chars.begin(), chars.end(), str[end]) == chars.end())
            break;

    return str.substr(0, end+1);
}



// private

// face methods

void FaceMesh::calc_face_centroids()
{
    for(int face_id=0; face_id<num_faces; face_id++)
    {
        Point centroid = calc_face_centroid(faces[face_id]);
        faces[face_id].set_centroid(centroid);
    }
}

// cell methods

void FaceMesh::calc_cell_centroids()
{
    for(int cell_id=0; cell_id<num_cells; cell_id++)
    {
        Point centroid = calc_cell_centroid(cells[cell_id]);
        cells[cell_id].set_centroid(centroid);
    }
}

// other

void FaceMesh::read_from_input()
{
    // declarations
    std::string points_path, faces_path, cells_path, patches_path;
    std::ifstream points_file, faces_file, cells_file, patches_file;

    // read points
    std::cout<<"Write file path for points:"<<std::endl;
    std::cin>>points_path;
    points_file.open(points_path);

    read_points(points_file);

    // read faces
    std::cout<<"Write file path for faces:"<<std::endl;
    std::cin>>faces_path;
    faces_file.open(faces_path);

    read_faces(faces_file);

    // read cells
    std::cout<<"Write file path for cells:"<<std::endl;
    std::cin>>cells_path;
    cells_file.open(cells_path);

    read_cells(cells_file);

    // read patches
    std::cout<<"Write file path for boundary patches:"<<std::endl;
    std::cin>>patches_path;
    patches_file.open(patches_path);

    read_patches(patches_file);
}

void FaceMesh::read_points(std::ifstream& points_file)
{
    // declarations
    std::string line;
    std::array<double, N_DIMS> coords;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t', '('},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of points
    std::getline(points_file, line);
    
    num_points = std::stoi(line);
    points.reserve(num_points);

    // check for open bracket to begin read
    std::getline(points_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_points; i++)
    {
        std::getline(points_file, line);

        // trim line from left and right
        line = trim(line, lchars, rchars);

        std::istringstream isstream(line);

        // read point coordinates
        for(Dimensions dim: ALL_DIMS)
            isstream >> coords[dim];

        // assert reaching end of line - it was fully read
        assert(isstream.peek() == -1);

        // add point to points vector
        Point p(coords);
        points.push_back(p);
    }

    // check for close bracket to end read
    std::getline(points_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close file
    points_file.close();
}

void FaceMesh::read_faces(std::ifstream& faces_file)
{
    // declarations
    std::string line;
    std::vector<int> face_point_ids;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t'},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of faces
    std::getline(faces_file, line);
    
    num_faces = std::stoi(line);
    faces.reserve(num_faces);

    // check for open bracket to begin read
    std::getline(faces_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_faces; i++)
    {
        std::getline(faces_file, line);

        // trim line from left and right
        line = trim(line, lchars, rchars);

        std::istringstream isstream(line);

        // read number of points in face
        int num_face_points;
        isstream >> num_face_points;
        face_point_ids.resize(num_face_points);

        // check next character is '('
        assert(isstream.get() == '(');

        // read point ids of face
        for(int j=0; j<num_face_points; j++)
            isstream >> face_point_ids[j];

        // assert reaching end of line - it was fully read
        assert(isstream.peek() == -1);

        // initialise face and add it to faces vector
        Face f(face_point_ids);
        faces.push_back(f);
    }

    // check for close bracket to end read
    std::getline(faces_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close file
    faces_file.close();
}

void FaceMesh::read_cells(std::ifstream& cells_file)
{
    // declarations
    std::string line;
    std::vector<int> cell_face_ids;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t'},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of cells
    std::getline(cells_file, line);
    
    num_cells = std::stoi(line);
    cells.reserve(num_cells);

    // check for open bracket to begin read
    std::getline(cells_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_cells; i++)
    {
        std::getline(cells_file, line);

        // trim line from left and right
        line = trim(line, lchars, rchars);

        std::istringstream isstream(line);

        // read number of faces in cell
        int num_cell_faces;
        isstream >> num_cell_faces;
        cell_face_ids.resize(num_cell_faces);

        // check next character is '('
        assert(isstream.get() == '(');

        // read face ids of cell
        for(int j=0; j<num_cell_faces; j++)
            isstream >> cell_face_ids[j];

        // assert reaching end of line - it was fully read
        assert(isstream.peek() == -1);

        // initialise cell and add it to cells vector
        Cell cell(cell_face_ids);
        cells.push_back(cell);
    }

    // check for close bracket to end read
    std::getline(cells_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close file
    cells_file.close();
}

void FaceMesh::read_patches(std::ifstream& patches_file)
{
    // declarations
    std::string line, patch_name;
    std::vector<int> patch_face_ids;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t'},
                        rchars{' ', '\t', '\r'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of patches
    std::getline(patches_file, line);
    
    num_patches = std::stoi(line);
    patches.reserve(num_patches);

    // check for open bracket to begin read
    std::getline(patches_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_patches; i++)
    {
        // read patch name
            std::getline(patches_file, line);

            // trim line from left and right
            line = trim(line, lchars, rchars);

            patch_name = line;

        // read number of boundary faces
            std::getline(patches_file, line);

            // trim line from left and right
            line = trim(line, lchars, rchars);
            int num_bdr_faces = std::stoi(line);
            patch_face_ids.resize(num_bdr_faces);

        // read next line, trim and check if it is "("
            std::getline(patches_file, line);

            // trim line from left and right
            line = trim(line, lchars, rchars);

            assert(line == "(");

        // read face ids at patch
            std::getline(patches_file, line);

            // trim line from left and right
            line = trim(line, lchars, rchars);

            std::istringstream isstream(line);
            
            for(int j=0; j<num_bdr_faces; j++)
                isstream >> patch_face_ids[j];

            // assert reaching end of line - it was fully read
            assert(isstream.peek() == -1);

            // initialise patch and add it to patches vector
            BoundaryPatch patch(patch_name, patch_face_ids);
            patches.push_back(patch);

        // read next line, trim and check if it is ")"
            std::getline(patches_file, line);

            // trim line from left and right
            line = trim(line, lchars, rchars);

            assert(line == ")");
    }

    // check for close bracket to end read
    std::getline(patches_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close file
    patches_file.close();
}



void FaceMesh::calc_parameters()
{
    // calling functions in this order is crucial as parameters from earlier calls
    // are used in later calls
    calc_owners_neighbours();
    calc_face_centroids();
    calc_cell_centroids();
    calc_face_area_vectors();
    calc_volumes();
}

void FaceMesh::calc_owners_neighbours()
{
    // loop through cells, assign owners and neighbours to faces
    for(int cell_id=0; cell_id<num_cells; cell_id++)
    {
        for(int face_id: cells[cell_id].get_face_ids())
        {
            if(faces[face_id].get_owner_id() == -1) // assign owner id if unassigned
                faces[face_id].set_owner_id(cell_id);
            else if(faces[face_id].get_neighbour_id() == -1)    // assign neighbour id if unassigned
                faces[face_id].set_neighbour_id(cell_id);
            else
            {
                std::cerr<<"FACE "<<face_id<<" IS THE FACE OF MORE THAN 2 CELLS!\n"
                            <<"CELL IDs: "<<faces[face_id].get_owner_id()<<", "
                            <<faces[face_id].get_neighbour_id()<<", "
                            <<cell_id<<std::endl;
                std::abort();
            }
        }
    }
}

void FaceMesh::calc_face_area_vectors()
{
    // loop through faces, set face area vectors
    for(int face_id=0; face_id<num_faces; face_id++)
    {
        Point face_area_vector = calc_face_area_vector(faces[face_id]);
        faces[face_id].set_area_vector(face_area_vector);
    }
}

void FaceMesh::calc_volumes()
{
    // loop through cells, set volumes
    for(int cell_id=0; cell_id<num_cells; cell_id++)
    {
        double volume = calc_cell_volume(cells[cell_id]);
        cells[cell_id].set_volume(volume);
    }
}



// ----------------------------- FaceMesh2D methods -----------------------------

// public

FaceMesh2D::FaceMesh2D(): FaceMesh()
{
    // check correct number of dimensions
    assert(N_DIMS == 2);
}



// protected

void FaceMesh2D::check_mesh()
{
    // initialisations
    bool abort_flag = false;

    // perform usual checks
    FaceMesh::check_mesh();

    // check every face has only 2 points defining it
    for(int face_id=0; face_id<num_faces; face_id++)
    {
        int num_face_points = faces[face_id].get_point_ids().size();
        if(num_face_points != NUM_FACE_POINTS_2D)
        {
            std::cerr<<"2D FACE "<<face_id<<" CONSISTS OF "<<num_face_points<<" != "<<NUM_FACE_POINTS_2D<<" POINTS!"<<std::endl;
            abort_flag = true;
        }
    }

    // perform 2D topological checks

    // topological closedness of cell:
    // within each 2D cell, check that each point is used exactly twice
    for(int cell_id=0; cell_id<num_cells; cell_id++)
    {
        std::vector<int> cell_face_ids = cells[cell_id].get_face_ids();
        int num_cell_faces = cell_face_ids.size();

        int total_points_in_cell = num_cell_faces;  // due to all faces being connected, the total number of points
                                                    // is the same as the total number of faces
        std::vector<int> cell_point_ids, cell_point_counts;
        cell_point_ids.reserve(total_points_in_cell);
        cell_point_counts.resize(total_points_in_cell, 0);    // initialised to 0 for counting

        for(int face_id: cell_face_ids)   // loop through faces of cell to populate cell_point_ids
        {
            for(int point_id: faces[face_id].get_point_ids())  // loop through points of face_id
            {
                bool find_flag = false;

                int num_cell_points = cell_point_ids.size();    // current num of points in cell

                for(int p=0; p<num_cell_points; p++)    // loop through cell points already written
                {
                    if(point_id == cell_point_ids[p])   // if found, raise flag and increase count
                    {
                        find_flag = true;
                        cell_point_counts[p]++;
                        break;
                    }
                }

                if(!find_flag)  // if not found, add point_id to cell_point_ids and increase its count
                {
                    cell_point_ids.push_back(point_id);
                    cell_point_counts[num_cell_points]++;
                }
            }
        }

        for(int p=0; p<total_points_in_cell; p++) // loop to check all point counts are exactly 2
        {
            if(cell_point_counts[p] != POINT_COUNT_VAL_2D)
            {
                std::cerr<<"CELL "<<cell_id<<" CONTAINS POINT "<<cell_point_ids[p]<<" A NUMBER OF "
                        <<cell_point_counts[p]<<" TIMES!"<<std::endl;
                abort_flag = true;
            }
        }
    }

    // topological closedness of mesh boundary:
    // within all boundary faces, check that each edge is used exactly twice

    // get total number of points in the boundary
    int total_points_in_boundary = 0;
    for(int patch_id=0; patch_id<num_patches; patch_id++)
        total_points_in_boundary += patches[patch_id].get_face_ids().size();    // due to all faces being connected,
                                                                                // the total number of points is the same
                                                                                // as the total number of faces

    std::vector<int> patch_point_ids, patch_point_counts;
    patch_point_ids.reserve(total_points_in_boundary);
    patch_point_counts.resize(total_points_in_boundary, 0);    // initialised to 0 for counting

    // loop through patches to populate patch_point_ids
    for(int patch_id=0; patch_id<num_patches; patch_id++)
    {
        std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();

        for(int face_id: patch_face_ids)   // loop through faces of patch to populate patch_point_ids
        {
            for(int point_id: faces[face_id].get_point_ids())  // loop through points of face_id
            {
                bool find_flag = false;

                int num_patch_points = patch_point_ids.size();    // current num of points in patch

                for(int p=0; p<num_patch_points; p++)    // loop through patch points already written
                {
                    if(point_id == patch_point_ids[p])   // if found, raise flag and increase count
                    {
                        find_flag = true;
                        patch_point_counts[p]++;
                        break;
                    }
                }

                if(!find_flag)  // if not found, add point_id to patch_point_ids and increase its count
                {
                    patch_point_ids.push_back(point_id);
                    patch_point_counts[num_patch_points]++;
                }
            }
        }
    }

    // loop to check all point counts are exactly 2
    for(int p=0; p<total_points_in_boundary; p++)
    {
        if(patch_point_counts[p] != POINT_COUNT_VAL_2D)
        {
            std::cerr<<"BOUNDARY CONTAINS POINT "<<patch_point_ids[p]<<" A NUMBER OF "
                    <<patch_point_counts[p]<<" TIMES!"<<std::endl;
            abort_flag = true;
        }
    }

    // --------------------------------------------

    // abort if flagged
    if(abort_flag)
        std::abort();
}



// private

// face methods

Point FaceMesh2D::calc_face_centroid(Face& face)
{
    /* Although the algorithm is independent of the number of faces, the 2-point check
     * was performed in check_mesh();
     * as such, the face centroid is a plain average, no adjustments needed
     */

    // initialisations
    std::vector<int> face_point_ids = face.get_point_ids();
    int num_face_points = face_point_ids.size();
    std::vector<Point> face_points;
    face_points.resize(num_face_points);
    Point centroid;

    for(int i=0; i<num_face_points; i++)
        face_points[i] = points[face_point_ids[i]];

    // sum up all points coordinates into centroid
    for(Point p: face_points)
        centroid.move_by(p);

    // divide by number of points to get centroid
    centroid.scale_by(1./num_face_points);

    return centroid;
}

Point FaceMesh2D::calc_face_area_vector(Face& face)
{
    // note: every 2D every face has 2 points defining it;
    //          check performed in check_mesh()

    // initialisations
    std::vector<int> face_point_ids = face.get_point_ids();
    std::array<double, N_DIMS> coords0 = points[face_point_ids[0]].get_coordinates();
    std::array<double, N_DIMS> coords1 = points[face_point_ids[1]].get_coordinates();
    std::array<double, N_DIMS> coords_owner_centroid = cells[face.get_owner_id()].get_centroid().get_coordinates();
    std::array<double, N_DIMS> face_plane_coords, area_vector_coords, dist_centroid_coords0;

    for(Dimensions dim: ALL_DIMS)
    {
        face_plane_coords[dim] = coords1[dim] - coords0[dim];
        dist_centroid_coords0[dim] = coords0[dim] - coords_owner_centroid[dim];
    }

    // get area_vector_coords by rotating counterclockwise;
    // assign each coordinate separately in case of dimension misalignment
    area_vector_coords[X] = - face_plane_coords[Y];
    area_vector_coords[Y] = face_plane_coords[X];

    Point area_vector(area_vector_coords);
    Point dist_vector(dist_centroid_coords0);

    // if area_vector_coords is pointing towards owner centroid,
    // multiply coordinates by -1
    if(area_vector.dot_product_with(dist_vector) < 0)
        area_vector.negate();


    return area_vector;
}

// cell methods

Point FaceMesh2D::calc_cell_centroid(Cell& cell)
{
    // initialisations
    std::vector<int> cell_face_ids = cell.get_face_ids();
    int num_cell_faces = cell_face_ids.size();
    Point init_centroid, centroid;
    double averaging_area = 0;

    // sum up all face centroids into init_centroid;
    // in 2D equivalent to summing all vertices
    for(int face_id: cell_face_ids)
    {
        Point face_centroid = faces[face_id].get_centroid();
        init_centroid.move_by(face_centroid);
    }

    // divide by number of vertices (=num_cell_faces)
    init_centroid.scale_by(1./num_cell_faces);

    // loop through faces, calculate centroid of triangle connecting face to init_centroid,
    // add weighted by area to centroid
    for(int face_id: cell_face_ids)
    {
        // initialisations
        Point centroid_i;
        std::vector<int> face_point_ids = faces[face_id].get_point_ids();  // must have NUM_FACE_POINTS_2D=2 points
        double area_i = 0;

        // populate centroid_i
        centroid_i.move_by(init_centroid);

        for(int point_id: face_point_ids)
            centroid_i.move_by(points[point_id]);

        centroid_i.scale_by(1./3);

        // calculate area of triangle made from init_centroid and face_points
        area_i = calc_triangle_area_2D(init_centroid,
                                        points[face_point_ids[0]],
                                        points[face_point_ids[1]]);

        // make updates to centroid and area
        averaging_area += area_i;

        centroid_i.scale_by(area_i);
        centroid.move_by(centroid_i);
    }

    // average centroid for final value
    centroid.scale_by(1./averaging_area);

    return centroid;
}

double FaceMesh2D::calc_cell_volume(Cell& cell)
{
    // initialisations
    std::vector<int> cell_face_ids = cell.get_face_ids();
    int num_faces = cell_face_ids.size();
    Point cell_centroid = cell.get_centroid();
    double volume = 0;

    // calculate volume of i-th pyramid (in the 2D case, area of i-th triangle)
    for(int i=0; i<num_faces; i++)
    {
        std::vector<int> face_i_point_ids = faces[cell_face_ids[i]].get_point_ids();
        
        // every 2D face only has 2 points at each face;
        // we calculate the area of the triangle in 2D
        volume += calc_triangle_area_2D(cell_centroid,
                                        points[face_i_point_ids[0]],
                                        points[face_i_point_ids[1]]);
    }

    return volume;
}

// other

double FaceMesh2D::calc_triangle_area_2D(Point& a, Point& b, Point& c)
{
    // extract coordinates
    std::array<double, N_DIMS> coords_a = a.get_coordinates();
    std::array<double, N_DIMS> coords_b = b.get_coordinates();
    std::array<double, N_DIMS> coords_c = c.get_coordinates();

    double area = 0.5 * fabs( coords_a[X] * (coords_b[Y] - coords_c[Y])
                            + coords_b[X] * (coords_c[Y] - coords_a[Y])
                            + coords_c[X] * (coords_a[Y] - coords_b[Y]) );

    return area;
}



// ----------------------------- FaceMesh3D methods -----------------------------

// public

FaceMesh3D::FaceMesh3D(): FaceMesh()
{
    // check correct number of dimensions
    assert(N_DIMS == 3);
}



// protected

void FaceMesh3D::check_mesh()
{
    // initialisations
    bool abort_flag = false;

    // perform usual checks
    FaceMesh::check_mesh();

    // perform 3D topological checks

    // topological closedness of cell:
    // within each cell, check that each edge is used exactly twice
    for(int cell_id=0; cell_id<num_cells; cell_id++)
    {
        std::vector<int> cell_face_ids = cells[cell_id].get_face_ids();

        // calculate number of edges in cell
        int total_edges_in_cell = calc_num_edges(cells[cell_id]);

        std::vector<std::array<int, NUM_EDGE_POINTS_3D>> cell_edges;
        std::vector<int> cell_edge_counts;
        cell_edges.reserve(total_edges_in_cell);
        cell_edge_counts.resize(total_edges_in_cell, 0);    // initialised to 0 for counting

        for(int face_id: cell_face_ids)   // loop through faces of cell to populate cell_edges
        {
            // get edges of face
            std::vector<std::array<int, NUM_EDGE_POINTS_3D>> face_edges = calc_edges(faces[face_id]);

            // get points of face
            std::vector<int> face_point_ids = faces[face_id].get_point_ids();
            int num_face_points = face_point_ids.size();

            for(int p=0; p<num_face_points; p++)  // loop through points of face_id
            {
                bool find_flag = false;

                int num_cell_edges = cell_edges.size();    // current num of edges in cell

                // edge to compare to when iterating through cell_edges and its reverse
                std::array<int, NUM_EDGE_POINTS_3D> ref_edge {
                    face_point_ids[p], face_point_ids[(p+1)%num_face_points]
                };
                std::array<int, NUM_EDGE_POINTS_3D> ref_edge_reverse {
                    ref_edge[1], ref_edge[0]
                };

                for(int e=0; e<num_cell_edges; e++)    // loop through cell edges already written
                {
                    // if found, raise flag and increase count
                    if(cell_edges[e] == ref_edge || cell_edges[e] == ref_edge_reverse)
                    {
                        find_flag = true;
                        cell_edge_counts[e]++;
                        break;
                    }
                }

                if(!find_flag)  // if not found, add ref_edge to cell_edges and increase its count
                {
                    cell_edges.push_back(ref_edge);
                    cell_edge_counts[num_cell_edges]++;
                }
            }
        }

        for(int e=0; e<total_edges_in_cell; e++) // loop to check all edge counts are exactly 2
        {
            if(cell_edge_counts[e] != EDGE_COUNT_VAL_3D)
            {
                std::cerr<<"CELL "<<cell_id<<" CONTAINS EDGE { ";
                for(int point_id: cell_edges[e])
                    std::cerr<<point_id<<" ";
                std::cerr<<" } A NUMBER OF "<<cell_edge_counts[e]<<" TIMES!"<<std::endl;
                abort_flag = true;
            }
        }
    }


    // topological closedness of mesh boundary:
    // within all boundary faces, check that each edge is used exactly twice

    int total_edges_in_boundary = calc_num_edges();
    std::vector<std::array<int, NUM_EDGE_POINTS_3D>> boundary_edges;
    std::vector<int> boundary_edge_counts;
    boundary_edges.reserve(total_edges_in_boundary);
    boundary_edge_counts.resize(total_edges_in_boundary, 0);    // initialised to 0 for counting

    for(int patch_id=0; patch_id<num_patches; patch_id++)
    {
        std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();

        for(int face_id: patch_face_ids)   // loop through faces of patch to populate boundary_edges
        {
            // get edges of face
            std::vector<std::array<int, NUM_EDGE_POINTS_3D>> face_edges = calc_edges(faces[face_id]);

            // get points of face
            std::vector<int> face_point_ids = faces[face_id].get_point_ids();
            int num_face_points = face_point_ids.size();

            for(int p=0; p<num_face_points; p++)  // loop through points of face_id
            {
                bool find_flag = false;

                int num_boundary_edges = boundary_edges.size();    // current num of edges in boundary

                // edge to compare to when iterating through boundary_edges and its reverse
                std::array<int, NUM_EDGE_POINTS_3D> ref_edge {
                    face_point_ids[p], face_point_ids[(p+1)%num_face_points]
                };
                std::array<int, NUM_EDGE_POINTS_3D> ref_edge_reverse {
                    ref_edge[1], ref_edge[0]
                };

                for(int e=0; e<num_boundary_edges; e++)    // loop through boundary edges already written
                {
                    // if found, raise flag and increase count
                    if(boundary_edges[e] == ref_edge || boundary_edges[e] == ref_edge_reverse)
                    {
                        find_flag = true;
                        boundary_edge_counts[e]++;
                        break;
                    }
                }

                if(!find_flag)  // if not found, add ref_edge to boundary_edges and increase its count
                {
                    boundary_edges.push_back(ref_edge);
                    boundary_edge_counts[num_boundary_edges]++;
                }
            }
        }
    }

    for(int e=0; e<total_edges_in_boundary; e++) // loop to check all edge counts are exactly 2
    {
        if(boundary_edge_counts[e] != EDGE_COUNT_VAL_3D)
        {
            std::cerr<<"BOUNDARY CONTAINS EDGE { ";
            for(int point_id: boundary_edges[e])
                std::cerr<<point_id<<" ";
            std::cerr<<" } A NUMBER OF "<<boundary_edge_counts[e]<<" TIMES!"<<std::endl;
            abort_flag = true;
        }
    }

    // --------------------------------------------

    // abort if flagged
    if(abort_flag)
        std::abort();
}



// private

// face methods

Point FaceMesh3D::calc_face_centroid(Face& face)
{
    // initialisations
    std::vector<int> face_point_ids = face.get_point_ids();
    int num_face_points = face_point_ids.size();
    std::vector<Point> face_points;
    face_points.resize(num_face_points);
    Point init_centroid, centroid;
    double averaging_area = 0;

    for(int i=0; i<num_face_points; i++)
        face_points[i] = points[face_point_ids[i]];

    // sum up all points coordinates into init_centroid
    for(Point p: face_points)
        init_centroid.move_by(p);

    // divide by number of points to get init_centroid
    init_centroid.scale_by(1./num_face_points);

    // create triangles from init_centroid to each edge, calculate its centroid,
    // add weighted by area to centroid
    for(int i=0; i<num_face_points; i++)
    {
        Point centroid_i;
        double area_i = 0;

        // populate centroid_i
        centroid_i.move_by(init_centroid);
        centroid_i.move_by(face_points[i]);
        centroid_i.move_by(face_points[(i+1)%num_face_points]);
        centroid_i.scale_by(1./3);

        // calculate area_i
        area_i = calc_triangle_area_3D(init_centroid,
                                        face_points[i],
                                        face_points[(i+1)%num_face_points]).calc_modulus();

        // perform update
        averaging_area += area_i;
        centroid_i.scale_by(area_i);
        centroid.move_by(centroid_i);
    }

    // average centroid for final value
    centroid.scale_by(1./averaging_area);

    return centroid;
}

Point FaceMesh3D::calc_face_area_vector(Face& face)
{
    // initialisations
    Point face_centroid = face.get_centroid();
    Point owner_centroid = cells[face.get_owner_id()].get_centroid();
    std::vector<int> face_point_ids = face.get_point_ids();
    int num_face_points = face_point_ids.size();

    // calculate distance vector from owner_centroid to face_centroid
    Point dist_vector = Point::distance(owner_centroid, face_centroid);

    // calculate total area vector of face
    Point total_area_vector;

    // create triangles from cell centroid to every face edge and sum areas of all of those triangles
    for(int i=0; i<num_face_points; i++)
    {
        Point area_vector_i = calc_triangle_area_3D(face_centroid,
                                                    points[face_point_ids[i]],
                                                    points[face_point_ids[(i+1)%num_face_points]]);
        total_area_vector.move_by(area_vector_i);
    }

    // if total_area_vector is pointing towards owner centroid,
    // multiply coordinates by -1
    if(total_area_vector.dot_product_with(dist_vector) < 0)
        total_area_vector.negate();

    return total_area_vector;
}

std::vector<std::array<int, NUM_EDGE_POINTS_3D>> FaceMesh3D::calc_edges(Face& face)
{
    // initialisations
    std::vector<int> face_point_ids = face.get_point_ids();
    int num_face_points = face_point_ids.size();

    std::vector<std::array<int, NUM_EDGE_POINTS_3D>> edges;
    edges.reserve(num_face_points);

    // populate edges vector
    for(int i=0; i<num_face_points; i++)
        edges.push_back(
            std::array<int, NUM_EDGE_POINTS_3D>{
                face_point_ids[i], face_point_ids[(i+1)%num_face_points]
            }
        );

    return edges;
}

// cell methods

Point FaceMesh3D::calc_cell_centroid(Cell& cell)
{
    // initialisations
    std::vector<int> cell_face_ids = cell.get_face_ids();
    int num_cell_faces = cell_face_ids.size();
    Point init_centroid, centroid;
    double averaging_volume = 0;

    // sum up all face centroids into init_centroid
    for(int i=0; i<num_cell_faces; i++)
    {
        Point centroid_i = faces[cell_face_ids[i]].get_centroid();
        init_centroid.move_by(centroid_i);
    }

    // divide by number of init_centroids (=num_cell_faces)
    init_centroid.scale_by(1./num_cell_faces);

    for(int face_id: cell_face_ids)
    {
        std::vector<int> face_point_ids = faces[face_id].get_point_ids();
        int num_face_points = face_point_ids.size();
        Point face_centroid = faces[face_id].get_centroid();

        for(int i=0; i<num_face_points; i++)
        {
            Point tet_centroid;
            double tet_volume = 0;

            // tet_centroid is the average of the 4 vertices
            tet_centroid.move_by(init_centroid);
            tet_centroid.move_by(face_centroid);
            tet_centroid.move_by(points[face_point_ids[i]]);
            tet_centroid.move_by(points[face_point_ids[(i+1)%num_face_points]]);
            tet_centroid.scale_by(1./4);

            // calculate the volume of the tetrahedron
            tet_volume = calc_tet_volume_3D(init_centroid,
                                            face_centroid,
                                            points[face_point_ids[i]],
                                            points[face_point_ids[(i+1)%num_face_points]]);

            // perform update
            averaging_volume += tet_volume;
            tet_centroid.scale_by(tet_volume);
            centroid.move_by(tet_centroid);
        }
    }

    // average centroid for final value
    centroid.scale_by(1./averaging_volume);

    return centroid;
}

double FaceMesh3D::calc_cell_volume(Cell& cell)
{
    // initialisations
    std::vector<int> cell_face_ids = cell.get_face_ids();
    Point cell_centroid = cell.get_centroid();
    double volume = 0;

    // calculate volume of pyramid with base face_id via tetrahedra and add it to volume
    for(int face_id: cell_face_ids)
    {
        Point face_centroid = faces[face_id].get_centroid();
        std::vector<int> face_point_ids = faces[face_id].get_point_ids();
        int num_face_points = face_point_ids.size();

        for(int i=0; i<num_face_points; i++)
        {
            volume += calc_tet_volume_3D(cell_centroid,
                                            face_centroid,
                                            points[face_point_ids[i]],
                                            points[face_point_ids[(i+1)%num_face_points]]);
        }
    }

    return volume;
}

int FaceMesh3D::calc_num_edges(Cell& cell)
{
    // calculate the number of edges in a cell via Euler's formula

    // extract all points in cell and calculate the size of that vector
    std::vector<int> cell_point_ids;
    cell_point_ids.reserve(num_points); // the real number is likely going to be less,
                                        // but this is a safe upper bound

    for(int face_id: cell.get_face_ids())  // loop through faces of cell to populate cell_point_ids
    {
        for(int point_id: faces[face_id].get_point_ids())   // loop through points of face_id to populate cell_point_ids
        {
            bool find_flag = false;
            int num_cell_points = cell_point_ids.size();    // current num of points in cell

            // look for point_id in cell_point_ids
            for(int p=0; p<num_cell_points; p++)
            {
                if(point_id == cell_point_ids[p])
                {
                    find_flag = true;
                    break;
                }
            }

            if(!find_flag)
                cell_point_ids.push_back(point_id); // if not found, add point_id to cell_point_ids
        }
    }
    
    // calculate number of points in cell
    int num_cell_points = cell_point_ids.size();

    // calculate the number of faces in the cell
    int num_cell_faces = cell.get_face_ids().size();

    // calculate total number of edges in cell via Euler's formula
    int total_edges_in_cell = num_cell_points + num_cell_faces - 2;

    return total_edges_in_cell;
}

// other

Point FaceMesh3D::calc_triangle_area_3D(Point& a, Point& b, Point& c)
{
    // get coordinates of i-th triangle
    std::array<double, N_DIMS> coords_a = a.get_coordinates();
    std::array<double, N_DIMS> coords_b = b.get_coordinates();
    std::array<double, N_DIMS> coords_c = c.get_coordinates();

    // calculate vectors of two triangle edges in 3D
    std::array<double, 3> vec1{0, 3}, vec2{0, 3};   // N_DIMS not used in case of dimension misalignment
    for(Dimensions dim: ALL_DIMS)
    {
        vec1[dim] = coords_c[dim] - coords_a[dim];
        vec2[dim] = coords_b[dim] - coords_a[dim];
    }

    // calculate area vector coordinates of i-th triangle;
    // N_DIMS not used in case of dimension misalignment
    std::array<double, 3> cross_product_coords = {(vec1[Y]*vec2[Z] - vec1[Z]*vec2[Y]),
                                                (vec1[Z]*vec2[X] - vec1[X]*vec2[Z]),
                                                (vec1[X]*vec2[Y] - vec1[Y]*vec2[X])};

    // populate a N_DIMS coordinates array used for output alignment
    std::array<double, N_DIMS> area_vector_coords;

    for(Dimensions dim: ALL_DIMS)
        area_vector_coords[dim] = 0.5 * cross_product_coords[dim];

    return Point(area_vector_coords);
}

double FaceMesh3D::calc_tet_volume_3D(Point& a, Point& b, Point& c, Point& d)
{
    // v1 = b-a
    Point v1 = Point::distance(a, b);

    Point acd_area_vector = calc_triangle_area_3D(a, c, d);

    double signed_volume = 1./3 * v1.dot_product_with(acd_area_vector);

    return fabs(signed_volume);
}

int FaceMesh3D::calc_num_edges()
{
    // calculate the number of edges in the mesh boundary via Euler's formula

    // extract all points in the boudnary and calculate the size of that vector
    std::vector<int> boundary_point_ids;
    boundary_point_ids.reserve(num_points); // the real number is likely going to be less,
                                            // but this is a safe upper bound
    for(int patch_id=0; patch_id<num_patches; patch_id++)
    {
        for(int face_id: patches[patch_id].get_face_ids())  // loop through faces of patch_id to populate boundary_point_ids
        {
            for(int point_id: faces[face_id].get_point_ids())   // loop through points of face_id to populate boundary_point_ids
            {
                bool find_flag = false;
                int num_boundary_points = boundary_point_ids.size();    // current num of points in boundary

                // look for point_id in boundary_point_ids
                for(int p=0; p<num_boundary_points; p++)
                {
                    if(point_id == boundary_point_ids[p])
                    {
                        find_flag = true;
                        break;
                    }
                }

                if(!find_flag)
                    boundary_point_ids.push_back(point_id); // if not found, add point_id to boundary_point_ids
            }
        }
    }
    
    // calculate number of points in boundary
    int num_boundary_points = boundary_point_ids.size();

    // calculate the number of faces in the boundary
    int num_boundary_faces = 0;
    for(int patch_id=0; patch_id<num_patches; patch_id++)
        num_boundary_faces += patches[patch_id].get_face_ids().size();

    // calculate total number of edges in boundary via Euler's formula
    int total_edges_in_boundary = num_boundary_points + num_boundary_faces - 2;

    return total_edges_in_boundary;
}