#include <array>
#include <string>
#include <vector>

#include "constants.hpp"
#include "point.hpp"
#include "face.hpp"
#include "cell.hpp"
#include "bdr_patch.hpp"

#ifndef FACE_MESH
#define FACE_MESH

class FaceMesh
{
    // attributes
    protected:
        int num_points;
        int num_faces;
        int num_cells;
        int num_patches;

        std::vector<Point> points;
        std::vector<Face> faces;
        std::vector<Cell> cells;
        std::vector<BoundaryPatch> patches;

    // methods
    public:
        FaceMesh();
        void initialise();

        int get_num_points();
        int get_num_faces();
        int get_num_cells();
        int get_num_patches();

        std::vector<double> get_cell_volumes();
        std::vector<Point> get_face_area_vectors();
        std::vector<Point> get_cell_centroids();
        std::vector<Point> get_face_centroids();

        std::vector<std::string> get_patch_names();
        std::vector<int> get_patch_face_ids(int patch_id);
        std::vector<int> get_patch_cell_ids(int patch_id);
        std::vector<int> get_neighbouring_cell_ids(int cell_id);

        std::vector<int> get_owner_ids();
        std::vector<int> get_neighbour_ids();

    protected:
        virtual void check_mesh();   // additional implementation dependent on N_DIMS
        void check_parameters();

        // methods below are protected in case child classes have read methods
        std::string trim(std::string& str, std::vector<char>& lchars, std::vector<char>& rchars);
        std::string ltrim(std::string& str, std::vector<char>& chars);
        std::string rtrim(std::string& str, std::vector<char>& chars);

    private:
        // face methods
            void calc_face_centroids();

            virtual Point calc_face_centroid(Face& face) = 0;       // purely virtual,
            virtual Point calc_face_area_vector(Face& face) = 0;    // implementation depends on N_DIMS

        // cell methods
            void calc_cell_centroids();

            virtual Point calc_cell_centroid(Cell& cell) = 0;       // purely virtual,
            virtual double calc_cell_volume(Cell& cell) = 0;        // implementation depends on N_DIMS

        // other
            void read_from_input();
            void read_points(std::ifstream& points_file);
            void read_faces(std::ifstream& faces_file);
            void read_cells(std::ifstream& cells_file);
            void read_patches(std::ifstream& patches_file);

            void calc_parameters();
            void calc_owners_neighbours();
            void calc_face_area_vectors();
            void calc_volumes();
};

class FaceMesh2D: public FaceMesh
{
    public:
        FaceMesh2D();

    protected:
        void check_mesh() override;

    private:
        // face methods
            Point calc_face_centroid(Face& face) override;
            Point calc_face_area_vector(Face& face) override;

        // cell methods
            Point calc_cell_centroid(Cell& cell) override;
            double calc_cell_volume(Cell& cell) override;

        // other
            double calc_triangle_area_2D(Point& a, Point& b, Point& c);
};

class FaceMesh3D: public FaceMesh
{
    public:
        FaceMesh3D();

    protected:
        void check_mesh() override;

    private:
        // face methods
            Point calc_face_centroid(Face& face) override;
            Point calc_face_area_vector(Face& face) override;
            std::vector<std::array<int, NUM_EDGE_POINTS_3D>> calc_edges(Face& face);

        // cell methods
            Point calc_cell_centroid(Cell& cell) override;
            double calc_cell_volume(Cell& cell) override;
            int calc_num_edges(Cell& cell);

        // other
            Point calc_triangle_area_3D(Point& a, Point& b, Point& c);
            double calc_tet_volume_3D(Point& a, Point& b, Point& c, Point& d);
            int calc_num_edges();   // calculate the number of edges of boundary
};

#endif