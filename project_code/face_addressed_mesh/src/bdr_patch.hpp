#include <string>
#include <vector>

#ifndef BDR_PATCH
#define BDR_PATCH

class BoundaryPatch
{
    std::string name;
    std::vector<int> face_ids;

    public:
        BoundaryPatch(std::string& name, std::vector<int>& face_ids);

        std::string get_name();
        std::vector<int> get_face_ids();
};

#endif