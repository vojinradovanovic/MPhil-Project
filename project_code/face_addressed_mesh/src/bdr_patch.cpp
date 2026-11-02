#include <string>
#include <vector>

#include "bdr_patch.hpp"

BoundaryPatch::BoundaryPatch(std::string& name, std::vector<int>& face_ids)
{
    this->name = name;
    this->face_ids = face_ids;
}



std::string BoundaryPatch::get_name()
{
    return name;
}

std::vector<int> BoundaryPatch::get_face_ids()
{
    return face_ids;
}