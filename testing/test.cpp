
#include "test.hpp"

#include <unordered_map>

struct FaceIndexHash {
    size_t operator()(const pair<int, int>& index) const {
        return std::hash<int>()(index.first) ^ std::hash<int>()(index.second);
    }
};


int main(){
    int n_cells = 6;

    Vector<Cell> cells;
    for (int i=0;i<n_cells;i++){
        cells.push_back(Cell());
    }

    Vector<Face> faces;
    faces.push_back({0,2});
    faces.push_back({1,3});
    faces.push_back({0,1});
    faces.push_back({1,4});
    faces.push_back({0,5});


}