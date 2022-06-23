//
// Created by agalex on 23/6/2022.
//

#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <functional>
#include "CPointCloud.h"

void readPointCloud(const std::string& fileName, std::vector<CPointEx3D>& points, std::vector<CPointEx3D>& normals)
{
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(fileName, std::ios::in);
    while (std::getline(myfile, line))
        ++number_of_lines;
    std::cout << "reading " << number_of_lines << " points." << std::endl;
    myfile.clear();
    myfile.seekg(0);
    float x, y, z, nx, ny, nz;
    while (myfile >> x >> y >> z >> nx >> ny >> nz) {
        points.emplace_back(x, y, z);
        normals.emplace_back(nx, ny, nz);
    }
}

int main(int argc, char** argv)
{
    std::vector<CPointEx3D> points;
    std::vector<CPointEx3D> normals;
    readPointCloud(argv[1], points, normals);
    std::random_device rd;
    std::mt19937 gen(rd());
    uniform_int_distribution<> distr(0, points.size() - 1); // define the range
    int idx1 = distr(gen);
    int idx2 = distr(gen);
    while (idx1 == idx2)
    {
        idx1 = distr(gen);
        idx2 = distr(gen);
    }

    LIBGEODESIC::CPointCloud gpc(points, normals);
    std::vector<std::vector<int>> EGG;
    gpc.calculateNeighbors(EGG, 40);

    std::vector<CPointEx3D> path;

    //definition of the distance metric also used as the heuristic distance measure.
    auto eucledianDistance = [&gpc](int i, int j) {
        return (gpc.getPnts()[i] - gpc.getPnts()[j]).norm();
    };
    gpc.GeodesicDijkstra(idx1, idx2, EGG, eucledianDistance, path);

    std::cout << "Dijkstra Geodesic" << std::endl;
    for (const auto& pnt : path)
    {
        std::cout << "(" << pnt.x << "," << pnt.y << "," << pnt.z << "), ";
    }
    std::cout << std::endl;

    gpc.GeodesicFlow(path);
    std::cout << "Geodesic correction" << std::endl;
    for (const auto& pnt : path)
    {
        std::cout << "(" << pnt.x << "," << pnt.y << "," << pnt.z << "), ";
    }
    std::cout << std::endl;

    path.clear();

    //definition of the distance metric

    gpc.GeodesicAstar(idx1, idx2, EGG, eucledianDistance, eucledianDistance, path);
    std::cout << "AStar Geodesic" << std::endl;
    for (const auto& pnt : path)
    {
        std::cout << "(" << pnt.x << "," << pnt.y << "," << pnt.z << "), ";
    }
    std::cout << std::endl;


    gpc.GeodesicFlow(path);
    std::cout << "Geodesic correction" << std::endl;
    for (const auto& pnt : path)
    {
        std::cout << "(" << pnt.x << "," << pnt.y << "," << pnt.z << "), ";
    }
    std::cout << std::endl;


}


