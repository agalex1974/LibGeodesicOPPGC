#include "pch.h"
#include <tbb/parallel_for.h>
#include <mkl.h>
#include <stack>
#include <functional>
#include <queue>
#include <vector>
#include <fstream>
#include <iostream>
#include <utility> // defines std::pair
#include <list>
#include <execution>
#include "CPointCloud.h"
#include "cmathutilities.h"

// CGAL HEADERS ///////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
//#include <CGAL/compute_average_spacing.h>
//include <CGAL/remove_outliers.h>
//#include <CGAL/jet_smooth_point_set.h>
//#include <CGAL/grid_simplify_point_set.h>
//#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/point_generators_3.h>

// CGAL TYPES ///////////////////////////////////
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Point_and_int = boost::tuple<Point_3, int>;
using Traits_base = CGAL::Search_traits_3<Kernel>;
using Traits = CGAL::Search_traits_adapter<Point_and_int,
	CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
	Traits_base>;
using Random_points_iterator = CGAL::Random_points_in_cube_3<Point_3>;
using N_Random_points_iterator = CGAL::Counting_iterator<Random_points_iterator>;
using Tree = CGAL::Kd_tree<Traits>;
using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;
using Fuzzy_iso_box = CGAL::Fuzzy_iso_box<Traits>;
using K_neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Distance = K_neighbor_search::Distance;
using Concurrency_tag = CGAL::Parallel_if_available_tag;
// Point_3 with normal vector stored in a std::pair.
using CGALPointVectorPair = std::pair<Point_3, Vector_3>;
using CGALPointList = std::vector<CGALPointVectorPair>;
////////////////////////////////////////////////////////////////////////


using namespace LIBGEODESIC;

CPointCloud::CPointCloud() {
	points = NULL;
	normals = NULL;
	pointCloudSpatialTree = NULL;
	mem_vecs.K = 0;
	mem_vecs.nStar = 0;
	mem_vecs.kG = 0;
	mem_vecs.nT = 0;
	mem_vecs.p = 0;
	mem_vecs.flow = 0;
	mem_vecs.DF = 0;
	mem_vecs.D2F = 0;
	mem_vecs.DF = 0;

}


CPointCloud::CPointCloud(const std::vector<CPointEx3D>& inputPoints, const std::vector<CPointEx3D>& inputNormals) {
	init(inputPoints, inputNormals);
	mem_vecs.K = 0;
	mem_vecs.nStar = 0;
	mem_vecs.kG = 0;
	mem_vecs.nT = 0;
	mem_vecs.p = 0;
	mem_vecs.flow = 0;
	mem_vecs.DF = 0;
	mem_vecs.D2F = 0;
	mem_vecs.DF = 0;
}

void CPointCloud::init(const std::vector<CPointEx3D>& inputPoints, const std::vector<CPointEx3D>& inputNormals) {
    clear(); // make sure we have no doublicate data
	points = new std::vector<CPointEx3D>;
	*points = inputPoints;
	normals = new std::vector<CPointEx3D>;
	*normals = inputNormals;
    std::vector<Point_3> pointsCGAL;
	std::vector<int> indicesCGAL;
	for (int i = 0; i < points->size(); i++)
	{
		pointsCGAL.push_back(Point_3((*points)[i][0], (*points)[i][1], (*points)[i][2]));
		indicesCGAL.push_back(i);
	}
	pointCloudSpatialTree = new Tree(boost::make_zip_iterator(boost::make_tuple(pointsCGAL.begin(), indicesCGAL.begin())),
		boost::make_zip_iterator(boost::make_tuple(pointsCGAL.end(), indicesCGAL.end())));

}

void CPointCloud::init_pnts(const std::vector<CPointEx3D>& inputPoints, bool bCalcNormals) {
	clear(); // make sure we have no doublicate data
	points = new std::vector<CPointEx3D>;
	*points = inputPoints;
	if (bCalcNormals){
		calculateNormals();
		return;
	}

	std::vector<Point_3> pointsCGAL;
	std::vector<int> indicesCGAL;
	for (int i = 0; i < points->size(); i++)
	{
		pointsCGAL.push_back(Point_3((*points)[i][0], (*points)[i][1], (*points)[i][2]));
		indicesCGAL.push_back(i);

	}
	pointCloudSpatialTree = new Tree(boost::make_zip_iterator(boost::make_tuple(pointsCGAL.begin(), indicesCGAL.begin())),
		boost::make_zip_iterator(boost::make_tuple(pointsCGAL.end(), indicesCGAL.end())));
}

CDoubleMatrix CPointCloud::CreateEGGLocalCoordinateSystem(const CPointEx3D& o, const CPointEx3D& d)
{
	CPointEx3D firstPerpendicular = GetNormalizedPerpendicularVectorToVector(d);
	CPointEx3D secondPerpendicular = cross(d, firstPerpendicular);
	secondPerpendicular.normalize();
	CDoubleMatrix transformationMatrix(4, 4);
	transformationMatrix.Identity();
	transformationMatrix(1, 0) = firstPerpendicular[0];
	transformationMatrix(1, 1) = firstPerpendicular[1];
	transformationMatrix(1, 2) = firstPerpendicular[2];
	transformationMatrix(2, 0) = secondPerpendicular[0];
	transformationMatrix(2, 1) = secondPerpendicular[1];
	transformationMatrix(2, 2) = secondPerpendicular[2];
	transformationMatrix(0, 0) = d[0];
	transformationMatrix(0, 1) = d[1];
	transformationMatrix(0, 2) = d[2];
	auto t = Mult_External(transformationMatrix, PointToVector(o));
	transformationMatrix(0, 3) = -t[0];
	transformationMatrix(1, 3) = -t[1];
	transformationMatrix(2, 3) = -t[2];
	return transformationMatrix;
}

void CPointCloud::CreateEGGLocalCoordinateSystem(Matrix4x4& transformationMatrix, const CPointEx3D& o, const CPointEx3D& d)
{
	CPointEx3D firstPerpendicular = GetNormalizedPerpendicularVectorToVector(d);
	CPointEx3D secondPerpendicular = cross(d, firstPerpendicular);
	secondPerpendicular.normalize();
	identity4x4(transformationMatrix);
	transformationMatrix[1][0] = firstPerpendicular[0];
	transformationMatrix[1][1] = firstPerpendicular[1];
	transformationMatrix[1][2] = firstPerpendicular[2];
	transformationMatrix[2][0] = secondPerpendicular[0];
	transformationMatrix[2][1] = secondPerpendicular[1];
	transformationMatrix[2][2] = secondPerpendicular[2];
	transformationMatrix[0][0] = d[0];
	transformationMatrix[0][1] = d[1];
	transformationMatrix[0][2] = d[2];
	auto t = mult4x4(transformationMatrix, o);
	transformationMatrix[0][3] = -t[0];
	transformationMatrix[1][3] = -t[1];
	transformationMatrix[2][3] = -t[2];
}

// modified to support static Matrix4x4
bool CPointCloud::isEllipicGabrielNeighbor(int i, const std::vector<int>& NNs, double a) const
{
	const CPointEx3D& p = (*points)[NNs[0]];
	const CPointEx3D& qi = (*points)[NNs[i]];
	const CPointEx3D origin = 0.5 * (p + qi);
	const double d = (qi - p).norm() / 2.0;
	CPointEx3D localXaxis = qi - p;
	localXaxis.normalize();
	Matrix4x4 transformationMatrix;
	CreateEGGLocalCoordinateSystem(transformationMatrix, origin, localXaxis);
	for (int j = 1; j < i; j++)
	{
		CPointEx3D pnt = mult4x4(transformationMatrix, (*points)[NNs[j]]);
		double x = pnt.x, y = pnt.y, z = pnt.z;
		double ellipsoidValue = x * x + y * y / (a * a) + z * z / (a * a);
		if (ellipsoidValue < d * d) return false;
	}
	return true;
}


void CPointCloud::calculateEGG(std::vector<std::vector<int>>& EGG, int knn, double ratio) const
{
	EGG.resize(points->size());
	int pointsCount = (int)EGG.size();
	parallel_for(tbb::blocked_range<int>(0, pointsCount),
		[&](tbb::blocked_range<int> r)
		{
			for (int i = r.begin(); i < r.end(); i++)
			{
				auto NNs = GetKNearestNeighborIndex((*points)[i], knn);
				for (int j = 1; j < NNs.size(); j++)
				{
					if (isEllipicGabrielNeighbor(j, NNs, ratio))
						EGG[i].push_back(NNs[j]);
				}
			}
		});
}

std::vector<int> CPointCloud::GetKNearestNeighborIndex(const CPointEx3D& pnt, int K) const
{
	std::vector<int> indexes;
	const Point_3 point(pnt.x, pnt.y, pnt.z);
	K_neighbor_search search(*((Tree*) pointCloudSpatialTree), point, K);
	for (K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
	{
		indexes.push_back(boost::get<1>(it->first));
	}
	return indexes;
}

CDoubleVector CPointCloud::PointToVector(const CPointEx3D& pnt, bool homogenous)
{
	CDoubleVector vec(4);
	vec[0] = pnt[0];
	vec[1] = pnt[1];
	vec[2] = pnt[2];
	if (homogenous) vec[3] = 1.0;
	else vec[3] = 0.0;
	return vec;
}

CPointEx3D CPointCloud::VectorToPoint(const CVector<double>& vec)
{
	return CPointEx3D(vec[0], vec[1], vec[2]);
}

CPointEx3D CPointCloud::GetNormalizedPerpendicularVectorToVector(const CPointEx3D& inVector)
{
	double max = fabs(inVector[0]);
	int cordIndex = 0;

	if (max < fabs(inVector[1]))
	{
		cordIndex = 1;
		max = fabs(inVector[1]);
	}

	if (max < fabs(inVector[2]))
	{
		cordIndex = 2;
	}

	CPointEx3D outVector;
	outVector[0] = 1.0;
	outVector[1] = 1.0;
	outVector[2] = 1.0;

	switch (cordIndex)
	{
	case 0:
		outVector[0] = (-inVector[1] * outVector[1] - inVector[2] * outVector[2]) / inVector[0];
		break;
	case 1:
		outVector[1] = (-inVector[0] * outVector[0] - inVector[2] * outVector[2]) / inVector[1];
		break;
	case 2:
		outVector[2] = (-inVector[0] * outVector[0] - inVector[1] * outVector[1]) / inVector[2];
		break;
	}
	outVector.normalize();
	return outVector;
}

bool CPointCloud::GeodesicAstar(const CPointEx3D& start, const CPointEx3D& end, const std::vector<std::vector<int>>& NN,
	std::function<double(int i, int j)> metricDistance, std::function<double(int i, int j)> metricHeuristic, std::vector<CPointEx3D>& path)
{
	int index_start = GetKNearestNeighborIndex(start, 1)[0];
	int index_end = GetKNearestNeighborIndex(end, 1)[0];
	return GeodesicAstar(index_start, index_end, NN, metricDistance, metricHeuristic, path);
}

bool CPointCloud::GeodesicAstar(const int& index_start, const int& index_end, const std::vector<std::vector<int>>& NN,
	std::function<double(int i, int j)> metricDistance, std::function<double(int i, int j)> metricHeuristic, std::vector<CPointEx3D>& path)
{
	struct Element {
		int index;
		std::shared_ptr<Element> parent;
		double geodistance;
		double heuristicDistance;
		double fdistance;
	};
	auto compSet = [](const Element& x, const Element& y) { return x.index < y.index; };
	auto closed = std::set<Element, decltype(compSet)>(compSet);
	Element start{ index_start, nullptr, 0.0, 0.0, 0.0 };
	std::map<int, double> openIndexes;
	auto cmp = [](const Element left, const Element right) { return left.fdistance > right.fdistance; };
	std::priority_queue<Element, std::vector<Element>, decltype(cmp)> open(cmp);
	open.push(start);
	Element el;
	while (true) {
		el = open.top();
		open.pop();
		int v = el.index;
		if (v == index_end) break;
		closed.insert(el);
		auto neighbors = NN[v];
		std::shared_ptr<Element> parent = std::make_shared<Element>();
		*parent = el;
		for (auto va : neighbors)
		{
			double distG = el.geodistance + metricDistance(v, va);
			double distH = metricHeuristic(va, index_end);
			double Css = el.parent ? (((*points)[v] - (*points)[va]) - ((*points)[el.parent->index] - (*points)[v])).norm() : 0.0;
			double distF = distG + distH + Css;
			Element newEl = { va, parent, distG, distH, distF };
			if (!closed.count(newEl) && !openIndexes.count(va)) {
				open.push(newEl);
				openIndexes[va] = distG;
			}
			else if (!closed.count(newEl) && openIndexes[va] > distG) {
				open.push(newEl);
				openIndexes[va] = distG;
			}
		}
	}
	std::stack<int> reversePath;
	while (el.parent) {
		reversePath.push(el.index);
		el = *el.parent;
	}
	reversePath.push(el.index);
	while (!reversePath.empty()) {
		int idx = reversePath.top();
		reversePath.pop();
		path.push_back((*points)[idx]);
	}
	return true;
}

bool CPointCloud::GeodesicDijkstra(const int& index_start, const int& index_end, const std::vector<std::vector<int>>& NN,
	std::function<double(int i, int j)> metricDistance, std::vector<CPointEx3D>& path)
{
	struct GeoElement {
		double geodistance;
		int index;
	};
	struct GeoDistance
	{
		double d;
		int parentIdx;
	};
	std::vector<GeoDistance> g;

	int pointsCount = (int)points->size();
	g.resize(pointsCount, { 1e10, -1 });
	g[index_start] = { metricDistance(index_start, index_start), -1 };
	GeoElement gi{ metricDistance(index_start, index_start), index_start };
	auto cmp = [](const GeoElement left, const GeoElement right) { return left.geodistance > right.geodistance; };
	std::priority_queue<GeoElement, std::vector<GeoElement>, decltype(cmp)> VLIST(cmp);
	VLIST.push(gi);
	while (!VLIST.empty()) {
		GeoElement ge = VLIST.top();
		VLIST.pop();
		int v = ge.index;
		auto neighbors = NN[v];
		for (auto va : neighbors)
		{
			double dist = g[v].d + metricDistance(v, va);
			if (dist < g[va].d) {
				GeoElement gel = { dist, va };
				VLIST.push(gel);
				g[va].d = dist;
				g[va].parentIdx = v;
			}
		}
	}
	std::stack<int> reversePath;
	auto currentIndex = index_end;
	do {
		reversePath.push(currentIndex);
		currentIndex = g[currentIndex].parentIdx;
	} while (currentIndex != -1);
	while (!reversePath.empty()) {
		int idx = reversePath.top();
		reversePath.pop();
		path.push_back((*points)[idx]);
	}
	return true;
}

bool CPointCloud::GeodesicDijkstra(const CPointEx3D& start, const CPointEx3D& end, const std::vector<std::vector<int>>& NN,
	std::function<double(int i, int j)> metricDistance, std::vector<CPointEx3D>& path)
{
	int index_start = GetKNearestNeighborIndex(start, 1)[0];
	int index_end = GetKNearestNeighborIndex(end, 1)[0];
	return GeodesicDijkstra(index_start, index_end, NN, metricDistance, path);
}


void CPointCloud::calculateGeodesicPositionAndNormals(const std::vector<CPointEx3D>& normalizedPoints,
	std::vector<CPointEx3D>& pathPoints, std::vector<CPointEx3D>& pathNormals, double min, double max, double alpha) const
{
    auto isEllipicGabrielNeighbor = [&normalizedPoints](const CPointEx3D& originPoint, int i, const std::vector<int>& NNs, double a)
    {
        const CPointEx3D& p = originPoint;
        const CPointEx3D& qi = normalizedPoints[NNs[i]];
        const CPointEx3D origin = 0.5 * (p + qi);
        const double d = (qi - p).norm() / 2.0;
        CPointEx3D localXaxis = qi - p;
        localXaxis.normalize();
        //		CDoubleMatrix transformationMatrix = CreateEGGLocalCoordinateSystem(origin, localXaxis);
        Matrix4x4 transformationMatrix;
        CreateEGGLocalCoordinateSystem(transformationMatrix, origin, localXaxis);
        for (int j = 0; j < i; j++)
        {
            //			CPointEx3D pnt = VectorToPoint(Mult(transformationMatrix, PointToVector(normalizedPoints[NNs[j]])));
            //			CPointEx3D pnt = VectorToPoint(Mult_External(transformationMatrix, PointToVector(normalizedPoints[NNs[j]])));			double x = pnt.x, y = pnt.y, z = pnt.z;
            CPointEx3D pnt = mult4x4(transformationMatrix, normalizedPoints[NNs[j]]);
            double x = pnt.x, y = pnt.y, z = pnt.z;
            double ellipsoidValue = x * x + y * y / (a * a) + z * z / (a * a);
            if (ellipsoidValue < d * d) return false;
        }
        return true;
    };

    int curvePointCount = pathPoints.size();
    //pathNormals.push_back(normals[]);
    tbb::parallel_for(tbb::blocked_range<int>(1, curvePointCount - 1),
                      [&](tbb::blocked_range<int> r)
                      {
                          for (int i = r.begin(); i < r.end(); i++)
                          {
                              CPointEx3D originPoint = pathPoints[i];
                              std::vector<int> curveGabrielNeighbors;
                              auto pointNNs = GetKNearestNeighborIndex(denormalizeCoordinate(originPoint, min, max), 30);
                              if (originPoint != normalizedPoints[pointNNs[0]]) {
                                  for (int j = 0; j < pointNNs.size(); j++) {
                                      //double alpha = 1e-12; when noise is high
                                      //double alpha = 0.5; when noise is normal
                                      if (isEllipicGabrielNeighbor(originPoint, j, pointNNs, alpha)) {
                                          curveGabrielNeighbors.push_back(pointNNs[j]);
                                      }
                                  }
                                  double c0 = 0.;
                                  CPointEx3D normal(0., 0., 0.);
                                  CPointEx3D c(0., 0., 0.);
                                  for (auto idx : curveGabrielNeighbors) {
                                      auto direction = normalizedPoints[idx] - originPoint;
                                      double dist = direction.norm();
                                      double w = exp(-dist * dist);
                                      normal += w * (*normals)[idx];
                                      c.x += w * normalizedPoints[idx].x;
                                      c.y += w * normalizedPoints[idx].y;
                                      c.z += w * normalizedPoints[idx].z;
                                      c0 += w;
                                  }
                                  normal /= c0;
                                  normal.normalize();
                                  double lambda = (c.x * normal.x + c.y * normal.y + c.z * normal.z) / c0;
                                  double t = lambda - dot(pathPoints[i], normal);
                                  pathPoints[i] += t * normal;
                                  pathNormals[i] = normal;
                              }
                              else
                              {
                                  pathNormals[i] = (*normals)[pointNNs[0]];
                              }
                          }
                      });
}

CPointEx3D CPointCloud::normalizeCoordinate(const CPointEx3D& p, double minNormalizedExtent, double maxNormalizedExtent) {
	return CPointEx3D(
		2.0 * (p.x - minNormalizedExtent) / (maxNormalizedExtent - minNormalizedExtent) - 1.0,
		2.0 * (p.y - minNormalizedExtent) / (maxNormalizedExtent - minNormalizedExtent) - 1.0,
		2.0 * (p.z - minNormalizedExtent) / (maxNormalizedExtent - minNormalizedExtent) - 1.0);
}

CPointEx3D CPointCloud::denormalizeCoordinate(const CPointEx3D& p, double minNormalizedExtent, double maxNormalizedExtent) {
	return CPointEx3D(
		0.5 * (p.x + 1.0) * (maxNormalizedExtent - minNormalizedExtent) + minNormalizedExtent,
		0.5 * (p.y + 1.0) * (maxNormalizedExtent - minNormalizedExtent) + minNormalizedExtent,
		0.5 * (p.z + 1.0) * (maxNormalizedExtent - minNormalizedExtent) + minNormalizedExtent);
}

double CPointCloud::sqr_length(const std::vector<CPointEx3D>& path) {
	double d2 = 0;
	for (int i = 0; i < path.size() - 1; i++)
		d2 += (path[i + 1] - path[i]).sqr_norm();
	return d2;
}

bool CPointCloud::GeodesicFlow(std::vector<CPointEx3D>& path, double lambda, double mu, double alpha)
{
    std::vector<CPointEx3D> pathNormals(path.size());
    getNormalWithGabrielNeighborWeighting(path, pathNormals);
    double maxx = -1e10;
    double minx = 1e10;
    double maxy = -1e10;
    double miny = 1e10;
    double maxz = -1e10;
    double minz = 1e10;
    for (int i = 0; i < points->size(); i++)
    {
        if ((*points)[i].x > maxx) maxx = (*points)[i].x;
        if ((*points)[i].x < minx) minx = (*points)[i].x;
        if ((*points)[i].y > maxy) maxy = (*points)[i].y;
        if ((*points)[i].y < miny) miny = (*points)[i].y;
        if ((*points)[i].z > maxz) maxz = (*points)[i].z;
        if ((*points)[i].z < minz) minz = (*points)[i].z;
    }

    double max = maxx;
    if (max < maxy) max = maxy;
    if (max < maxz) max = maxz;

    double min = minx;
    if (min > miny) min = miny;
    if (min > minz) min = minz;
    auto EucledianDistance = [this](int i, int j)
    {
        return ((*points)[i] - (*points)[j]).norm();
    };

    std::for_each(path.begin(), path.end(),
                  [min, max](CPointEx3D& p) {p = normalizeCoordinate(p, min, max); });

    std::vector<CPointEx3D> normalizedPoints = *points;
    std::for_each(normalizedPoints.begin(), normalizedPoints.end(),
                  [min, max](CPointEx3D& p) {p = normalizeCoordinate(p, min, max); });

    std::vector<CPointEx3D> correctedCurve = path;
    int m = path.size();
    double length1 = sqr_length(path);
    double length2 = -1e12;
    const double TOL = 1e-6 * 1e-6; // convergence criterion
    int k = 0;
    auto start = std::chrono::steady_clock::now();
    _mem_vecs mem_vectors;
    mem_vectors.Init(m);
    FILE* file = fopen("Flows.txt", "w");
    while (k++ < 2453 && fabs(length2 - length1)>TOL) {
        double flow = getGeodesicFlowMKL(path, pathNormals, &mem_vectors);
        CMathUtilities::ConjugateGradientMKL(mem_vectors.D2F, mem_vectors.DF, mem_vectors.DF, 3 * m);
        for (int i = 1; i < m - 1; i++) {
            CPointEx3D pk(mem_vectors.DF[3 * i + 0], mem_vectors.DF[3 * i + 1], mem_vectors.DF[3 * i + 2]);
            auto p = path[i] + pk;
            correctedCurve[i] = p;
        }
        calculateGeodesicPositionAndNormals(normalizedPoints, correctedCurve, pathNormals, min, max, alpha);
        path = correctedCurve;
        std::vector<CPointEx3D> pathSmooth = path;

        for (int l = 0; l < 4; l++)
        {
            for (int i = 1; i < path.size() - 1; i++)
            {
                auto Laplacian = 0.5 * (path[i + 1] + path[i - 1]) - path[i];
                pathSmooth[i] = path[i] + lambda * Laplacian;
            }
            swap(path, pathSmooth);
            for (int i = 1; i < path.size() - 1; i++)
            {
                auto Laplacian = 0.5 * (path[i + 1] + path[i - 1]) - path[i];
                pathSmooth[i] = path[i] + mu * Laplacian;
            }
            swap(path, pathSmooth);
        }

        if (length2<0)
            length2 = sqr_length(path);
        else
        {
            length1 = length2;
            length2 = sqr_length(path);
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in milliseconds: "
              << std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count()
              << " ms" << std::endl;
    std::cout << "iterations taken:" << k << std::endl;
    std::for_each(path.begin(), path.end(),
                  [min, max](CPointEx3D& p) {p = denormalizeCoordinate(p, min, max); });
    mem_vectors.Destroy();
    return true;
}

void CPointCloud::calculateNeighbors(std::vector<std::vector<int>>& NNs, int neighborCount) const
{
	int pointsCount = points->size();
	NNs.resize(pointsCount);
	tbb::parallel_for(tbb::blocked_range<int>(0, pointsCount),
		[&](tbb::blocked_range<int> r)
		{
			for (int i = r.begin(); i < r.end(); i++)
			{
				NNs[i].resize(neighborCount);
				auto pointNNs = GetKNearestNeighborIndex((*points)[i], neighborCount + 1);
				for (int j = 1; j < pointNNs.size(); j++) {
					NNs[i][j - 1] = pointNNs[j];
				}
			}
		});
}

void CPointCloud::_mem_vecs::Init(int m){
	K = (double*)mkl_malloc(9 * m * m * sizeof(double), 64);
	nStar = (double*)mkl_malloc(9 * m * m * sizeof(double), 64);
	kG = (double*)mkl_malloc(9 * m * m * sizeof(double), 64);
	nT = (double*)mkl_malloc(3 * m * m * sizeof(double), 64);
	p = (double*)mkl_malloc(3 * m * sizeof(double), 64);
	flow = (double*)mkl_malloc(3 * m * sizeof(double), 64);
	D2F = (double*)mkl_malloc(9 * m * m * sizeof(double), 64);
	DF = (double*)mkl_malloc(3 * m * sizeof(double), 64);
}

void CPointCloud::_mem_vecs::Destroy() {
	if(K) mkl_free(K); K = 0;
	if(nStar) mkl_free(nStar); nStar = 0;
	if(kG) mkl_free(kG); kG = 0;
	if(nT) mkl_free(nT); nT = 0;
	if(p) mkl_free(p); p = 0;
	if (flow) mkl_free(flow); flow = 0;
	if (D2F) mkl_free(D2F); D2F = 0;
	if (DF) mkl_free(DF); DF = 0;
}


double CPointCloud::getGeodesicFlowMKL(const std::vector<CPointEx3D>& curvePoints, const std::vector<CPointEx3D>& curveNormals, _mem_vecs* mv) {
    int m = curvePoints.size();
    bool memory_already_allocated = false;
    if (!mv) {
        mv = new _mem_vecs;
        mv->Init(m);
    }
    else
    {
        memory_already_allocated = true;
    }

    memset(mv->K, 0, 9 * m * m * sizeof(double));

    int j = 0;
    for (int i = 3; i < 3 * m - 3; i++) {
        mv->K[3 * m * i + j] = -1;
        mv->K[3 * m * i + j + 3] = 2;
        mv->K[3 * m * i + j + 6] = -1;
        j++;
    }

    cblas_dcopy(9 * m * m, mv->K, 1, mv->D2F, 1);

    memset(mv->nT, 0, 3 * m * m * sizeof(double));
    j = 0;
    for (int i = 0; i < m; i++) {
        mv->nT[3 * m * i + j] = curveNormals[i].x;
        mv->nT[3 * m * i + j + 1] = curveNormals[i].y;
        mv->nT[3 * m * i + j + 2] = curveNormals[i].z;
        j += 3;
    }

    int mm = 3 * m, kk = m, nn = 3 * m;
    double alpha = 1.0, beta = 0.0;

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, mm, nn, kk, alpha, mv->nT, mm, mv->nT, mm, beta, mv->nStar, nn);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mm, mm, mm, -1.0, mv->nStar, mm, mv->K, mm, 1.0, mv->D2F, mm);
    cblas_dcopy(9 * m * m, mv->D2F, 1, mv->kG, 1);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, mm, mm, mm, 1.0, mv->D2F, mm, mv->D2F, mm, 0.0, mv->K, mm);
    cblas_dcopy(9 * m * m, mv->K, 1, mv->D2F, 1);

    for (int i = 0; i < m; i++) {
        mv->p[3 * i + 0] = curvePoints[i].x;
        mv->p[3 * i + 1] = curvePoints[i].y;
        mv->p[3 * i + 2] = curvePoints[i].z;
    }

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3 * m, 3 * m, -1.0, mv->D2F, 3 * m, mv->p, 1, 0.0, mv->DF, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3 * m, 3 * m, 1.0, mv->kG, 3 * m, mv->p, 1, 0.0, mv->flow, 1);

    double flow_norm = cblas_dnrm2(3 * m, mv->flow, 1);
    if (!memory_already_allocated)
    {
        mv->Destroy();
        mv = nullptr;
    }
    return 1.0/sqrt(m) * flow_norm;
}

CPointCloud::~CPointCloud()
{
	clear();
}

void CPointCloud::clear() {
	if (points != NULL) {
		points->clear();
		delete points;
		points = NULL;
	}
	if (normals != NULL) {
		normals->clear();
		delete normals;
		normals = NULL;
	}
	if (pointCloudSpatialTree) {
		((Tree*)pointCloudSpatialTree)->clear();
		delete pointCloudSpatialTree;
		pointCloudSpatialTree = NULL;
	}
}


std::vector<CPointEx3D>& CPointCloud::getPnts()
{
	return *points;
}
std::vector<CPointEx3D>& CPointCloud::getNorms()
{
	return *normals;
}
std::vector<CPointEx3D>* CPointCloud::getPntsPtr()
{
	return points;
}
std::vector<CPointEx3D>* CPointCloud::getNormsPtr()
{
	return normals;
}


void CPointCloud::updateSpatialTree()
{
	std::vector<Point_3> CGALpoints(points->size());
	std::vector<int> indices(points->size());
#pragma omp parallel for
	for (int i = 0; i < points->size(); i++)
	{
		CGALpoints[i] = Point_3((*points)[i][0], (*points)[i][1], (*points)[i][2]);
		indices[i] = i;
	}
	delete pointCloudSpatialTree;
	pointCloudSpatialTree = new Tree(boost::make_zip_iterator(boost::make_tuple(CGALpoints.begin(), indices.begin())),
		boost::make_zip_iterator(boost::make_tuple(CGALpoints.end(), indices.end())));
}

void CPointCloud::calculateNormals(bool bOriented) {
    if (!points || points->size() == 0)
        return; // nothing to do

    // Reads a .xyz point set file in points[].
    std::list<CGALPointVectorPair> points_n_normals;
    for (int i = 0; i < points->size(); i++)
    {
        Point_3 pnt((*points)[i].x, (*points)[i].y, (*points)[i].z);
        Vector_3 normal(0, 0, 0);
        points_n_normals.push_back(std::make_pair(pnt, normal));
    }

    const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
    // Estimates normals direction.
    // Note: pca_estimate_normals() requires an iterator over points
    // as well as property maps to access each point's position and normal.
    CGAL::pca_estimate_normals<Concurrency_tag>
            (points_n_normals, nb_neighbors,
             CGAL::parameters::point_map(CGAL::First_of_pair_property_map<CGALPointVectorPair>()).
                     normal_map(CGAL::Second_of_pair_property_map<CGALPointVectorPair>()));

/*
		// First compute a spacing using the K parameter
		double spacing
			= CGAL::compute_average_spacing<Concurrency_tag>
			(points_n_normals, nb_neighbors,
				CGAL::parameters::point_map(CGAL::First_of_pair_property_map<CGALPointVectorPair>()));
		// Then, estimate normals with a fixed radius
		CGAL::pca_estimate_normals<Concurrency_tag>
			(points_n_normals,
				0, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
				CGAL::parameters::point_map(CGAL::First_of_pair_property_map<CGALPointVectorPair>())
				.normal_map(CGAL::Second_of_pair_property_map<CGALPointVectorPair>())
				.neighbor_radius(2. * spacing)); // use 2*spacing as neighborhood radius
*/

    if (bOriented)
    {
        // Orients normals.
        // Note: mst_orient_normals() requires a range of points
        // as well as property maps to access each point's position and normal.
        std::list<CGALPointVectorPair>::iterator unoriented_points_begin =
                CGAL::mst_orient_normals(points_n_normals, nb_neighbors,
                                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<CGALPointVectorPair>())
                                                 .normal_map(CGAL::Second_of_pair_property_map<CGALPointVectorPair>()));
        // Optional: delete points with an unoriented normal
        // if you plan to call a reconstruction algorithm that expects oriented normals.
        points_n_normals.erase(unoriented_points_begin, points_n_normals.end());

        // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
        std::list<CGALPointVectorPair>(points_n_normals).swap(points_n_normals);
    }

    // copy points back to points and normal struct

    if (points->size() != points_n_normals.size()) // points have been changed; reinitialize data struct
    {
        points->resize(points_n_normals.size());
    }

    // clear existing data if any
    if (normals != NULL) {
        normals->clear();
        delete normals;
        normals = NULL;
    }
    if (pointCloudSpatialTree) {
        ((Tree*)pointCloudSpatialTree)->clear();
        delete pointCloudSpatialTree;
        pointCloudSpatialTree = NULL;
    }

    normals = new std::vector<CPointEx3D>(points->size());


    std::list<CGALPointVectorPair>::iterator it = points_n_normals.begin();
    for (int i = 0; i < points->size(); i++)
    {
        (*points)[i].x = it->first.x();
        (*points)[i].y = it->first.y();
        (*points)[i].z = it->first.z();
        (*normals)[i].x = it->second.x();
        (*normals)[i].y = it->second.y();
        (*normals)[i].z = it->second.z();
        it++;
    }

    std::vector<Point_3> pointsCGAL;
    std::vector<int> indicesCGAL;
    for (int i = 0; i < points->size(); i++)
    {
        pointsCGAL.push_back(Point_3((*points)[i][0], (*points)[i][1], (*points)[i][2]));
        indicesCGAL.push_back(i);

    }
    pointCloudSpatialTree = new Tree(boost::make_zip_iterator(boost::make_tuple(pointsCGAL.begin(), indicesCGAL.begin())),
                                     boost::make_zip_iterator(boost::make_tuple(pointsCGAL.end(), indicesCGAL.end())));
}


void CPointCloud::getNormalWithGabrielNeighborWeighting(const std::vector<CPointEx3D>& pointPath,
                                                        std::vector<CPointEx3D>& normalPath) const
{
    std::vector<CPointEx3D> normalizedPoints = *points;
    auto isEllipicGabrielNeighbor = [&normalizedPoints](const CPointEx3D& originPoint, int i, const std::vector<int>& NNs, double a)
    {
        const CPointEx3D& p = originPoint;
        const CPointEx3D& qi = normalizedPoints[NNs[i]];
        const CPointEx3D origin = 0.5 * (p + qi);
        const double d = (qi - p).norm() / 2.0;
        CPointEx3D localXaxis = qi - p;
        localXaxis.normalize();
        //		CDoubleMatrix transformationMatrix = CreateEGGLocalCoordinateSystem(origin, localXaxis);
        Matrix4x4 transformationMatrix;
        CreateEGGLocalCoordinateSystem(transformationMatrix, origin, localXaxis);
        for (int j = 0; j < i; j++)
        {
            //			CPointEx3D pnt = VectorToPoint(Mult(transformationMatrix, PointToVector(normalizedPoints[NNs[j]])));
            //			CPointEx3D pnt = VectorToPoint(Mult_External(transformationMatrix, PointToVector(normalizedPoints[NNs[j]])));			double x = pnt.x, y = pnt.y, z = pnt.z;
            CPointEx3D pnt = mult4x4(transformationMatrix, normalizedPoints[NNs[j]]);
            double x = pnt.x, y = pnt.y, z = pnt.z;
            double ellipsoidValue = x * x + y * y / (a * a) + z * z / (a * a);
            if (ellipsoidValue < d * d) return false;
        }
        return true;
    };
    double maxx = -1e10;
    double minx = 1e10;
    double maxy = -1e10;
    double miny = 1e10;
    double maxz = -1e10;
    double minz = 1e10;
    for (int i = 0; i < points->size(); i++)
    {
        if ((*points)[i].x > maxx) maxx = (*points)[i].x;
        if ((*points)[i].x < minx) minx = (*points)[i].x;
        if ((*points)[i].y > maxy) maxy = (*points)[i].y;
        if ((*points)[i].y < miny) miny = (*points)[i].y;
        if ((*points)[i].z > maxz) maxz = (*points)[i].z;
        if ((*points)[i].z < minz) minz = (*points)[i].z;
    }

    double max = maxx;
    if (max < maxy) max = maxy;
    if (max < maxz) max = maxz;

    double min = minx;
    if (min > miny) min = miny;
    if (min > minz) min = minz;

    std::vector<CPointEx3D> pointsPathNorm = pointPath;
    std::for_each(std::execution::par, pointsPathNorm.begin(), pointsPathNorm.end(),
                  [min, max](CPointEx3D& p) {p = normalizeCoordinate(p, min, max); });


    std::for_each(std::execution::par, normalizedPoints.begin(), normalizedPoints.end(),
                  [min, max](CPointEx3D& p) {p = normalizeCoordinate(p, min, max); });
    normalPath.resize(pointsPathNorm.size());

    int curvePointCount = pointsPathNorm.size();
    tbb::parallel_for(tbb::blocked_range<int>(0, curvePointCount),
                      [&](tbb::blocked_range<int> r)
                      {
                          for (int i = r.begin(); i < r.end(); i++)
                          {
                              CPointEx3D originPoint = pointsPathNorm[i];

                              std::vector<int> curveGabrielNeighbors;
                              auto pointNNs = GetKNearestNeighborIndex(denormalizeCoordinate(originPoint, min, max), 30);
                              if (originPoint != normalizedPoints[pointNNs[0]]) {
                                  for (int j = 0; j < pointNNs.size(); j++) {
                                      if (isEllipicGabrielNeighbor(originPoint, j, pointNNs, 0.5)) {
                                          curveGabrielNeighbors.push_back(pointNNs[j]);
                                      }
                                  }
                                  double c0 = 0.;
                                  CPointEx3D normal(0.0, 0.0, 0.0);
                                  CPointEx3D c(0., 0., 0.);
                                  for (auto idx : curveGabrielNeighbors) {
                                      auto direction = normalizedPoints[idx] - originPoint;
                                      double dist = direction.norm();
                                      double w = exp(-dist * dist);
                                      normal += w * (*normals)[idx];
                                      c0 += w;
                                  }
                                  normal /= c0;
                                  normal.normalize();
                                  normalPath[i] = normal;
                              }
                              else
                              {
                                  normalPath[i] = (*normals)[pointNNs[0]];
                              }
                          }
                      });
}
