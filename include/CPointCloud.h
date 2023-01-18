#pragma once

#include <functional>
#include <vector>
#ifdef BUILD_DLL_EXPORTS
#include "./NurbLibSt/MCLSEXST.h"
//#else
//#include "MCLSEXST.h"
#endif
#ifndef DECLSPEC
# if defined (BUILD_DLL_EXPORTS)
#     define DECLSPEC __declspec( dllexport ) 
# else
#     define DECLSPEC __declspec( dllimport ) 
# endif
#endif

namespace LIBGEODESIC {

	typedef double Matrix4x4[4][4];

	class DECLSPEC CPointCloud
	{
	private:
		// this struct is used by GeodesicFlow
		struct _mem_vecs {
			double* K;
			double* nStar;
			double* kG;
			double* nT;
			double* p;
			double* flow;
			double* DF;
			double* D2F;
			void Init(int m);
			void Destroy();
		} mem_vecs;  

	public:
		CPointCloud();
		/**
		 * Initialize the class with point data and normals
		 *
		 * @param inputPoints Vector containing the point cloud
		 * @param inputNormals Vector containing the normals
		 */
		CPointCloud(const std::vector<CPointEx3D>& inputPoints, const std::vector<CPointEx3D>& inputNormals);

		/**
		 * init Initialize the class with point data and normals
		 *
		 * @param inputPoints Vector containing the point cloud
		 * @param inputNormals Vector containing the normals
		 */
		void init(const std::vector<CPointEx3D>& inputPoints, const std::vector<CPointEx3D>& inputNormals);
		
		/**
		 * init_pnts Initialize the class with point data only
		 *
		 * @param inputPoints Vector containing the point cloud
		 * @param bCalcNormals Will compute the oriented normal vectors is true
		 */
		void init_pnts(const std::vector<CPointEx3D>& inputPoints, bool bCalcNormals = true);

		bool isInit() {
			return (points != NULL && normals != NULL && points->size() == normals->size() && pointCloudSpatialTree != NULL);
		}

		void clear();

		~CPointCloud();
		/**
		 * Ellipitic Gabriel Neighbors calculation
		 *
		 * @param EGG vector containing the gabriel neighbors of each point
		 * @param knn The kNN to search the neighbors from
		 * @param ratio The alpha value of the elliptic gabriel graph
		 */
		void calculateEGG(std::vector<std::vector<int>>& EGG, int knn = 20, double ratio = 0.75) const;
		/**
		 * Nearest Neighbors calculation
		 *
		 * @param NNs vector containing the knn neighbors of each point
		 * @param neighborCount The k number of kNN
		 * @param ratio The alpha value of the elliptic gabriel graph
		 */
		void calculateNeighbors(std::vector<std::vector<int>>& NNs, int neighborCount) const;

        /**
		 * calculateNormals Calculate normal vectors to the point data
		 *
		 * @param bOriented If bOriented == true then the final point cloud can be different than the original one
		 *					because points with unoriented normals will be removed
		 */
        void calculateNormals(bool bOriented = true);

		// Access to internal data; use with causious
		std::vector<CPointEx3D>& getPnts();
		std::vector<CPointEx3D>* getPntsPtr();
		std::vector<CPointEx3D>& getNorms();
		std::vector<CPointEx3D>* getNormsPtr();

		/**
		 * Geodesic Astar algorithm from index to index
		 * The user needs to provide a metric to measure the distance from
		 * point with index i to point with index j
		 * The user needs to provide a heuristic metric to measure the distance from
		 * point with index i to termination point with index j
		 *
		 * @param index_start The index of the starting point
		 * @param index_end The index of the termination point
		 * @param NN The connectivity of the points (for example kNN)
		 * @param metricDistance The distance metric between two neighbor points with indexes i and j
		 * @param metricHeuristic The heuristic metric distance between the point of index i and finish point of index j
		 * @param path The output path
		 * @return The path return success
		 */
		bool GeodesicAstar(const int& index_start, const int& index_end, const std::vector<std::vector<int>>& NN,
			std::function<double(int i, int j)> metricDistance, std::function<double(int i, int j)> metricHeuristic, std::vector<CPointEx3D>& path);
		/**
		 * Geodesic Astar algorithm from point to point
		 * The user needs to provide a metric to measure the distance from
		 * point with index i to point with index j
		 * The user needs to provide a heuristic metric to measure the distance from
		 * point with index i to termination point with index j
		 *
		 * @param start The starting point
		 * @param end The termination point
		 * @param NN The connectivity of the points (for example kNN)
		 * @param metricDistance The distance metric between two neighbor points with indexes i and j
		 * @param metricHeuristic The heuristic metric distance between the point of index i and finish point of index j
		 * @param path The output path
		 * @return The path return success
		 */
		bool GeodesicAstar(const CPointEx3D& start, const CPointEx3D& end, const std::vector<std::vector<int>>& NN,
			std::function<double(int i, int j)> metricDistance, std::function<double(int i, int j)> metricHeuristic, std::vector<CPointEx3D>& path);
		/**
		 * Geodesic Dijkstra algorithm from index to index
		 * The user needs to provide a metric to measure the distance from
		 * point with index i to point with index j
		 * point with index i to termination point with index j
		 *
		 * @param index_start The index of the starting point
		 * @param index_end The index of the termination point
		 * @param NN The connectivity of the points (for example kNN)
		 * @param metricDistance The distance metric between two neighbor points with indexes i and j
		 * @param path The output path
		 * @return The path return success
		 */
		bool GeodesicDijkstra(const int& index_start, const int& index_end, const std::vector<std::vector<int>>& NN,
			std::function<double(int i, int j)> metricDistance, std::vector<CPointEx3D>& path);
		/**
		 * Geodesic Dijkstra algorithm from point to point
		 * The user needs to provide a metric to measure the distance from
		 * point with index i to point with index j
		 * point with index i to termination point with index j
		 *
		 * @param start The starting point
		 * @param end The termination point
		 * @param NN The connectivity of the points (for example kNN)
		 * @param metricDistance The distance metric between two neighbor points with indexes i and j
		 * @param path The output path
		 * @return The path return success
		 */
		bool GeodesicDijkstra(const CPointEx3D& start, const CPointEx3D& end, const std::vector<std::vector<int>>& NN,
			std::function<double(int i, int j)> metricDistance, std::vector<CPointEx3D>& path);
		/**
		 * Geodesic path correction with geodesic curvature flow
		 *
		 * @param path The path initial input and the corrected path as output
		 * return patch correction success
		 */
		bool GeodesicFlow(std::vector<CPointEx3D>& path, double lambda = 0.05, double mu = -0.06, double alpha = 0.5);
        double getGeodesicFlowMKL(const std::vector<CPointEx3D>& curvePoints, const std::vector<CPointEx3D>& curveNormals,
                                  _mem_vecs* memory_vecs = nullptr);
        void getNormalWithGabrielNeighborWeighting(const std::vector<CPointEx3D>& pointPath,
                                                   std::vector<CPointEx3D>& normalPath) const;
	private:
		void updateSpatialTree();
		static CDoubleMatrix CreateEGGLocalCoordinateSystem(const CPointEx3D& o, const CPointEx3D& d);
		static void CreateEGGLocalCoordinateSystem(Matrix4x4& transformationMatrix, const CPointEx3D& o, const CPointEx3D& d);

		bool isEllipicGabrielNeighbor(int i, const std::vector<int>& NNs, double a) const;
		static CPointEx3D GetNormalizedPerpendicularVectorToVector(const CPointEx3D& inVector);
		static CDoubleVector PointToVector(const CPointEx3D& pnt, bool homogenous = true);
		static CPointEx3D VectorToPoint(const CVector<double>& pnt);
		void calculateGeodesicPositionAndNormals(const std::vector<CPointEx3D>& normalizedPoints, std::vector<CPointEx3D>& pathPoints,
			std::vector<CPointEx3D>& pathNormals, double min, double max, double alpha = 0.5) const;
		std::vector<int> GetKNearestNeighborIndex(const CPointEx3D& pnt, int K) const;
		static CPointEx3D normalizeCoordinate(const CPointEx3D& p, double minNormalizedExtent, double maxNormalizedExtent);
		static CPointEx3D denormalizeCoordinate(const CPointEx3D& p, double minNormalizedExtent, double maxNormalizedExtent);
		static double sqr_length(const std::vector<CPointEx3D>& path);
		std::vector<CPointEx3D>* points = nullptr;
		std::vector<CPointEx3D>* normals = nullptr;
		void* pointCloudSpatialTree = nullptr; // Tree*
	};

	inline void copyData(const CPointEx3DVector& inputData, std::vector<CPointEx3D>& outData) {
		outData.resize(inputData.GetSize());
#pragma omp parallel for
		for (int i = 0; i < inputData.GetSize(); i++)
		{
			outData[i] = inputData[i];
		}
	}


	inline void copyData(const std::vector<CPointEx3D>& inputData, CPointEx3DVector& outData) {
		outData.ReSize(inputData.size());
#pragma omp parallel for
		for (int i = 0; i < inputData.size(); i++)
		{
			outData[i] = inputData[i];
		}
	}

	inline CVector<double> Mult_External(const CMatrix<double>& M, const CVector<double>& v)
	{
		int m = M.GetRows();
		int n = M.GetCols();
		CVector<double> res(m);
		for (int i = 0; i < m; i++)
		{
			double data = 0.0;
			for (int j = 0; j < n; j++)
				data += (M(i, j) * v[j]);
			res[i] = data;
		}
		return res;
	}

	inline CPointEx3D mult4x4(const Matrix4x4& M, const CPointEx3D& p)
	{
		return CPointEx3D(
			p.x * M[0][0] + p.y * M[0][1] + p.z * M[0][2] + /*p.w **/ M[0][3],
			p.x * M[1][0] + p.y * M[1][1] + p.z * M[1][2] + /*p.w **/ M[1][3],
			p.x * M[2][0] + p.y * M[2][1] + p.z * M[2][2] + /*p.w **/ M[2][3]);
	}

	inline void identity4x4(Matrix4x4& M) {
		M[0][0] = 1.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = 0.0;
		M[1][0] = 0.0; M[1][1] = 1.0; M[1][2] = 0.0; M[1][3] = 0.0;
		M[2][0] = 0.0; M[2][1] = 0.0; M[2][2] = 1.0; M[2][3] = 0.0;
		M[3][0] = 0.0; M[3][1] = 0.0; M[3][2] = 0.0; M[3][3] = 1.0;
	}


} // NAMESPACE

