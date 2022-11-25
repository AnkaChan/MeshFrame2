#ifndef _TMESHLIB_SURFACE_MESH_H_
#define _TMESHLIB_SURFACE_MESH_H_

//#include "../TIterators.h"

namespace MF
{
	namespace TetMesh
	{
		// VertexType mush be child class of SurfaceVertexStatic (defined in SurfaceVertex.h)
		template <typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
		class CSurfaceMeshStatic : public TriMesh::CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType> {
		public:
			template <typename CTetMeshStatic>
			void initialize(CTetMeshStatic* pTetMesh) {
				typedef MF::TetMesh::TIterators<CTetMeshStatic> TItParent;

				std::vector<typename CTetMeshStatic::VPtr> _surfaceVerts;
				std::vector<typename CTetMeshStatic::HFPtr> _surfaceFaces;

				pTetMesh->tetMeshSurfaceMesh(_surfaceVerts, _surfaceFaces);

				for (int iV = 0; iV < _surfaceVerts.size(); iV++)
				{
					VertexType* currentVertex = createVertexWithId(_surfaceVerts[iV]->id());
					currentVertex->setTetMeshVertPtr(_surfaceVerts[iV]);
				}

				for (int iF = 0; iF < _surfaceFaces.size(); iF++)
				{
					VertexType* currentVs[3];
					/*Assume each face is a triangle*/
					//int fV = 0;
					//for (typename CTetMeshStatic::VPtr pFV : TItParent::HF_VIterator(_surfaceFaces[iF])) {
					//	int vId = pFV->id();
					//	currentVs[fV] = mVMap[vId];
					//	fV++;
					//}

	 				typename CTetMeshStatic::HEPtr pHE = (CTetMeshStatic::HEPtr)_surfaceFaces[iF]->half_edge();
					for (int fV = 0; fV < 3; ++fV) {
						int vId = pHE->target()->id();
						currentVs[fV] = mVMap[vId];

						pHE = (CTetMeshStatic::HEPtr)pHE->next();
					}

					FaceType* currentFace = createFace(currentVs, iF);
					mFIdMap.insert(FIdMapPair(iF, currentFace));
				}

				/*Label boundary edges*/
				for (int i = 0; i < mEContainer.getCurrentIndex(); ++i)
				{
					EdgeType* currentE = mEContainer.getPointer(i);
					HalfEdgeType* currentHE0 = edgeHalfedge(currentE, 0);
					HalfEdgeType* currentHE1 = edgeHalfedge(currentE, 1);
					if (currentHE1 == NULL)
					{
						currentHE0->source()->boundary() = true;
						currentHE0->target()->boundary() = true;
					}
				}
				/*Remove isolated vertex*/
				int  numIsolatedVertices = 0;
				
				std::vector<VertexType*> isolatedVertexs;
				//#pragma omp parallel for
				for (int i = 0; i < mVContainer.getCurrentIndex(); ++i)
				{
					VertexType* currentV = mVContainer.getPointer(i);
					if (mVContainer.hasBeenDeleted(currentV->index()) == false)
					{
						if (currentV->halfedge() != NULL) continue;
						isolatedVertexs.push_back(currentV);
					}
				}
				
				assert(isolatedVertexs.size()== 0);

				/*
				*	Arrange the boundary half_edge of boundary vertices, to make its halfedge
				*	to be the most ccw in half_edge
				*/
				for (int i = 0; i < mVContainer.getCurrentIndex(); ++i)
				{
					VertexType* currentV = mVContainer.getPointer(i);
					if (mVContainer.hasBeenDeleted(currentV->index()) == false)
					{
						if (!currentV->boundary()) continue;
						HalfEdgeType* currentHE = vertexMostCcwInHalfEdge(currentV);
						currentV->halfedge() = currentHE;
					}
				}
				mVMap.clear();
			}

		private:
			void readVFList(const std::vector<std::array<double, 3>>* verts,
				const std::vector<std::array<int, 3>>* faces, const std::vector<int>* vIds = nullptr, bool removeIsolatedVerts = true) = delete;

			void read_obj(const char* filename, bool removeIsolatedVertices = true) = delete;
			
			void read_ply(const char* input) = delete;
			void read_off(const char* input) = delete;

			void copy() = delete;

		};
	}
}

#endif