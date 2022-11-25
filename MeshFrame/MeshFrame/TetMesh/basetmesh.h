/*!
*      \file BaseTetMesh.h
*      \brief Base TetMesh Class for all types of Tetrahedron Mesh Classes
*
*		This is the fundamental class for tetrahedral meshes
*	   \author Anka Chen
*      \date 10/10/2022
*
*/

#ifndef _MESHFRAME_BASE_TET_MESH_H_
#define _MESHFRAME_BASE_TET_MESH_H_

#include <assert.h>
#include <list>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <iterator>
#include <memory>

#include "../Parser/StrUtil_fast.h"
#include "../Parser/IOFuncDef.h"
#include "../Memory/MemoryPool.h"
#include "../Memory/Array.h"

#include "TetMeshTypeDefs.h"
#include "../Types/TypeDefs.h"

#include "TProps.h"


#ifndef MAX_LINE 
#define MAX_LINE 2048
#endif

#define TMESH_ARRAY_PRE_ALLOC_SIZE 32
namespace MF
{
	namespace TetMesh
	{
		class CVertexBase;

		/*!
		* \brief CBaseTMesh, base class for all types of tet-mesh classes; 
		* 
		*  This is the fundamental class for tet-meshes. All the geometric objects are connected by pointers,
		*  vertex, edge, face, tet are connected by halffaces. The mesh class has file IO functionalities,
		*  supporting .tet file formats. It offers Euler operators, each geometric primative
		*  can access its neighbors freely.
		*  Note that this base class only contains topological informations. 
		*
		* \tparam VertexType   vertex   class, derived from MF::TetMesh::VertexType     class
		* \tparam TVertexType  tetrahedron vertex   class, derived from MF::TetMesh::TVertexType   class
		* \tparam HalfEdgeType halfedge class, derived from MF::TetMesh::HalfEdgeType class
		* \tparam TEdgeType	tetrahedron edge class, derived from MF::TetMesh::TEdgeType class
		* \tparam EdgeType     edge     class, derived from MeshLib::EdgeType     class
		* \tparam FaceType     face     class, derived from MF::TetMesh::FaceType     class
		* \tparam HalfFaceType half face     class, derived from MF::TetMesh::HalfFaceType     class
		* \tparam TetType      tetrahedron class, derived from MF::TetMesh::TetType class
		*/

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		class CTMeshBase
		{
		public:
			typedef DType     DType;
			typedef TVertexType     TVType;
			typedef VertexType      VType;
			typedef HalfEdgeType    HEType;
			typedef TEdgeType       TEType;
			typedef TetType         TType;
			typedef EdgeType        EType;
			typedef HalfFaceType    HFType;
			typedef FaceType        FType;

			typedef TVertexType* TVPtr;
			typedef VertexType* VPtr;
			typedef HalfEdgeType* HEPtr;
			typedef TEdgeType* TEPtr;
			typedef TetType* TPtr;
			typedef EdgeType* EPtr;
			typedef HalfFaceType* HFPtr;
			typedef FaceType* FPtr;

			typedef TVec3<DType> Vec3;

			typedef  std::pair<int, VertexType*>		VMapPair;
			typedef  std::pair<int, TetType*>			TMapPair;


			typedef MemoryPool<VertexType>				VContainer;
			typedef MemoryPool<TVertexType>				TVContainer;
			typedef MemoryPool<FaceType>				FContainer;
			typedef MemoryPool<HalfFaceType>			HFContainer;
			typedef MemoryPool<EdgeType>				EContainer;
			typedef MemoryPool<TEdgeType>				TEContainer;
			typedef MemoryPool<HalfEdgeType>			HEContainer;
			typedef MemoryPool<TetType>					TContainer;

			typedef CTMeshBase<DType, TVType, VType, HEType, TEType, EType, HFType, FType, TType>* Ptr;
			typedef std::shared_ptr<CTMeshBase<DType, TVType, VType, HEType, TEType, EType, HFType, FType, TType>> SharedPtr;

			typedef CPArray<HFPtr, TMESH_ARRAY_PRE_ALLOC_SIZE> HFArray;
			typedef CPArray<TEPtr, TMESH_ARRAY_PRE_ALLOC_SIZE> TEArray;

			/*!
				CTMeshBase constructor
				*/
			CTMeshBase() { std::ios::sync_with_stdio(false); };
			/*!
				CTMeshBase desctructor
				*/
			~CTMeshBase() { _clear(); };
		
			/*!
			Load tet mesh from a ".vtk" file
			*/
			void _load_vtk(const char* input, bool checkOrientation = false);

			/*!
			Load tet mesh from a ".t" file
			*/
			void _load_vtArray(const std::vector<std::array<double, 3>>& verts, const std::vector<std::array<int, 4>>& tetVIds, bool checkOrientation = false);

			/*!
				Write tet mesh to a file
				*/
			void _write(const char *);

			/*!
			Write tet mesh to a .t file
			*/
			void _write_t(const char * filename, bool highPrecision = false);


			/*!
			Write tet mesh to a .vtk file
			*/
			void _write_vtk(const char* filename, bool highPrecision = false);

			/*!
			Write selected tets as a simple tet mesh to a .vtk file
			*/
			void _write_tet_list_to_vtk(const char* filename, std::vector<TPtr> tets, bool highPrecision = false);

			/*!
				access the list of half faces
				*/
			HFContainer & half_faces() { return mHFContainer; };
			/*!
				access the list of edges
				*/
			EContainer       & edges() { return mEContainer; };
			/*!
				access list of faces
				*/
			FContainer        & faces() { return mFContainer; };

			HEContainer& half_edges() { return mHEContainer; };

			TEContainer& tedges() { return mTEContainer; };

			/*!
			access list of vertices
			*/
			VContainer & vertices() { return mVContainer; };
			TVContainer & tvertices() { return mTVContainer; };

			/*! number of tets */
			int numTets() { return m_nTets; };

			/*! number of edges */
			int numEdges() { return mEContainer.size(); };

			/*! number of vertices */
			int numVertices() { return m_nVertices; };

			/*! max vertex id*/
			int maxVertexId() { return m_maxVertexId; };

			/*! Access the array of tets */
			TContainer& tets() { return mTContainer; };
			const TContainer & tets() const { return mTContainer; };

			/*! access the vertex with ID */
			virtual VertexType * idVertex(int id) { return m_map_Vertices[id]; };

			/*! access the tet with ID */
			virtual TetType      * idTet(int id) { return m_map_Tets[id]; };


			//Access Vertex data members
			/*! Vertex->Edge List */
			static std::vector<EdgeType*>* VertexEdgeList(VertexType* pVertex);
			/*! Vertex->TEdge List */
			TEArray& VertexTEdgeList(VertexType* pVertex);
			/*! Vertex->HalfFace List */
			static std::list<HalfFaceType*>* VertexHalfFaceList(VertexType* pVertex);
			/*! Vertex->TVertex List */
			static std::vector<TVertexType*>* VertexTVertexList(VertexType* pVertex);

			/*! Vertex->Edge */
			static EdgeType* VertexEdge(VertexType* v1, VertexType* v2);

			//Access TVertex data memebers
			static VertexType* TVertexVertex(TVertexType* pTVertex);
			static TetType* TVertexTet(TVertexType* pTVertex);
			static HalfEdgeType* TVertexHalfEdge(TVertexType* pTVertex);
			static HalfFaceType* TVertexOppositeHalfFace(TVertexType* pTVertex);

			//Access TEdge data memebers
			static HalfEdgeType* TEdgeLeftHalfEdge(TEdgeType* pTEdge);
			static HalfEdgeType* TEdgeRightHalfEdge(TEdgeType* pTEdge);
			static EdgeType* TEdgeEdge(TEdgeType* pTEdge);
			static TEdgeType* TEdgeDualTEdge(TEdgeType* pTEdge);
			static TetType* TEdgeTet(TEdgeType* pTEdge);
			static HalfFaceType* TEdgeLeftHalfFace(TEdgeType* pTEdge);
			static HalfFaceType* TEdgeRightHalfFace(TEdgeType* pTEdge);

			//Access HalfEdge data members
			/*! HalfEdge->source vertex */
			static VertexType* HalfEdgeSource(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->target vertex */
			static VertexType* HalfEdgeTarget(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->source tvertex */
			static TVertexType* HalfEdgeTSource(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->target tvertex */
			static TVertexType* HalfEdgeTTarget(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->tet */
			static TetType* HalfEdgeTet(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->dual halfedge */
			static HalfEdgeType* HalfEdgeDual(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->next HalfEdge */
			static HalfEdgeType* HalfEdgeNext(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->prev HalfEdge */
			static HalfEdgeType* HalfEdgePrev(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->Edge Edge */
			static EdgeType* HalfEdgeEdge(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->TEdge TEdge */
			static TEdgeType* HalfEdgeTEdge(HalfEdgeType* pHalfEdge);
			/*! HalfEdge->HalfFace */
			static HalfFaceType* HalfEdgeHalfFace(HalfEdgeType* pHalfEdge);
			/*! Turn halfedge into vector in CPoint */
			template <typename EigenDerived3x1>
			static void  HalfEdgeVec(HalfEdgeType* pHalfEdge, EigenDerived3x1& v);

			//Access Edge data members
			/*! TEdge list of the edge */
			static std::vector<TEdgeType*>* EdgeTEdgeList(EdgeType* pEdge);
			/*! Edge->Vertex1 */
			static VertexType* EdgeVertex1(EdgeType* pEdge);
			/*! Edge->Vertex2 */
			static VertexType* EdgeVertex2(EdgeType* pEdge);
			/*! length of the edge*/
			static double EdgeLength(EdgeType* pEdge);
			/*! squared length of the edge*/
			static double EdgeLengthSquare(EdgeType* pEdge);

			//Access HalfFace data memebers
			/*! HalfFace->HalfEdge */
			static HalfEdgeType* HalfFaceHalfEdge(HalfFaceType* pHalfFace);
			/*! HalfFace->face */
			static FaceType* HalfFaceFace(HalfFaceType* pHalfFace);
			/*! HalfFace->Tet */
			static TetType* HalfFaceTet(HalfFaceType* pHalfFace);
			/*! HalfFace->dual half face */
			static HalfFaceType* HalfFaceDual(HalfFaceType* pHalfFace);
			/*! HalfFace's opposite tvertex, i.e, the tvertex not contain in the halfface */
			static TVertexType* HalfFaceOppositeTVertex(HalfFaceType* pHalfFace);

			template <typename EigenDerived3x3>
			void HalfFace3Points(HalfFaceType* pHF, EigenDerived3x3& vs);

			template <typename EigenDerived3x1>
			static void HalfFaceNormal(HalfFaceType* pHF, EigenDerived3x1& normal);

			template <typename EigenDerived3x1>
			static void HalfFaceOrientedArea(HalfFaceType* pHF, EigenDerived3x1& orientedArea);

			template <typename EigenDerived3x1>
			static bool PointInTet(TPtr pT, const HalfFaceType& p);

			//Operations of new a object
			/*! New a vertex in tmesh */
			VertexType* newVertex();

			/*! Create a vertex with input id, used in the process of reading .t file, id is user defined
			\return pointer to the new vertex
			*/
			VertexType* createVertexWithId(int id);

			/*! Create a vertex with automatically assigned index, faster than createVertexWithId
			* the index is its index in its containter, its id will be set to =index
						\return pointer to the new vertex
						*/
			VertexType* createVertexWithIndex();

			TVertexType* createTVertex();

			/*! Create a face with automatically assigned index, 
			* the index is the index in its containter, its id will be set to =index
			*
			\return pointer to the new face
			*/
			FaceType* createFace();


			/*! Create a halfface with automatically assigned index,
			* the index is its index in its containter, its id will be set to =index
			\return pointer to the new halfface
			*/
			HalfFaceType* createHalfFaceWithIndex();

			/*! Create a halfedge with automatically assigned index,
			* the index is its index in its containter, its id will be set to =index
			\return pointer to the new halfedge
			*/
			HalfEdgeType* createHalfEdgeWithIndex();

			/*! Create a edge with automatically assigned index,
			* the index is its index in its containter, its id will be set to =index
			\return pointer to the new edge
			*/
			EdgeType* createEdgeWithIndex();

			/*! Create a tedge with automatically assigned index,
			* the index is its index in its containter, its id will be set to =index
			\return pointer to the new tedge
			*/
			TEdgeType* createTEdgeWithIndex();

			/*! Create a tedge with automatically assigned index,
			* the index is its index in its containter, its id will be set to =index
				\return pointer to the new tedge
				*/
			TetType* createTetWithIndex();


			/*! Create a vertex with input id, used in the process of reading .t file, id is user defined
				\return pointer to the new tedge
				*/
			TetType* createTetWithId(int id);

			//Face
			/*! access the left half face of a face */
			static HalfFaceType* FaceLeftHalfFace(FaceType* pFace);
			/*! access the right half face of a face */
			static HalfFaceType* FaceRightHalfFace(FaceType* pFace);

			//Tetrahedron

			/*! access the j-th half edge of a tet */
			static HalfFaceType* TetHalfFace(TetType* pT, int j);
			/*! access the j-th tvertex of a tet */
			static TVertexType* TetTVertex(TetType* pT, int j);
			/*! access the j-th vertex of a tet */
			static VertexType* TetVertex(TetType* pT, int j);

			static double TetOrientedVolume(TetType* pT);

			template <typename EigenDerived3x1>
			static void TetCentroid(TetType* pT, EigenDerived3x1& centroid);


			void reinitializeVIds();

			void tetMeshSurfaceMesh(std::vector<VertexType *>& verts, std::vector<HalfFaceType *>& faces);
		protected:

			/*!
			construct tetrahedron
			\tparam v array of vertex ids
			\tparam pT retulting tetrahedron
			*/

			void  _construct_tet(TetType* pT, int tID, int * v);
			void  _construct_tet_orientation(TetType* pT, int tId, int  v[4]);
			/*! construct faces */
			void  _construct_faces();
			/*! construct edges */
			void  _construct_edges();
			/*!  construct half face
			\tparam array of vertex pointers
			*/
			HalfFaceType*   _construct_half_face(TVertexType **);

			/*! release all the memory allocations */
			void _clear();



		protected:
			/*!
			list of faces
			*/
			FContainer        mFContainer;
			/*!
			list of half faces
			*/
			HFContainer	 mHFContainer;
			/*!
			list of half edges
			*/
			HEContainer	 mHEContainer;
			/*!
			list of edges
			*/
			EContainer       mEContainer;
			/*!
			list of tetrahedra
			*/
			TEContainer		 mTEContainer;

			/*!
			 array of vertices
			 */
			VContainer		 mVContainer;

			/*!
			 array of vertices
			 */
			TVContainer		 mTVContainer;
			//VertexType *				 mVContainer;

			/*!
			map of VertexType id and pointer
			*/
			std::map<int, VertexType *> m_map_Vertices;

			/*!
			array of tets
			*/
			TContainer		 mTContainer;
			//TetType*                    mTContainer;

			std::map<int, TetType*>     m_map_Tets;

			/*! number of vertices */
			int m_nVertices;

			/*! number of tets */
			int m_nTets;

			/*! number of edges */
			int m_nEdges;

			/*! max vertex id */
			int m_maxVertexId;


			MAKE_PROP_OF(V);
			MAKE_PROP_OF(E);
			MAKE_PROP_OF(F);
			MAKE_PROP_OF(T);
			MAKE_PROP_OF(HF);
			MAKE_PROP_OF(HE);

			VPropHandle<HFArray> mVHFArrayHandle;
			VPropHandle<TEArray> mVTEArrayHandle;

		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_clear()
		{

		};


		

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_load_vtArray(
			const std::vector<std::array<double, 3>>& verts, const std::vector<std::array<int, 4>>& tetVIds, bool checkOrientation)
		{
			addVProp(mVHFArrayHandle);
			addVProp(mVTEArrayHandle);

			m_maxVertexId = verts.size()-1;

			m_nVertices = verts.size();
			m_nTets = tetVIds.size();
			m_nEdges = 0;

			//read in the vertices
			for (int i = 0; i < m_nVertices; i++)
			{
				int vid = i;

				Vec3 p(verts[i][0], verts[i][1], verts[i][2]);

				VertexType* v = createVertexWithId(vid);
				v->position() = p;
			}


			//read in tets 
			for (int id = 0; id < m_nTets; id++)
			{
				int tid = id;
				int vIds[4] = { tetVIds[id][0], tetVIds[id][1], tetVIds[id][2], tetVIds[id][3]};

				TetType* pT = createTetWithId(tid);

				if (checkOrientation) {
					_construct_tet_orientation(pT, tid, vIds);
				}
				else {
					_construct_tet(pT, tid, vIds);
				}

			}

			_construct_faces();
			_construct_edges();

			m_nEdges = (int)mEContainer.size();

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				if (pV->id() > m_maxVertexId)
				{
					m_maxVertexId = pV->id();
				}
			}

			// label the boundary for faces and vertices
			for (auto fIter = mFContainer.begin(); fIter != mFContainer.end(); ++fIter)
			{
				FPtr pF = *fIter;
				if (this->FaceLeftHalfFace(pF) == NULL || this->FaceRightHalfFace(pF) == NULL)
				{
					pF->boundary() = true;
					HalfFaceType* pH =
						FaceLeftHalfFace(pF) == NULL ? FaceRightHalfFace(pF) : FaceLeftHalfFace(pF);
					//added by Anka, mark edge as boundary

					for (int i = 0; i < 3; ++i)
					{
						int vid = pH->key(i);
						VertexType* v = idVertex(vid);
					}
				}
			}

			// read in traits
			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				pV->edges()->shrink_to_fit();
				pV->tvertices()->shrink_to_fit();
			}

			for (auto tIter = mTContainer.begin(); tIter != mTContainer.end(); tIter++)
			{
				TetType* pT = *tIter;
			}

			for (auto eIter = mEContainer.begin(); eIter != mEContainer.end(); eIter++)
			{
				EdgeType* pE = *eIter;
			}

			removeVProp(mVTEArrayHandle);
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>

		HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_construct_half_face(TVertexType ** pTV)
		{
			HalfFaceType* pHF = createHalfFaceWithIndex();

			VertexType * pV[3];

			for (int i = 0; i < 3; i++)
			{
				pV[i] = TVertexVertex(pTV[i]);
			}

			HalfEdgeType * pH[3];
			for (int i = 0; i < 3; i++)
			{
				pH[i] = createHalfEdgeWithIndex();

				pH[i]->SetHalfFace(pHF);
				//pH[i]->
				
				//SetSource(pTV[i]);
				pH[i]->SetTarget(pTV[(i + 1) % 3]);
				pTV[i]->set_halfedge(pH[i]);
			}

			for (int i = 0; i < 3; i++)
			{
				pH[i]->SetNext(pH[(i + 1) % 3]);
				pH[i]->SetPrev(pH[(i + 2) % 3]);
			}

			pHF->SetHalfEdge(pH[0]);

			for (int i = 0; i < 3; i++)
			{
				pHF->key(i) = pV[i]->id();
			}

			//bubble

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2 - i; j++)
				{
					if (pHF->key(j) > pHF->key(j + 1))
					{
						int tmp = pHF->key(j);
						pHF->key(j) = pHF->key(j + 1);
						pHF->key(j + 1) = tmp;
					}
				}
			}

			assert(pHF->key(0) < pHF->key(1));
			assert(pHF->key(1) < pHF->key(2));

			VertexType * pv = m_map_Vertices[pHF->key(0)];
			pv->HalfFaces()->push_back(pHF);

			return pHF;
		};

		//construct faces
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_construct_faces()
		{
			VertexType * pV = NULL;

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				pV = *vIter;

				std::list<HalfFaceType*> * pL = VertexHalfFaceList(pV);

				while (!pL->empty())
				{
					HalfFaceType * pF = pL->front();
					pL->pop_front();
					FaceType* f = createFace();
					f->SetLeft(pF);
					pF->SetFace(f);

					for (typename std::list<HalfFaceType*>::iterator it = pL->begin(); it != pL->end(); it++)
					{
						HalfFaceType * pH = *it;

						if (*pH == *pF)
						{
							pH->SetDual(pF);
							pF->SetDual(pH);
							f->SetRight(pH);
							pH->SetFace(f);
							break;
						}
					}

					if (pF->dual() != NULL)
					{
						pL->remove(HalfFaceDual(pF));
					}
				}
			}
		};

		//construct edges
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_construct_edges()
		{
			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType * pV = *vIter;
				TEArray& vteArr = VertexTEdgeList(pV);

				while (vteArr.size() != 0)
				{
					TEdgeType * pTE = vteArr.back();
					vteArr.pop_back();

					EdgeType * e = createEdgeWithIndex();
					assert(e != NULL);

					int id1 = pTE->key(0);
					VertexType * v1 = m_map_Vertices[id1];
					e->SetVertex1(v1);

					int id2 = pTE->key(1);
					VertexType * v2 = m_map_Vertices[id2];
					e->SetVertex2(v2);

					e->edges()->push_back(pTE);
					pTE->SetEdge(e);

					// find tedges that share the same edge
					// adding them to the correponding edge, and delete it from the array
					TEArray tmp_tedges;
					size_t arrSize = vteArr.size();
					for (size_t arrIdx =0; arrIdx < arrSize; arrIdx++)
					{
						TEPtr pTEMatch = vteArr[arrIdx];
						if (*pTEMatch == *pTE)
						{
							pTEMatch->SetEdge(e);
							tmp_tedges.push_back(pTEMatch);

							e->edges()->push_back(pTEMatch);
							vteArr.erase(arrIdx);
							--arrIdx;
							--arrSize;
							// if ereasing element at arrIdx we do not increase arrIdx
						}

					}

					//// adding those 

					//for (std::list<TEdgeType*>::iterator it = tmp_edges.begin(); it != tmp_edges.end(); it++)
					//{
					//	TEdgeType * pH = *it;
					//	pL->remove(pH);
					//	e->edges()->push_back(pH);
					//}

				}

			}


			for (auto it = mEContainer.begin(); it != mEContainer.end(); it++)
			{
				EdgeType * pE = *it;
				VertexType * v1 = EdgeVertex1(pE);
				VertexType * v2 = EdgeVertex2(pE);
				v1->edges()->push_back(pE);
				v2->edges()->push_back(pE);
			}
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_construct_tet(TetType* pT, int tId, int * v)
		{
			//set the tet->id

			pT->id() = tId;

			//set TVertices of the Tet

			for (int k = 0; k < 4; k++)
			{
				TVertexType * pTV = createTVertex();
				pT->setTVertex(pTV, k);
				pTV->id() = k;

				VertexType * pV = m_map_Vertices[v[k]];
				pTV->set_vert((CVertexBase*)pV);
				pV->tvertices()->push_back(pTV);

				pTV->set_tet(pT);
			}

			//set half faces

			int order[4][3] = { { 1, 2, 3 },{ 2, 0, 3 },{ 0, 1, 3 },{ 1, 0, 2 } };

			TVertexType   * pTV[3];
			HalfFaceType * pHF[4];

			for (int i = 0; i < 4; i++)
			{
				for (int k = 0; k < 3; k++)
				{
					pTV[k] = TetTVertex(pT, order[i][k]);
				}
				pT->setHalfFace(_construct_half_face(pTV), i);
				pHF[i] = TetHalfFace(pT, i);
			}

			// connect the four half faces

			for (int i = 0; i < 4; i++)
			{
				pHF[i]->SetTet(pT);
			}

			//Seting the dual half edges

			for (int i = 0; i < 3; i++)
			{
				HalfEdgeType * pH0 = HalfFaceHalfEdge(pHF[i]);
				pH0 = HalfEdgeNext(pH0);
				HalfEdgeType * pH1 = HalfFaceHalfEdge(pHF[(i + 1) % 3]);
				pH1 = HalfEdgePrev(pH1);

				pH0->SetDual(pH1);
				pH1->SetDual(pH0);

				TEdgeType * pTE = createTEdgeWithIndex();
				assert(pTE != NULL);
				pTE->SetTet(pT);
				pH0->SetTEdge(pTE);
				pH1->SetTEdge(pTE);

				if (pH0->source()->id() < pH0->target()->id())
				{
					pTE->SetLeft(pH0);
					pTE->SetRight(pH1);
				}
				else
				{
					pTE->SetLeft(pH1);
					pTE->SetRight(pH0);
				}

				pTE->key(0) = pTE->left()->source()->id();
				pTE->key(1) = pTE->left()->target()->id();

				VertexType * v = m_map_Vertices[pTE->key(0)];
				VertexTEdgeList(v).push_back(pTE);
			}

			HalfEdgeType * pH0 = HalfFaceHalfEdge(pHF[3]);
			for (int i = 0; i < 3; i++)
			{
				HalfEdgeType * pH1 = HalfFaceHalfEdge(pHF[2 - i]);
				pH0->SetDual(pH1);
				pH1->SetDual(pH0);

				TEdgeType * pTE = createTEdgeWithIndex();
				assert(pTE != NULL);
				//set TEdge->Tet
				pTE->SetTet(pT);
				//set HalfEdge->TEdge
				pH0->SetTEdge(pTE);
				pH1->SetTEdge(pTE);

				if (pH0->source()->id() < pH0->target()->id())
				{
					pTE->SetLeft(pH0);
					pTE->SetRight(pH1);
				}
				else
				{
					pTE->SetLeft(pH1);
					pTE->SetRight(pH0);
				}
				pTE->key(0) = pTE->left()->source()->id();
				pTE->key(1) = pTE->left()->target()->id();

				VertexType * v = m_map_Vertices[pTE->key(0)];
				VertexTEdgeList(v).push_back(pTE);

				pH0 = HalfEdgeNext(pH0);
			}
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_construct_tet_orientation(TetType* pT, int tId, int  v[4])
		{
			//orient the tet 
			Vec3 A = m_map_Vertices[v[0]]->position();
			Vec3 B = m_map_Vertices[v[1]]->position();
			Vec3 C = m_map_Vertices[v[2]]->position();
			Vec3 D = m_map_Vertices[v[3]]->position();
			Vec3 AB = B - A;
			Vec3 AC = C - A;
			Vec3 AD = D - A;

			double orientation_product = AB.dot(AC.cross(AD));
			if (orientation_product < 0) {
				int temp = v[2];
				v[2] = v[3];
				v[3] = temp;
			}
			//set the tet->id

			pT->id() = tId;

			//set TVertices of the Tet

			for (int k = 0; k < 4; k++)
			{
				TVertexType * pTV = createTVertex();
				pT->setTVertex(pTV, k);
				pTV->id() = k;

				VertexType * pV = m_map_Vertices[v[k]];
				pTV->set_vert(pV);
				pV->tvertices()->push_back(pTV);

				pTV->set_tet(pT);
			}

			//set half faces

			int order[4][3] = { { 1, 2, 3 }, { 2, 0, 3 }, { 0, 1, 3 }, { 1, 0, 2 } };

			TVertexType   * pTV[3];
			HalfFaceType * pHF[4];

			for (int i = 0; i < 4; i++)
			{
				for (int k = 0; k < 3; k++)
				{
					pTV[k] = TetTVertex(pT, order[i][k]);
				}
				pT->setHalfFace(_construct_half_face(pTV), i);
				pHF[i] = TetHalfFace(pT, i);
			}

			// connect the four half faces

			for (int i = 0; i < 4; i++)
			{
				pHF[i]->SetTet(pT);
			}

			//Seting the dual half edges

			for (int i = 0; i < 3; i++)
			{
				HalfEdgeType * pH0 = HalfFaceHalfEdge(pHF[i]);
				pH0 = HalfEdgeNext(pH0);
				HalfEdgeType * pH1 = HalfFaceHalfEdge(pHF[(i + 1) % 3]);
				pH1 = HalfEdgePrev(pH1);

				pH0->SetDual(pH1);
				pH1->SetDual(pH0);

				TEdgeType * pTE = createTEdgeWithIndex();
				pTE->SetTet(pT);
				pH0->SetTEdge(pTE);
				pH1->SetTEdge(pTE);

				if (pH0->source()->id() < pH0->target()->id())
				{
					pTE->SetLeft(pH0);
					pTE->SetRight(pH1);
				}
				else
				{
					pTE->SetLeft(pH1);
					pTE->SetRight(pH0);
				}

				pTE->key(0) = pTE->left()->source()->id();
				pTE->key(1) = pTE->left()->target()->id();

				VertexType * v = m_map_Vertices[pTE->key(0)];
				VertexTEdgeList(v).push_back(pTE);
			}

			HalfEdgeType * pH0 = HalfFaceHalfEdge(pHF[3]);
			for (int i = 0; i < 3; i++)
			{
				HalfEdgeType * pH1 = HalfFaceHalfEdge(pHF[2 - i]);
				pH0->SetDual(pH1);
				pH1->SetDual(pH0);

				TEdgeType * pTE = createTEdgeWithIndex();
				//set TEdge->Tet
				pTE->SetTet(pT);
				//set HalfEdge->TEdge
				pH0->SetTEdge(pTE);
				pH1->SetTEdge(pTE);

				if (pH0->source()->id() < pH0->target()->id())
				{
					pTE->SetLeft(pH0);
					pTE->SetRight(pH1);
				}
				else
				{
					pTE->SetLeft(pH1);
					pTE->SetRight(pH0);
				}
				pTE->key(0) = pTE->left()->source()->id();
				pTE->key(1) = pTE->left()->target()->id();

				VertexType * v = m_map_Vertices[pTE->key(0)];
				VertexTEdgeList(v).push_back(pTE);

				pH0 = HalfEdgeNext(pH0);
			}
		};

		//write tet mesh to the file

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_write(const char * output)
		{

			std::fstream _os(output, std::fstream::out);
			if (_os.fail())
			{
				fprintf(stderr, "Error is opening file %s\n", output);
				return;
			}
			_os << m_nVertices << " vertices" << std::endl;
			_os << m_nTets << " tets" << std::endl;

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType * pV = *vIter;
				Vec3  p = pV->position();
				for (int k = 0; k < 3; k++)
				{
					_os << " " << p[k];
				}
				if (pV->string().size() > 0)
				{
					_os << " " << "{" << pV->string() << "}";
				}
				_os << std::endl;
			}

			for (int i = 0; i < m_nTets; i++)
			{
				TetType * pt = m_map_Tets[i];
				_os << 4;
				for (int k = 0; k < 4; k++)
				{
					_os << " " << pt->tvertex(k)->vert()->id();
				}
				//if( pt->string().size() > 0 )
				//{
				//	_os << " " << "{"<< pt->string() << "}";
				//}
				_os << std::endl;
			}

			for (auto eIter = mEContainer.begin(); eIter != mEContainer.end(); eIter++)
			{
				EdgeType * pE = *eIter;
				if (pE->string().size() > 0)
				{
					_os << "Edge " << EdgeVertex1(pE)->id() << " " << EdgeVertex2(pE)->id() << " ";
					_os << "{" << pE->string() << "}" << std::endl;
				}
			}

			_os.close();
		};

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_write_t(const char * output, bool highPrecision)
		{
			//write traits to string, add by Wei Chen, 11/23/2015
			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType * pV = *vIter;
				pV->_to_string();
			}

			for (auto tIter = mTContainer.begin(); tIter != mTContainer.end(); tIter++)
			{
				TetType * pT = *tIter;
				pT->_to_string();
			}

			for (auto eIter = mEContainer.begin(); eIter != mEContainer.end(); eIter++)
			{
				EdgeType * pE = *eIter;
				pE->_to_string();
			}
			//write traits end

			std::fstream _os(output, std::fstream::out);
			if (_os.fail())
			{
				fprintf(stderr, "Error while opening file %s\n", output);
				return;
			}

			if (highPrecision) {
				_os << std::setiosflags(std::ios::fixed);
			}

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType * pV = *vIter;
				Vec3 p = pV->position();
				_os << "Vertex " << pV->id();
				for (int k = 0; k < 3; k++)
				{
					if (highPrecision) {
						_os << " " << std::setprecision(16) << p[k];
					}
					else {
						_os << " " << p[k];
					}
				}
				//if (pV->string().size() > 0)
				//{
				//	_os << " " << "{" << pV->string() << "}";
				//}
				_os << std::endl;
			}

			for (auto tIter = mTContainer.begin(); tIter != mTContainer.end(); tIter++)
			{
				TetType * pT = *tIter;
				_os << "Tet " << pT->id();
				for (int k = 0; k < 4; k++)
				{
					_os << " " << pT->tvertex(k)->vert()->id();
				}
				//if (pT->string().size() > 0)
				//{
				//	_os << " " << "{" << pT->string() << "}";
				//}
				_os << std::endl;
			}

			for (auto eIter = mEContainer.begin(); eIter != mEContainer.end(); eIter++)
			{
				EdgeType * pE = *eIter;
				//if (pE->string().size() > 0)
				//{
				//	_os << "Edge " << EdgeVertex1(pE)->id() << " " << EdgeVertex2(pE)->id() << " ";
				//	_os << "{" << pE->string() << "}" << std::endl;
				//}
			}

			_os.close();
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_write_vtk(const char* filename, bool highPrecision)
		{
			reinitializeVIds();

			std::fstream _os(filename, std::fstream::out);
			if (_os.fail())
			{
				fprintf(stderr, "Error while opening file %s\n", filename);
				return;
			}

			if (highPrecision) {
				_os << std::setiosflags(std::ios::fixed);
			}

			_os << "# vtk DataFile Version 2.0\n";
			_os << "Unstructured Grid\n";
			_os << "ASCII\n";
			_os << "DATASET UNSTRUCTURED_GRID\n";
			_os << "POINTS " << mVContainer.size() << " double\n";

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				Vec3 p = pV->position();
				for (int k = 0; k < 3; k++)
				{
					if (highPrecision) {
						_os << " " << std::setprecision(16) << p[k];
					}
					else {
						_os << " " << p[k];
					}
				}

				_os << std::endl;
			}
			_os << "\n";

			_os << "CELLS " << mTContainer.size() << ' ' << mTContainer.size()*5 << '\n';


			for (auto tIter = mTContainer.begin(); tIter != mTContainer.end(); tIter++)
			{
				TetType* pT = *tIter;
				_os << "4 ";
				for (int k = 0; k < 4; k++)
				{
					_os << " " << pT->tvertex(k)->vert()->id();
				}
				_os << std::endl;
			}
			_os << "\n";

			_os << "CELL_TYPES " << mTContainer.size()  << '\n';
			for (auto tIter = mTContainer.begin(); tIter != mTContainer.end(); tIter++)
			{
				_os << "10\n";
			}

			_os.close();
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::_write_tet_list_to_vtk(const char* filename, std::vector<TPtr> tets, bool highPrecision)
		{
			VPropHandle<int> vNewIdsHdl;
			addVProp(vNewIdsHdl);

			for (auto vIter = mVContainer.begin(); vIter != mVContainer.end(); vIter++)
			{
				VertexType* pV = *vIter;
				gVP(vNewIdsHdl, pV) = -1;
			}

			int currentVIds = 0;
			std::vector<VPtr> appearingVerts;
			for (auto tIter = tets.begin(); tIter != tets.end(); tIter++)
			{
				TPtr pT = *tIter;
				for (int k = 0; k < 4; k++)
				{
					VertexType* pV = (VertexType*)pT->tvertex(k)->vert();

					int & vId = gVP(vNewIdsHdl, pV);
					if (vId == -1) {
						vId = currentVIds;
						++currentVIds;
						appearingVerts.push_back(pV);
					}

				}
			}

			std::fstream _os(filename, std::fstream::out);
			if (_os.fail())
			{
				fprintf(stderr, "Error while opening file %s\n", filename);
				return;
			}

			if (highPrecision) {
				_os << std::setiosflags(std::ios::fixed);
			}

			_os << "# vtk DataFile Version 2.0\n";
			_os << "Unstructured Grid\n";
			_os << "ASCII\n";
			_os << "DATASET UNSTRUCTURED_GRID\n";
			_os << "POINTS " << appearingVerts.size() << " double\n";

			for (auto vIter = appearingVerts.begin(); vIter != appearingVerts.end(); vIter++)
			{
				VertexType* pV = *vIter;
				Vec3 p = pV->position();
				for (int k = 0; k < 3; k++)
				{
					if (highPrecision) {
						_os << " " << std::setprecision(16) << p[k];
					}
					else {
						_os << " " << p[k];
					}
				}

				_os << std::endl;
			}
			_os << "\n";

			_os << "CELLS " << tets.size() << ' ' << tets.size() * 5 << '\n';


			for (auto tIter = tets.begin(); tIter != tets.end(); tIter++)
			{
				TetType* pT = *tIter;
				_os << "4 ";
				for (int k = 0; k < 4; k++)
				{
					_os << " " << gVP(vNewIdsHdl, (VertexType*)pT->tvertex(k)->vert());
				}
				_os << std::endl;
			}
			_os << "\n";

			_os << "CELL_TYPES " << tets.size() << '\n';
			for (auto tIter = tets.begin(); tIter != tets.end(); tIter++)
			{
				_os << "10\n";
			}

			removeVProp(vNewIdsHdl);

			_os.close();
		}

		/*------------------------------------------------------------------------------------------------
		Access Vertex data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline std::vector<EdgeType*>* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::VertexEdgeList(VertexType* pVertex)
		{
			return (std::vector<EdgeType*>*) pVertex->edges();
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline CPArray<TEdgeType*, TMESH_ARRAY_PRE_ALLOC_SIZE> & CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::VertexTEdgeList(VertexType* pVertex)
		{
			return gVP(mVTEArrayHandle, pVertex);
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline std::list<HalfFaceType*>* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::VertexHalfFaceList(VertexType* pVertex)
		{
			return (std::list<HalfFaceType*>*) pVertex->HalfFaces();
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline std::vector<TVertexType*>* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::VertexTVertexList(VertexType* pVertex)
		{
			return (std::vector<TVertexType*>*) pVertex->tvertices();
		};

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline EdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::VertexEdge(VertexType* v1, VertexType* v2)
		{
			std::vector<EdgeType*>* vEdgeList = VertexEdgeList(v1);

			for (typename std::vector<EdgeType*>::iterator titer = (*vEdgeList).begin(); titer != (*vEdgeList).end(); titer++)
			{
				EdgeType* pE = *titer;

				VertexType* w1 = EdgeVertex1(pE);
				VertexType* w2 = EdgeVertex2(pE);

				if (w1 == v1 && w2 == v2)
				{
					return pE;
				}
				if (w1 == v2 && w2 == v1)
				{
					return pE;
				}
			}
			return NULL;
		}
		/*------------------------------------------------------------------------------------------------
		Access TVertex data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TVertexVertex(TVertexType* pTVertex)
		{
			return (VertexType*)pTVertex->vert();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TVertexTet(TVertexType* pTVertex)
		{
			return (TetType*)pTVertex->tet();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TVertexHalfEdge(TVertexType* pTVertex)
		{
			return (HalfEdgeType*)pTVertex->halfedge();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TVertexOppositeHalfFace(TVertexType* pTVertex)
		{
			return (HalfFaceType*)pTVertex->halfedge()->next()->dual()->half_face();
		}
		/*------------------------------------------------------------------------------------------------
		Access TEdge data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeLeftHalfEdge(TEdgeType* pTEdge)
		{
			return (HalfEdgeType*)pTEdge->left();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeRightHalfEdge(TEdgeType* pTEdge)
		{
			return (HalfEdgeType*)pTEdge->right();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline EdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeEdge(TEdgeType* pTEdge)
		{
			return (EdgeType*)pTEdge->edge();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeDualTEdge(TEdgeType* pTEdge)
		{
			return (TEdgeType*)pTEdge->dual();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeTet(TEdgeType* pTEdge)
		{
			return (TetType*)pTEdge->tet();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeLeftHalfFace(TEdgeType* pTEdge)
		{
			return HalfEdgeHalfFace(TEdgeLeftHalfEdge(pTEdge));
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TEdgeRightHalfFace(TEdgeType* pTEdge)
		{
			return HalfEdgeHalfFace(TEdgeRightHalfEdge(pTEdge));
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template<typename EigenDerived3x1>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeVec(HalfEdgeType* pHalfEdge, EigenDerived3x1& v)
		{
			v = HalfEdgeTarget(pHalfEdge)->position() - HalfEdgeSource(pHalfEdge)->position();
			
		}

		/*------------------------------------------------------------------------------------------------
		Access HalfEdge data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeSource(HalfEdgeType* pHalfEdge)
		{
			return (VertexType*)pHalfEdge->source();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeTarget(HalfEdgeType* pHalfEdge)
		{
			return (VertexType*)pHalfEdge->target();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TVertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeTSource(HalfEdgeType* pHalfEdge)
		{
			return (TVertexType*)pHalfEdge->tSource();
		}
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TVertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeTTarget(HalfEdgeType* pHalfEdge)
		{
			return (TVertexType*)pHalfEdge->tTarget();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeTet(HalfEdgeType* pHalfEdge)
		{
			return (TetType*)HalfEdgeTTarget(pHalfEdge)->tet();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeDual(HalfEdgeType* pHalfEdge)
		{
			return (HalfEdgeType*)pHalfEdge->dual();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeNext(HalfEdgeType* pHalfEdge)
		{
			return (HalfEdgeType*)pHalfEdge->next();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgePrev(HalfEdgeType* pHalfEdge)
		{
			return (HalfEdgeType*)pHalfEdge->prev();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline EdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeEdge(HalfEdgeType* pHalfEdge)
		{
			return (EdgeType*)pHalfEdge->tedge()->edge();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeTEdge(HalfEdgeType* pHalfEdge)
		{
			return (TEdgeType*)pHalfEdge->tedge();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfEdgeHalfFace(HalfEdgeType* pHalfEdge)
		{
			return (HalfFaceType*)pHalfEdge->half_face();
		}


		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline std::vector<TEdgeType*>* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::EdgeTEdgeList(EdgeType* pEdge)
		{
			return (std::vector<TEdgeType*>*) pEdge->edges();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::EdgeVertex1(EdgeType* pEdge)
		{
			return (VertexType*)pEdge->vertex1();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::EdgeVertex2(EdgeType* pEdge)
		{
			return (VertexType*)pEdge->vertex2();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline double CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::EdgeLength(EdgeType* pEdge)
		{
			VertexType* pV1 = EdgeVertex1(pEdge);
			VertexType* pV2 = EdgeVertex2(pEdge);

			return (pV1->position() - pV2->position()).norm();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline double CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::EdgeLengthSquare(EdgeType* pEdge)
		{
			return (EdgeVertex1(pEdge)->position() - EdgeVertex2(pEdge)->position()).norm2();
		}

		/*------------------------------------------------------------------------------------------------
		Access HalfFace data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceHalfEdge(HalfFaceType* pHalfFace)
		{
			return (HalfEdgeType*)pHalfFace->half_edge();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline FaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceFace(HalfFaceType* pHalfFace)
		{
			return (FaceType*)pHalfFace->face();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceTet(HalfFaceType* pHalfFace)
		{
			return (TetType*)pHalfFace->tet();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceDual(HalfFaceType* pHalfFace)
		{
			return (HalfFaceType*)pHalfFace->dual();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TVertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceOppositeTVertex(HalfFaceType* pHalfFace)
		{
			HalfEdgeType* pHE = HalfFaceHalfEdge(pHalfFace);
			HalfEdgeType* pHEDualNext = HalfEdgeNext(HalfEdgeDual(pHE));

			return HalfEdgeTTarget(pHEDualNext);
		}
		//template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		//inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFace3Points(HalfFaceType* pHF, CPoint* v)
		//{
		//	HalfEdgeType* pHE = HalfFaceHalfEdge(pHF);
		//	v[0] = HalfEdgeSource(pHE)->position();
		//	v[1] = HalfEdgeTarget(pHE)->position();
		//	v[2] = HalfEdgeTarget(HalfEdgeNext(pHE))->position();
		//}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template<typename EigenDerived3x3>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFace3Points(HalfFaceType* pHF, EigenDerived3x3& vs)
		{
			HalfEdgeType* pHE = HalfFaceHalfEdge(pHF);
			vs.block<3, 1>(0, 0) = HalfEdgeSource(pHE)->position();
			vs.block<3, 1>(0, 1) = HalfEdgeTarget(pHE)->position();
			vs.block<3, 1>(0, 2) = HalfEdgeTarget(HalfEdgeNext(pHE))->position();
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template<typename EigenDerived3x1>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceNormal(HalfFaceType* pHF, EigenDerived3x1& normal)
		{
			HalfEdgeType* pHE1 = HalfFaceHalfEdge(pHF);;
			HalfEdgeType* pHE2 = HalfEdgeNext(pHE1);;

			TVec3<DType> v1, v2;

			HalfEdgeVec(pHE1, v1);
			HalfEdgeVec(pHE1, v2);

			normal = v1.cross(v2);
			normal /= normal.norm();
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template<typename EigenDerived3x1>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::HalfFaceOrientedArea(HalfFaceType* pHF, EigenDerived3x1& orientedArea)
		{
			HalfEdgeType* pHE1 = HalfFaceHalfEdge(pHF);;
			HalfEdgeType* pHE2 = HalfEdgeNext(pHE1);;

			TVec3<DType> v1, v2;

			HalfEdgeVec(pHE1, v1);
			HalfEdgeVec(pHE1, v2);

			orientedArea = v1.cross(v2);
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template<typename EigenDerived3x1>
		inline bool CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::PointInTet(TPtr pT, const HalfFaceType& p)
		{
			register Vec3 vs4[4] = {
				pT->vertex(0)->position(),
				pT->vertex(1)->position(),
				pT->vertex(2)->position(),
				pT->vertex(3)->position()
			};

			Vec3 AB = vs4[1] - vs4[0];
			Vec3 AC = vs4[2] - vs4[0];
			Vec3 AD = vs4[3] - vs4[0];

			double tetOrientedVol = AB.dot(AC.cross(AD));

			const int order[4][3] = { { 1, 2, 3 },{ 2, 0, 3 },{ 0, 1, 3 },{ 1, 0, 2 } };

			for (int i = 0; i < 4; ++i) {

				Vec3 v1 = vs4[order[i][1]] - vs4[order[i][0]]; // HalfEdgeVec(pHE1);
				Vec3 v2 = vs4[order[i][2]] - vs4[order[i][1]];  // HalfEdgeVec(pHE2);
				Vec3 orientedAreaF = v1.cross(v2);
				if ((p - vs4[order[i][0]]).dot(orientedAreaF) * tetOrientedVol > 0)
				{
					return false;
				}
			}

			return true;
		}

		

		/*------------------------------------------------------------------------------------------------
		Access Face data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::FaceLeftHalfFace(FaceType* pFace)
		{
			return (HalfFaceType*)pFace->left();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::FaceRightHalfFace(FaceType* pFace)
		{
			return (HalfFaceType*)pFace->right();
		}

		/*------------------------------------------------------------------------------------------------
		Access Tetrahedron data members
		--------------------------------------------------------------------------------------------------*/
		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TetHalfFace(TetType* pT, int j)
		{
			return (HalfFaceType*)pT->half_face(j);
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TVertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TetTVertex(TetType* pT, int j)
		{
			return (TVertexType*)pT->tvertex(j);
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TetVertex(TetType* pT, int j)
		{
			return (VertexType*)pT->tvertex(j)->vert();
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline double CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TetOrientedVolume(TetType* pT)
		{
			auto& A = pT->vertex(0)->position();
			auto& B = pT->vertex(1)->position();
			auto& C = pT->vertex(2)->position();
			auto& D = pT->vertex(3)->position();
			Vec3 AB = B - A;
			Vec3 AC = C - A;
			Vec3 AD = D - A;

			double orientation_product = AB.dot(AC.cross(AD));
			return orientation_product;;
		}

		template <typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		template <typename EigenDerived3x1>
		inline void  CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::TetCentroid(TetType* pT, EigenDerived3x1& centroid)
		{
			centroid = TVec3<DType>::Zero();
			for (int i = 0; i < 4; ++i) {
				centroid += TetVertex(pT, i)->position();
			}
			centroid /= 4;
			return centroid;
		}

		//Operations of new a object
		/*! New a vertex in tmesh */
		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::newVertex()
		{
			size_t index;
			VertexType* pV = mVContainer.newMember(index);
			assert(pV != NULL);

			//for (int i = 0; i < mVProps.size(); i++) {
			//	BasicPropHandle* pProp = mVProps[i];
			//	pProp->initializePropMember(index);
			//}
			pV->index() = index;
			return pV;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createVertexWithId(int id)
		{
			VertexType* pV;
			pV = newVertex();
			m_map_Vertices.insert(VMapPair(id, pV));
			pV->id() = id;
			return pV;;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline VertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createVertexWithIndex()
		{
			VertexType* pV;
			pV = newVertex();
			pV->id() = (int)pV->index();
			m_map_Vertices.insert(VMapPair(pV->id(), pV));
			return pV;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TVertexType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createTVertex()
		{
			size_t index;
			TVertexType* pTV = mTVContainer.newMember(index);
			assert(pTV != NULL);
			pTV->index() = index;

			return pTV;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline FaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createFace()
		{
			size_t index;
			FaceType* pF = mFContainer.newMember(index);
			assert(pF != NULL);
			pF->index() = index;

			return pF;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfFaceType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createHalfFaceWithIndex()
		{
			size_t index;
			HalfFaceType* pHF = mHFContainer.newMember(index);
			assert(pHF != NULL);
			pHF->index() = index;

			return pHF;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline HalfEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createHalfEdgeWithIndex()
		{
			size_t index;
			HalfEdgeType* pHE = mHEContainer.newMember(index);
			assert(pHE != NULL);
			pHE->index() = index;

			return pHE;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline EdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createEdgeWithIndex()
		{
			size_t index;
			EdgeType* pE = mEContainer.newMember(index);
			assert(pE != NULL);
			pE->index() = index;

			return pE;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TEdgeType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createTEdgeWithIndex()
		{
			size_t index;
			TEdgeType* pTE = mTEContainer.newMember(index);
			assert(pTE != NULL);
			pTE->index() = index;

			return pTE;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createTetWithIndex()
		{
			size_t index;
			TetType* pT = mTContainer.newMember(index);
			assert(pT != NULL);
			m_map_Tets.insert(TMapPair(index, pT));
			pT->index() = index;

			return pT;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline TetType* CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::createTetWithId(int id)
		{
			size_t index;
			TetType* pT = mTContainer.newMember(index);
			assert(pT != NULL);
			m_map_Tets.insert(TMapPair(id, pT));
			pT->index() = index;

			return pT;
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::reinitializeVIds()
		{
			int currentId = 0;
			m_map_Vertices.clear();
			for (VPtr pV : mVContainer)
			{
				m_map_Vertices.insert(VMapPair(currentId, pV));
				pV->id() = currentId;
				++currentId;
			}
		}

		template<typename DType, typename TVertexType, typename VertexType, typename HalfEdgeType, typename TEdgeType, typename EdgeType, typename HalfFaceType, typename FaceType, typename TetType>
		inline void CTMeshBase<DType, TVertexType, VertexType, HalfEdgeType, TEdgeType, EdgeType, HalfFaceType, FaceType, TetType>::tetMeshSurfaceMesh(std::vector<VertexType*>& verts, std::vector<HalfFaceType*>& faces)
		{
			verts.clear();
			faces.clear();

			std::vector<HFPtr> surfaceHFList;
			auto pVless = [](VPtr pVa, VPtr pVb) {
				if (pVa->id() < pVb->id())
					return true;
				else
					return false;
			};
			std::set<VPtr, decltype(pVless)> vSet(pVless);
			for (HFPtr pHF : mHFContainer) {
				if (mHFContainer.hasBeenDeleted(pHF->index()))
				{
					continue;
				}
				if (HalfFaceDual(pHF) == NULL) {
					surfaceHFList.push_back(pHF);
					//for (auto pV : TIt::HF_VIterator(pHF)) {
					//	vSet.insert(pV);
					//}

					HalfEdgeType* pHE = (HalfEdgeType*)pHF->half_edge();

					for (int i = 0; i < 3; ++i)
					{
						int vid = pHF->key(i);
						VertexType* pV = idVertex(vid);
						vSet.insert(pV);
					}

					//vSet.insert((VertexType*)pHE->target());
					//pHE = (HalfEdgeType*)pHE->next();
					//vSet.insert((VertexType*)pHE->target());
					//pHE = (HalfEdgeType*)pHE->next();
					//vSet.insert((VertexType*)pHE->target());
				}
			}


			std::copy(vSet.begin(), vSet.end(), std::back_inserter(verts));
			std::copy(surfaceHFList.begin(), surfaceHFList.end(), std::back_inserter(faces));
		}

	};


};
#endif _MESHLIB_BASE_TET_MESH_H