/*!
*      \file BaseMesh.h
*      \brief Base Mesh Class for all types of Mesh Classes
*
*		This is the fundamental class for meshes
*
*	   \Version 2.0
*	   \Update 31/05/2018
*/

#ifndef _MESHLIB_BASE_MESH_H_
#define _MESHLIB_BASE_MESH_H_

#ifdef _MSC_VER
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 
#endif 
#endif

#include <math.h>
#include <assert.h>
#include <list>
#include <vector>
#include <map>
#include <omp.h>
#include <ctime>

#include "../Parser/strutil.h"
#include "../Parser/IOFuncDef.h"
#include "../Memory/MemoryPool.h"
#include "../FileIO/PlyFile.h"
#include "HalfEdge.h"
#include "Props.h"

#include "../Types/TypeDefs.h"


namespace MF {
	namespace TriMesh{
	/*!
	* \brief CMeshBase, base class for all types of mesh classes
	*
	*  This is the fundamental class for meshes. It includes a list of vertices,
	*  a list of edges, a list of faces. All the geometric objects are connected by pointers,
	*  vertex, edge, face are connected by halfedges. The mesh class has file IO functionalities,
	*  supporting .obj, .m and .off file formats. It offers Euler operators, each geometric primative
	*  can access its neighbors freely.
	*
	* \tparam VertexType   vertex   class, derived from MeshLib::CVertex   class
	* \tparam EdgeType     edge     class, derived from MeshLib::CEdge     class
	* \tparam FaceType     face     class, derived from MeshLib::CFace     class
	* \tparam HalfEdgeType halfedge class, derived from MeshLib::CHalfEdge class
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	class CMeshBase
	{
	public:
		typedef std::map<int, VertexType*>			VMap;
		typedef std::map<int, FaceType*>			FIdMap;
		typedef std::pair<int, VertexType*>			VMapPair;
		typedef std::pair<int, FaceType*>			FIdMapPair;
		typedef MemoryPool<VertexType>				VContainer;
		typedef MemoryPool<FaceType>				FContainer;
		typedef MemoryPool<EdgeType>				EContainer;
		typedef MemoryPool<HalfEdgeType>			HEContainer;

		typedef TVec3<DType> Vec3;
		typedef TVec2<DType> Vec2;

		// pointer to Vertices, Halfedges, Edges, Face and Solid
		typedef CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>* Ptr;
		typedef std::shared_ptr<CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>> SharedPtr;

		typedef DType DType;
		typedef VertexType   * VPtr;
		typedef EdgeType	 * EPtr;
		typedef FaceType	 * FPtr;
		typedef HalfEdgeType * HEPtr;

		typedef VertexType     VType;
		typedef EdgeType	   EType;
		typedef FaceType       FType;
		typedef HalfEdgeType   HEType;

		//constructor and destructor
		/*!
		CMeshBase constructor.
		*/
		CMeshBase() {};
		/*!
		CMeshBase destructor
		*/
		~CMeshBase();

		//copy operator
		/*!
		Copy operator
		*/
		void			copy(CMeshBase & rMesh);

		//file io
		/*!
		Read an .obj file.
		\param filename the input .obj file name
		*/
		void			read_obj(const char * filename, bool removeIsolatedVertices = true);

		void            readVFList(const std::vector<std::array<double, 3>>* verts, const std::vector<std::array<int, 3>>* faces, const std::vector<int>* vIds=nullptr, bool removeIsolatedVerts=true);

		/*!
		Write an .obj file.
		\param output the output .obj file name
		*/
		void			write_obj(const char * output);
		/*!
		Read an .m file.
		\param input the input .m file name
		*/
		void			read_m(const char * input);
		/*!
		Write an .m file.
		\param output the output .m file name
		*/
		void			write_m(const char * output, bool highPrecisionFloats = false);
		/*!
		Read an .ply file.
		\param input the input .ply file name
		*/
		/*Mesh type:: vertex = 1, face = 2, edge = 3, helfedge = 4, other edge = 5,*/
		enum			MeshType { PLY_V_Type = 1, PLY_F_Type = 2, PLY_E_Type = 3, PLY_HE_Type = 4, PLY_OTHER_Type = 5};
		void			read_ply(const char * input);
		void			mesh_ply_put_element(PlyFile* plyFile, void* voidPtr, MeshType meshType, PlyFileReader* plyFileReader);
		/*Get current element information from ply file*/
		void*			mesh_ply_get_element(PlyFile* plyFile, PlyElement* currentEle, PlyOtherElem* currentOtherEle, MeshType meshType, PlyFileReader* plyFileReader);

		/* reinitialize id() to make sure all the undeleted vertices' id() arrange tightly from 0 to N, may need to call it after deleting vertices */
		void            reinitializeVId();
		/* reinitialize id() to make sure all the undeleted vertices' id() arrange tightly from 0 to N, may need to call it after deleting faces */
		void            reinitializeFId();
		/*!
		Write an .ply file.
		\param output the output .ply file name
		\param ply's fileType
		\\fileType:
		\\ASCII = 1
		\\BINARY_BE = 2(binary PLY file, big endian)
		\\BINARY_LE = 3(binary PLY file, little endian)
		\\PLY_BINARY_NATIVE = 4(binary PLY file, same endianness as current architecture))
		\param comment vector
		\param comment's number
		*/
		bool			write_ply(const char * fileName, int fileType = 1, char** comments = NULL, const int& commentNum = 0);
		/*!
		Write an .g file.
		\param output the output .g file name
		*/
		void			write_g(const char * output);
		/*!
		Read an .off file
		\param input the input .off filename
		*/
		void			read_off(const char * input);
		/*!
		Write an .off file.
		\param output the output .off file name
		*/
		void			write_off(const char * output);

		//number of vertices, faces, edges
		/*! number of vertices */
		int				numVertices();
		/*! number of edges */
		int				numEdges();
		/*! number of faces map */
		int				numFaces();
		/*! number of halfedges */
		int				numHalfEdges();

		MemoryPool<VertexType>	 & vertices() { return mVContainer; };
		MemoryPool<FaceType>	 & faces() { return mFContainer; };
		MemoryPool<EdgeType>	 & edges() { return mEContainer; };
		MemoryPool<HalfEdgeType> & halfedges() { return mHEContainer;; };

		//is boundary
		/*! whether a vertex is on the boundary
		\param v the pointer to the vertex
		*/
		static bool		isBoundary(VPtr	pV);
		/*! whether an edge is on the boundary
		\param e the pointer to the edge
		*/
		static bool		isBoundary(EPtr	pE);
		/*! whether a halfedge is on the boundary
		\param he the pointer to the halfedge
		*/
		static bool		isBoundary(HEPtr pHE);

		//acess vertex - id
		/*!
		Access a vertex by its id
		\param id the vertex id
		\return the vertex, whose ID equals to id. NULL, if there is no such a vertex.
		*/
		//VPtr			idVertex(int id);
		/*!
		The vertex id
		\param v the input vertex
		\return the vertex id.
		*/
		//static int		vertexId(VPtr pV);

		//access face - id
		/*!
		Access a face by its id
		\param id the face id
		\return the face, whose ID equals to id. NULL, if there is no such a face.
		*/
		//FPtr			idFace(int id);
		/*!
		Access a face by its fKey
		\param fKey the face fKey
		\return the face, whose fKey equals to fKey. NULL, if there is no such a face.
		*/
		//FPtr			idFace(KeyVec & fKey);
		/*!
		The face id
		\param f the input face
		\return the face id.
		*/
		//static int		faceId(FPtr pF);

		//access edge - eKey
		/*!
		Access a edge by its eKey
		\param eKey the edge eKey
		\return the edge, whose eKey equals to eKey. NULL, if there is no such a edge.
		*/
		//EPtr			idEdge(KeyVec & eKey);
		/*!
		The edge id
		\param e the input edge
		\return the edge id.
		*/
		//static int		edgeId(EPtr pE);

		//access halfedge - heKey
		/*!
		Access a halfedge by its heKey
		\param heKey the halfedge heKey
		\return the halfedge, whose heKey equals to heKey. NULL, if there is no such a halfedge.
		*/
		//HEPtr			idHalfEdge(KeyVec & heKey);

		//access edge - edge key, vertex
		/*!
		Access an edge by its two end vertices
		\param v0 one vertex of the edge
		\param v1 the other vertex of the edge
		\return the edge connecting both v0 and v1, NULL if no such edge exists.
		*/
		static EPtr		vertexEdge(VPtr pV0, VPtr pV1);

		//access halfedge - halfedge key, vertex
		/*!
		Access a halfedge by its two end vertices
		\param v0 one vertex of the halfedge
		\param v1 the other vertex of the halfedge
		\return the halfedge connecting both v0 and v1, NULL if no such edge exists.
		*/
		static HEPtr	vertexHalfedge(VPtr pV0, VPtr pV1);

		//access corner(halfedge) by attaced face and vertex
		/*!
		Access a halfedge by its target vertex, and attaching face.
		\param v target vertex
		\param f attaching face
		\return halfedge, whose target is v, attaching face is f. NULL if no such an halfedge exists.
		*/
		static HEPtr	corner(VPtr pV, FPtr pF);

		//halfedge->face
		/*!
		The face a halfedge attaching to.
		\param he the input halfedge
		\return the face he attaches
		*/
		static FPtr		halfedgeFace(HEPtr pHE);

		//halfedge->vertex
		/*!
		The target vertex of a halfedge.
		\param he the input halfedge.
		\return the target vertex of he.
		*/
		static VPtr		halfedgeVertex(HEPtr pHE);

		//halfedge->vertex
		/*!
		The target vertex of a halfedge.
		\param he the input halfedge.
		\return the target vertex of he.
		*/
		static VPtr		halfedgeTarget(HEPtr pHE);

		//halfedge->vertex
		/*!
		The source vertex of a halfedge.
		\param he the input halfedge.
		\return the source vertex of he.
		*/
		static VPtr		halfedgeSource(HEPtr pHE);

		//halfedge->next
		/*!
		The next halfedge of a halfedge.
		\param he the input halfedge.
		\return the next halfedge of he.
		*/
		static HEPtr	halfedgeNext(HEPtr pHE);

		//halfedge->prev
		/*!
		The previous halfedge of a halfedge.
		\param he the input halfedge.
		\return the next halfedge of he.
		*/
		static HEPtr	halfedgePrev(HEPtr pHE);

		//halfedge->sym
		/*!
		The dual halfedge of a halfedge.
		\param he the input halfedge.
		\return the dual halfedge of he.
		*/
		static HEPtr	halfedgeSym(HEPtr pHE);

		//halfedge->edge
		/*!
		The edge of a halfedge.
		\param he the input halfedge.
		\return the edge of he.
		*/
		static EPtr		halfedgeEdge(HEPtr pHE);

		//halfedge->vector
		/*!
		Return the vector \in R^3 of a halfedge.
		\param he the input halfedge.
		*/
		template<typename EigenDerived3x1>
		static void halfedgeVec(HEPtr he, EigenDerived3x1& vec);

		//v->halfedge
		/*!
		The halfedge targeting at a vertex.
		\param v the input vertex.
		\return the halfedge targeting at v, which is the most ccw in halfedge of v.
		*/
		static HEPtr	vertexHalfedge(VPtr pV);

		//v->edges
		/*!
		The edge list attaching to the vertex v, such that v is the first vertex of the edge
		\param v the input vertex.
		\return the reference to the edge list
		*/
		//static std::list<EPtr> & vertexEdges(VPtr pV);

		//edge->vertex
		/*!
		The first vertex of an edge.
		\param e the input edge.
		\return the first vertex of e.
		*/
		static VPtr		edgeVertex1(EPtr pE);
		/*!
		The second vertex of an edge.
		\param e the input edge.
		\return the second vertex of e.
		*/
		static VPtr		edgeVertex2(EPtr pE);

		//edge->face
		/*!
		The first face attaching to an edge.
		\param e the input edge.
		\return the first face attaching to e.
		*/
		static FPtr		edgeFace1(EPtr pE);
		/*!
		The second face attaching to an edge.
		\param e the input edge.
		\return the second face attaching to e.
		*/
		static FPtr		edgeFace2(EPtr pE);

		//edge->halfedge
		/*!
		The halfedge attaching to an edge.
		\param e the input edge.
		\param id the index of the halfedge, either 0 or 1
		\return the halfedge[i] attaching to edge e.
		*/
		static HEPtr	edgeHalfedge(EPtr pE, int id);

		//edge->halfedge
		/*!
		The halfedge attaching to an edge.
		\param e the input edge.
		\return the halfedge attaching to edge e.
		*/
		static HEPtr edgeHalfedge(EPtr  e);

		//face->halfedge
		/*!
		The first halfedge attaching to a face f.
		\param f the input face.
		\return the first halfedge attaching to f.
		*/
		static HEPtr	faceHalfedge(FPtr pF);

		//Euler operations
		/*!
		The most Clw Out HalfEdge of a vertex
		\param v the input vertex.
		\return the most Clw Out HalfEdge of v.
		*/
		static HEPtr	vertexMostClwOutHalfEdge(VPtr pV);
		/*!
		The next Ccw Out HalfEdge
		\param he the input halfedge .
		\return the next Ccw Out HalfEdge, sharing the same source of he.
		*/
		static HEPtr	vertexNextCcwOutHalfEdge(HEPtr pHE);
		/*!
		The most Ccw Out HalfEdge of a vertex
		\param v the input vertex.
		\return the most Ccw Out HalfEdge of v.
		*/
		static HEPtr	vertexMostCcwOutHalfEdge(VPtr pV);
		/*!
		The next Clw Out HalfEdge
		\param he the input halfedge .
		\return the next Clw Out HalfEdge, sharing the same source of he.
		*/
		static HEPtr	vertexNextClwOutHalfEdge(HEPtr pHE);
		/*!
		The most Clw In HalfEdge of a vertex
		\param v the input vertex.
		\return the most Clw In HalfEdge of v.
		*/
		static HEPtr	vertexMostClwInHalfEdge(VPtr pV);
		/*!
		The next Ccw In HalfEdge
		\param he the input halfedge .
		\return the next Ccw In HalfEdge, sharing the same target of he.
		*/
		static HEPtr	vertexNextCcwInHalfEdge(HEPtr pHE);
		/*!
		The most Ccw In HalfEdge of a vertex
		\param v the input vertex.
		\return the most Ccw In HalfEdge of v.
		*/
		static HEPtr	vertexMostCcwInHalfEdge(VPtr pV);
		/*!
		The next Clw In HalfEdge
		\param he the input halfedge .
		\return the next Clw In HalfEdge, sharing the same target of he.
		*/
		static HEPtr	vertexNextClwInHalfEdge(HEPtr pHE);
		/*!
		The most Clw HalfEdge of a face
		\param face the input face.
		\return the most Clw HalfEdge of f.
		*/
		static HEPtr	faceMostClwHalfEdge(FPtr pF);
		/*!
		The most Ccw HalfEdge of a face
		\param face the input face.
		\return the most Ccw HalfEdge of f.
		*/
		static HEPtr	faceMostCcwHalfEdge(FPtr pF);
		/*!
		The next Ccw HalfEdge of a halfedge in a face
		\param he the input halfedge.
		\return the next Ccw HalfEdge of he in a face.
		*/
		static HEPtr	faceNextCcwHalfEdge(HEPtr pHE);
		/*!
		The next Clw HalfEdge of a halfedge in a face
		\param he the input halfedge.
		\return the next Clw HalfEdge of he in a face.
		*/
		static HEPtr	faceNextClwHalfEdge(HEPtr pHE);
		/*!
		Edge length
		\param e the input edge
		\return the length of the edge e
		*/
		static double	edgeLength(EPtr pE);

		template<typename EigenDerived3x1>
		static void faceNormal(FPtr pF, EigenDerived3x1& normal);

		template<typename EigenDerived3x1>
		static void vertexNormal(VPtr pV, EigenDerived3x1& normal);

		template<typename EigenDerived3x1>
		static void faceOrientedArea(FPtr pF, EigenDerived3x1& area);

		static DType faceArea(FPtr pF);

		//Function to access maps
		/*!
		Map of the vertices of the mesh.
		*/
		//typename VMap &		getVMap() { return mVMap; };
		/*!
		Map1 of the faces of the mesh.
		*/
		//typename FIdMap &	getFIdMap() { return mFIdMap; };
		/*!
		Map2 of the faces of the mesh.
		*/
		//typename FKeyMap &	getFKeyMap() { return mFKeyMap; };
		/*!
		Map of the edges of the mesh.
		*/
		//typename EMap &		getEMap() { return mEMap; };
		/*!
		Map of the halfedges of the mesh.
		*/
		//typename HEMap &	getHEMap() { return mHEMap; };

		//Function to access container
		/*!
		Container of the vertices of the mesh.
		*/
		VContainer &		getVContainer() { return mVContainer; };
		/*!
		Container of the faces of the mesh.
		*/
		FContainer &	getFContainer() { return mFContainer; };
		/*!
		Container of the edges of the mesh.
		*/
		EContainer &		getEContainer() { return mEContainer; };
		/*!
		Container of the halfedges of the mesh.
		*/
		HEContainer &	getHEContainer() { return mHEContainer; };

	protected:
		//Maps
		/*! Map of vertices */
		VMap				mVMap;
		/*! Map1 of faces */
		FIdMap				mFIdMap;
		/*! Map2 of faces */
		//FKeyMap				mFKeyMap;
		/*! Map of edges */
		//EMap				mEMap;
		/*! Map of halfedge */
		//HEMap				mHEMap;

		//Container
		/*Container of vertices*/
		VContainer			mVContainer;
		/*Container of faces*/
		FContainer			mFContainer;
		/*Container of edges*/
		EContainer			mEContainer;
		/*Container of halfedge*/
		HEContainer			mHEContainer;

		//Element container
		/*! Element container of vertices */
		//StringsContainer	mVStrings;
		/*! Element container of edges */
		//StringsContainer	mEStrings;
		/*! Element container of faces */
		//StringsContainer	mFStrings;
		/*! Element container of halfedge */
		//StringsContainer    mHEStrings;
		
		MAKE_PROP_OF(V);
		MAKE_PROP_OF(E);
		MAKE_PROP_OF(F);
		MAKE_PROP_OF(HE);
	public:
		//Operations of new a object
		/*! New a vertex in mesh */
		VPtr			newVertex();
		/*! New a edge in mesh */
		EPtr			newEdge();
		/*! New a face in mesh */
		FPtr			newFace();
		/*! New a halfedge in mesh */
		HEPtr			newHalfEdge();

		//Operations of delete a object
		/*! Delete a vertex in mesh */
		bool			deleteVertex(VPtr pV);
		/*! Delete a edge in mesh */
		bool			deleteEdge(EPtr pE);
		/*! Delete a face in mesh */
		bool			deleteFace(FPtr pF);
		/*! Delete a halfedge in mesh*/
		bool			deleteHalfEdge(HEPtr pHE);

	public:
		//Operations of modifying mesh's data
		/*! Create a vertex with id, used in the process of reading .m file
		\return pointer to the new vertex
		*/
		VPtr			createVertexWithId(int id);
		//Operations of modifying mesh's data
		/*! Create a vertex with index, used in the process of reading .obj file
		* faster than createVertexWithId
		\return pointer to the new vertex
		*/
		VPtr			createVertexWithIndex();
		/*! Create an edge
		\param v1 end vertex of the edge
		\param v2 end vertex of the edge
		\return pointer to the new edge
		*/
		EPtr			createEdge(VPtr pV0, VPtr pV1);
		/*! Create an halfedge
		\param pSource, source vertex of the halfedge
		\param pTarget, target vertex of the halfedge
		\return pointer to the new halfedge
		*/
		HEPtr			createHalfEdge(VPtr pSource, VPtr pTarget);
		/*! Create a face
		\param v an array of vertices
		\param id face id
		\return pointer to the new face
		*/
		FPtr			createFace(VPtr pVs[], int id);	//create a triangle, param v[n] (n = 3)
		FPtr			createFace(VPtr pVs[]);	//create a triangle, param v[n] (n = 3)

		/*! Create a face
		\param v an array of vertices
		\param id face id
		\return pointer to the new face
		*/
		FPtr			createFace(std::vector<VPtr> & pVs, int id);	//create a triangle, param v[n] (n = 3)
		/*! remove one face
		\param pFace the face to be removed
		*/
		//void			removeFace(FPtr pF);

		//Other operations
		/*! label boundary vertices, edges, faces */
		void			labelBoundary(void);

		//Marks of mesh's character
		/*! whether the vertex is with texture coordinates */
		bool			m_with_texture;
		/*! whether the mesh is with normal */
		bool			m_with_normal;

	public:
		/*!
		*   the input traits of the mesh, there are 64 bits in total
		*/
		static unsigned long long m_input_traits;
		/*!
		*	 the output triats of the mesh, there are 64 bits in total
		*/
		static unsigned long long m_output_traits;
	};

	//access e->vertex
	/*!
	The first vertex of an edge.
	\param e the input edge.
	\return the first vertex of e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType *	CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeVertex1(EPtr pE)
	{
		return (VertexType*)pE->halfedge()->source(); 
	};
	/*!
	The second vertex of an edge.
	\param e the input edge.
	\return the first vertex of e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType *	CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeVertex2(EPtr pE)
	{
		return (VertexType*)pE->halfedge()->target(); 
	};

	//access e->f
	/*!
	The first face attaching to an edge.
	\param e the input edge.
	\return the first face attaching to e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeFace1(EPtr pE)
	{
		assert(pE->halfedge() != NULL);
		return (FaceType*)pE->halfedge()->face();
	};
	/*!
	The second face attaching to an edge.
	\param e the input edge.
	\return the first face attaching to e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeFace2(EPtr pE)
	{
		assert(pE->halfedge() != NULL);
		if (pE->halfedge()->he_sym() == NULL) {
			return NULL;
		}
		else
		{
			return (FaceType*)pE->halfedge()->he_sym()->face();
		}
	};

	//access e->he
	/*!
	The halfedge attaching to an edge.
	\param e the input edge.
	\param id the index of the halfedge, either 0 or 1
	\return the halfedge[i] attaching to edge e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeHalfedge(EPtr pE, int id)
	{
		return (HalfEdgeType*)(id == 0 ? pE->halfedge() : pE->halfedge()->he_sym());
	}

	//edge->halfedge
	/*!
	The halfedge attaching to an edge.
	\param e the input edge.
	\return the halfedge attaching to edge e.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeHalfedge(EPtr pE)
	{
		return (HalfEdgeType*)pE->halfedge();
	};

	//access he->f
	/*!
	The face a halfedge attaching to.
	\param he the input halfedge
	\return the face he attaches
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeFace(HEPtr pHE)
	{
		return (FaceType*)pHE->face();
	};

	//access f->he
	/*!
	The first halfedge attaching to a face f.
	\param f the input face.
	\return the first halfedge attaching to f.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceHalfedge(FPtr pF)
	{
		return (HalfEdgeType*)pF->halfedge();
	};

	//access he->next
	/*!
	The next halfedge of a halfedge.
	\param he the input halfedge.
	\return the next halfedge of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeNext(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->he_next();
	};

	//access he->prev
	/*!
	The previous halfedge of a halfedge.
	\param he the input halfedge.
	\return the next halfedge of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgePrev(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->he_prev();
	};

	//access he->sym
	/*!
	The dual halfedge of a halfedge.
	\param he the input halfedge.
	\return the dual halfedge of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeSym(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->he_sym();
	};

	//access he->edge
	/*!
	The edge of a halfedge.
	\param he the input halfedge.
	\return the edge of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	EdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeEdge(HEPtr pHE)
	{
		return (EdgeType*)pHE->edge();
	}

	//turn the halfedge to a vector 
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	template<typename EigenDerived3x1>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeVec(HEPtr he, EigenDerived3x1& vec)
	{
		vec = he->target()->position() - he->source()->position();
	};

	//access he->v
	/*!
	The target vertex of a halfedge.
	\param he the input halfedge.
	\return the target vertex of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeVertex(HEPtr pHE)
	{
		return (VertexType*)pHE->vertex();
	};

	//access he->v
	/*!
	The target vertex of a halfedge.
	\param he the input halfedge.
	\return the target vertex of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeTarget(HEPtr pHE)
	{
		return (VertexType*)pHE->vertex();
	};

	//access he->v
	/*!
	The source vertex of a halfedge.
	\param he the input halfedge.
	\return the source vertex of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::halfedgeSource(HEPtr pHE)
	{
		return (VertexType*)pHE->he_prev()->vertex();
	};

	//access (v,f)->corner
	/*!
	Access a halfedge by its target vertex, and attaching face.
	\param v target vertex
	\param f attaching face
	\return halfedge, whose target is v, attaching face is f. NULL if no such an halfedge exists.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::corner(VPtr pV, FPtr pF)
	{
		HalfEdgeType * pHE = faceMostCcwHalfEdge(pF);
		do {
			if (pHE->vertex() == pV)
				return (HalfEdgeType*)pHE;
			pHE = faceNextCcwHalfEdge(pHE);
		} while (pHE != faceMostCcwHalfEdge(pF));
		return NULL;
	};

	//Judge if on boundary
	/*! whether a vertex is on the boundary
	\param v the pointer to the vertex
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::isBoundary(VPtr pV)
	{
		return pV->boundary();
	};
	/*! whether an edge is on the boundary
	\param e the pointer to the edge
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::isBoundary(EPtr pE)
	{
		return (pE->halfedge()->he_sym() == NULL);
	};
	/*! whether a halfedge is on the boundary
	\param he the pointer to the halfedge
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::isBoundary(HEPtr pHE)
	{
		return (pHE->he_sym() == NULL);
	};

	//Function to get number
	/*! Number of vertices of the mesh
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::numVertices()
	{
		return (int)this->mVContainer.size();
	};
	/*! Number of edges of the mesh
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::numEdges()
	{
		return (int)this->mEContainer.size();
	};
	/*! Number of faces map1 of the mesh
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::numFaces()
	{
		return (int)this->mFContainer.size();
	};
	/*! Number of halfedge of the mesh
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::numHalfEdges()
	{
		return (int)this->mHEContainer.size();
	};

	//Euler operation
	/*!
	The most Clw Out HalfEdge of a vertex
	\param v the input vertex.
	\return the most Clw Out HalfEdge of v.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexMostClwOutHalfEdge(VPtr pV)
	{
		return (HalfEdgeType*)pV->most_clw_out_halfedge();
	};
	/*!
	The most Ccw Out HalfEdge of a vertex
	\param v the input vertex.
	\return the most Ccw Out HalfEdge of v.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexMostCcwOutHalfEdge(VPtr pV)
	{
		return (HalfEdgeType*)pV->most_ccw_out_halfedge();
	};
	/*!
	The next Ccw Out HalfEdge
	\param he the input halfedge .
	\return the next Ccw Out HalfEdge, sharing the same source of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexNextCcwOutHalfEdge(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->ccw_rotate_about_source();
	};
	/*!
	The next Clw Out HalfEdge
	\param he the input halfedge .
	\return the next Clw Out HalfEdge, sharing the same source of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexNextClwOutHalfEdge(HEPtr pHE)
	{
		if (pHE->he_sym() != NULL)
			return NULL;
		return (HalfEdgeType*)pHE->clw_rotate_about_source();
	};
	/*!
	The most Clw In HalfEdge of a vertex
	\param v the input vertex.
	\return the most Clw In HalfEdge of v.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexMostClwInHalfEdge(VPtr pV)
	{
		return (HalfEdgeType*)pV->most_clw_in_halfedge();
	};
	/*!
	The most Ccw In HalfEdge of a vertex
	\param v the input vertex.
	\return the most Clw In HalfEdge of v.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexMostCcwInHalfEdge(VPtr pV)
	{
		return (HalfEdgeType*)pV->most_ccw_in_halfedge();
	};
	/*!
	The next Ccw In HalfEdge
	\param he the input halfedge .
	\return the next Ccw In HalfEdge, sharing the same target of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexNextCcwInHalfEdge(HEPtr pHE)
	{
		assert(pHE->he_sym() != NULL);
		return (HalfEdgeType*)pHE->ccw_rotate_about_target();
	};
	/*!
	The next Clw In HalfEdge
	\param he the input halfedge .
	\return the next Clw In HalfEdge, sharing the same target of he.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexNextClwInHalfEdge(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->clw_rotate_about_target();
	};
	/*!
	The next Clw HalfEdge of a halfedge in a face
	\param he the input halfedge.
	\return the next Clw HalfEdge of he in a face.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceNextClwHalfEdge(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->he_prev();
	};
	/*!
	The next Ccw HalfEdge of a halfedge in a face
	\param he the input halfedge.
	\return the next Ccw HalfEdge of he in a face.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceNextCcwHalfEdge(HEPtr pHE)
	{
		return (HalfEdgeType*)pHE->he_next();
	};
	/*!
	The most Ccw HalfEdge of a face
	\param face the input face.
	\return the most Ccw HalfEdge of f.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceMostCcwHalfEdge(FPtr pF)
	{
		return (HalfEdgeType*)pF->halfedge();
	};
	/*!
	The most Clw HalfEdge in a face
	\param face the input face.
	\return the most Clw HalfEdge in a face.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceMostClwHalfEdge(FPtr pF)
	{
		return (HalfEdgeType*)pF->halfedge()->he_next();
	};

	/*!
	CMeshBase destructor
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::~CMeshBase()
	{
		if (mVProps.size() != 0) {
			for (int i = 0; i < mVProps.size(); ++i)
			{
				BasicPropHandle * pVProp = mVProps[i];
				pVProp->destructProp();
			}
		}

		if (mFProps.size() != 0) {
			for (int i = 0; i < mFProps.size(); ++i)
			{
				BasicPropHandle * pFProp = mFProps[i];
				pFProp->destructProp();
			}
		}

		if (mEProps.size() != 0) {
			for (int i = 0; i < mEProps.size(); ++i)
			{
				BasicPropHandle * pEProp = mEProps[i];
				pEProp->destructProp();
			}
		}

		if (mHEProps.size() != 0) {
			for (int i = 0; i < mHEProps.size(); ++i)
			{
				BasicPropHandle * pHEProp = mHEProps[i];
				pHEProp->destructProp();
			}
		}
	};

	/*!
	Edge length
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	double CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeLength(EPtr pE)
	{
		VertexType * pV1 = edgeVertex1(pE);
		VertexType * pV2 = edgeVertex2(pE);
		return (pV1->position() - pV2->position()).norm();
	}
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	template<typename EigenDerived3x1>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceNormal(FPtr pF, EigenDerived3x1 & fNormal)
	{
		fNormal = faceOrientedArea(pF);
	}

	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	template<typename EigenDerived3x1>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexNormal(VPtr pV, EigenDerived3x1& vNormal)
	{
		vNormal = Vec3::Zero();
		for (CHalfEdge * pHE : pV->outHEs())
		{
			FaceType * pF = (FaceType *)pHE->face();
			Vec3 fn;
			faceOrientedArea(pF, fn);
			vNormal += fn;
		}
		vNormal /= vNormal.norm();
	}

	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	template<typename EigenDerived3x1>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceOrientedArea(FPtr pF, EigenDerived3x1 & area)
	{
		HalfEdgeType * pHE = faceHalfedge(pF);
		HalfEdgeType * pHENext = halfedgeNext(pHE);
		Vec3 v1, v2;
		halfedgeVec(pHE, v1);
		v1 = -v1;
		halfedgeVec(pHENext, v2);
		area = 0.5 * v2.cross(v1);
	}
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline DType CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceArea(FPtr pF)
	{
		return faceOrientedArea(pF).norm();
	}
	;

	//access id->v
	/*!
	Access a vertex by its id
	\param id the vertex id
	\return the vertex, whose ID equals to id. NULL, if there is no such a vertex.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::idVertex(int id)
	//{
	//	auto pV = mVMap.find(id);
	//	if (pV != mVMap.end()) { return pV->second; }
	//	else { return NULL; }
	//};

	//access v->id
	/*!
	The vertex id
	\param v the input vertex
	\return the vertex id.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexId(VPtr pV)
	//{
	//	return pV->id();
	//};

	//access id->f
	/*!
	Access a face by its id
	\param id the face id
	\return the face, whose ID equals to id. NULL, if there is no such a face.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::idFace(int id)
	//{
	//	auto pF = mFIdMap.find(id);
	//	if (pF != mFIdMap.end()) { return pF->second; }
	//	else { return NULL; }
	//};

	//access fKey->f
	/*!
	Access a face by its fKey
	\param fKey the face fKey
	\return the face, whose fKey equals to fKey. NULL, if there is no such a face.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::idFace(KeyVec & fKey)
	//{
	//	auto pF = mFKeyMap.find(fKey);
	//	if (pF != mFKeyMap.end()) { return pF->second; }
	//	else { return NULL; }
	//};

	//acess f->id
	/*!
	The face id
	\param f the input face
	\return the face id.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::faceId(FPtr pF)
	//{
	//	return pF->id();
	//};

	//access eKey->e
	/*!
	Access a edge by its eKey
	\param eKey the edge eKey
	\return the edge, whose eKey equals to eKey. NULL, if there is no such a edge.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//EdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::idEdge(KeyVec & eKey)
	//{
	//	auto pE = mEMap.find(eKey);
	//	if (pE != mEMap.end()) { return pE->second; }
	//	else { return NULL; }
	//};

	//access e->id
	/*!
	The edge id
	\param e the input edge
	\return the edge id.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//int CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::edgeId(EPtr pE)
	//{
	//	return pE->id();
	//};

	//access heKey->he
	/*!
	Access a halfedge by its heKey
	\param heKey the halfedge heKey
	\return the halfedge, whose heKey equals to heKey. NULL, if there is no such a halfedge.
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::idHalfEdge(KeyVec & heKey)
	//{
	//	auto pHE = mHEMap.find(heKey);
	//	if (pHE != mHEMap.end()) { return pHE->second; }
	//	else { return NULL; }
	//};

	//access vertex->edge
	/*!
	Access an edge by its two end vertices
	\param v0 one vertex of the edge
	\param v1 the other vertex of the edge
	\return the edge connecting both v0 and v1, NULL if no such edge exists.
	*/
	//use the edge list associated with each vertex to locate the edge
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	EdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexEdge(VPtr pV0, VPtr pV1)
	{
		//VertexType * pV = (v0->id() < v1->id() )? v0: v1;
		//CPArray<HEPtr, 8> & lhes0 = (CPArray<HEPtr, 8>&)pV0->outHEs();//vertexEdges(v0);
		//CPArray<HEPtr, 8> & lhes1 = (CPArray<HEPtr, 8>&)pV1->outHEs();//vertexEdges(v1);

		for (CHalfEdge* pH : pV0->outHEs())
		{
			if (pH->source() == pV0 && pH->target() == pV1) return (EPtr)pH->edge();
		}
		for (CHalfEdge* pH : pV1->outHEs())
		{
			if (pH->source() == pV0 && pH->target() == pV1) return  (EPtr)pH->edge();
		}
		return NULL;
	};

	//access vertex->halfedge
	/*!
	Access a halfedge by its two end vertices
	\param v0 one vertex of the halfedge
	\param v1 the other vertex of the halfedge
	\return the halfedge connecting both v0 and v1, NULL if no such edge exists.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexHalfedge(VPtr pV0, VPtr pV1)
	{
		//std::vector<HEPtr> & lhes0 = (std::vector<HEPtr>&)();
		for (CHalfEdge* pH : pV0->outHEs())
		{
			if (pH->source() == pV0 && pH->target() == pV1) return (HEPtr)pH;
		}
		return NULL;
	};

	//access vertex->edges
	/*!
	Access the edge list of a vertex, {e} such that e->vertex1() == v
	\param v vertex
	\return the list of adjacent edges
	*/
	//template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	//std::list<EdgeType*> & CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexEdges(VPtr pV)
	//{
	//	return (std::list<EdgeType*> &)pV->edges();
	//}

	//access vertex->halfedge
	/*!
	The halfedge targeting at a vertex.
	\param v the input vertex.
	\return the halfedge targeting at v, which is the most ccw in halfedge of v.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::vertexHalfedge(VPtr pV)
	{
		return (HalfEdgeType*)pV->halfedge();
	}

	//Operations of new a object
	/*! New a vertex in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::newVertex()
	{
		size_t index;
		VertexType * pV = mVContainer.newMember(index);
		assert(pV != NULL);

		for (int i = 0; i < mVProps.size(); i++) {
			BasicPropHandle * pProp = mVProps[i];
			pProp->initializePropMember(index);
		}
		pV->index() = index;
		return pV;
	}
	/*! New a edge in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline EdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::newEdge()
	{
		size_t index;
		EdgeType * pE = mEContainer.newMember(index);
		assert(pE != NULL);
		pE->index() = index;
		pE->id() = index;
		return pE;
	}
	/*! New a face in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::newFace()
	{
		size_t index;
		FaceType * pF = mFContainer.newMember(index);
		assert(pF != NULL);
		pF->index() = index;
		return pF;
	}
	/*! New a halfedge in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::newHalfEdge()
	{
		size_t index;
		HalfEdgeType * pHE = mHEContainer.newMember(index);
		assert(pHE != NULL);
		pHE->index() = index;
		return pHE;
	}

	//Operations of delete a object
	/*! Delete a vertex in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::deleteVertex(VPtr pV)
	{
		//remove current vertex from pool
		return mVContainer.deleteMember(pV->index());
	}
	/*! Delete a edge in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::deleteEdge(EPtr pE)
	{
		//remove current edge from pool
		return mEContainer.deleteMember(pE->index());
	}
	/*! Delete a face in mesh */
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::deleteFace(FPtr pF)
	{
		//remove current face from pool
		return mFContainer.deleteMember(pF->index());
	}
	/*! Delete a halfedge in mesh*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline bool CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::deleteHalfEdge(HEPtr pHE)
	{
		//remove current halfedge from pool
		return mHEContainer.deleteMember(pHE->index());
	}

	//create new gemetric simplexes
	/*! Create a vertex
	\return pointer to the new vertex
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createVertexWithId(int id)
	{
		VertexType * pV;
		pV = newVertex();
		mVMap.insert(VMapPair(id, pV));
		pV->id() = id;
		return pV;

	}
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline VertexType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createVertexWithIndex()
	{
		VertexType * pV;
		pV = newVertex();
		pV->id() = (int)pV->index();
		mVMap.insert(VMapPair(pV->id(), pV));
		return pV;
	}
	/*! Create an edge
	\param v1 end vertex of the edge
	\param v2 end vertex of the edge
	\return pointer to the new edge
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	EdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createEdge(VPtr pV0, VPtr pV1)
	{
		VertexType * pV = (pV0->id() < pV1->id()) ? pV0 : pV1;
		std::list<EdgeType*> & ledges = vertexEdges(pV);
		for (typename std::list<EdgeType*>::iterator eiter = ledges.begin(); eiter != ledges.end(); eiter++)
		{
			EdgeType * pE = *eiter;
			HalfEdgeType * pHE = edgeHalfedge(pE, 0);
			if (pHE->source() == pV0 && pHE->target() == pV1) return pE;
			if (pHE->source() == pV1 && pHE->target() == pV0) return pE;
		}
		EdgeType * pE = newEdge();
		pE->id() = (int)pE->index();
		ledges.push_back(pE);
		//int idV0 = (pV0->id() < pV1->id()) ? pV0->id() : pV1->id();
		//int idV1 = (pV0->id() > pV1->id()) ? pV0->id() : pV1->id();
		//KeyVec eKey = { idV0, idV1 };
		//mEMap.insert(EMapPair(eKey, pE));
		return pE;
	}
	/*! Create an halfedge
	\param pSource, source vertex of the halfedge
	\param pTarget, target vertex of the halfedge
	\return pointer to the new halfedge
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	HalfEdgeType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createHalfEdge(VPtr pSource, VPtr pTarget)
	{
		HalfEdgeType * pHE;
		pHE = newHalfEdge();
		//KeyVec heKey = { pSource->id(), pTarget->id() };
		//mHEMap.insert(HEMapPair(heKey, pHE));
		
		pTarget->halfedge() = pHE;
		pHE->vertex() = pTarget;
		pSource->outHEs().push_back(pHE);

		return pHE;
	}
	/*! Create a face
	\param v an array of vertices
	\param id face id
	\return pointer to the new face
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createFace(VPtr pVs[], int id)
	{
		FaceType * pF;
		pF = newFace();
		//mFIdMap.insert(FIdMapPair(id, pF));
		pF->id() = id;
		//create halfedges
		HalfEdgeType * pHEs[3];
		for (int i = 0; i < 3; i++)
		{
			VertexType * pTarget = pVs[i];
			VertexType * pSource = pVs[(i + 2) % 3];
			pHEs[i] = createHalfEdge(pSource, pTarget);
			pHEs[i]->face() = pF;
			//Search for the symmetrical halfedge
			HalfEdgeType * pHESym = NULL;
			for (int i = 0; i < pTarget->outHEs().size(); ++i)
			{
				HalfEdgeType * pHE = (HalfEdgeType *)pTarget->outHEs()[i];
				if ((VertexType *)pHE->vertex() == pSource) {
					pHESym = pHE;
					break;
				}
			}
			// if symmetrical halfedge already exist, link current halfedge to it 
			if (pHESym) {
				//assertx(!hes->he_sym());
				pHESym->he_sym() = pHEs[i];
				pHEs[i]->he_sym() = pHESym;
				pHEs[i]->edge() = (EdgeType *)pHESym->edge();
			}
			// if symmetrical halfedge does not exist, create a new edge and link to it
			else {                                          
				EdgeType * pE = newEdge();
				pE->halfedge() = pHEs[i];
				pHEs[i]->edge() = pE;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			pHEs[i]->he_next() = pHEs[(i + 1) % 3];
			pHEs[i]->he_prev() = pHEs[(i + 2) % 3];
		}

		pF->halfedge() = pHEs[0];
		//connecting with edge
		return pF;
	}
	/*! Create a face
	\param v an array of vertices
	\param id face id
	\return pointer to the new face
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	FaceType * CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::createFace(VPtr pVs[])
	{
		FaceType * pF;
		pF = newFace();
		//mFIdMap.insert(FIdMapPair(id, pF));
		pF->id() = (int)pF->index();
		//create halfedges
		HalfEdgeType * pHEs[3];
		for (int i = 0; i < 3; i++)
		{
			VertexType * pTarget = pVs[i];
			VertexType * pSource = pVs[(i + 2) % 3];
			pHEs[i] = createHalfEdge(pSource, pTarget);
			pHEs[i]->face() = pF;
			//Search for the symmetrical halfedge
			HalfEdgeType * pHESym = NULL;
			for (int i = 0; i < pTarget->outHEs().size(); ++i)
			{
				HalfEdgeType * pHE = (HalfEdgeType *)pTarget->outHEs()[i];
				if ((VertexType *)pHE->vertex() == pSource) {
					pHESym = pHE;
					break;
				}
			}
			// if symmetrical halfedge already exist, link current halfedge to it 
			if (pHESym) {
				//assertx(!hes->he_sym());
				pHESym->he_sym() = pHEs[i];
				pHEs[i]->he_sym() = pHESym;
				pHEs[i]->edge() = (EdgeType *)pHESym->edge();
			}
			// if symmetrical halfedge does not exist, create a new edge and link to it
			else {
				EdgeType * pE = newEdge();
				pE->halfedge() = pHEs[i];
				pHEs[i]->edge() = pE;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			pHEs[i]->he_next() = pHEs[(i + 1) % 3];
			pHEs[i]->he_prev() = pHEs[(i + 2) % 3];
		}

		pF->halfedge() = pHEs[0];
		//connecting with edge
		return pF;
	}

	//Boundary
	/*!
	Label boundary edges, vertices
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::labelBoundary(void)
	{
		/*Label boundary edges*/
		for (auto edge : mEContainer)
		{
			EdgeType* currentE = edge;
			HalfEdgeType* currentHE[2];
			currentHE[0] = (HalfEdgeType*)currentE->halfedge(0);
			currentHE[1] = (HalfEdgeType*)currentE->halfedge(1);
			assert(currentHE[0] != NULL);
			if (currentHE[1] != NULL)
			{
				assert(currentHE[0]->target() == currentHE[1]->source() && currentHE[0]->source() == currentHE[1]->target());
				if (currentHE[0]->target()->id() < currentHE[0]->source()->id())
				{
					currentE->halfedge(0) = currentHE[1];
					currentE->halfedge(1) = currentHE[0];
				}
				assert(edgeVertex1(currentE)->id() < edgeVertex2(currentE)->id());
			}
			else
			{
				currentHE[0]->vertex()->boundary() = true;
				currentHE[0]->he_prev()->vertex()->boundary() = true;
			}
		}
		/*Remove isolated vertex*/
		std::vector<VertexType*> isolatedVertexs;
		for (auto pV : mVContainer)
		{
			if (pV->halfedge() != NULL) continue;
			isolatedVertexs.push_back(pV);
		}
		for (auto pV : isolatedVertexs)
		{
			mVContainer.deleteMember(pV->index());
		}
		/*
		*	Arrange the boundary half_edge of boundary vertices, to make its halfedge
		*	to be the most ccw in half_edge
		*/
		for (auto pV : mVContainer)
		{
			
			if (!pV->boundary()) continue;
			HalfEdgeType* currentHE = (HalfEdgeType*)pV->halfedge();
			while (currentHE->he_sym() != NULL)
			{
				currentHE = (HalfEdgeType*)currentHE->ccw_rotate_about_target();
			}
			pV->halfedge() = currentHE;
		}
	}

	//file io
	/*
		Write .ply file
	*/
	/*!
	\ Put element info into .ply file.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::mesh_ply_put_element(PlyFile* plyFile, void* voidPtr, MeshType meshType, PlyFileReader* plyFileReader)
	{
		FILE* fp = plyFile->fp;
		int fileType = plyFile->file_type;
		VertexType* VTypePtr = NULL;
		EdgeType* ETypePtr = NULL;
		FaceType* FTypePtr = NULL;
		HalfEdgeType* HETypePtr = NULL;
		int intValue = 0;
		unsigned int uintValue = 0;
		double doubleValue = 0;
		if (meshType == MeshType::PLY_V_Type)
			VTypePtr = (VertexType*)voidPtr;
		else if (meshType == MeshType::PLY_F_Type)
			FTypePtr = (FaceType*)voidPtr;
		else if (meshType == MeshType::PLY_E_Type)
			ETypePtr = (EdgeType*)voidPtr;
		else if (meshType == MeshType::PLY_HE_Type)
			HETypePtr = (HalfEdgeType*)voidPtr;
		PlyElement* currentElem = plyFile->which_elem;
		/* Write out either to an ascii or binary file */
		/* Write out each property of the element */
		for (int i = 0; i < currentElem->propNum; i++)
		{
			PlyProperty* currentProp = currentElem->propList[i];
			int nameMark = currentProp->nameMark;
			if (!currentProp->isList)		/*Not a list*/
			{
				if (meshType == MeshType::PLY_V_Type)
				{
					if (nameMark == PLY_NAME_MARK::PLY_X)
						doubleValue = VTypePtr->position()[0];
					else if (nameMark == PLY_NAME_MARK::PLY_Y)
						doubleValue = VTypePtr->position()[1];
					else if (nameMark == PLY_NAME_MARK::PLY_Z)
						doubleValue = VTypePtr->position()[2];
					//else if (nameMark == PLY_NAME_MARK::PLY_RED)
					//	doubleValue = VTypePtr->color().r;
					//else if (nameMark == PLY_NAME_MARK::PLY_GREEN)
					//	doubleValue = VTypePtr->color().g;
					//else if (nameMark == PLY_NAME_MARK::PLY_BLUE)
					//	doubleValue = VTypePtr->color().b;
					//else if (nameMark == PLY_NAME_MARK::PLY_NX)
					//	doubleValue = VTypePtr->normal()[0];
					//else if (nameMark == PLY_NAME_MARK::PLY_NY)
					//	doubleValue = VTypePtr->normal()[1];
					//else if (nameMark == PLY_NAME_MARK::PLY_NZ)
					//	doubleValue = VTypePtr->normal()[2];
					//else if (nameMark == PLY_NAME_MARK::PLY_U)
					//	doubleValue = VTypePtr->uv()[0];
					//else if (nameMark == PLY_NAME_MARK::PLY_V)
					//	doubleValue = VTypePtr->uv()[1];
					else
						;
					//Some operation......
				}
				if (fileType == PLY_FILE_TYPE::PLY_ASCII)
				{
					plyFileReader->write_ascii_item(fp, intValue, uintValue, doubleValue, currentProp->externalType);
				}
					
				else
					plyFileReader->write_binary_item(fp, fileType, intValue, uintValue, doubleValue, currentProp->externalType);
			}
			else	/*Is a list*/
			{
				if (meshType == MeshType::PLY_F_Type)
				{
					int vertexId[3];
					if (fileType == PLY_FILE_TYPE::PLY_ASCII)
						plyFileReader->write_ascii_item(fp, intValue, 3, doubleValue, currentProp->countType);
					else
						plyFileReader->write_binary_item(fp, fileType, intValue, 3, doubleValue, currentProp->countType);
					HalfEdgeType* currentHE = FTypePtr->halfedge();
					vertexId[0] = currentHE->vertex()->id();
					currentHE = currentHE->he_next();
					vertexId[1] = currentHE->vertex()->id();
					currentHE = currentHE->he_next();
					vertexId[2] = currentHE->vertex()->id();
					for (int j = 0; j < 3; j++) {
						if (fileType == PLY_FILE_TYPE::PLY_ASCII)
							plyFileReader->write_ascii_item(fp, vertexId[j], uintValue, doubleValue, currentProp->externalType);
						else
							plyFileReader->write_binary_item(fp, fileType, vertexId[j], uintValue, doubleValue, currentProp->externalType);
					}
				}
				else
				{
					//Some other operation......
				}
			}
		}
		if (fileType == PLY_FILE_TYPE::PLY_ASCII)
			fprintf(fp, "\n");
	}
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::reinitializeVId()
	{	
		mVMap.clear();

		int currentId = 0;
		for (auto pV : mVContainer)
		{
			pV->id() = currentId;
			mVMap.insert(VMapPair(currentId, pV));
			++currentId;

		}
	}
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::reinitializeFId()
	{
		mFIdMap.clear();
		int currentId = 0;
		for (auto pF : mFContainer)
		{
			pF->id() = currentId;
			mFIdMap.insert(FIdMapPair(currentId, pF));
			++currentId;
		}
	}
	/*!
		Write an .ply file.
		\param input the input .ply file name and fileType(ASCII = 1	BINARY_BE = 2	BINARY_LE = 3	PLY_BINARY_NATIVE = 4) 
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline bool	CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::write_ply(const char * fileName, int fileType, char** comments, const int& commentNum)
	{
		/*VertexType vExample;
		bool hasColor = (vExample.hasColor()) ? true : false;
		bool hasNormal = (vExample.hasNormal()) ? true : false;
		bool hasUV = (vExample.hasUV()) ? true : false;*/
		int numVertexs = mVContainer.size();
		int numFaces = mFContainer.size();
		const char *elem_names[] = { "vertex" , "face" };
		float version;
		PlyFileReader plyFileReader(false, false, false);    //used to write .ply file
		PlyFile* plyFile = plyFileReader.ply_open_for_writing(fileName, 2, elem_names, fileType, &version);
		if (!plyFile) {
			printf("Can't create plyFile object!\n");
			return false;
		}
		int vId = 0;
		/*Reset all vertexs' id*/
		for (int i = 0; i < mVContainer.getCurrentIndex(); i++)
		{
			VertexType* currentVertex = mVContainer.getPointer(i);
			if (mVContainer.hasBeenDeleted(currentVertex->index()) == false)
			{
				currentVertex->id() = vId;
				vId++;
			}
		}
		plyFileReader.ply_element_count(plyFile, "vertex", numVertexs);
		for (int i = 0; i < 11; i++)
		{
			if (i >= 0 && i <= 2)
				plyFileReader.ply_describe_property(plyFile, "vertex", i);
	/*		else if (i >= 3 && i <= 5 && vExample.hasColor())
				plyFileReader.ply_describe_property(plyFile, "vertex", i);
			else if (i >= 6 && i <= 8 && vExample.hasNormal())
				plyFileReader.ply_describe_property(plyFile, "vertex", i);
			else if (i >= 9 && i <= 10 && vExample.hasUV())
				plyFileReader.ply_describe_property(plyFile, "vertex", i);*/
		}
		plyFileReader.ply_element_count(plyFile, "face", numFaces);
		plyFileReader.ply_describe_property(plyFile, "face", 11);
		/*Write comments*/
		if (comments && commentNum)
			for (int i = 0; i < commentNum; i++)
				plyFileReader.ply_put_comment(plyFile, comments[i]);
		plyFileReader.ply_header_complete(plyFile);
		/*Write vertexs' info*/
		plyFileReader.ply_put_element_setup(plyFile, "vertex");
		for (int i = 0; i < mVContainer.getCurrentIndex(); i++)
		{
			VertexType* currentVertex = mVContainer.getPointer(i);
			if (mVContainer.hasBeenDeleted(currentVertex->index()) == false)
			{
				mesh_ply_put_element(plyFile, (void *)currentVertex, MeshType::PLY_V_Type, &plyFileReader);
			}
		}
		/*Write faces' info*/
		plyFileReader.ply_put_element_setup(plyFile, "face");
		for (int i = 0; i < mFContainer.getCurrentIndex(); i++)
		{
			FaceType* currentFace = mFContainer.getPointer(i);
			if (mFContainer.hasBeenDeleted(currentFace->index()) == false)
			{
				mesh_ply_put_element(plyFile, (void*)currentFace, MeshType::PLY_F_Type, &plyFileReader);
			}
		}
		/*Free ply object's memory*/
		plyFileReader.free_ply_memory(plyFile, false);
		return true;
	}

	/*
		Read .ply file
	*/
	/*!
	\ Get element info from .ply file.
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void* CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::
		mesh_ply_get_element(PlyFile* plyFile, PlyElement* currentEle, PlyOtherElem* currentOtherEle,MeshType meshType, PlyFileReader* plyFileReader)
	{
		FILE* fp = plyFile->fp;
		int fileType = plyFile->file_type;
		int numWords;
		char *origLine;
		int currentWord;
		int intValue;
		unsigned int uintValue;
		double doubleValue;
		/*Create a void pointer*/
		void* voidPtr = NULL;
		VertexType* VTypePtr = NULL;
		EdgeType* ETypePtr = NULL;
		FaceType* FTypePtr = NULL;
		HalfEdgeType* HETypePtr = NULL;
		VertexType* currentVs[3];
		/*Tokenizer*/
		strutil::Tokenizer stokenizer;
		const char * token;
		char lineBuffer[MAX_LINE_SIZE];
		if (meshType == MeshType::PLY_V_Type) 
		{
			VTypePtr = createVertexWithIndex();
			voidPtr = VTypePtr;
		}
		else if (meshType == MeshType::PLY_E_Type)
		{
			voidPtr = ETypePtr;
			//Some operation......
		}
		else if (meshType == MeshType::PLY_HE_Type)
		{
			voidPtr = HETypePtr;
			//Some operation......
		}
		if (fileType == PLY_FILE_TYPE::PLY_ASCII)
		{
			/*If get nothing, free buffer, return!*/
			if (!fgets(lineBuffer, MAX_LINE_SIZE, plyFile->fp))
			{
				return NULL;
			}
			stokenizer.initToken(lineBuffer, " ");
		}
		/*Get Properties number in current element*/
		int numProp = currentEle->propNum;
		/*Reset plyfile's propindex*/
		int propIndex = 0;
		/*read in a number of elements*/
		for (int i = 0; i < numProp; i++)
		{
			PlyProperty* currentProp = currentEle->propList[i];
			char* propName = currentProp->propName;
			int nameMark = currentProp->nameMark;
			if (!currentProp->isList)		/* Not a list */
			{
				if (fileType == PLY_FILE_TYPE::PLY_ASCII)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					plyFileReader->get_ascii_item(token, currentProp->externalType, &intValue, &uintValue, &doubleValue);
				}
				else
					plyFileReader->get_binary_item(fp, fileType, currentProp->externalType, &intValue, &uintValue, &doubleValue);
				if (meshType == MeshType::PLY_V_Type)
				{
					if (nameMark == PLY_NAME_MARK::PLY_X)
						VTypePtr->position()[0] = doubleValue;
					else if (nameMark == PLY_NAME_MARK::PLY_Y)
						VTypePtr->position()[1] = doubleValue;
					else if (nameMark == PLY_NAME_MARK::PLY_Z)
						VTypePtr->position()[2] = doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_RED)
					//	VTypePtr->color().r = (float)doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_GREEN)
					//	VTypePtr->color().g = (float)doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_BLUE)
					//	VTypePtr->color().b = (float)doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_NX)
					//	VTypePtr->normal()[0] = doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_NY)
					//	VTypePtr->normal()[1] = doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_NZ)
					//	VTypePtr->normal()[2] = doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_U)
					//	VTypePtr->uv()[0] = doubleValue;
					//else if (nameMark == PLY_NAME_MARK::PLY_V)
					//	VTypePtr->uv()[1] = doubleValue;
					else {
						PlyOtherProperty* currentOtherProp = currentOtherEle->propList[propIndex];
						propIndex++;
						if (!currentOtherProp->ifMalloc) {
							plyFileReader->malloc_other_property(currentOtherEle, currentOtherProp, currentOtherProp->isList);
						}
						plyFileReader->store_other_item(currentOtherProp, intValue, uintValue, doubleValue);
					}
				}
				else
				{
					PlyOtherProperty* currentOtherProp = currentOtherEle->propList[propIndex];
					propIndex++;
					if (!currentOtherProp->ifMalloc) {
						plyFileReader->malloc_other_property(currentOtherEle, currentOtherProp, currentOtherProp->isList);
					}
					plyFileReader->store_other_item(currentOtherProp, intValue, uintValue, doubleValue);
				}
			}
			else	/* Is a list */
			{
				/* get and store the number of items in the list */
				if (fileType == PLY_FILE_TYPE::PLY_ASCII)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					plyFileReader->get_ascii_item(token, currentProp->countType, &intValue, &uintValue, &doubleValue);
				}
				else
					plyFileReader->get_binary_item(fp, fileType, currentProp->countType, &intValue, &uintValue, &doubleValue);
				int listCount = intValue;
				if (meshType == MeshType::PLY_F_Type)
				{
					if (nameMark == PLY_NAME_MARK::PLY_VERTEXS)
					{
						for (int j = 0; j < 3; j++)
						{
							if (fileType == PLY_FILE_TYPE::PLY_ASCII)
							{
								stokenizer.nextToken();
								token = stokenizer.getToken();
								plyFileReader->get_ascii_item(token, currentProp->externalType, &intValue, &uintValue, &doubleValue);
							}
							else
								plyFileReader->get_binary_item(fp, fileType, currentProp->externalType, &intValue, &uintValue, &doubleValue);
							int vId = intValue;
							currentVs[j] = mVContainer.getPointer(vId);
						}
					}
					else
					{
						PlyOtherProperty* currentOtherProp = currentOtherEle->propList[propIndex];
						propIndex++;
						if (!currentOtherProp->ifMalloc) {
							currentOtherProp->listCount = listCount;
							plyFileReader->malloc_other_property(currentOtherEle, currentOtherProp, currentOtherProp->isList);
						}
						for (int j = 0; j < listCount; j++)
						{
							if (fileType == PLY_FILE_TYPE::PLY_ASCII)
							{
								stokenizer.nextToken();
								token = stokenizer.getToken();
								plyFileReader->get_ascii_item(token, currentProp->externalType, &intValue, &uintValue, &doubleValue);
							}
							else
								plyFileReader->get_binary_item(fp, fileType, currentProp->externalType, &intValue, &uintValue, &doubleValue);
							plyFileReader->store_other_item(currentOtherProp, intValue, uintValue, doubleValue);
						}
					}
				}
				else
				{
					PlyOtherProperty* currentOtherProp = currentOtherEle->propList[propIndex];
					propIndex++;
					if (!currentOtherProp->ifMalloc) {
						currentOtherProp->listCount = listCount;
						plyFileReader->malloc_other_property(currentOtherEle, currentOtherProp, currentOtherProp->isList);
					}
					for (int j = 0; j < listCount; j++)
					{
						if (fileType == PLY_FILE_TYPE::PLY_ASCII)
						{
							stokenizer.nextToken();
							token = stokenizer.getToken();
							plyFileReader->get_ascii_item(token, currentProp->externalType, &intValue, &uintValue, &doubleValue);
						}
						else
							plyFileReader->get_binary_item(fp, fileType, currentProp->externalType, &intValue, &uintValue, &doubleValue);
						plyFileReader->store_other_item(currentOtherProp, intValue, uintValue, doubleValue);
					}
				}
			}
		}
		if (meshType == MeshType::PLY_F_Type)
		{
			FTypePtr = createFace(currentVs);
			voidPtr = FTypePtr;
		}
		return voidPtr;
	}
	/*!
	Read an .ply file.
	\param input the input .ply file name
	*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void	CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::read_ply(const char * fileName)
	{
		//clock_t  clockBegin, clockEnd;
		//clockBegin = clock();
		VertexType vExample;
		//bool hasColor = (vExample.hasColor()) ? true : false;
		//bool hasNormal = (vExample.hasNormal()) ? true : false;
		//bool hasUV = (vExample.hasUV()) ? true : false;
		int numElements;				//num of element in object    
		char* currentEleName;			//current element's name
		int currentEleNum;              //current element's number
		PlyProperty** propList;         //properties that current elment contains.
		PlyFile* plyFile;				//current .ply file object
		// PlyFileReader plyFileReader(hasColor, hasNormal, hasUV);    //used to read .ply file
		PlyFileReader plyFileReader(false, false, false);    //used to read .ply file
		/*Open a polygon file for reading*/
		plyFile = plyFileReader.ply_open_for_reading(fileName);
		if (!plyFile) {
			printf("Can't create plyFile object!\n");
			return;
		}
		/*Iterator all element tpye in object*/
		numElements = plyFile->nelems;
		for (int i = 0; i < numElements; i++)
		{
			PlyElement* currentElement = plyFile->elems[i];
			PlyOtherElem* currentOtherElement = plyFile->other_elems->elemList[i];
			currentEleName = currentElement->elemName;
			/*Get current element's propertys list*/
			propList = plyFileReader.ply_get_element_description(plyFile, currentEleName, &currentEleNum);
			/*If current element doesn't contain property, free some objects' memory and return*/
			if (!propList)
			{
				printf("Error: The element( %s )doesn't contain any properties!\n", currentEleName);
				plyFileReader.free_ply_memory(plyFile, true);
				return;
			}
			if (plyFileReader.equal_strings("vertex", currentEleName))
			{
				mVContainer.reserve(currentEleNum);
				for (int j = 0; j < currentEleNum; j++)
					mesh_ply_get_element(plyFile, currentElement, currentOtherElement, MeshType::PLY_V_Type, &plyFileReader);
			}
			else if (plyFileReader.equal_strings("face", currentEleName))
			{
				mFContainer.reserve(currentEleNum);
				for (int j = 0; j < currentEleNum; j++)
					mesh_ply_get_element(plyFile, currentElement, currentOtherElement, MeshType::PLY_F_Type, &plyFileReader);
			}
			else if (plyFileReader.equal_strings("edge", currentEleName))
			{
				mEContainer.reserve(currentEleNum);
				for (int j = 0; j < currentEleNum; j++)
					mesh_ply_get_element(plyFile, currentElement, currentOtherElement, MeshType::PLY_E_Type, &plyFileReader);
			}
			else if (plyFileReader.equal_strings("halfedge", currentEleName))
			{
				mHEContainer.reserve(currentEleNum);
				for (int j = 0; j < currentEleNum; j++)
					mesh_ply_get_element(plyFile, currentElement, currentOtherElement, MeshType::PLY_HE_Type, &plyFileReader);
			}
			else
			{
				for (int j = 0; j < currentEleNum; j++)
					mesh_ply_get_element(plyFile, currentElement, currentOtherElement, MeshType::PLY_OTHER_Type, &plyFileReader);
			}
		}
		MyPlyOtherData myPOD(plyFile->other_elems);
		//myPOD.ShowData();
		/*Free memory*/
		plyFileReader.free_ply_memory(plyFile, false);
		/*Label boundary edges*/
#pragma omp parallel
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
		std::vector<VertexType*> isolatedVertexs;
		for (int i = 0; i < mVContainer.getCurrentIndex(); i++)
		{
			VertexType* currentVertex = mVContainer.getPointer(i);
			if (mVContainer.hasBeenDeleted(currentVertex->index()) == false)
			{
				if (currentVertex->halfedge() != NULL) {
					continue;
				}
				else {
					isolatedVertexs.push_back(currentVertex);
				}
			}
		}
		for (auto vertex : isolatedVertexs)
		{
			VertexType* currentV = vertex;
			mVContainer.deleteMember(currentV->index());
		}
		/*
		*	Arrange the boundary half_edge of boundary vertices, to make its halfedge
		*	to be the most ccw in half_edge
		*/
		for (int i = 0; i < mVContainer.getCurrentIndex(); i++)
		{
			VertexType* currentVertex = mVContainer.getPointer(i);
			if (mVContainer.hasBeenDeleted(currentVertex->index()) == false)
			{
				if (!currentVertex->boundary()) continue;
				HalfEdgeType* currentHE = vertexMostCcwInHalfEdge(currentVertex);
				while (currentHE->he_sym() != NULL)
				{
					currentHE = vertexNextCcwInHalfEdge(currentHE);
				}
				currentVertex->halfedge() = currentHE;
			}
		}
		mVMap.clear();
	}

	///*!
	//Read an .obj file.
	//\param input the input obj file name
	//*/
	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::read_obj(const char * fileName, bool removeIsolatedVertices)
	{
		FILE * pFile;
		/*Open file*/
		if ((pFile = fopen(fileName, "r")) == NULL){
			printf("Error in opening file: %s!\n", fileName);
			return;
		}
		char lineBuffer[MAX_LINE_SIZE];

		bool with_uv = false;
		bool with_normal = false;

		std::vector<Vec2> uvs;
		std::vector<Vec3> normals;

		VertexType vExample;

		while (true){
			/*If get nothing, free buffer, break!*/
			if (!fgets(lineBuffer, MAX_LINE_SIZE, pFile)){
				break;
			}

			strutil::Tokenizer stokenizer(lineBuffer, " \r\n");
			stokenizer.nextToken();
			const char * token;

			if (strcmp(stokenizer.getToken(), "v") == 0){
				VertexType * v = createVertexWithIndex();
				for (int i = 0; i < 3; i++){
					stokenizer.nextToken();
					token  = stokenizer.getToken();
					 v->position()[i] = strutil::parseStringToFloat(token);
				}

				//if (v->hasColor() && stokenizer.nextToken() != strutil::StringOver) {
				//	token = stokenizer.getToken();
				//	v->color().r = strutil::parseStringToFloat(token);
				//	stokenizer.nextToken();
				//	token = stokenizer.getToken();
				//	v->color().g = strutil::parseStringToFloat(token);
				//	token = stokenizer.getToken();
				//	v->color().b = strutil::parseStringToFloat(token);
				//}
				continue;
			}
			/*else if (strcmp(stokenizer.getToken(), "vt") == 0 && vExample.hasUV()) {
				with_uv = true;
				Vec2 uv;
				for (int i = 0; i < 2; i++){
					stokenizer.nextToken();
					token = stokenizer.getToken();
					uv[i] = strutil::parseStringToFloat(token);
				}
				uvs.push_back(uv);
				continue;
			}
			else if (strcmp(stokenizer.getToken(), "vn") == 0 && vExample.hasNormal()){
				with_normal = true;
				Vec3 n;
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();
					n[i] = strutil::parseStringToFloat(token);
				}
				normals.push_back(n);
				continue;
			}*/
			else if (strcmp(stokenizer.getToken(), "f") == 0){
				VertexType* v[3];
				for (int i = 0; i < 3; i++)
				{
					stokenizer.nextToken();
					token = stokenizer.getToken();

					strutil::Tokenizer tokenizer(token, "/\\");

					int indices[3];
					int k = 0;
					while (tokenizer.nextToken() != strutil::StringOver)
					{
						token = tokenizer.getToken();
						indices[k] = strutil::parseStringToInt(token);
						k++;
					}

					v[i] = mVContainer.getPointer(indices[0]-1);
					//if (with_uv && vExample.hasUV())
					//	v[i]->uv() = uvs[indices[1] - 1];
					//if (with_normal && vExample.hasNormal())
					//	v[i]->normal() = normals[indices[2] - 1];
				}
				createFace(v);
			}
		}
		fclose(pFile);
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
		if (removeIsolatedVertices)
		{
			/*Remove isolated vertex*/
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
			for (auto vertex : isolatedVertexs)
			{
				VertexType* currentV = vertex;
				mVContainer.deleteMember(currentV->index());
			}
		}
		
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
	}

	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::write_obj(const char * output)
	{
		VPropHandle<int> vObjIdHdl;
		addVProp(vObjIdHdl);

		FILE *fp = fopen(output, "wb");
		if (!fp) {
			printf("Fail to open output file: %s\n", output);
			return;
		}

		int objId = 1;
		for (int i = 0; i < mVContainer.getCurrentIndex(); ++i)
		{
			
			if (mVContainer.hasBeenDeleted(i)) {
				continue;
			}
			VertexType * pV = mVContainer.getPointer(i);
			gVP(vObjIdHdl, pV) = objId;
			++objId;
		}

		for (int i = 0; i < mVContainer.getCurrentIndex(); ++i)
		{
			if (mVContainer.hasBeenDeleted(i)) {
				continue;
			}
			VertexType * pV = mVContainer.getPointer(i);
			fprintf(fp, "v %f %f %f\n", pV->position()[0], pV->position()[1], pV->position()[2]);
		}

		for (int i = 0; i < mFContainer.getCurrentIndex(); ++i)
		{
			if (mFContainer.hasBeenDeleted(i)) {
				continue;
			}

			FaceType * pF = mFContainer.getPointer(i);

			HalfEdgeType * pHBegin = faceHalfedge(pF);
			HalfEdgeType * pH = pHBegin;
			fprintf(fp, "f");

			do {
				fprintf(fp, " %d", gVP(vObjIdHdl, halfedgeTarget(pH)));
				pH = faceNextCcwHalfEdge(pH);
			} while (pH != pHBegin);

			fprintf(fp, "\n");
		}

		removeVProp(vObjIdHdl);
		fclose(fp);
	}

	template<typename DType, typename VertexType, typename EdgeType, typename FaceType, typename HalfEdgeType>
	inline void CMeshBase<DType, VertexType, EdgeType, FaceType, HalfEdgeType>::readVFList(const std::vector<std::array<double, 3>>* verts, const std::vector<std::array<int, 3>>* faces, 
		const std::vector<int>* vIds, bool removeIsolatedVerts)
	{

		for (int iV = 0; iV < verts->size(); iV++)
		{
			int id = vIds != nullptr ? (*vIds)[iV] : iV;
			VertexType* currentVertex = createVertexWithId(id);
			currentVertex->position()[0] = (*verts)[iV][0];
			currentVertex->position()[1] = (*verts)[iV][1];
			currentVertex->position()[2] = (*verts)[iV][2];

		}

		for (int iF = 0; iF < faces->size(); iF++)
		{
			VertexType* currentVs[3];
			/*Assume each face is a triangle*/
			for (int j = 0; j < 3; j++)
			{
				int vId = (*faces)[iF][j];
				currentVs[j] = mVMap[vId];
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
		if (removeIsolatedVerts) {
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
			for (auto vertex : isolatedVertexs)
			{
				++numIsolatedVertices;

				VertexType* currentV = vertex;
				mVContainer.deleteMember(currentV->index());
			}

			if (numIsolatedVertices)
			{
				std::cout << "Removed " << numIsolatedVertices << " isolated vertices.";
			}
		}
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

	}//name space TriMesh

}//name space MF



#endif //_MESHLIB_BASE_MESH_H_ defined
