#ifndef _MESHFRAME_FACE_H_
#define _MESHFRAME_FACE_H_
/*!
*      \file Face.h
*      \brief Base class of face
*
*	   \Version 1.0
*	   \Update 11/03/2022
*/
#include <assert.h>

#define STRINGSIZE 256

namespace MF{
	namespace TriMesh {


		class CHalfEdge;

		/*!
			\brief CFace base class of all kinds of face classes
		*/
		class CFace
		{
		public:
			/*!
			CFace constructor
			*/
			CFace() { m_halfedge = NULL; };
			/*!
			CFace destructor
			*/
			~CFace() {};
			/*!
				One of the halfedges attaching to the current face.
			*/
			CHalfEdge*& halfedge() { return m_halfedge; };
			/*!
				index of the current face
			*/
			size_t& index() { return m_index; };
			/*!
				The reference to the current face id
			*/
			int& id() { return m_id; };
			/*!
				The value of the current face id.
			*/
			const int             id() const { return m_id; };

			//virtual bool hasNormal() {
			//	return false;
			//}
			//virtual CPoint & normal() {
			//	static CPoint _normal;
			//	printf("Face does not have normal!\n");
			//	assert(false);
			//	system("pause");
			//	return _normal;
			//}

			//virtual bool hasColor() {
			//	return false;
			//}
			//virtual ColorUnion & color() {
			//	printf("Face does not have color!\n");
			//	static ColorUnion _color;
			//	assert(false);
			//	system("pause");
			//	return _color;
			//}
		protected:
			/*!
				index of the current face
			*/
			size_t			   m_index;
			/*!
				id of the current face
			*/
			int			       m_id;
			/*!
				One halfedge  attaching to the current face.
			*/
			CHalfEdge* m_halfedge;
		};

	}//name space TriMesh
}//name space MF

#endif 