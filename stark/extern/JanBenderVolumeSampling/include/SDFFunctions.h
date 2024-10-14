#ifndef SDFFunctions_H
#define SDFFunctions_H

#include <Eigen/Eigen>
#include <Discregrid/All>

namespace Utilities
{
	/** \brief Functions for generating and querying an SDF. 
	*/
	class SDFFunctions
	{
	public:
		/** Generate SDF from mesh.
		*/
		static Discregrid::CubicLagrangeDiscreteGrid* generateSDF(const unsigned int numVertices, 
			const Eigen::Vector3d *vertices, const unsigned int numFaces, const unsigned int *faces,
			const Eigen::AlignedBox<double, 3> &bbox, const std::array<unsigned int, 3> &resolution, const bool invert=false);
	
		/** Compute the bounding box of a mesh. 
		 */
		static Eigen::AlignedBox<double, 3> computeBoundingBox(const unsigned int numVertices, const Eigen::Vector3d *vertices);

		/** Determine distance of a point x to the surface represented by the SDF and corresponding surface normal and 
		* next point on the surface.
		*/
		static double distance(Discregrid::CubicLagrangeDiscreteGrid* sdf, const Eigen::Vector3d &x,
			const double thickness, Eigen::Vector3d &normal, Eigen::Vector3d &nextSurfacePoint);

		/** Determine distance of a point x to the surface represented by the SDF. 
		 */
		static double distance(Discregrid::CubicLagrangeDiscreteGrid* sdf,
			const Eigen::Vector3d &x, const double thickness);

	};
}

#endif // SDFFunctions_H
