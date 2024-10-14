#include "SDFFunctions.h"
#include "JanBenderUtilities/OBJLoader.h"

using namespace Eigen;
using namespace std;
using namespace Utilities;

Eigen::AlignedBox<double, 3> SDFFunctions::computeBoundingBox(const unsigned int numVertices, const Eigen::Vector3d *vertices)
{
	Eigen::AlignedBox<double, 3> box;

	// compute bounding box	 
	box.min() = vertices[0];
	box.max() = box.min();
	box.setEmpty();
	for (unsigned int i = 1; i < numVertices; ++i)
	{
		const Eigen::Vector3d& p = vertices[i];
		box.extend(p);
	}
	return box;
}

double SDFFunctions::distance(Discregrid::CubicLagrangeDiscreteGrid* sdf, const Eigen::Vector3d &x,
	const double thickness, Eigen::Vector3d &normal, Eigen::Vector3d &nextSurfacePoint)
{
	Eigen::Vector3d n;
	const double dist = sdf->interpolate(0, x.template cast<double>(), &n);
	if (dist == std::numeric_limits<double>::max())
		return dist;
	normal.normalize();
	normal = n.template cast<double>();

	nextSurfacePoint = (x - dist * normal);

	return dist - thickness;
}

double SDFFunctions::distance(Discregrid::CubicLagrangeDiscreteGrid* sdf, const Eigen::Vector3d &x,
	const double thickness)
{
	const double dist = sdf->interpolate(0, x.template cast<double>());
	if (dist == std::numeric_limits<double>::max())
		return dist;
	return dist - thickness;
}

Discregrid::CubicLagrangeDiscreteGrid* SDFFunctions::generateSDF(const unsigned int numVertices,
	const Eigen::Vector3d *vertices, const unsigned int numFaces, const unsigned int *faces,
	const Eigen::AlignedBox<double, 3> &bbox, const std::array<unsigned int, 3> &resolution, const bool invert)
{
	//START_TIMING("SDF Generation");
	//////////////////////////////////////////////////////////////////////////
	// Generate distance field of object using Discregrid
	//////////////////////////////////////////////////////////////////////////
	Discregrid::TriangleMesh sdfMesh(&vertices[0][0], faces, numVertices, numFaces);

	Discregrid::MeshDistance md(sdfMesh);
	Eigen::AlignedBox3d domain;
	domain.extend(bbox.min().cast<double>());
	domain.extend(bbox.max().cast<double>());
	domain.max() += 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();
	domain.min() -= 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();

	Discregrid::CubicLagrangeDiscreteGrid *distanceField = new Discregrid::CubicLagrangeDiscreteGrid(domain, resolution);
	auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
	double factor = 1.0;
	if (invert)
		factor = -1.0;
	func = [&md,&factor](Eigen::Vector3d const& xi) {return factor * md.signedDistanceCached(xi); };

	distanceField->addFunction(func, false);
	//STOP_TIMING_PRINT;

	return distanceField;
}
