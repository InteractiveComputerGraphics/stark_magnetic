/** Created by jfernandez / José Antonio Fernández-Fernández
 *  Modified by Lukas Westhofen / BjoernBinaer to implement the magnetic behaviour
 */

#pragma once
#include "RigidBodyHandler.h"
#include "../types.h"

namespace stark
{
	enum class MagneticMaterial { Permanent = 0, Linear };
	inline std::ostream& operator<<(std::ostream& os, const MagneticMaterial& m)
	{
		switch (m)
		{
			case MagneticMaterial::Permanent:	return (os << "MagneticMaterial::Permanent");
			case MagneticMaterial::Linear:		return (os << "MagneticMaterial::Linear");
		}
		return os << "MagneticMaterial::Error";
	}

	class EnergyRigidBodyMagnetic
	{
	public:
        /* Types */
        struct Params
        {
            STARK_PARAM_NON_NEGATIVE(double, magnetic_permeability, 0.0)
			STARK_PARAM_NO_VALIDATION(stark::MagneticMaterial, magnetic_material, stark::MagneticMaterial::Permanent)
			STARK_PARAM_NO_VALIDATION(bool, enabled, true)
        };
        struct Handler { STARK_COMMON_HANDLER_CONTENTS(EnergyRigidBodyMagnetic, Params) };

    private:
		const double MU_0 = 1.2566370614e-6;

		MagneticMethod magnetic_method = MagneticMethod::DipoleMoment;

        /* Fields */
        spRigidBodyDynamics rb;
		std::shared_ptr<EnergyRigidBodyInertia> rb_inertia;
        symx::LabelledConnectivity<6> conn_implicit{ { "rb_a", "rb_b", "mrb_a", "mrb_b", "sample_i", "sample_j" }};
		symx::LabelledConnectivity<8> conn_fully_implicit{ { "rb_a", "rb_b", "mrb_a", "mrb_b", "type_a", "type_b", "sample_i", "sample_j" }};

		/* Mapping from magnetic rb to rb */
		std::unordered_map<int, int> magnetized_rb_to_global_rb_map;

		/* Per magnetic body fields */
		std::vector<MagneticMaterial> magnetic_material;
		std::vector<double> magnetic_permeability;

		std::vector<int> magnetic_rb_id;
		std::vector<std::array<int, 2>> magnetic_sample_range;

		std::vector<bool> enabled;

		/* Per sample fields */
		std::vector<Eigen::Vector3d> magnetic_sample_pos_local;
		std::vector<Eigen::Vector3d> magnetic_dipole_moment_local;
		std::vector<double> magnetic_sample_radius;
		std::vector<double> magnetic_sample_volume;

		/* Per sample debug fields */
		std::vector<std::vector<bool>> is_strongly_coupled;
		std::vector<std::vector<Eigen::Vector3d>> force_ij;
		std::vector<std::vector<Eigen::Vector3d>> torque_ij;

		double average_sample_radius;
		double distance_threshold;

		// Magnetic dipole moment in global space as a degree of freedom for linear materials
		symx::DoF dof_m;

	public:
        /* Methods */
        EnergyRigidBodyMagnetic(core::Stark& stark, spRigidBodyDynamics rb, std::shared_ptr<EnergyRigidBodyInertia> rb_inert);
        Handler add(const RigidBodyHandler& rb, const Params& params);
        Params get_params(const Handler& handler) const;
        void set_params(const Handler& handler, const Params& params);

		void add_magnetic_samples(const Handler& handler, const std::vector<Eigen::Vector3d>& loc_pos, const double radius, const Eigen::Vector3d& glob_mag_dip0 = Eigen::Vector3d::Zero());
		void add_magnetic_samples(const Handler& handler, const std::vector<Eigen::Vector3d>& loc_pos, const double radius, const std::vector<Eigen::Vector3d>& glob_mag_dip0_vec);

		void set_distance_threshold(double threshold) { distance_threshold = threshold; }
		double get_distance_threshold() { return distance_threshold; }

    private:
        // Stark callbacks
        void _before_time_step(core::Stark& stark);
        bool _is_converged_state_valid(core::Stark& stark);
		void _on_time_step_accepted(core::Stark& stark);
		void _after_time_step(core::Stark& stark);
		void _write_frame(core::Stark& stark);
	};
}
