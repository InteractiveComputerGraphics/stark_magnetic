/** Created by jfernandez / José Antonio Fernández-Fernández
 *  Modified by Lukas Westhofen / BjoernBinaer to implement the magnetic behaviour
 */

#include "EnergyRigidBodyMagnetic.h"

#include "rigidbody_transformations.h"
#include "../time_integration.h"

#include "../../utils/mesh_generators.h"
#include "../../utils/mesh_utils.h"
#include <set>

stark::EnergyRigidBodyMagnetic::EnergyRigidBodyMagnetic(core::Stark& stark, spRigidBodyDynamics rb, std::shared_ptr<EnergyRigidBodyInertia> rb_inert)
	: rb(rb), rb_inertia(rb_inert), magnetic_method(stark.settings.simulation.magnetic_method)
{
	const MagneticMethod method = stark.settings.simulation.magnetic_method;

	stark.callbacks.add_before_time_step([&]() { this->_before_time_step(stark); });
	stark.callbacks.add_is_converged_state_valid([&]() { return this->_is_converged_state_valid(stark); });
	stark.callbacks.add_after_time_step([&]() { return this->_after_time_step(stark); });
	stark.callbacks.add_write_frame([&]() { return this->_write_frame(stark); });

	// Ensure that this callback is run
	stark.callbacks.on_time_step_accepted.insert(stark.callbacks.on_time_step_accepted.begin(), [&]() { return this->_on_time_step_accepted(stark); });

	average_sample_radius = 0.0;
	distance_threshold = std::numeric_limits<double>::max();

	if (method == MagneticMethod::DipoleMoment)
	{
		stark.global_energy.add_energy("EnergyRigidBodyMagnetic_linear_dipole", this->conn_implicit,
			[&](symx::Energy &energy, symx::Element &conn)
			{
				const std::vector<symx::Index> &bodyIDs = conn.slice(0, 2);
				const std::vector<symx::Index> &magneticIDs = conn.slice(2, 4);
				const std::vector<symx::Index> &sampleIDs = conn.slice(4, 6);

				std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->rb->dof_v,
																	   this->rb->v1,
																	   bodyIDs);
				std::vector<symx::Vector> w1 = energy.make_dof_vectors(this->rb->dof_w,
																	   this->rb->w1,
																	   bodyIDs);
				std::vector<symx::Vector> t0 = energy.make_vectors(this->rb->t0, bodyIDs);
				std::vector<symx::Vector> w0 = energy.make_vectors(this->rb->w0, bodyIDs);
				std::vector<symx::Vector> q0 = energy.make_vectors(this->rb->q0_, bodyIDs);

				symx::Scalar dt = energy.make_scalar(stark.dt);
				std::vector<symx::Vector> magn_samp_offset_loc = energy.make_vectors(
					this->magnetic_sample_pos_local, sampleIDs);

				std::vector<symx::Vector> t1 = {
					time_integration(t0[0], v1[0], dt),
					time_integration(t0[1], v1[1], dt)
				};
				std::vector<symx::Matrix> R0 = {
					quat_to_rotation(q0[0]),
					quat_to_rotation(q0[1])
				};
				std::vector<symx::Matrix> R1 = {
					quat_time_integration_as_rotation_matrix(q0[0], w1[0], dt),
					quat_time_integration_as_rotation_matrix(q0[1], w1[1], dt)
				};
				std::vector<symx::Vector> m0 = energy.make_vectors(magnetic_dipole_moment_local, sampleIDs);
				std::vector<symx::Vector> m1 = {
					local_to_global_direction(m0[0], R1[0]),
					local_to_global_direction(m0[1], R1[1]),
				};

				std::vector<symx::Vector> global_offset_from_cog = {
					local_to_global_direction(magn_samp_offset_loc[0], R1[0]),
					local_to_global_direction(magn_samp_offset_loc[1], R1[1])
				};

				std::vector<symx::Vector> pos = {
					t1[0] + global_offset_from_cog[0],
					t1[1] + global_offset_from_cog[1]
				};

				symx::Vector dist_vec = pos[1] - pos[0];
				symx::Scalar dist = dist_vec.norm();
				symx::Vector norm_dist = dist_vec / dist;

				symx::Vector B_ext = (3.0 * norm_dist * norm_dist.dot(m1[1]) - m1[1]) * MU_0 / (4.0 * M_PI * dist * dist * dist);

				symx::Scalar E = -0.5 * m1[0].dot(B_ext);
				energy.set(E);
			}
		);
	}
	else if (method == MagneticMethod::DipoleMoment_DoF)
	{
		stark.global_energy.add_energy("EnergyRigidBodyMagnetic_linear_dipole", this->conn_fully_implicit,
		    [&](symx::Energy &energy, symx::Element &conn)
			{
				const std::vector<symx::Index> &bodyIDs = conn.slice(0, 2);
				const std::vector<symx::Index> &magneticIDs = conn.slice(2, 4);
				const std::vector<symx::Index> &magneticMaterials = conn.slice(4, 6);
				const std::vector<symx::Index> &sampleIDs = conn.slice(6, 8);

				std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->rb->dof_v,
																	   this->rb->v1,
																	   bodyIDs);
				std::vector<symx::Vector> w1 = energy.make_dof_vectors(this->rb->dof_w,
																	   this->rb->w1,
																	   bodyIDs);


				std::vector<symx::Vector> t0 = energy.make_vectors(this->rb->t0, bodyIDs);
				std::vector<symx::Vector> w0 = energy.make_vectors(this->rb->w0, bodyIDs);
				std::vector<symx::Vector> q0 = energy.make_vectors(this->rb->q0_, bodyIDs);

				symx::Scalar dt = energy.make_scalar(stark.dt);
				std::vector<symx::Vector> magn_samp_offset_loc = energy.make_vectors(
					this->magnetic_sample_pos_local, sampleIDs);

				std::vector<symx::Vector> t1 = {
					time_integration(t0[0], v1[0], dt),
					time_integration(t0[1], v1[1], dt)
				};
				std::vector<symx::Matrix> R0 = {
					quat_to_rotation(q0[0]),
					quat_to_rotation(q0[1])
				};
				std::vector<symx::Matrix> R1 = {
					quat_time_integration_as_rotation_matrix(q0[0], w1[0], dt),
					quat_time_integration_as_rotation_matrix(q0[1], w1[1], dt)
				};

				std::vector<symx::Vector> m1_glob;
				for (int i = 0; i < 2; i++)
				{
					if (static_cast<MagneticMaterial>(magneticMaterials[i].idx) == MagneticMaterial::Permanent)
					{
						m1_glob.push_back(energy.make_vector(this->magnetic_dipole_moment_local, magneticIDs[i]));
					}
					else if (static_cast<MagneticMaterial>(magneticMaterials[i].idx) == MagneticMaterial::Linear)
					{
						m1_glob.push_back(energy.make_dof_vector(this->dof_m, this->magnetic_dipole_moment_local, magneticIDs[i]));
					}
				}

				std::vector<symx::Vector> global_offset_from_cog = {
					local_to_global_direction(magn_samp_offset_loc[0], R1[0]),
					local_to_global_direction(magn_samp_offset_loc[1], R1[1])
				};

				std::vector<symx::Vector> pos = {
					t1[0] + global_offset_from_cog[0],
					t1[1] + global_offset_from_cog[1]
				};

				symx::Vector dist_vec = pos[1] - pos[0];
				symx::Scalar dist = dist_vec.norm();
				symx::Vector norm_dist = dist_vec / dist;

				symx::Vector B_ext = MU_0 * (3.0 * norm_dist * norm_dist.dot(m1_glob[1]) - m1_glob[1]) /
									 (4.0 * M_PI * dist * dist * dist);

				symx::Scalar E = -m1_glob[0].dot(B_ext);
				energy.set(E);
			}
		);
	}
	else if (method == MagneticMethod::Jackson_old)
	{
		stark.global_energy.add_energy("EnergyRigidBodyMagnetic_linear_Jackson", this->conn_implicit,
            [&](symx::Energy &energy, symx::Element &conn)
			{
				const std::vector<symx::Index> &bodyIDs = conn.slice(0, 2);
				const std::vector<symx::Index> &magneticIDs = conn.slice(2, 4);
				const std::vector<symx::Index> &sampleIDs = conn.slice(4, 6);

				std::vector<symx::Vector> v1 = energy.make_dof_vectors(this->rb->dof_v,
																	   this->rb->v1,
																	   bodyIDs);
				std::vector<symx::Vector> w1 = energy.make_dof_vectors(this->rb->dof_w,
																	   this->rb->w1,
																	   bodyIDs);
				std::vector<symx::Vector> t0 = energy.make_vectors(this->rb->t0, bodyIDs);
				std::vector<symx::Vector> w0 = energy.make_vectors(this->rb->w0, bodyIDs);
				std::vector<symx::Vector> q0 = energy.make_vectors(this->rb->q0_, bodyIDs);

				symx::Scalar dt = energy.make_scalar(stark.dt);
				std::vector<symx::Vector> magn_samp_offset_loc = energy.make_vectors(
					this->magnetic_sample_pos_local, sampleIDs);

				std::vector<symx::Vector> t1 = {
					time_integration(t0[0], v1[0], dt),
					time_integration(t0[1], v1[1], dt)
				};
				std::vector<symx::Matrix> R0 = {
					quat_to_rotation(q0[0]),
					quat_to_rotation(q0[1])
				};
				std::vector<symx::Matrix> R1 = {
					quat_time_integration_as_rotation_matrix(q0[0], w1[0], dt),
					quat_time_integration_as_rotation_matrix(q0[1], w1[1], dt)
				};

				std::vector<symx::Vector> magn_dip_mom_glob = energy.make_vectors(this->magnetic_dipole_moment_local, sampleIDs);

				std::vector<symx::Vector> global_offset_from_cog = {
					local_to_global_direction(magn_samp_offset_loc[0], R1[0]),
					local_to_global_direction(magn_samp_offset_loc[1], R1[1])
				};

				std::vector<symx::Vector> pos = {
					t1[0] + global_offset_from_cog[0],
					t1[1] + global_offset_from_cog[1]
				};

				symx::Vector dist_vec = pos[1] - pos[0];
				symx::Scalar dist = dist_vec.norm();
				symx::Vector norm_dist = dist_vec / dist;
				std::vector<symx::Scalar> volume = {
					energy.make_scalar(this->magnetic_sample_volume, magneticIDs[0]),
					energy.make_scalar(this->magnetic_sample_volume, magneticIDs[1])
				};

				symx::Vector H_ext = (3.0 * norm_dist * norm_dist.dot(magn_dip_mom_glob[1]) -
									  magn_dip_mom_glob[1]) /
									 (4.0 * M_PI * dist * dist * dist);
				symx::Vector H_int = -magn_dip_mom_glob[0] / (3.0 * volume[0]);
				symx::Vector H = H_ext + H_int;

				symx::Scalar magn_permeab_0 = energy.make_scalar(this->magnetic_permeability,
																 magneticIDs[0]);

				symx::Scalar E = 0.5 * (magn_permeab_0 - MU_0) * volume[0] * H.dot(H_ext);
				energy.set(E);
			}
		);
	}
}

stark::EnergyRigidBodyMagnetic::Handler stark::EnergyRigidBodyMagnetic::add(const RigidBodyHandler& rb, const Params& params)
{
	rb.exit_if_not_valid("EnergyRigidBodyMagnetic::add()");
	const int magnetized_rb_idx = (int)this->magnetic_rb_id.size();
	this->magnetic_rb_id.push_back(rb.get_idx());

	this->magnetic_material.push_back(params.magnetic_material);
	this->magnetic_permeability.push_back(params.magnetic_permeability);

	this->enabled.push_back(params.enabled);

	this->magnetized_rb_to_global_rb_map[magnetized_rb_idx] = rb.get_idx();
	this->magnetic_sample_range.push_back({0, 0});

	return Handler(this, magnetized_rb_idx);
}
stark::EnergyRigidBodyMagnetic::Params stark::EnergyRigidBodyMagnetic::get_params(const Handler& handler) const
{
	handler.exit_if_not_valid("EnergyRigidBodyMagnetic::get_params()");
	stark::EnergyRigidBodyMagnetic::Params params;
	params.magnetic_material = this->magnetic_material[handler.get_idx()];
	params.magnetic_permeability = this->magnetic_permeability[handler.get_idx()];
	params.enabled = this->enabled[handler.get_idx()];

	return params;
}
void stark::EnergyRigidBodyMagnetic::set_params(const Handler& handler, const Params& params)
{
	handler.exit_if_not_valid("EnergyRigidBodyMagnetic::set_params()");
	this->magnetic_material[handler.get_idx()] = params.magnetic_material;
	this->magnetic_permeability[handler.get_idx()] = params.magnetic_permeability;
	this->enabled[handler.get_idx()] = params.enabled;
}

void stark::EnergyRigidBodyMagnetic::add_magnetic_samples(const stark::EnergyRigidBodyMagnetic::Handler &handler,
														  const std::vector<Eigen::Vector3d> &loc_pos,
														  const double radius,
														  const Eigen::Vector3d &glob_mag_dip0)
{
	handler.exit_if_not_valid("EnergyRigidBodyMagnetic::add_magnetic_samples()");
	const unsigned int num_samples = loc_pos.size();
	const int sample_offset = static_cast<int>(this->magnetic_sample_pos_local.size());

	this->magnetic_sample_pos_local.insert(this->magnetic_sample_pos_local.end(), loc_pos.begin(), loc_pos.end());
	this->magnetic_dipole_moment_local.insert(this->magnetic_dipole_moment_local.end(), num_samples, global_to_local_direction(glob_mag_dip0, this->rb->R0[this->magnetic_rb_id[handler.get_idx()]]));
	this->magnetic_sample_radius.insert(this->magnetic_sample_radius.end(), num_samples, radius);
	this->magnetic_sample_volume.insert(this->magnetic_sample_volume.end(), num_samples, 4.0 / 3.0 * M_PI * pow(radius, 3.0));
	this->is_strongly_coupled.insert(this->is_strongly_coupled.end(), num_samples, {});
	this->force_ij.insert(this->force_ij.end(), num_samples, {});
	this->torque_ij.insert(this->torque_ij.end(), num_samples, {});

	const int new_sample_offset = static_cast<int>(this->magnetic_sample_pos_local.size());
	this->magnetic_sample_range[handler.get_idx()] = { sample_offset, new_sample_offset };

	// Small hack to calculate the cumulative average. Breaks if the sample offset array does not correspond to the amount of different radii in the simulation!
	this->average_sample_radius += (radius - this->average_sample_radius) / magnetic_sample_range.size();
}

void stark::EnergyRigidBodyMagnetic::add_magnetic_samples(const stark::EnergyRigidBodyMagnetic::Handler &handler,
														  const std::vector<Eigen::Vector3d> &loc_pos,
														  const double radius,
														  const std::vector<Eigen::Vector3d>& glob_mag_dip0_vec)
{
	handler.exit_if_not_valid("EnergyRigidBodyMagnetic::add_magnetic_samples()");
	const unsigned int num_samples = loc_pos.size();
	const int sample_offset = static_cast<int>(this->magnetic_sample_pos_local.size());

	if (loc_pos.size() != glob_mag_dip0_vec.size()){
		std::cout << "EnergyRigidBodyMagnetic::add_magnetic_samples(): Magnetic moment vector does not match positions in size!" << std::endl;
		return;
	}

	this->magnetic_sample_pos_local.insert(this->magnetic_sample_pos_local.end(), loc_pos.begin(), loc_pos.end());

	for (const auto& glob_mag_dip0 : glob_mag_dip0_vec){
		this->magnetic_dipole_moment_local.push_back(global_to_local_direction(glob_mag_dip0, this->rb->R0[this->magnetic_rb_id[handler.get_idx()]]));
	}
	this->magnetic_sample_radius.insert(this->magnetic_sample_radius.end(), num_samples, radius);
	this->magnetic_sample_volume.insert(this->magnetic_sample_volume.end(), num_samples, 4.0 / 3.0 * M_PI * pow(radius, 3.0));
	this->is_strongly_coupled.insert(this->is_strongly_coupled.end(), num_samples, {});
	this->force_ij.insert(this->force_ij.end(), num_samples, {});
	this->torque_ij.insert(this->torque_ij.end(), num_samples, {});

	const int new_sample_offset = static_cast<int>(this->magnetic_sample_pos_local.size());
	this->magnetic_sample_range[handler.get_idx()] = { sample_offset, new_sample_offset };

	// Small hack to calculate the cumulative average. Breaks if the sample offset array does not correspond to the amount of different radii in the simulation!
	this->average_sample_radius += (radius - this->average_sample_radius) / magnetic_sample_range.size();
}

void stark::EnergyRigidBodyMagnetic::_before_time_step(core::Stark& stark)
{
	stark.logger.start_timing("magnetic_before_time_step");

	this->conn_implicit.clear();
	this->conn_fully_implicit.clear();
	const int n = (int)this->magnetic_rb_id.size();

	std::set<std::pair<unsigned int, unsigned int>> magnetic_hessian_pairs;

	std::vector<Eigen::Vector3d> new_magnetic_dipole_moments_local(this->magnetic_dipole_moment_local.size(), Eigen::Vector3d::Zero());

	// Determine dipole moment of linear magnets
	#pragma omp parallel for schedule(static) default(shared)
	for (int mrb_a = 0; mrb_a < n; mrb_a++)
	{
		if (!this->enabled[mrb_a])
			continue;

		const int rb_a = this->magnetic_rb_id[mrb_a];
		const std::array<int, 2> &magn_sample_range_i = this->magnetic_sample_range[mrb_a];

		for (int sampleID_i = magn_sample_range_i[0]; sampleID_i < magn_sample_range_i[1]; sampleID_i++)
		{
			// Compute magnetization and magnetic dipole moment for linear magnets
			if (this->magnetic_material[mrb_a] == MagneticMaterial::Linear)
			{
				Eigen::Vector3d global_magn_dip_mom_i = Eigen::Vector3d::Zero();

				const Eigen::Vector3d sample_pos_global_i =
					local_to_global_direction(this->magnetic_sample_pos_local[sampleID_i],
											  this->rb->R0[rb_a]) +
					this->rb->t0[rb_a];

				for (int mrb_b = 0; mrb_b < n; mrb_b++)
				{
					const int rb_b = this->magnetic_rb_id[mrb_b];
					const std::array<int, 2> &magn_sample_range_j = this->magnetic_sample_range[mrb_b];

					if (!this->enabled[mrb_b] || rb_a == rb_b)
						continue;

					for (int sampleID_j = magn_sample_range_j[0]; sampleID_j < magn_sample_range_j[1]; sampleID_j++)
					{
						const Eigen::Vector3d sample_pos_global_j =
							local_to_global_direction(this->magnetic_sample_pos_local[sampleID_j],
													  this->rb->R0[rb_b]) +
							this->rb->t0[rb_b];
						const Eigen::Vector3d dist_vec = sample_pos_global_i - sample_pos_global_j;

						const Eigen::Vector3d mag_dip_mom_j = local_to_global_direction(this->magnetic_dipole_moment_local[sampleID_j], this->rb->R0[rb_b]);
						const Eigen::Vector3d norm_dist_vec = dist_vec.normalized();

						const Eigen::Vector3d H_ext =
							(3.0 * norm_dist_vec * norm_dist_vec.dot(mag_dip_mom_j) - mag_dip_mom_j) /
							(4.0 * M_PI * pow(dist_vec.norm(), 3.0));
						global_magn_dip_mom_i += H_ext;
					}
				}

				const double magn_permeab_i = this->magnetic_permeability[mrb_a];
				const double magn_susc_i = magn_permeab_i / MU_0 - 1.0;
				global_magn_dip_mom_i *=
					this->magnetic_sample_volume[mrb_a] * (magn_susc_i / (1.0 + magn_susc_i / 3.0));

				new_magnetic_dipole_moments_local[sampleID_i] = global_to_local_direction(global_magn_dip_mom_i, this->rb->R0[rb_a]);
			}
		}
	}

	std::vector<std::vector<std::array<int32_t, 6>>> conn_implicit_per_body(n);

	#pragma omp parallel for schedule(static) default(shared)
	for (int mrb_a = 0; mrb_a < n; mrb_a++)
	{
		if (!this->enabled[mrb_a])
			continue;

		const int rb_a = this->magnetic_rb_id[mrb_a];
		const std::array<int, 2> &magn_sample_range_i = this->magnetic_sample_range[mrb_a];

		for (int mrb_b = 0; mrb_b < n; mrb_b++)
		{
			const int rb_b = this->magnetic_rb_id[mrb_b];
			const std::array<int, 2> &magn_sample_range_j = this->magnetic_sample_range[mrb_b];

			if (!this->enabled[mrb_b] || rb_a == rb_b)
				continue;

			for (int sampleID_i = magn_sample_range_i[0]; sampleID_i < magn_sample_range_i[1]; sampleID_i++)
			{
				this->is_strongly_coupled[sampleID_i] = std::vector<bool>(this->magnetic_dipole_moment_local.size(), false);
				this->force_ij[sampleID_i] = std::vector<Eigen::Vector3d>(this->magnetic_dipole_moment_local.size(), Eigen::Vector3d::Zero());
				this->torque_ij[sampleID_i] = std::vector<Eigen::Vector3d>(this->magnetic_dipole_moment_local.size(), Eigen::Vector3d::Zero());

				if (this->magnetic_material[mrb_a] == MagneticMaterial::Linear)
					this->magnetic_dipole_moment_local[sampleID_i] = new_magnetic_dipole_moments_local[sampleID_i];

				// Find magnetic pairs
				for (int sampleID_j = magn_sample_range_j[0]; sampleID_j < magn_sample_range_j[1]; sampleID_j++)
				{
					const Eigen::Vector3d sample_pos_global_i =
						local_to_global_direction(this->magnetic_sample_pos_local[sampleID_i], this->rb->R0[rb_a]) +
						this->rb->t0[rb_a];
					const Eigen::Vector3d sample_pos_global_j =
						local_to_global_direction(this->magnetic_sample_pos_local[sampleID_j], this->rb->R0[rb_b]) +
						this->rb->t0[rb_b];
					const Eigen::Vector3d r_ij = sample_pos_global_i - sample_pos_global_j;

					if (r_ij.norm() < distance_threshold)
					{
						this->is_strongly_coupled[sampleID_i][sampleID_j] = true;
						conn_implicit_per_body[mrb_a].push_back({rb_a, rb_b, mrb_a, mrb_b, sampleID_i, sampleID_j});
						#pragma omp critical
						{
							magnetic_hessian_pairs.insert({rb_a, rb_b});
						};
					}
					else
					{
						// Note: Updating to R1 and t1 does not yield any improvements!

						// Magnetic values in world coordinates
						const Eigen::Vector3d mag_dip_mom_i = local_to_global_direction(this->magnetic_dipole_moment_local[sampleID_i], this->rb->R0[rb_a]);
						const Eigen::Vector3d mag_dip_mom_j = local_to_global_direction(this->magnetic_dipole_moment_local[sampleID_j], this->rb->R0[rb_b]);
						const Eigen::Vector3d r_ij_hat = r_ij.normalized();

						// Magnetic force
						const Eigen::Vector3d force_ij = MU_0 / (8.0 * M_PI * pow(r_ij.norm(), 4.0)) *
														 (-15.0 * r_ij_hat * (mag_dip_mom_i.dot(r_ij_hat)) *
														  (mag_dip_mom_j.dot(r_ij_hat))
														  + 3.0 * r_ij_hat * (mag_dip_mom_j.dot(mag_dip_mom_i))
														  + 3.0 * mag_dip_mom_i * (mag_dip_mom_j.dot(r_ij_hat))
														  + 3.0 * mag_dip_mom_j * (mag_dip_mom_i.dot(r_ij_hat)));

						const Eigen::Vector3d torque_ij = MU_0 / (8.0 * M_PI * pow(r_ij.norm(), 3.0)) *
														  (3.0 * r_ij_hat.dot(mag_dip_mom_j) *
														   mag_dip_mom_i.cross(r_ij_hat)
														   - mag_dip_mom_i.cross(mag_dip_mom_j));

						// Magnetic and mechanical torque
						this->rb->force[rb_a] += force_ij;
						this->rb->torque[rb_a] += (sample_pos_global_i - this->rb->t0[rb_a]).cross(force_ij);

						// Magnetic torque
						this->rb->torque[rb_a] += torque_ij;

						// Debug fields
						this->force_ij[sampleID_i][sampleID_j] = force_ij;
						this->torque_ij[sampleID_i][sampleID_j] = (sample_pos_global_i - this->rb->t0[rb_a]).cross(force_ij) + torque_ij;
					}
				}
			}
		}
	}

	for (const auto& rb_conns : conn_implicit_per_body)
		for (const auto& rb_conn : rb_conns)
			conn_implicit.push_back(rb_conn);

	std::cout << "Implicit pairs: " << conn_implicit.size() << std::endl;

	stark.logger.stop_timing_add("magnetic_before_time_step");
}
bool stark::EnergyRigidBodyMagnetic::_is_converged_state_valid(core::Stark& stark)
{
	// TODO
	return true;
}
void stark::EnergyRigidBodyMagnetic::_on_time_step_accepted(stark::core::Stark &stark)
{

}
void stark::EnergyRigidBodyMagnetic::_after_time_step(core::Stark &stark)
{
	for (std::vector<bool>& interaction_vector : this->is_strongly_coupled)
		interaction_vector.clear();
}
void stark::EnergyRigidBodyMagnetic::_write_frame(core::Stark &stark)
{
	Mesh<3> cylinder_arrow = make_cylinder(0.1 * average_sample_radius, 0.5 * average_sample_radius);
	Mesh<3> cone_arrow = make_cone(0.25 * average_sample_radius, 0.25 * average_sample_radius);

	std::for_each(cone_arrow.conn.begin(), cone_arrow.conn.end(), [&](std::array<int, 3>& elem) {
		elem[0] += cylinder_arrow.vertices.size();
		elem[1] += cylinder_arrow.vertices.size();
		elem[2] += cylinder_arrow.vertices.size();
	});

	move(cylinder_arrow.vertices, -0.125 * average_sample_radius * Eigen::Vector3d::UnitZ());
	move(cone_arrow.vertices, 0.125 * average_sample_radius * Eigen::Vector3d::UnitZ());

	// Note that arrow is aligned with the z-axis!
	Mesh<3> arrow = Mesh<3>();
	arrow.vertices.insert(arrow.vertices.end(), cylinder_arrow.vertices.begin(), cylinder_arrow.vertices.end());
	arrow.vertices.insert(arrow.vertices.end(), cone_arrow.vertices.begin(), cone_arrow.vertices.end());
	arrow.conn.insert(arrow.conn.end(), cylinder_arrow.conn.begin(), cylinder_arrow.conn.end());
	arrow.conn.insert(arrow.conn.end(), cone_arrow.conn.begin(), cone_arrow.conn.end());

	std::vector<Eigen::Vector3d> vtk_arrow_vertices;
	std::vector<std::array<int, 3>> vtk_conn_buffer;

	std::vector<Eigen::Vector3d> vtk_rb_cog_vertices;
	std::vector<Eigen::Vector3d> vtk_rb_directions;

	const int n = (int)this->magnetic_rb_id.size();
	for (int mrb_id = 0; mrb_id < n; mrb_id++)
	{
		const int rb_id = this->magnetic_rb_id[mrb_id];
		const std::array<int, 2> &magn_sample_range_i = this->magnetic_sample_range[mrb_id];

		const Eigen::Vector3d& pos = this->rb->t0[rb_id];
		const Eigen::Matrix3d& R0 = this->rb->R0[rb_id];

		vtk_rb_cog_vertices.push_back(pos);
		vtk_rb_directions.push_back(R0 * Eigen::Vector3d::UnitZ());

		for (int sampleID_i = magn_sample_range_i[0]; sampleID_i < magn_sample_range_i[1]; sampleID_i++)
		{
			const Eigen::Vector3d& sample_pos_global = local_to_global_point(this->magnetic_sample_pos_local[sampleID_i], this->rb->R0[rb_id], pos);
			const Eigen::Vector3d& global_magn_dipole = local_to_global_direction(this->magnetic_dipole_moment_local[sampleID_i], this->rb->R0[rb_id]);
			const Eigen::Vector3d& global_magn_dip_dir = global_magn_dipole.normalized();

			// Axis for rotation of x-axis to µ
			const Eigen::Vector3d viz_rotation_axis = Eigen::Vector3d::UnitZ().cross(global_magn_dip_dir).normalized();
			// Angle of rotation
			const float viz_rotation_angle_rad = acos(Eigen::Vector3d::UnitZ().dot(global_magn_dip_dir));

			Mesh<3> local_arrow = arrow;
			rotate_deg(local_arrow.vertices, rad2deg(viz_rotation_angle_rad), viz_rotation_axis);

			const int out_vertex_offset = (int)vtk_arrow_vertices.size();
			for (const auto& vert: local_arrow.vertices)
			{
				vtk_arrow_vertices.emplace_back(vert + sample_pos_global);
			}
			for (const auto& face: local_arrow.conn)
			{
				vtk_conn_buffer.push_back({
					face[0] + out_vertex_offset,
					face[1] + out_vertex_offset,
					face[2] + out_vertex_offset
				});
			}
		}
	}

	stark::write_VTK(stark.get_frame_path("magnetic_arrows") + ".vtk", vtk_arrow_vertices, vtk_conn_buffer);
	//stark::write_VTK(stark.get_frame_path("magnetic_dir") + ".vtk", vtk_rb_cog_vertices, vtk_rb_directions);
	//stark::write_VTK(stark.get_frame_path("magnetic_arrows") + ".vtk", arrow.vertices, arrow.conn);
}
