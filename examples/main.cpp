#include <iostream>

#include <stark>
#include <filesystem>
#include "paths.h"
#include "rb_constraint_test_scenes.h"

#include <VolumeSampling.h>

using MagneticMaterial = stark::MagneticMaterial;

float mu_neodymium = 1.32e-6;
float mu_nickel = 1.26e-4;
float mu_iron = 6.3e-3;


void sample_magnetic_object(const std::vector<Eigen::Vector3d>& vertices, const std::vector<unsigned int>& faces, const double radius, std::vector<Eigen::Vector3d>& magnetic_pos_samples)
{
	Utilities::VolumeSampling::sampleMesh(vertices.size(), vertices.data(), faces.size() / 3, faces.data(), radius,
										  nullptr, {50, 50, 50}, false, 0, magnetic_pos_samples);

	std::cout << "Created N=" << magnetic_pos_samples.size() << " samples for the magnetic needles!" << std::endl;
}

void sample_magnetic_object(const stark::Mesh<3> mesh, const double radius, std::vector<Eigen::Vector3d>& magnetic_pos_samples)
{
	std::vector<unsigned int> casted_faces;
	for (const auto& face : mesh.conn)
	{
		casted_faces.push_back(face[0]);
		casted_faces.push_back(face[1]);
		casted_faces.push_back(face[2]);
	}

	sample_magnetic_object(mesh.vertices, casted_faces, radius, magnetic_pos_samples);
}

void sample_magnetic_object(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int,3>>& faces, const double radius, std::vector<Eigen::Vector3d>& magnetic_pos_samples)
{
	std::vector<unsigned int> casted_faces;
	for (const auto& face : faces)
	{
		casted_faces.push_back(face[0]);
		casted_faces.push_back(face[1]);
		casted_faces.push_back(face[2]);
	}

	sample_magnetic_object(vertices, casted_faces, radius, magnetic_pos_samples);
}

void magnetic()
{
	// Parameters
	double mass = 0.1;
	double size = 0.1;
	double contact_thickness = 0.0005;


	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_implicit";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 15.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	stark::Simulation simulation(settings);


	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
		.set_default_contact_thickness(contact_thickness)
		//.set_friction_enabled(false)
		.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
		.set_min_contact_stiffness(1e12)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

	int n = 2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 1; j++) {
			for (int k = 0; k < 1; k++) {
				//auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_sphere("spheres", mass, 0.5*size, 3);
				auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass, size);
				handler.rigidbody.add_translation({ 1.5*i*size, 1.5*j*size, 1.5*k*size });
				auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody, stark::EnergyRigidBodyMagnetic::Params()
					.set_magnetic_permeability(mu_neodymium)
					.set_magnetic_material(MagneticMaterial::Permanent)
				);
				simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, { Eigen::Vector3d::Zero() }, 0.5*size, (double)((2*i)-1) * 10.0 * Eigen::Vector3d::UnitY());
			}
		}
	}

	// Run
	simulation.run();
}

void magnetic_distance_criterion()
{
	// Parameters
	double mass = 0.25;
	double size = 0.1;
	double contact_thickness = 0.0005;
	double magnetic_strength = 1000.0;

	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_distance_criterion";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_distance_criterion";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 10.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.max_time_step_size = 0.001;
	settings.newton.max_line_search_iterations = 100;
	settings.newton.max_newton_iterations = 1000;
	settings.simulation.use_adaptive_time_step = true;

	std::vector<double> thresholds = {
		10000.0 * size,
		10.0 * size,
		5.0 * size,
		2.0 * size,
		1.0 * size
	};

	for (auto thresh : thresholds)
	{
		stark::Simulation simulation(settings);


		simulation.interactions->contact->set_global_params(
			stark::EnergyFrictionalContact::GlobalParams()
				.set_default_contact_thickness(contact_thickness)
					//.set_friction_enabled(false)
				.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
				.set_min_contact_stiffness(
					1e8)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
		);


		int n = 2;
		for (int i = 0; i < n; i++)
		{
			//auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_sphere("spheres", mass, 0.5*size, 3);
			auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box_thresh_" + std::to_string(thresh), mass, size);
			handler.rigidbody.add_translation({2.0 * i - 1.0, 0.0, 0.0});
			auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																		  stark::EnergyRigidBodyMagnetic::Params()
																			  .set_magnetic_permeability(mu_neodymium)
																			  .set_magnetic_material(
																				  MagneticMaterial::Permanent)
			);
			simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, {Eigen::Vector3d::Zero()},
																   0.5 * size,
																   magnetic_strength * Eigen::Vector3d::UnitX());
			simulation.rigidbodies->magnetic->set_distance_threshold(thresh);
		}

		// Run
		simulation.run();
	}
}

void magnetic_coupling()
{
	// Parameters
	double mass = 0.1;
	double size = 0.1;
	double contact_thickness = 0.00005;


	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_coupling";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_coupling";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 1.5;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.max_time_step_size = 0.001;
	settings.simulation.use_adaptive_time_step = false;
	settings.newton.max_line_search_iterations = 100;
	settings.newton.max_newton_iterations = 1000;
	settings.output.fps = 240;

	stark::Simulation simulation(settings);


	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
			.set_friction_enabled(false)
			.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e14)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

    const double magnetic_strength = 14.2;

	int n = 2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 1; j++) {
			for (int k = 0; k < 1; k++) {
				//auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_sphere("spheres", mass, 0.5*size, 3);
				auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass, size);
				handler.rigidbody.add_translation({ 2.5*i*size, 0.0, 0.0 });
				auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody, stark::EnergyRigidBodyMagnetic::Params()
					.set_magnetic_permeability(mu_neodymium)
					.set_magnetic_material(MagneticMaterial::Permanent)
				);
				simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, { Eigen::Vector3d::Zero() }, 0.5*size, magnetic_strength*Eigen::Vector3d::UnitX());
				simulation.rigidbodies->magnetic->set_distance_threshold(1.01*size);
			}
		}
	}

	// Run
	simulation.run();
}

void magnetic_coupling_adaptive()
{
	// Parameters
	double mass = 0.25;
	double size = 0.1;
	double contact_thickness = 0.0005;
	double magnetic_strength = 1420.0;


	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_coupling_adaptive_time";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_coupling_adaptive_time";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 10.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.max_time_step_size = 0.001;
	settings.newton.max_line_search_iterations = 100;
	settings.newton.max_newton_iterations = 1000;
	settings.simulation.use_adaptive_time_step = false;
	stark::Simulation simulation(settings);


	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
			.set_friction_enabled(false)
			//.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e10)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);


	int n = 2;
	for (int i = 0; i < n; i++)
	{
		//auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_sphere("spheres", mass, 0.5*size, 3);
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass, size);
		handler.rigidbody.add_translation({2.0 * i - 1.0, 0.0, 0.0});
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, {Eigen::Vector3d::Zero()}, 0.5 * size,
															   magnetic_strength * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->magnetic->set_distance_threshold(1.1*size);
	}

	// Run
	simulation.run();
}

void magnetic_field_torque_test()
{
	// Parameters
	double mass = 0.1;
	double size = 0.1;
	double contact_thickness = 0.0005;
	double angular_damping = 0.99;


	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_torque_test";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_torque_test";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 5.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	stark::Simulation simulation(settings);


	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
				//.set_friction_enabled(false)
			.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e12)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

	const Eigen::Vector3d needle_size = {2.0 * size, 0.5 * size, 0.5 * size};
	const double radius = 0.25 * size;

	stark::Mesh<3> box = stark::make_box(needle_size);

	std::vector<unsigned int> casted_box_faces;
	for (const auto& face : box.conn)
	{
		casted_box_faces.push_back(face[0]);
		casted_box_faces.push_back(face[1]);
		casted_box_faces.push_back(face[2]);
	}

	std::vector<Eigen::Vector3d> magnetic_samples;
	Utilities::VolumeSampling::sampleMesh(box.vertices.size(), box.vertices.data(), box.conn.size(), casted_box_faces.data(), radius,
										  nullptr, {50, 50, 50}, false, 0, magnetic_samples);

	std::cout << "Created N=" << magnetic_samples.size() << " samples for the magnetic needles!" << std::endl;

	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(0.5 * Eigen::Vector3d::UnitX() - 0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   5 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(-0.5 * Eigen::Vector3d::UnitX() - 0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(-0.5 * Eigen::Vector3d::UnitX() + 0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(0.5 * Eigen::Vector3d::UnitX() + 0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(- 0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.add_translation(0.5 * Eigen::Vector3d::UnitZ());
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
															handler.rigidbody.get_translation());
		simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
	}
	{
		auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																					   needle_size);
		handler.rigidbody.set_angular_damping(angular_damping);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples/*{Eigen::Vector3d::Zero()}*/, radius,
															   1000000 * Eigen::Vector3d::UnitX());
		simulation.rigidbodies->add_constraint_fix(handler.rigidbody);
	}

	// Run
	simulation.run();
}

void magnetic_field_orientation()
{
	// Parameters
	double mass = 0.1;
	double size = 0.1;
	double contact_thickness = 0.0005;
	double angular_damping = 0.75;


	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_orientation";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_orientation";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 15.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.gravity = { 0.0, 0.0, 0.0 };
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	stark::Simulation simulation(settings);


	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
				//.set_friction_enabled(false)
			.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e12)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

	const Eigen::Vector3d needle_half_size = {0.025, 0.01, 0.05};
	const double radius = 0.01;

	auto needle_mesh = stark::load_obj(MODELS_PATH + "/CompassNeedle.obj");
	std::vector<Eigen::Vector3d> verts;
	std::vector<std::array<int, 3>> tris;
	stark::clean_triangle_mesh(verts, tris, needle_mesh[0].vertices, needle_mesh[0].conn);
	const auto needle_inertia = stark::inertia_tensor_box(mass, 2.0 * needle_half_size);

	std::vector<Eigen::Vector3d> needle_samples = {
		0.03 * Eigen::Vector3d::UnitZ(),
		0.01 * Eigen::Vector3d::UnitZ(),
		-0.01 * Eigen::Vector3d::UnitZ(),
		-0.03 * Eigen::Vector3d::UnitZ(),
	};

    const double magnetic_strength = 0.142;

	int n = 11;
	for (int xx = 0; xx < n; xx++)
	{
		for (int zz = 0; zz < n; zz++)
		{
			if (xx == std::floor(0.5*n) && zz == std::floor(0.5*n))
			{
				const Eigen::Vector3d bar_size = {0.25 ,0.02, 1.5};
				auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_box("box", mass,
																							   bar_size);
				auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																			  stark::EnergyRigidBodyMagnetic::Params()
																				  .set_magnetic_permeability(
																					  mu_neodymium)
																				  .set_magnetic_material(
																					  MagneticMaterial::Permanent)
				);
				std::vector<Eigen::Vector3d> samples;
				sample_magnetic_object(vertices, triangles, 0.99 * radius, samples);

				simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, samples, radius,
																	   50 * 1.42 *  Eigen::Vector3d::UnitZ());
				simulation.rigidbodies->add_constraint_fix(handler.rigidbody);
			}
			else if (!(xx == std::floor(0.5*n) && (zz == std::floor(0.5*n) || zz == std::floor(0.5*n) - 1  || zz == std::floor(0.5*n) + 1)))
			{
				const Eigen::Vector3d translation = 0.25 * Eigen::Vector3d{ 2.0 * xx - (n - 1), 0.0, 2.0 * zz - (n - 1)};

				//auto [vertices, triangles, handler] = simulation.presets->rigidbodies->add_sphere("spheres", mass, 0.5*size, 3);
				auto handler = simulation.presets->rigidbodies->add("needles", 1.0, needle_inertia, verts, tris);
				handler.rigidbody.add_translation(translation);
				handler.rigidbody.set_angular_damping(angular_damping);
				auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																			  stark::EnergyRigidBodyMagnetic::Params()
																				  .set_magnetic_permeability(
																					  mu_neodymium)
																				  .set_magnetic_material(
																					  MagneticMaterial::Permanent)
				);
				simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, needle_samples, radius,
																	   magnetic_strength * Eigen::Vector3d::UnitZ());
				simulation.rigidbodies->add_constraint_global_point(handler.rigidbody,
																	handler.rigidbody.get_translation());
				simulation.rigidbodies->add_constraint_global_direction(handler.rigidbody, Eigen::Vector3d::UnitY());
			}
		}
	}

	//simulation.rigidbodies->magnetic->set_distance_threshold(0.0);

	// Run
	simulation.run();
}
void magnetic_perpetuum_mobile()
{
	const double contact_thickness = 0.0005;

	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_perpetuum_mobile";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_perpetuum_mobile";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 30.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.gravity = {0.0, -9.81, 0.0};
	settings.simulation.max_time_step_size = 0.001;
	stark::Simulation simulation(settings);

	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
			.set_friction_enabled(true)
			.set_friction_stick_slide_threshold(0.001)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e7)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

	stark::RigidBody::Handler ball_handler;
	stark::RigidBody::Handler base_handler;
	stark::RigidBody::Handler track_handler;
	stark::EnergyRigidBodyMagnetic::Handler electro_magnet_handler;

	// Ballin'
	{
		auto handler = simulation.presets->rigidbodies->add_sphere("ball", 0.5, 0.05, 2);
		//handler.handler.rigidbody.set_translation(1.0 * Eigen::Vector3d::UnitY());
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_nickel)
																		  .set_magnetic_material(
																			  MagneticMaterial::Linear)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler,
															   {handler.handler.rigidbody.get_translation()}, 0.05,
															   1 * Eigen::Vector3d::UnitX());

		ball_handler = handler.handler;
	}
	// Base
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Perpetuum_Mobile/Base.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_box(1.0, 1.0);
		auto handler = simulation.presets->rigidbodies->add("fix", 1.0, dummy_inertia, verts, tris);
		simulation.rigidbodies->add_constraint_fix(handler.rigidbody);

		base_handler = handler;
	}
	// Tracks
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Perpetuum_Mobile/Tracks.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_box(1.0, 1.0);
		auto handler = simulation.presets->rigidbodies->add("fix", 1.0, dummy_inertia, verts, tris);
		simulation.rigidbodies->add_constraint_fix(handler.rigidbody);

		track_handler = handler;
	}
	// Electromagnet
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Perpetuum_Mobile/ElectroMagnet.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_box(1.0, 1.0);
		auto handler = simulation.presets->rigidbodies->add("fix", 1.0, dummy_inertia, verts, tris);
		simulation.rigidbodies->add_constraint_fix(handler.rigidbody);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_neodymium)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
		);

		std::vector<Eigen::Vector3d> samples;
		const double magnetic_sample_radius = 0.01;
		sample_magnetic_object(verts, tris, magnetic_sample_radius, samples);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, samples, magnetic_sample_radius,
															   3.25 * Eigen::Vector3d::UnitY());

		electro_magnet_handler = magnetic_handler;
	}

	// Friction
	simulation.interactions->contact->set_friction(ball_handler.contact, base_handler.contact, 0.1);
	simulation.interactions->contact->set_friction(ball_handler.contact, track_handler.contact, 0.001);

	simulation.run([&]() {
		const Eigen::Vector3d min_point(-0.1, -0.5, -0.2);
		const Eigen::Vector3d max_point(0.4, -0.1, 0.2);
		Eigen::AlignedBox3d aabb(min_point, max_point);

		if (aabb.contains(ball_handler.rigidbody.get_translation()))
		{
			electro_magnet_handler.set_params(stark::EnergyRigidBodyMagnetic::Params()
												  .set_magnetic_permeability(
													  mu_neodymium)
												  .set_magnetic_material(
													  MagneticMaterial::Permanent)
												  .set_enabled(true));
		}
		else
		{
			electro_magnet_handler.set_params(stark::EnergyRigidBodyMagnetic::Params()
												  .set_magnetic_permeability(
													  mu_neodymium)
												  .set_magnetic_material(
													  MagneticMaterial::Permanent)
												  .set_enabled(false));
		}

	});
}

void magnetic_sphere_pile()
{
	const double contact_thickness = 0.0005;

	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_sphere_pile";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_sphere_pile";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 25.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.gravity = {0.0, 0.0, -9.81};
	settings.simulation.max_time_step_size = 0.001;
	stark::Simulation simulation(settings);

	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
			.set_default_contact_thickness(contact_thickness)
			.set_friction_enabled(false)
			.set_friction_stick_slide_threshold(0.01)  // Warning: very low quality friction. Better value 0.001
			.set_min_contact_stiffness(1e10)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
	);

	stark::ContactHandler boundary_contact_handler;
	std::vector<stark::ContactHandler> object_contact_handler;

	const double damping = 0.99;

	struct MagneticObject {
		stark::Mesh<3> mesh;
		Eigen::Matrix3d inertia_tensor;
		std::vector<Eigen::Vector3d> magnetic_samples;
		double magnetic_sphere_radius;
	} bunny, armadillo, dragon, star;

	const double mass = 0.2;
	const double bounding_sphere_radius = 0.08;
	const double magnetic_sphere_radius = 0.03;

	// Bunny
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Bunny.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_sphere(mass, bounding_sphere_radius);

		bunny.mesh = {verts, tris};
		bunny.inertia_tensor = dummy_inertia;
		bunny.magnetic_samples = { Eigen::Vector3d::Zero() };
		bunny.magnetic_sphere_radius = 0.029;
	}

	// Armadillo
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Armadillo.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_cylinder(mass, bounding_sphere_radius, 0.15);

		armadillo.mesh = {verts, tris};
		armadillo.inertia_tensor = dummy_inertia;
		armadillo.magnetic_samples = { Eigen::Vector3d::Zero() };
		armadillo.magnetic_sphere_radius = 0.015;
	}

	// Star
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/JackStar2.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_sphere(mass, bounding_sphere_radius);

		star.mesh = {verts, tris};
		star.inertia_tensor = dummy_inertia;
		star.magnetic_samples = { Eigen::Vector3d::Zero() };
		star.magnetic_sphere_radius = 0.008;
	}

	// Dragon
	{
		const auto mesh = stark::load_obj(MODELS_PATH + "/Dragon.obj");
		std::vector<Eigen::Vector3d> verts;
		std::vector<std::array<int, 3>> tris;
		stark::clean_triangle_mesh(verts, tris, mesh[0].vertices, mesh[0].conn);
		const auto dummy_inertia = stark::inertia_tensor_sphere(mass, bounding_sphere_radius);

		dragon.mesh = {verts, tris};
		dragon.inertia_tensor = dummy_inertia;
		dragon.magnetic_samples = { Eigen::Vector3d::Zero() };
		dragon.magnetic_sphere_radius = 0.008;
	}


	const unsigned int N = 3;
	const double distance = 0.2;

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				Eigen::Vector3d position = Eigen::Vector3d(distance * i - std::floor(N / 2.) * distance, distance * j - std::floor(N / 2.) * distance, distance * k - std::floor(N / 2.) * distance);

				stark::RigidBody::Handler handler;

				int randomNumber = rand() % 2;
				if (randomNumber == 0)
				{
					handler = simulation.presets->rigidbodies->add("objects", mass, bunny.inertia_tensor, bunny.mesh.vertices, bunny.mesh.conn);
					auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																				  stark::EnergyRigidBodyMagnetic::Params()
																					  .set_magnetic_permeability(
																						  mu_iron)
																					  .set_magnetic_material(
																						  MagneticMaterial::Linear)
					);
					simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler,
																		   {Eigen::Vector3d::Zero()}, magnetic_sphere_radius);
				}
				else if (randomNumber == 1)
				{
					handler = simulation.presets->rigidbodies->add("objects", mass, armadillo.inertia_tensor, armadillo.mesh.vertices, armadillo.mesh.conn);
					auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																				  stark::EnergyRigidBodyMagnetic::Params()
																					  .set_magnetic_permeability(
																						  mu_iron)
																					  .set_magnetic_material(
																						  MagneticMaterial::Linear)
					);
					simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler,
																		   {Eigen::Vector3d::Zero()}, magnetic_sphere_radius);
				}
				/*else if (randomNumber == 2)
				{
					handler = simulation.presets->rigidbodies->add("objects", mass, dragon.inertia_tensor, dragon.mesh.vertices, dragon.mesh.conn);
					auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.rigidbody,
																				  stark::EnergyRigidBodyMagnetic::Params()
																					  .set_magnetic_permeability(
																						  mu_nickel)
																					  .set_magnetic_material(
																						  MagneticMaterial::Linear)
					);
					simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler,
																		   {Eigen::Vector3d::Zero()}, magnetic_sphere_radius);
				}*/

				handler.rigidbody.set_translation(position);
				handler.rigidbody.set_rotation(Eigen::Quaterniond::UnitRandom());

				handler.rigidbody.set_linear_damping(damping);
				handler.rigidbody.set_angular_damping(damping);

				object_contact_handler.push_back(handler.contact);
			}
		}
	}

	{
		auto handler = simulation.presets->rigidbodies->add_sphere("boundary", 100., 0.8, 3);
		simulation.rigidbodies->add_constraint_fix(handler.handler.rigidbody);

		boundary_contact_handler = handler.handler.contact;
	}

	stark::RBCFixHandler magnet_translation_handler;
	stark::EnergyRigidBodyMagnetic::Handler magnet_activation_handler;
	const Eigen::Vector3d magnet_pos = -1.25 * Eigen::Vector3d::UnitZ();
	const double magnetic_strength = 12500.;

	{
		auto handler = simulation.presets->rigidbodies->add_cylinder("magnets", 100., 0.1, 0.8);
		handler.handler.rigidbody.set_translation(magnet_pos);
		auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.handler.rigidbody,
																	  stark::EnergyRigidBodyMagnetic::Params()
																		  .set_magnetic_permeability(
																			  mu_nickel)
																		  .set_magnetic_material(
																			  MagneticMaterial::Permanent)
																		  .set_enabled(false)
		);
		simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, {-0.3*Eigen::Vector3d::UnitZ(), -0.1 * Eigen::Vector3d::UnitZ(), 0.1 * Eigen::Vector3d::UnitZ(), 0.3*Eigen::Vector3d::UnitZ()}, 0.1, magnetic_strength * Eigen::Vector3d::UnitZ());

		magnet_activation_handler = magnetic_handler;
		magnet_translation_handler = simulation.rigidbodies->add_constraint_fix(handler.handler.rigidbody);
	}

	/*const double friction = 2.0;
	for (auto obj : object_contact_handler) {
		simulation.interactions->contact->set_friction(boundary_contact_handler, obj, friction);
	}*/

	simulation.rigidbodies->magnetic->set_distance_threshold(0.5);

	// Events
	std::vector<double> times = {
		2.0,
		6.0,
		10.0,
		14.0,
		16.0,
		20.0
	};

	simulation.add_time_event(times[0], times[1], [&](const double t){
		const double blend_val = stark::blend(0.0, 90.0, times[0], times[1], t, stark::BlendType::EaseInOut);
	 	const Eigen::Vector3d new_pos = 1.25 * Eigen::Vector3d(0.0,sin(M_PI / 180.0 * blend_val),-cos(M_PI / 180.0 * blend_val));
		magnet_translation_handler.set_transformation(new_pos, blend_val, Eigen::Vector3d::UnitX());
		magnet_activation_handler.set_params(stark::EnergyRigidBodyMagnetic::Params().set_magnetic_permeability(mu_nickel).set_magnetic_material(MagneticMaterial::Permanent).set_enabled(true));
	});
	simulation.add_time_event(times[1], times[2], [&](const double t){
		const double blend_val = stark::blend(0.0, 180.0, times[1], times[2], t, stark::BlendType::EaseInOut);
		const Eigen::Vector3d new_pos = 1.25 * Eigen::Vector3d(-sin(M_PI / 180.0 * blend_val),cos(M_PI / 180.0 * blend_val),0.0);
		const Eigen::Matrix3d rot = Eigen::AngleAxisd(M_PI / 180.0 * blend_val, Eigen::Vector3d::UnitZ()).toRotationMatrix() * Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d::UnitX()).toRotationMatrix();
		magnet_translation_handler.set_transformation(new_pos, rot);
	});
	simulation.add_time_event(times[2], times[3], [&](const double t){
		const double blend_val = stark::blend(0.0, 90.0, times[2], times[3], t, stark::BlendType::EaseInOut);
		const Eigen::Vector3d new_pos = 1.25 * Eigen::Vector3d(0.0,-cos(M_PI / 180.0 * blend_val),sin(M_PI / 180.0 * blend_val));
		const Eigen::Matrix3d rot = Eigen::AngleAxisd(-M_PI / 180.0 * blend_val, Eigen::Vector3d::UnitX()).toRotationMatrix() * Eigen::AngleAxisd(-M_PI / 2.0, Eigen::Vector3d::UnitX()).toRotationMatrix();
		magnet_translation_handler.set_transformation(new_pos, rot);
	});
	simulation.add_time_event(times[4], times[5], [&](const double t){
		magnet_activation_handler.set_params(stark::EnergyRigidBodyMagnetic::Params().set_magnetic_permeability(mu_nickel).set_magnetic_material(MagneticMaterial::Permanent).set_enabled(false));
	});

	simulation.run();
}



std::vector<Eigen::Vector3d> generatePointsInCylinder(float r, float h, float d) {
	// Calculate the number of layers and points per layer
	int layers = std::ceil(h / d);
	int pointsPerLayer = std::ceil(2 * M_PI * r / d);

	// Create a matrix to store the points
	std::vector<Eigen::Vector3d> points(layers * pointsPerLayer);
	float dz = h / (layers - 1);  // Distance between layers

	// Loop over each layer
	for (int i = 0; i < layers; ++i) {
		float z = -h / 2 + i * dz;  // Current layer's z-coordinate
		// Loop over each point in the layer
		for (int j = 0; j < pointsPerLayer; ++j) {
			float angle = 2 * M_PI * j / pointsPerLayer;  // Angle for the current point
			float x = r * cos(angle);  // x-coordinate
			float y = r * sin(angle);  // y-coordinate
			// Store the point in the matrix
			points[i * pointsPerLayer + j][0] = x;
			points[i * pointsPerLayer + j][1] = y;
			points[i * pointsPerLayer + j][2] = z;
		}
	}

	return points;
}

void magnetic_crane()
{
	bool script_only = false;

	stark::Settings settings = stark::Settings();
	settings.output.simulation_name = "magnetic_crane";
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_crane";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.execution.end_simulation_time = 22.0;
	settings.debug.symx_check_for_NaNs = true;

	//settings.simulation.max_time_step_size = 0.002;
	if (script_only) {
		settings.simulation.init_frictional_contact = false;
	}

	stark::Simulation simulation(settings);
	double damping = 2.0;

	// Contact
	simulation.interactions->contact->set_global_params(
		stark::EnergyFrictionalContact::GlobalParams()
		//.set_default_contact_thickness(0.00025)
		.set_default_contact_thickness(0.001)
		.set_min_contact_stiffness(1e7)
		//.set_friction_stick_slide_threshold(0.0025)
		.set_friction_stick_slide_threshold(0.1)
		//.set_friction_enabled(false)
	);

	// Add floor
	auto [floor_vertices, floor_triangles, floor] = simulation.presets->rigidbodies->add_box("floor", 1.0, { 10.0, 10.0, 0.01 });
	floor.rigidbody.add_translation({ 0.0, 0.0, -0.01 });
	simulation.rigidbodies->add_constraint_fix(floor.rigidbody);

	// Add containers
	auto container_mesh_raw = stark::load_obj(MODELS_PATH + "/box_no_lid.obj")[0];
	auto [container_vertices, container_triangles] = stark::clean_triangle_mesh(container_mesh_raw.vertices, container_mesh_raw.conn);
	stark::scale(container_vertices, { 1.0, 1.0, 0.5 });

	// Starting container
	auto container = simulation.presets->rigidbodies->add("container", 1.0, stark::inertia_tensor_box(1.0, { 1.0, 1.0, 0.5 }), container_vertices, container_triangles);
	container.rigidbody.add_translation({ 0.0, 0.0, 0.25 });
	simulation.rigidbodies->add_constraint_fix(container.rigidbody);

	// Final container
	auto container_final = simulation.presets->rigidbodies->add("container", 1.0, stark::inertia_tensor_box(1.0, { 1.0, 1.0, 0.5 }), container_vertices, container_triangles);
	container_final.rigidbody.add_translation({ 0.0, 0.0, 0.25 });
	container_final.rigidbody.add_rotation(60.0, Eigen::Vector3d::UnitZ(), { -2.8549, 0.0, 0.5 });  // Crane position
	simulation.rigidbodies->add_constraint_fix(container_final.rigidbody);

	// Add objects
	std::vector<stark::RigidBody::Handler> rbs;
	double obj_mass = 0.02;
	if (!script_only) {
		double size = 0.1;
		std::array<int, 3> grid = { 5, 5, 5 };
		double spacing = size * 1.8;
		double height = 0.5;
		Eigen::Vector3d center = { 0.5 * (grid[0] - 1) * spacing, 0.5 * (grid[1] - 1) * spacing, 0.5 * (grid[2] - 1) * spacing };
		for (int i = 0; i < grid[0]; i++) {
			for (int j = 0; j < grid[1]; j++) {
				for (int k = 0; k < grid[2]; k++) {
					auto [V, C, obj] = simulation.presets->rigidbodies->add_box("object", obj_mass, size);
					auto mag = simulation.rigidbodies->magnetic->add(obj.rigidbody, stark::EnergyRigidBodyMagnetic::Params()
						.set_magnetic_permeability(mu_nickel)
						.set_magnetic_material(MagneticMaterial::Linear));
					simulation.rigidbodies->magnetic->add_magnetic_samples(mag, {Eigen::Vector3d::Zero()}, 0.5 * size, Eigen::Vector3d::Zero());
					obj.rigidbody.add_rotation(Eigen::Vector3d::Random().x() * 90.0, Eigen::Vector3d::Random());
					obj.rigidbody.add_translation({ i * spacing - center.x(), j * spacing - center.y(), k * spacing - center.z() + height });
					obj.rigidbody.set_linear_damping(damping);
					obj.rigidbody.set_angular_damping(damping);
					rbs.push_back(obj);
				}
			}
		}
	}

	// Crane
	auto [crane_vertices, crane_triangles, crane] = simulation.presets->rigidbodies->add_box("crane", 1.0, 1.0);
	Eigen::Vector3d crane_position = { -2.8549, 0.0, 0.5 };
	crane.rigidbody.add_translation(crane_position);
	auto crane_fix = simulation.rigidbodies->add_constraint_fix(crane.rigidbody);

	// Magnet
	double magnet_strength = 1.42 * 500.0;
	Eigen::Vector3d magnet_center = { -2.8649, -2.85459, 1.273 };

	auto magnet_mesh = stark::load_obj(MODELS_PATH + "/CylinderMagnet.obj");
	std::vector<Eigen::Vector3d> verts;
	std::vector<std::array<int, 3>> tris;
	stark::clean_triangle_mesh(verts, tris, magnet_mesh[0].vertices, magnet_mesh[0].conn);
	stark::rotate_deg(verts, 90.0, Eigen::Vector3d::UnitX());

	Eigen::Matrix3d magnet_inertia = stark::inertia_tensor_cylinder(1.0, 0.25, 0.2);

	auto magnet = simulation.presets->rigidbodies->add("magnet", 1.0, magnet_inertia, verts, tris);
	std::vector<Eigen::Vector3d> magnetic_samples;
	sample_magnetic_object(verts, tris, 0.05, magnetic_samples);
	//std::vector<Eigen::Vector3d> magnetic_samples = generatePointsInCylinder(0.25, 0.2, 0.1);

	magnet.rigidbody.add_translation(magnet_center);
	magnet.rigidbody.set_linear_damping(damping);
	magnet.rigidbody.set_angular_damping(damping);
	auto crane_magnetic_handler = simulation.rigidbodies->magnetic->add(magnet.rigidbody, stark::EnergyRigidBodyMagnetic::Params()
		.set_magnetic_permeability(mu_nickel)
		.set_magnetic_material(MagneticMaterial::Permanent)
		.set_enabled(false));
	simulation.rigidbodies->magnetic->add_magnetic_samples(crane_magnetic_handler, magnetic_samples, 0.05, magnet_strength * Eigen::Vector3d::UnitZ());


	// Cable
	auto cable_material = stark::Line::Params::Elastic_Rubberband();
	cable_material.strain.elasticity_only = false;
	cable_material.strain.youngs_modulus = 1e5;
	cable_material.strain.section_radius = 0.01;
	cable_material.strain.strain_limit = 0.01;
	cable_material.strain.strain_limit_stiffness = 1e8;

	int n = 50;
	auto [cable_vertices, cable_segments, cable] = simulation.presets->deformables->add_line_as_segments("cable", magnet_center + (0.1 + 0.96)*Eigen::Vector3d::UnitZ(), magnet_center + 0.1*Eigen::Vector3d::UnitZ(), n, cable_material);
	simulation.interactions->attachments->add(crane.rigidbody, cable.point_set, { 0 }, stark::EnergyAttachments::Params().set_stiffness(1e7));
	simulation.interactions->attachments->add(magnet.rigidbody, cable.point_set, { n }, stark::EnergyAttachments::Params().set_stiffness(1e7));
	cable.contact.disable_collision(crane.contact);
	cable.contact.disable_collision(magnet.contact);

	// Friction
	double friction = 0.5;
	for (int i = 0; i < (int)rbs.size(); i++) {
		simulation.interactions->contact->set_friction(rbs[i].contact, container.contact, friction);
		simulation.interactions->contact->set_friction(rbs[i].contact, container_final.contact, friction);
		simulation.interactions->contact->set_friction(rbs[i].contact, magnet.contact, friction);
		simulation.interactions->contact->set_friction(rbs[i].contact, floor.contact, friction);
		for (int j = i + 1; j < (int)rbs.size(); j++) {
			simulation.interactions->contact->set_friction(rbs[i].contact, rbs[j].contact, friction);
		}
	}

	// Script
	//// Times
	std::vector<double> durations =
	{
		2.0,  // [0] Objects fall
		3.0,  // [1] Crane moves
		2.0,  // [2] Wait
		1.0,  // [3] Cable extends
		2.0,  // [4] Magnet is on and wait
		1.0,  // [5] Cable retracts
		2.0,  // [6] Wait
		3.0,  // [7] Crane moves
		2.0,  // [8] Wait
		3.0,  // [9] Magnet is off and end
	};
	std::vector<double> times = { 0.0 };
	for (double duration : durations) {
		times.push_back(times.back() + duration);
	}
	double magnet_starts = times[4];
	double magnet_ends = times[9];

	//// Events
	simulation.add_time_event(times[1], times[2], [&](double t) { crane_fix.set_transformation(crane_position, stark::blend(0.0, 90.0, times[1], times[2], t, stark::BlendType::Linear), Eigen::Vector3d::UnitZ()); });
	simulation.add_time_event(times[3], times[4], [&](double t) { cable.strain.set_params(cable_material.strain.set_strain_limit(stark::blend(0.01, 0.75, times[3], times[4], t, stark::BlendType::Linear))); });
	simulation.add_time_event(times[5], times[6], [&](double t) { cable.strain.set_params(cable_material.strain.set_strain_limit(stark::blend(0.75, 0.01, times[5], times[6], t, stark::BlendType::Linear))); });
	simulation.add_time_event(times[7], times[8], [&](double t) { crane_fix.set_transformation(crane_position, stark::blend(90.0, 150.0, times[7], times[8], t, stark::BlendType::Linear), Eigen::Vector3d::UnitZ()); });
	
	simulation.add_time_event(magnet_starts, magnet_ends,
		[&](double t) {
			// LUKAS: The magnets is turned on in this function
			crane_magnetic_handler.set_params(stark::EnergyRigidBodyMagnetic::Params()
												  .set_magnetic_permeability(
													  mu_neodymium)
												  .set_magnetic_material(
													  MagneticMaterial::Permanent)
												  .set_enabled(true));
		}
	);
	simulation.add_time_event(magnet_ends, times.back(),
		[&](double t) {
			// LUKAS: The magnets is turned off in this function
			crane_magnetic_handler.set_params(stark::EnergyRigidBodyMagnetic::Params()
												  .set_magnetic_permeability(
													  mu_neodymium)
												  .set_magnetic_material(
													  MagneticMaterial::Permanent)
												  .set_enabled(false));
		}
	);

	// Run
	simulation.rigidbodies->magnetic->set_distance_threshold(0.4);
	simulation.run(times.back());
}

void magnetic_optimization_timings_old()
{
	const double contact_thickness = 0.0005;
	omp_set_num_threads(1);

	stark::Settings settings = stark::Settings();
	settings.output.output_directory = OUTPUT_PATH + "/magnetic_optimization_timings";
	settings.output.codegen_directory = COMPILE_PATH;
	settings.output.fps = 1;
	settings.execution.end_simulation_time = 1.0;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.gravity = {0.0, 0.0, -9.81};
	settings.simulation.max_time_step_size = 0.001;
	settings.execution.n_threads = 1;

	const unsigned int num_stacked_samples = 1;

	for (unsigned int i = 0; i < 100; i++)
	{
		const double threshold = 10.0 * (100 - i) / 100. + 0.1;
		settings.output.simulation_name = "magnetic_optimization_timings_" + std::to_string(threshold);
		stark::Simulation simulation(settings);

		simulation.interactions->contact->set_global_params(
			stark::EnergyFrictionalContact::GlobalParams()
				.set_default_contact_thickness(contact_thickness)
				.set_friction_enabled(false)
				.set_friction_stick_slide_threshold(0.01)  // Warning: very low quality friction. Better value 0.001
				.set_min_contact_stiffness(1e10)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
		);

		// Central permanent magnet
		{
			auto handler = simulation.presets->rigidbodies->add_box("permanent", 1.0, 0.1);
			simulation.rigidbodies->add_constraint_fix(handler.handler.rigidbody);
			auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.handler.rigidbody,
																		  stark::EnergyRigidBodyMagnetic::Params()
																			  .set_magnetic_permeability(
																				  mu_neodymium)
																			  .set_magnetic_material(
																				  MagneticMaterial::Permanent)
																			  .set_enabled(true));
			std::vector<Eigen::Vector3d> magnetic_samples;
			for (unsigned int k = 0; k < num_stacked_samples; k++)
				magnetic_samples.push_back(Eigen::Vector3d::Zero());
			simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples, 0.05,
																   Eigen::Vector3d::UnitX());
		}

		// Magnetizable bar
		{
			auto handler = simulation.presets->rigidbodies->add_box("permanent", 1.0, {10.0, 0.1, 0.1});
			handler.handler.rigidbody.set_translation(Eigen::Vector3d::UnitX() * (5.0 + 0.1));
			simulation.rigidbodies->add_constraint_fix(handler.handler.rigidbody);
			auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.handler.rigidbody,
																		  stark::EnergyRigidBodyMagnetic::Params()
																			  .set_magnetic_permeability(
																				  mu_neodymium)
																			  .set_magnetic_material(
																				  MagneticMaterial::Permanent)
																			  .set_enabled(true));
			std::vector<Eigen::Vector3d> magnetic_samples;
			for (unsigned int j = 0; j < 100; j++)
				for (unsigned int k = 0; k < num_stacked_samples; k++)
					magnetic_samples.push_back({
						-4.95 + 0.1 * j, 0.0, 0.0
					});

			simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, magnetic_samples, 0.05,
																   Eigen::Vector3d::UnitX());
		}

		simulation.rigidbodies->magnetic->set_distance_threshold(threshold);
		simulation.run();
	}
}

void magnetic_optimization_timings()
{
	const double contact_thickness = 0.0005;
	//omp_set_num_threads(1);

	stark::Settings settings = stark::Settings();
	settings.output.codegen_directory = COMPILE_PATH;
	settings.output.fps = 1;
	settings.execution.end_simulation_time = 0.01;
	settings.debug.symx_check_for_NaNs = true;
	settings.simulation.magnetic_method = stark::MagneticMethod::DipoleMoment;
	settings.simulation.gravity = {0.0, -9.81, 0.0};
	settings.simulation.max_time_step_size = 0.001;
	//settings.execution.n_threads = 1;

	// Bowl mesh
	const auto bowl_mesh = stark::load_obj(MODELS_PATH + "/Bowl.obj");
	std::vector<Eigen::Vector3d> bowl_verts;
	std::vector<std::array<int, 3>> bowl_tris;
	stark::clean_triangle_mesh(bowl_verts, bowl_tris, bowl_mesh[0].vertices, bowl_mesh[0].conn);
	const auto bowl_dummy_inertia = stark::inertia_tensor_sphere(100.0, 1.0);

	const auto sphere_pos_mesh = stark::load_obj(MODELS_PATH + "/SpherePattern2.obj");
	std::vector<Eigen::Vector3d> sphere_pos = sphere_pos_mesh[0].vertices;

	for (unsigned int i = 0; i < 125; i++)
	{
		const double threshold = (125. - i) / 100. + 0.1;
		settings.output.output_directory = OUTPUT_PATH + "/magnetic_optimization_timings_" + std::to_string(threshold);
		settings.output.simulation_name = "magnetic_optimization_timings_" + std::to_string(threshold);
		stark::Simulation simulation(settings);

		simulation.interactions->contact->set_global_params(
			stark::EnergyFrictionalContact::GlobalParams()
				.set_default_contact_thickness(contact_thickness)
				.set_friction_enabled(false)
				.set_friction_stick_slide_threshold(0.01)  // Warning: very low quality friction. Better value 0.001
				.set_min_contact_stiffness(1e10)  // Note: In this branch, there is no adaptive contact stiffness. Just use a high value all the time.
		);

		// Bowl
		{
			auto handler = simulation.presets->rigidbodies->add("bowl", 100., bowl_dummy_inertia, bowl_verts, bowl_tris);
			simulation.rigidbodies->add_constraint_fix(handler.rigidbody);
		}

		// Balls
		for (unsigned int j = 0; j < sphere_pos.size(); j++)
		{
			auto handler = simulation.presets->rigidbodies->add_sphere("sphere", 1.0, 0.09, 3);
			handler.handler.rigidbody.add_translation(sphere_pos[j]);
			auto magnetic_handler = simulation.rigidbodies->magnetic->add(handler.handler.rigidbody,
																		  stark::EnergyRigidBodyMagnetic::Params()
																			  .set_magnetic_permeability(
																				  mu_neodymium)
																			  .set_magnetic_material(
																				  MagneticMaterial::Permanent)
																			  .set_enabled(true));
			simulation.rigidbodies->magnetic->add_magnetic_samples(magnetic_handler, {Eigen::Vector3d::Zero()}, 0.05,
																   0.01 * Eigen::Vector3d::UnitY());
		}

		simulation.rigidbodies->magnetic->set_distance_threshold(threshold);
		simulation.run();
	}
}



int main()
{
    // This will run all scenes of the paper consecutively. Comment out unwanted scenes ;)
	magnetic_distance_criterion();
	magnetic_coupling();
	magnetic_coupling_adaptive();
	magnetic_perpetuum_mobile();
	magnetic_field_orientation();
	magnetic_sphere_pile();
	magnetic_crane();
	magnetic_optimization_timings();
}
