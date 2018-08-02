#include "casmutils/stage.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
#include <casm/app/AppIO.hh>
#include <casm/casm_io/Log.hh>
#include <casm/casm_io/jsonParser.hh>
#include <casm/clex/ConfigDoF.hh>
#include <casm/clex/ConfigEnumInterpolation.hh>
#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/Configuration.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/clex/Supercell.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/strain/StrainConverter.hh>
#include <fstream>
#include <set>

namespace
{
CASM::ConfigDoF _sym_break_projection(const CASM::Supercell& scel, const CASM::ConfigDoF& mapped_dof)
{
    CASM::ConfigDoF ret_dof = mapped_dof;
    CASM::ConfigDoF ideal_dof = mapped_dof;
    Eigen::Matrix3d tmp_def = CASM::StrainConverter::right_stretch_tensor(mapped_dof.deformation());
    ideal_dof.set_deformation(tmp_def);
    Eigen::Matrix3d sym_strain = Eigen::Matrix3d::Zero();
    Eigen::MatrixXd sym_disp =
        Eigen::MatrixXd::Zero(mapped_dof.displacement().rows(), mapped_dof.displacement().cols());
    auto it = scel.permute_begin();
    for (auto it = scel.permute_begin(); it != scel.permute_end(); ++it)
    {
        CASM::ConfigDoF tmp_dof = copy_apply(it, ideal_dof);
        sym_strain = sym_strain + tmp_dof.deformation() / std::distance(scel.permute_begin(), scel.permute_end());
        sym_disp = sym_disp + tmp_dof.displacement() / std::distance(scel.permute_begin(), scel.permute_end());
    }
    ret_dof.set_deformation(sym_strain);
    ret_dof.set_displacement(sym_disp);
    CASM::Configuration tconfig(scel, CASM::jsonParser(), ret_dof);
    CASM::Supercell prim_scel(&(scel.primclex()), scel.primclex().prim().lattice());
    tconfig = sub_configuration(prim_scel, tconfig);
    Eigen::Matrix3d tiny_sym_strain = Eigen::Matrix3d::Zero();
    Eigen::MatrixXd tiny_sym_disp = Eigen::MatrixXd::Zero(tconfig.displacement().rows(), tconfig.displacement().cols());
    for (auto it = prim_scel.permute_begin(); it != prim_scel.permute_end(); ++it)
    {
        CASM::ConfigDoF tmp_dof = copy_apply(it, tconfig.configdof());
        tiny_sym_strain =
            tiny_sym_strain + tmp_dof.deformation() / std::distance(prim_scel.permute_begin(), prim_scel.permute_end());
        tiny_sym_disp =
            tiny_sym_disp + tmp_dof.displacement() / std::distance(prim_scel.permute_begin(), prim_scel.permute_end());
    }
    tconfig.set_deformation(tiny_sym_strain);
    tconfig.set_displacement(tiny_sym_disp);
    CASM::ConfigDoF sym_dof = tconfig.fill_supercell(scel).configdof();
    Eigen::Matrix3d sym_hencky = CASM::StrainConverter::hencky(sym_dof.deformation());
    Eigen::Matrix3d total_hencky = CASM::StrainConverter::hencky(ideal_dof.deformation());
    Eigen::Matrix3d sym_break_hencky = total_hencky - sym_hencky;
    ret_dof.set_deformation(CASM::StrainConverter::hencky_to_F(sym_break_hencky));
    ret_dof.set_displacement(mapped_dof.displacement() - sym_dof.displacement());
    return ret_dof;
}

Rewrap::Structure _apply_configdof(const CASM::Supercell& scel, const CASM::ConfigDoF& dof)
{
    return CASM::make_deformed_struc(CASM::Configuration(scel, CASM::jsonParser(), dof));
}

void _rigid_rotate(Rewrap::Structure* struc, Eigen::Matrix3d& rot_mat)
{
    struc->set_lattice(CASM::Lattice(rot_mat * struc->lattice().lat_column_mat() * rot_mat.transpose()), CASM::FRAC);
    return;
}



} // namespace
std::pair<Rewrap::Structure, Eigen::Matrix3i> _minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                                             const Rewrap::Structure& deformed_struc)
{
    CASM::PrimClex* pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, 0.5);
    CASM::ConfigDoF mapped_dof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_dof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_pair(_apply_configdof(scel, _sym_break_projection(scel, mapped_dof)), scel.transf_mat());
}

Rewrap::Structure minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                const Rewrap::Structure& deformed_struc)
{
    return _minimally_distorted_structure(ref_struc, deformed_struc).first;
}

std::pair<Rewrap::Structure, Eigen::Matrix3i> _distorted_structure(const Rewrap::Structure& ref_struc,
                                                                   const Rewrap::Structure& deformed_struc)
{
    CASM::PrimClex* pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, 0.5);
    CASM::ConfigDoF mapped_configdof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_configdof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_pair(_apply_configdof(scel, mapped_configdof), scel.transf_mat());
}
Rewrap::Structure distorted_structure(const Rewrap::Structure& ref_struc, const Rewrap::Structure& deformed_struc)
{
    return _distorted_structure(ref_struc, deformed_struc).first;
}
std::vector<Rewrap::Structure> interpolate(const Rewrap::Structure& init_struc, const Rewrap::Structure& final_struc,
                                           int n_images)
{
    std::cout << "initial lattice " << init_struc.lattice().lat_column_mat() << std::endl;
    std::cout << "final lattice " << final_struc.lattice().lat_column_mat() << std::endl;
    std::vector<Rewrap::Structure> images;
    CASM::PrimClex* pclex = new CASM::PrimClex(init_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, 0.5);
    CASM::ConfigDoF mapped_configdof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(final_struc, mapped_configdof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    CASM::Configuration init_config(scel);
    init_config.init_occupation();
    init_config.init_deformation();
    init_config.init_displacement();
    CASM::Configuration final_config(scel, CASM::jsonParser(), mapped_configdof);
    // std::cout << "init_config lattice" << init_config.supercell().lattice().lat_column_mat() << std::endl;
    // std::cout << "final_config lattice" << final_config.supercell().lattice().lat_column_mat() << std::endl;
    // std::cout << "final_config deformation" << final_config.deformation() << std::endl;
    auto new_F = CASM::StrainConverter::right_stretch_tensor(final_config.deformation());
    final_config.set_deformation(new_F);
    // std::cout << "final_config deformation new" << new_F << std::endl;
    // std::cout << "i_occ,i_def,i_disp,f_occ,f_def,f_disp" << init_config.has_occupation() <<
    // init_config.has_deformation() << init_config.has_displacement() << final_config.has_occupation() <<
    // final_config.has_deformation() << final_config.has_displacement() << std::endl; std::cout << "interpolator
    // attempt to initialize" << std::endl;
    CASM::ConfigEnumInterpolation interpolator(init_config, final_config, n_images);
    // std::cout << "interpolator initialized" << std::endl;
    // int count = 0;
    for (auto& image : interpolator)
    {
        //    std::cout << "image no " << count << std::endl;
        images.push_back(make_deformed_struc(image));
        // count++;
    }
    return images;
}

std::vector<Rewrap::Structure> deformation_pathway(const Rewrap::Structure& init_struc,
                                                   const Rewrap::Structure& final_struc, int n_images, bool minimal)
{
    Rewrap::Structure mapped_struc = final_struc;
    Eigen::Matrix3i transfmat = Eigen::Matrix3i::Identity();
    if (minimal)
    {
        auto pair = _minimally_distorted_structure(init_struc, final_struc);
        mapped_struc = pair.first;
        transfmat = pair.second;
    }
    else
    {
        auto pair = _distorted_structure(init_struc, final_struc);
        mapped_struc = pair.first;
        transfmat = pair.second;
    }
    CASM::PrimClex* pclex = new CASM::PrimClex(init_struc, CASM::null_log());
    CASM::Supercell scel(pclex, transfmat);
    CASM::Configuration config(scel);
    config.init_occupation();
    config.init_deformation();
    config.init_displacement();
    Rewrap::Structure first_image = make_deformed_struc(config);
    return interpolate(first_image, mapped_struc, n_images);
}

std::pair<double, bool> gus_entry(const Rewrap::Structure& host_struc, const Rewrap::Structure& test_struc,
                                  bool sym_break_only)
{
    double score = 1e9;
    bool grp_sbgrp = false;
    CASM::PrimClex pclex(host_struc, CASM::null_log());
    CASM::ConfigMapper mapper(pclex, 0.5,0.2);
    mapper.restricted();
    //double divisor = 1.0 * test_struc.basis.size() / host_struc.basis.size();
    CASM::ConfigDoF mapped_dof;
    CASM::Lattice mapped_lat;
    if (/*floor(divisor) == divisor &&*/ mapper.struc_to_configdof(test_struc, mapped_dof, mapped_lat))
    {
        CASM::Supercell scel(&pclex, mapped_lat);
        if (sym_break_only)
        {
            mapped_dof = _sym_break_projection(scel, mapped_dof);
        }
        CASM::Configuration config(scel, CASM::jsonParser(), mapped_dof);
        double sc = CASM::ConfigMapping::strain_cost(test_struc.lattice(), mapped_dof, test_struc.basis.size());
        double bc = CASM::ConfigMapping::basis_cost(mapped_dof, test_struc.basis.size());
        score = (mapper.lattice_weight() * sc + (1.0 - mapper.lattice_weight()) * bc);
        double ratio =
            1.0 * host_struc.factor_group().size() / config.multiplicity() / test_struc.factor_group().size();
		std::vector<int> sbgrp_sizes;
		auto tree_pair=host_struc.factor_group().unique_subgroup_tree();
		std::string curr_name=test_struc.factor_group().get_name();
		auto subgroup_it=std::find_if(tree_pair.first.begin(),tree_pair.first.end(),[&curr_name](const CASM::SymGroup&g){return g.get_name()==curr_name;});
		auto subgroup_index=std::distance(tree_pair.first.begin(),subgroup_it);
		if (subgroup_index < tree_pair.first.size()){
			grp_sbgrp= tree_pair.second[subgroup_index].size()==2;
		}
		std::cout << "DEBUGGING: host_struc.title " << host_struc.title  << std::endl;
		std::cout << "DEBUGGING: test_struc.title " << test_struc.title  << std::endl;
		std::cout << "DEBUGGING: host_struc.factor_group().name " << host_struc.factor_group().get_name()  << std::endl;
		std::cout << "DEBUGGING: test_struc.factor_group().name " << test_struc.factor_group().get_name()  << std::endl;
		std::cout << "DEBUGGING: host_struc.factor_group().size() / config.multiplicity()" << host_struc.factor_group().size() / config.multiplicity() << std::endl;
		std::set<int> fl_sizes;
		for ( int i=0;i< tree_pair.second.size();++i){
			if (tree_pair.second[i].size()==2){
				fl_sizes.insert(tree_pair.first[i].size());
			}
		}
		for (auto size : fl_sizes){
			std::cout << "first level subgroup size: " << size << std::endl;
		}	
        //grp_sbgrp = (floor(ratio) == ratio && std::find(fl_sizes.begin(),fl_sizes.end(),host_struc.factor_group().size()/config.multiplicity())!=fl_sizes.end());
    }
    return std::make_pair(score, grp_sbgrp);
}

std::vector<Rewrap::Structure> read_and_rename_json(const Rewrap::fs::path& struc_folder)
{
    std::vector<Rewrap::Structure> prim_list;
    for (auto& struc_path : Rewrap::fs::directory_iterator(struc_folder))
    {
        if (Rewrap::fs::is_regular(struc_path))
        {
            Rewrap::Structure prim = CASM::Structure(CASM::read_prim(CASM::jsonParser(struc_path), 1e-4));
            prim.title = struc_path.path().filename().string();
            prim_list.push_back(prim);
        }
    }
    return prim_list;
}

std::vector<Rewrap::Structure> read_and_rename_poscar(const Rewrap::fs::path& struc_folder)
{
    std::vector<Rewrap::Structure> struc_list;
    for (const auto& struc_path : Rewrap::fs::directory_iterator(struc_folder))
    {
        if (Rewrap::fs::is_regular(struc_path))
        {
            Rewrap::Structure struc = Rewrap::Structure::from_poscar(struc_path);
            struc.title = struc_path.path().filename().string();
            struc_list.push_back(struc);
        }
    }
    return struc_list;
}

Rewrap::Structure reassign_all_occs(const Rewrap::Structure& original, const std::set<std::string>& occ_list)
{
    Rewrap::Structure copy = original;
    std::vector<CASM::Molecule> dof_list;
    for (auto it = occ_list.begin(); it != occ_list.end(); ++it)
    {
        dof_list.emplace_back(CASM::Molecule::make_atom(*it));
    }
    for (auto& site : copy.basis)
    {
        site.set_allowed_species(dof_list);
    }
    return copy;
}

Rewrap::Structure symmetrize(const Rewrap::Structure& struc, double tol)
{
	std::set<int> point_group_sizes {1,2,3,4,6,8,12,16,24,48};
	Rewrap::Structure sym_struc=struc;
    int biggest = struc.factor_group().size();
    for (double f = 0.1; f < 1.1; f += 0.1)
    {
		CASM::Structure tmp = struc;
        // a) symmetrize the lattice vectors
		CASM::Lattice lat = tmp.lattice();
        lat.symmetrize(tol * f);
        lat.set_tol(tol * f);
		CASM::SymGroup pg;
        lat.generate_point_group(pg);
        if (point_group_sizes.find(pg.size()) == point_group_sizes.end())
        {
            continue;
        }
        tmp.set_lattice(lat, CASM::FRAC);

        tmp.factor_group();
        // b) find factor group with same tolerance
        tmp.fg_converge(tol * f);
        if (point_group_sizes.find(tmp.point_group().size()) == point_group_sizes.end())
        {
            continue;
        }

        // c) symmetrize the basis sites
        tmp.symmetrize(tmp.factor_group());
        if (tmp.factor_group().is_group(f * tol) && (tmp.factor_group().size() > biggest))
        {
            sym_struc=tmp;
            biggest = tmp.factor_group().size();
        }
    }
	return sym_struc;
}
