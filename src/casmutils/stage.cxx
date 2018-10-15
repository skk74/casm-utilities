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

Rewrap::Structure add_va(const Rewrap::Structure &struc){
	auto cpy = struc;
    for (auto& site : cpy.basis)
    {
		std::vector<CASM::Molecule> dof_list;
		for (auto occ : site.allowed_occupants())
		{
			dof_list.emplace_back(CASM::Molecule::make_atom(occ));
		}
		dof_list.emplace_back(CASM::Molecule::make_vacancy());
        site.set_allowed_species(dof_list);
    }
	return cpy;

}

double va_concentration(const Rewrap::Structure &struc){
	double conc=0.0;
	for (auto &site : struc.basis){
		if (site.is_vacant()){
			conc+=1.0/struc.basis.size();
		}
	
	}
	return conc;
}
/// This function alters the coordinates of the given struc to have fractional
/// coordinates between 0 and 1 only
void bring_coords_within(Rewrap::Structure* struc)
{
    for (auto& site : struc->basis)
    {
        site.within();
    }
    return;
}

bool is_planar(const Eigen::Matrix3d &cart_op){
	//std::cout << cart_op << std::endl;
	//std::cout << cart_op(2,0)==0 << cart_op(2,1)==0 << cart_op(2,2)==1 << cart_op(1,2) << cart_op(0,2))==0 << std::endl; 
	return ( cart_op(2,0)==0 && cart_op(2,1)==0 && cart_op(2,2)==1 && cart_op(1,2) ==0 && cart_op(0,2)==0);
	
}

bool struc_equal(const Rewrap::Structure &lhs, const Rewrap::Structure &rhs){

}

} // namespace


std::tuple<Rewrap::Structure, Eigen::Matrix3i, std::vector<int>> _minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                                             const Rewrap::Structure& deformed_struc, bool allow_va,double lattice_weight)
{
	CASM::PrimClex* pclex;
	if (allow_va){
     pclex = new CASM::PrimClex(add_va(ref_struc), CASM::null_log());
	}
	else {
     pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
	}
	CASM::ConfigMapper mapper(*pclex, lattice_weight);
	mapper.restricted();
	mapper.set_max_va_frac(0.5);
    CASM::ConfigDoF mapped_dof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_dof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_tuple(_apply_configdof(scel, _sym_break_projection(scel, mapped_dof)), scel.transf_mat(),mapped_dof.occupation());
}

Rewrap::Structure minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                const Rewrap::Structure& deformed_struc, bool allow_va,double lattice_weight)
{
    return std::get<0>(_minimally_distorted_structure(ref_struc, deformed_struc,allow_va,lattice_weight));
}

std::tuple<Rewrap::Structure, Eigen::Matrix3i,std::vector<int>> _distorted_structure(const Rewrap::Structure& ref_struc,
                                                                   const Rewrap::Structure& deformed_struc, bool allow_va,double lattice_weight)
{
    CASM::PrimClex* pclex;
	if (allow_va){
     pclex = new CASM::PrimClex(add_va(ref_struc), CASM::null_log());
	}
	else {
     pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
	}
    CASM::ConfigMapper mapper(*pclex, lattice_weight);
	mapper.restricted();
	mapper.set_max_va_frac(0.5);
    CASM::ConfigDoF mapped_configdof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_configdof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_tuple(_apply_configdof(scel, mapped_configdof), scel.transf_mat(),mapped_configdof.occupation());
}
Rewrap::Structure distorted_structure(const Rewrap::Structure& ref_struc, const Rewrap::Structure& deformed_struc, bool allow_va, double lattice_weight)
{
    return std::get<0>(_distorted_structure(ref_struc, deformed_struc,allow_va, lattice_weight));
}
std::vector<Rewrap::Structure> interpolate(const Rewrap::Structure& init_struc, const Rewrap::Structure& final_struc,
                                           int n_images,double lattice_weight)
{
    std::vector<Rewrap::Structure> images;
    CASM::PrimClex* pclex = new CASM::PrimClex(init_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, lattice_weight);
    CASM::ConfigDoF mapped_configdof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(final_struc, mapped_configdof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    CASM::Configuration init_config(scel);
    init_config.init_occupation();
	init_config.set_occupation(mapped_configdof.occupation());
    init_config.init_deformation();
    init_config.init_displacement();
    CASM::Configuration final_config(scel, CASM::jsonParser(), mapped_configdof);
	if (!final_config.has_deformation()){
		final_config.init_deformation();
	}
	if (!final_config.has_displacement()){
		final_config.init_displacement();
	}
    auto new_F = CASM::StrainConverter::right_stretch_tensor(final_config.deformation());
    final_config.set_deformation(new_F);
    CASM::ConfigEnumInterpolation interpolator(init_config, final_config, n_images);
    for (auto& image : interpolator)
    {
        images.push_back(make_deformed_struc(image));
    }
    return images;
}

std::vector<Rewrap::Structure> deformation_pathway(const Rewrap::Structure& init_struc,
                                                   const Rewrap::Structure& final_struc, int n_images, bool minimal, bool allow_va, double lattice_weight)
{
    Rewrap::Structure mapped_struc = final_struc;

    Eigen::Matrix3i transfmat = Eigen::Matrix3i::Identity();
	std::vector<int> occ_vec;
    if (minimal)
    {
        auto tuple = _minimally_distorted_structure(init_struc, final_struc,allow_va,lattice_weight);
        mapped_struc = std::get<0>(tuple);
        transfmat = std::get<1>(tuple);
		occ_vec= std::get<2>(tuple);
    }
    else
    {
        auto tuple = _distorted_structure(init_struc, final_struc,allow_va,lattice_weight);
        mapped_struc = std::get<0>(tuple);
        transfmat = std::get<1>(tuple);
		occ_vec= std::get<2>(tuple);
    }
	CASM::PrimClex* pclex; 
	if (allow_va){
    pclex =new CASM::PrimClex(add_va(init_struc), CASM::null_log());
	}
	else{ 
    pclex = new CASM::PrimClex(init_struc, CASM::null_log());
	}
	CASM::Supercell scel(pclex, transfmat);
    CASM::Configuration config(scel);
    config.init_occupation();
	config.set_occupation(occ_vec);
    config.init_deformation();
    config.init_displacement();
    Rewrap::Structure first_image = make_deformed_struc(config);
	std::cout << "Va concentration is " << va_concentration(first_image) << std::endl;
    return interpolate(first_image, mapped_struc, n_images,lattice_weight);
}

std::tuple<double,double,double,bool> gus_entry(const Rewrap::Structure& host_struc, const Rewrap::Structure& test_struc,
                                  bool sym_break_only,double lattice_weight)
{
    double score = 1e9;
    bool grp_sbgrp = false;
    CASM::PrimClex pclex(host_struc, CASM::null_log());
    CASM::ConfigMapper mapper(pclex, lattice_weight);
	mapper.set_max_va_frac(0.5);
    mapper.restricted();
    //double divisor = 1.0 * test_struc.basis.size() / host_struc.basis.size();
    CASM::ConfigDoF mapped_dof;
    CASM::Lattice mapped_lat;
	double sc = 1e9;
	double bc = 1e9;
    if (/*floor(divisor) == divisor &&*/ mapper.struc_to_configdof(test_struc, mapped_dof, mapped_lat))
    {
        CASM::Supercell scel(&pclex, mapped_lat);
        if (sym_break_only)
        {
            mapped_dof = _sym_break_projection(scel, mapped_dof);
        }
        CASM::Configuration config(scel, CASM::jsonParser(), mapped_dof);
        sc = CASM::ConfigMapping::strain_cost(test_struc.lattice(), mapped_dof, test_struc.basis.size());
        bc = CASM::ConfigMapping::basis_cost(mapped_dof, test_struc.basis.size());
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
		//std::cout << "DEBUGGING: host_struc.title " << host_struc.title  << std::endl;
		//std::cout << "DEBUGGING: test_struc.title " << test_struc.title  << std::endl;
		//std::cout << "DEBUGGING: host_struc.factor_group().name " << host_struc.factor_group().get_name()  << std::endl;
		//std::cout << "DEBUGGING: test_struc.factor_group().name " << test_struc.factor_group().get_name()  << std::endl;
		//std::cout << "DEBUGGING: host_struc.factor_group().size() / config.multiplicity()" << host_struc.factor_group().size() / config.multiplicity() << std::endl;
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
    return std::make_tuple(score, sc, bc, grp_sbgrp);
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

///This function takes a 2D hexagonal layer and returns all the equivalents (may have some repeats)
///due to operation of the lattice point group + basis shifts [1/3, 2/3, 0] and [2/3, 1/3, 0]
std::vector<Rewrap::Structure> enumerate_layer_equivalents(const Rewrap::Structure &original){
	std::vector<Rewrap::Structure> equivs;
	CASM::SymGroup g;
	original.lattice().generate_point_group(g);
	for( auto &op : g){
		if (!is_planar(op.matrix())){
			continue;
		}
		auto cpy = original;
		for (auto &site : cpy.basis){
			site = CASM::copy_apply(op,site);
		}
		bring_coords_within(&cpy);	
		equivs.push_back(cpy);
		Frankenstein::shift_coords_by(&cpy,Eigen::Vector3d(1.0/3.0,2.0/3.0,0));
		bring_coords_within(&cpy);	
		equivs.push_back(cpy);
		Frankenstein::shift_coords_by(&cpy,Eigen::Vector3d(1.0/3.0,2.0/3.0,0));
		bring_coords_within(&cpy);	
		equivs.push_back(cpy);
	}
	return equivs;
}

#include <fstream>
#include <set>
namespace Viewpoint
{
struct NeighborInfo
{
    int index;
    double distance;
};

Eigen::Vector3d find_layer_direction(const Rewrap::Structure& struc, int num_samples, double nn_bubble_radius,
                                     std::set<std::string>& layering_species)
{
    Eigen::Vector3d total_hkl(0.0, 0.0, 0.0);
    // repeat finding normal direction from random sample points
    for (size_t i = 0; i < num_samples; i += 1)
    {
        int rand_idx = rand() % (struc.basis.size() - 1);
        auto& pivot = struc.basis[rand_idx];
        // first find all neighbors in range
        std::vector<NeighborInfo> neighbors;
        for (size_t j = 0; j < struc.basis.size(); j++)
        {
            const auto& neighbor = struc.basis[j];
            double dist = (neighbor.const_cart() - pivot.const_cart()).norm();
            auto end = layering_species.end();
            if (0 < dist && dist < nn_bubble_radius &&
                (!layering_species.size() ||
                 (layering_species.find(pivot.occ_name()) != end && layering_species.find(neighbor.occ_name()) != end)))
            {
                neighbors.push_back(NeighborInfo{j, dist});
            }
        }
        // sort the neighbors
        std::sort(neighbors.begin(), neighbors.end(),
                  [](const NeighborInfo& a, const NeighborInfo& b) { return a.distance < b.distance; });
        // cut the number of neighbors
        if (neighbors.size() < 2)
        {
            continue;
        }
        // get normals from the sample point and two neighbors
        std::vector<Eigen::Vector3d> normals;
        for (size_t j = 0; j < neighbors.size() - 1; j++)
        {
            const auto& neighbor0 = struc.basis[neighbors[j].index];
            const auto& neighbor1 = struc.basis[neighbors[j + 1].index];
            Eigen::Vector3d vec0 = neighbor0.const_frac() - pivot.const_frac();
            Eigen::Vector3d vec1 = neighbor1.const_frac() - pivot.const_frac();
            auto n = (vec0.cross(vec1)).normalized();
            if (n(2) < 0.0)
            {
                n = -1.0 * n;
            }
            else if (n(2) == 0.0)
            {
                if (n(1) < 0.0)
                {
                    n = -1.0 * n;
                }
            }
            Eigen::Vector3d vec2 = neighbor0.const_cart() - pivot.const_cart();
            Eigen::Vector3d vec3 = neighbor1.const_cart() - pivot.const_cart();
            if (std::abs(vec2.dot(vec3) / vec2.norm() / vec3.norm()) > 0.99)
            {
                // This sample's test vectors are essentially parallel. Not using because it will ruin average.
                continue;
            }
            normals.push_back(n);
        }
        // average the normals
        Eigen::Vector3d averaged(0.0, 0.0, 0.0);
        for (const auto& n : normals)
        {
            averaged += n;
        }
        averaged /= double(normals.size());
        averaged.normalize();

        total_hkl += averaged;
    }

    // average the normals from each sample points
    total_hkl /= double(num_samples);
    total_hkl.normalize();
    // This is the layering plane in terms of h k l
    Eigen::Vector3d ret_val = total_hkl;
    // This is the layering plane in terms of u v w
    ret_val = struc.lattice().inv_lat_column_mat() * struc.lattice().reciprocal().lat_column_mat() * ret_val;
    return ret_val;
}

Rewrap::Structure reoriented_struc(const Rewrap::Structure &struc,Eigen::Vector3i millers){
	CASM::Lattice transf_lat = struc.lattice().lattice_in_plane(millers);
	return struc.create_superstruc(transf_lat);
}

} // namespace Viewpoint
