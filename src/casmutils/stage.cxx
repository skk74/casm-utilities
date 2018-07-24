#include "casmutils/stage.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
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
std::pair<Rewrap::Structure, Eigen::Matrix3i> _minimally_distorted_structure(Rewrap::Structure& ref_struc,
                                                                             Rewrap::Structure& deformed_struc)
{
    CASM::PrimClex* pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, 0.5);
    CASM::ConfigDoF mapped_dof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_dof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_pair(_apply_configdof(scel, _sym_break_projection(scel, mapped_dof)), scel.transf_mat());
}

Rewrap::Structure minimally_distorted_structure(Rewrap::Structure& ref_struc, Rewrap::Structure& deformed_struc)
{
    return _minimally_distorted_structure(ref_struc, deformed_struc).first;
}

std::pair<Rewrap::Structure, Eigen::Matrix3i> _distorted_structure(Rewrap::Structure& ref_struc,
                                                                   Rewrap::Structure& deformed_struc)
{
    CASM::PrimClex* pclex = new CASM::PrimClex(ref_struc, CASM::null_log());
    CASM::ConfigMapper mapper(*pclex, 0.5);
    CASM::ConfigDoF mapped_configdof;
    CASM::Lattice mapped_lat;
    mapper.struc_to_configdof(deformed_struc, mapped_configdof, mapped_lat);
    CASM::Supercell scel(pclex, mapped_lat);
    return std::make_pair(_apply_configdof(scel, mapped_configdof), scel.transf_mat());
}
Rewrap::Structure distorted_structure(Rewrap::Structure& ref_struc, Rewrap::Structure& deformed_struc)
{
    return _distorted_structure(ref_struc, deformed_struc).first;
}
std::vector<Rewrap::Structure> interpolate(Rewrap::Structure& init_struc, Rewrap::Structure& final_struc, int n_images)
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
    std::cout << "init_config lattice" << init_config.supercell().lattice().lat_column_mat() << std::endl;
    std::cout << "final_config lattice" << final_config.supercell().lattice().lat_column_mat() << std::endl;
    std::cout << "final_config deformation" << final_config.deformation() << std::endl;
    auto new_F = CASM::StrainConverter::right_stretch_tensor(final_config.deformation());
    final_config.set_deformation(new_F);
    std::cout << "final_config deformation new" << new_F << std::endl;

    CASM::ConfigEnumInterpolation interpolator(init_config, final_config, n_images);
    for (auto& image : interpolator)
    {
        images.push_back(make_deformed_struc(image));
    }
    return images;
}

std::vector<Rewrap::Structure> deformation_pathway(Rewrap::Structure& init_struc, Rewrap::Structure& final_struc,
                                                   int n_images, bool minimal)
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
