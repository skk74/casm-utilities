#include "casmutils/stage.hpp"
#include <casm/crystallography/SymTools.hh>
#include <casm/crystallography/StrucMapping.hh>
#include <casm/crystallography/SimpleStrucMapCalculator.hh>
std::tuple<Eigen::Matrix3d, Eigen::Matrix3d, Eigen::MatrixXd> structure_map(
    CASM::xtal::SimpleStructure &_parent,
    std::unordered_set<std::string> &allowed_occupants, CASM::xtal::SimpleStructure &_test,
    double max_va_concentration, double lattice_weight) {
    std::vector<std::unordered_set<std::string>> occs(_parent.n_atom(),
						      allowed_occupants);
    CASM::xtal::SimpleStrucMapCalculator calculator(
	_parent,
	CASM::xtal::make_point_group(CASM::xtal::Lattice(_parent.lat_column_mat)),
	CASM::xtal::SimpleStructure::SpeciesMode::ATOM, occs);

	double max_vol_change=0.5;
	int opts=0;
	double tol=1e-5;
	double min_va_frac=0;
	CASM::xtal::StrucMapper mapper(calculator,lattice_weight,max_vol_change,opts,tol,min_va_frac,max_va_concentration);
	std::set<CASM::xtal::MappingNode> nodeset= mapper.map_deformed_struc(_test);
	auto bestnode=*nodeset.begin();
	for (const auto&node : nodeset){
		if (node < bestnode){
			bestnode=node;
		}
	}
	Eigen::Matrix3d superlattice, left_stretch;
	Eigen::MatrixXd disp;
	left_stretch=bestnode.stretch();
	disp=bestnode.displacement;
	superlattice=bestnode.lat_node.parent.superlattice().lat_column_mat();
	return std::make_tuple(superlattice,left_stretch, disp);
}
