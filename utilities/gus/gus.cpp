#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"
#include "sqlite3.h"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void gus_initializer(po::options_description& gus_desc)
{
    UtilityProgramOptions::add_help_suboption(gus_desc);
    UtilityProgramOptions::add_desc_suboption(gus_desc);
    UtilityProgramOptions::add_output_suboption(gus_desc);
    gus_desc.add_options()("library,l", po::value<fs::path>()->required(),
                           "prim.json like files that represent the library of structures for prefiltering maps");
    gus_desc.add_options()("structure-folder,s", po::value<fs::path>()->required(),
                           "POS.vasp like files that represent the structures that you wish to categorize");
    gus_desc.add_options()("sym-break-only",
                           "if this flag is specified the symmetry preserving part of mapping score is removed");
    gus_desc.add_options()("weight,w", po::value<double>()->required(),
                           "dictates the lattice-basis tradeoff that determines the best mapping. 0 is only basis score, 1 is only lattice score");
	
    return;
}
} // namespace Utilities
static int callback(void* NotUsed, int argc, char** argv, char** azColName)
{
    int i;
    for (i = 0; i < argc; i++)
    {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
}

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler gus_launch(argc, argv, gus_initializer);

    if (gus_launch.count("help"))
    {
        std::cout << gus_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        gus_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto library_path = gus_launch.fetch<fs::path>("library");
    auto struc_path = gus_launch.fetch<fs::path>("structure-folder");
    auto lattice_weight = gus_launch.fetch<double>("weight");
    std::cout << "paths loaded" << std::endl;
    auto lib_list = read_and_rename_json(library_path);
    std::cout << "lib list ready" << std::endl;
    auto struc_list = read_and_rename_poscar(struc_path);
    std::cout << "poscar list ready" << std::endl;
    sqlite3* db;
    sqlite3_stmt* stmt;
    char* zErrMsg = 0;
    std::cout << "am i thread safe? " << sqlite3_threadsafe() << std::endl;
    sqlite3_open(("prefilter_"+ std::to_string(lattice_weight)+".db").c_str(), &db);
    sqlite3_exec(
        db, "create table if not exists prefilter (name text, host text, struc text, score real, lattice_score real, basis_score real, grp_sbgrp integer); ",
        callback, 0, &zErrMsg);
#pragma omp parallel for
    for (int i = 0; i < lib_list.size(); i++)
    {
#pragma omp parallel for
        for (int j = 0; j < struc_list.size(); j++)
        {
            int already_exists = 0;
            sqlite3_stmt* stmt;
            std::string name = struc_list[j].title + "_onto_" + lib_list[i].title;
            std::cout << name << std::endl;
            sqlite3_prepare(db, ("SELECT EXISTS(SELECT 1 FROM prefilter WHERE name='" + name + "');").c_str(), -1,
                            &stmt, NULL);
            sqlite3_step(stmt);
            auto res = sqlite3_column_text(stmt, 0);
            already_exists = res[0] - '0';
            sqlite3_finalize(stmt);
            if (!already_exists)
            {
                const auto info_tuple = gus_entry(lib_list[i], struc_list[j], gus_launch.count("sym-break-only"), lattice_weight);
                int rc = sqlite3_exec(db,
                                      ("INSERT into prefilter (name, host, struc, score,lattice_score, basis_score, grp_sbgrp) VALUES ('" + name +
                                       "'," + "'" + lib_list[i].title + "'," + "'" + struc_list[j].title + "'," +
                                       std::to_string(std::get<0>(info_tuple)) + "," + std::to_string(std::get<1>(info_tuple)) + 
									   "," + std::to_string(std::get<2>(info_tuple)) + "," + std::to_string(std::get<3>(info_tuple)) + ");")
                                          .c_str(),
                                      callback, 0, &zErrMsg);
                if (rc != SQLITE_OK)
                {
                    std::cout << "SQL error: " << zErrMsg << std::endl;
                    sqlite3_free(zErrMsg);
                }
                else
                {
                    std::cout << "insert success" << std::endl;
                }
            }
            else
            {
                std::cout << "entry has already been written moving on" << std::endl;
            }
        }
    }

    sqlite3_prepare(db, "SELECT struc from prefilter WHERE score<=0.07 GROUP BY struc;", -1, &stmt,
                    NULL);
    sqlite3_step(stmt);
    auto res = sqlite3_column_text(stmt, 0);
    std::set<std::string> prefiltered;
    while (res)
    {
        prefiltered.insert(std::string(reinterpret_cast<const char*>(res)));
        sqlite3_step(stmt);
        res = sqlite3_column_text(stmt, 0);
    }
    sqlite3_finalize(stmt);
    std::vector<Rewrap::Structure> non_prefiltered;
    std::remove_copy_if(struc_list.begin(), struc_list.end(), std::back_inserter(non_prefiltered),
                        [&prefiltered](Rewrap::Structure& struc) {
                            return std::find(prefiltered.begin(), prefiltered.end(), struc.title) != prefiltered.end();
                        });
    std::vector<Rewrap::Structure> prim_list = non_prefiltered;
    for (auto& struc : prim_list)
    {
        std::string name = struc.title;
        struc = reassign_all_occs(struc, {"A", "B", "X", "Z", "L", "M"});
        struc = Simplicity::make_primitive(struc);
        struc.title = name;
    }
    sqlite3_exec(
        db, "create table if not exists first_pass (name text, host text, struc text, score real, lattice_score real, basis_score real, grp_sbgrp integer); ",
        callback, 0, &zErrMsg);
#pragma omp parallel for
    for (int i = 0; i < prim_list.size(); i++)
    {
#pragma omp parallel for
        for (int j = 0; j < non_prefiltered.size(); j++)
        {
            int already_exists = 0;
            sqlite3_stmt* stmt;
            std::string name = non_prefiltered[j].title + "_onto_" + prim_list[i].title;
            std::cout << name << std::endl;
            sqlite3_prepare(db, ("SELECT EXISTS(SELECT 1 FROM first_pass WHERE name='" + name + "');").c_str(), -1,
                            &stmt, NULL);
            sqlite3_step(stmt);
            auto res = sqlite3_column_text(stmt, 0);
            already_exists = res[0] - '0';
            sqlite3_finalize(stmt);
            if (!already_exists)
            {
                const auto info_tuple = gus_entry(prim_list[i], non_prefiltered[j], gus_launch.count("sym-break-only"),lattice_weight);
                int rc = sqlite3_exec(db,
                                      ("INSERT into first_pass (name, host, struc, score, lattice_score, basis_score, grp_sbgrp) VALUES ('" + name +
                                       "'," + "'" + prim_list[i].title + "'," + "'" + non_prefiltered[j].title + "'," +
                                       std::to_string(std::get<0>(info_tuple)) + "," + std::to_string(std::get<1>(info_tuple)) + 
									   "," + std::to_string(std::get<2>(info_tuple)) + "," + std::to_string(std::get<3>(info_tuple)) + ");")
                                          .c_str(),
                                      callback, 0, &zErrMsg);
                if (rc != SQLITE_OK)
                {
                    std::cout << "SQL error: " << zErrMsg << std::endl;
                    sqlite3_free(zErrMsg);
                }
                else
                {
                    std::cout << "insert success" << std::endl;
                }
            }
            else
            {
                std::cout << "entry has already been written moving on" << std::endl;
            }
        }
    }
    sqlite3_prepare(
        db, "SELECT struc,MIN(score) from first_pass WHERE (NOT first_pass.host=first_pass.struc) GROUP BY struc;", -1,
        &stmt, NULL);
    sqlite3_step(stmt);
    res = sqlite3_column_text(stmt, 0);
    std::map<std::string, double> upper_sym_limits;
    while (res)
    {
        upper_sym_limits.insert(
            std::make_pair(std::string(reinterpret_cast<const char*>(res)),
                           std::stod(std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1))))));
        sqlite3_step(stmt);
        res = sqlite3_column_text(stmt, 0);
    }
    sqlite3_finalize(stmt);
    for (auto& item : upper_sym_limits)
    {
        item.second = sqrt(item.second);
        item.second = std::min(item.second, 0.1);
    }
    std::vector<Rewrap::Structure> sym_strucs;
    for (const auto& struc : prim_list)
    {
        std::string name = struc.title + "symmed";
        Rewrap::Structure symmed = symmetrize(struc, upper_sym_limits[struc.title]);
        if (symmed.factor_group().size() > struc.factor_group().size())
        {
            auto primmed_symmed = symmetrize(Simplicity::make_primitive(symmed), upper_sym_limits[struc.title]);
            symmed.title = name;
            sym_strucs.push_back(primmed_symmed);
        }
    }
    sqlite3_exec(
        db,
        "create table if not exists second_pass (name text, host text, struc text, score real, lattice_score real, basis_score real, grp_sbgrp integer); ",
        callback, 0, &zErrMsg);
#pragma omp parallel for
    for (int i = 0; i < sym_strucs.size(); i++)
    {
#pragma omp parallel for
        for (int j = 0; j < non_prefiltered.size(); j++)
        {
            int already_exists = 0;
            sqlite3_stmt* stmt;
            std::string name = non_prefiltered[j].title + "_onto_" + sym_strucs[i].title;
            std::cout << name << std::endl;
            sqlite3_prepare(db, ("SELECT EXISTS(SELECT 1 FROM second_pass WHERE name='" + name + "');").c_str(), -1,
                            &stmt, NULL);
            sqlite3_step(stmt);
            auto res = sqlite3_column_text(stmt, 0);
            already_exists = res[0] - '0';
            sqlite3_finalize(stmt);
            if (!already_exists)
            {
                const auto info_tuple = gus_entry(sym_strucs[i], non_prefiltered[j], gus_launch.count("sym-break-only"),lattice_weight);
                int rc =
                    sqlite3_exec(db,
                                 ("INSERT into second_pass (name, host, struc, score, lattice_score, basis_score, grp_sbgrp) VALUES ('" + name +
                                  "'," + "'" + sym_strucs[i].title + "'," + "'" + non_prefiltered[j].title + "'," +
                                  std::to_string(std::get<0>(info_tuple)) + "," + std::to_string(std::get<1>(info_tuple)) + 
									   "," + std::to_string(std::get<2>(info_tuple)) + "," + std::to_string(std::get<3>(info_tuple)) + ");")
                                     .c_str(),
                                 callback, 0, &zErrMsg);
                if (rc != SQLITE_OK)
                {
                    std::cout << "SQL error: " << zErrMsg << std::endl;
                    sqlite3_free(zErrMsg);
                }
                else
                {
                    std::cout << "insert success" << std::endl;
                }
            }
            else
            {
                std::cout << "entry has already been written moving on" << std::endl;
            }
        }
    }
    sqlite3_close(db);
    // if (gus_launch.vm().count("output"))
    //{
    //    auto out_path = gus_launch.fetch<fs::path>("output");
    //    int count = 0;
    //    for (auto& item : strucs)
    //    {
    //        // TODO: what if directory doesn't exist?
    //        std::ostringstream ostr;
    //        ostr << std::setfill('0') << std::setw(2) << count;
    //        Simplicity::write_poscar(item, out_path / Rewrap::fs::path("image" + ostr.str() + "POSCAR"));
    //        count++;
    //    }
    //    return 0;
    //}
    // int count = 0;
    // for (auto& item : strucs)
    //{
    //    std::cout << "image " << count << std::endl;
    //    Simplicity::print_poscar(item, std::cout);
    //    count++;
    //}
    return 0;
}
