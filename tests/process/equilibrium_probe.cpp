// Native evaluator used by the literature data comparison script.
#include <fstream>
#include <iostream>
#include "cycle/input.hpp"
#include "cycle/equilibrium.hpp"
using namespace cycle;
int main(int argc,char** argv) try {
 if(argc!=3)throw std::invalid_argument("Usage: equilibrium_probe simulation.yaml points.json");
 auto cfg=load_config(argv[1]);Json points;std::ifstream(argv[2])>>points;
 const auto names=cfg.at("comp_names").get<std::vector<std::string>>();
 MaterialEquilibrium model(cfg,cfg.at("R"),names.size(),1);
 Json out=Json::array();
 for(const auto& point:points){
  const double p=point.at("P_Pa"),t=point.at("T_K");
  auto y=point.at("y").get<Vec>();
  if(y.size()!=names.size())throw std::invalid_argument("Composition size mismatch");
  for(auto& v:y)v*=p;
  out.push_back({{"q_absolute_mol_kg",model.loading(y,{t})}});
 }
 std::cout<<out.dump(2)<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
