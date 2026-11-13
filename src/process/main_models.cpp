#include "cycle/main_models.hpp"
#include "inputreader.h"
#include "mixture_prediction.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <optional>
#include <set>

namespace cycle {
namespace {
std::string low(std::string s){for(auto& c:s)c=char(std::tolower(static_cast<unsigned char>(c)));return s;}
Json get(const Json& j,const std::string& key,Json fallback=nullptr){for(const auto& [k,v]:j.items())if(low(k)==low(key))return v;return fallback;}
Json clean_site(const Json& site){Json r;for(const auto& [k,v]:site.items())if(std::set<std::string>{"type","parameters","rateequation","rateequationparameters","heatofadsorption"}.contains(low(k)))r[k]=v;return r;}
Json parser_input(const Json& raw,const std::string& material){
  // Parse material laws without imposing legacy breakthrough carrier/boundary requirements.
  Json d={{"SimulationType","MixturePrediction"},{"MixturePredictionMethod","SPI"},{"BoundaryCondition","FixedVelocity"},
    {"BreakthroughIntegrator","RungeKutta3"},{"ConvectionScheme","Upwind"},{"Temperature",get(raw,"Temperature",298.15)},
    {"ColumnLength",.1},{"ColumnEntranceVelocity",.01},{"InletPressure",1e5},{"OutletPressure",1e5},{"NumberOfGridPoints",4},
    {"ParticleDensity",get(raw,"ParticleDensity",1000.)},{"FluidPhase",get(raw,"FluidPhase","Gas")},
    {"Geometry",{{"Type","PackedBed"},{"InternalDiameter",.01},{"OuterDiameter",.012},{"ColumnVoidFraction",.4},{"ParticleDiameter",.001}}}};
  if (low(get(raw,"MixturePredictionMethod","IAST").get<std::string>()) == "mpd") {
    d["MixturePredictionMethod"] = "MPD";
    if (!get(raw,"MPDSettings").is_null()) d["MPDSettings"] = get(raw,"MPDSettings");
  }
  const auto g=get(raw,"Geometry",Json::object());
  const std::set<std::string> geometry_keys={"type","internaldiameter","outerdiameter","columnvoidfraction","particlediameter","particlesphericity","channelshape","internalchanneldimension","numberofchannels","washcoatthickness","washcoatvolumeperchannelvolume","forchheimercoefficient"};
  for(const auto& [k,v]:g.items())if(geometry_keys.contains(low(k)))d["Geometry"][k]=v;
  for(auto key:{"ChannelDiameter","ChannelSize"})if(!get(g,key).is_null())d["Geometry"]["InternalChannelDimension"]=get(g,key);
  for(auto key:{"pHMode","pHValue","pKw","pHComponent","LiquidDensity"})if(!get(raw,key).is_null())d[key]=get(raw,key);
  d["Components"]=Json::array();
  const std::set<std::string> component_keys={"name","carriergas","gasphasemolfraction","liquidphaseconcentration","initialliquidphaseconcentration","masstransfercoefficient","axialdispersioncoefficient","molecularweight","heatofadsorption","chemisorptionsites"};
  for(const auto& c:get(raw,"Components",Json::array())){
    Json x=Json::object();for(const auto& [k,v]:c.items())if(component_keys.contains(low(k)))x[k]=v;
    x["PhysisorptionSites"]=Json::array();
    for(const auto& site:get(c,"PhysisorptionSites",Json::array()))x["PhysisorptionSites"].push_back(clean_site(site));
    d["Components"].push_back(x);
  }
  d["Reactions"]=Json::array();
  for(auto r:get(raw,"Reactions",Json::array())){
    const auto allowed=get(r,"Adsorbents");
    if(allowed.is_array()&&!material.empty()&&std::find(allowed.begin(),allowed.end(),Json(material))==allowed.end())continue;
    for(auto it=r.begin();it!=r.end();)if(low(it.key())=="adsorbents")it=r.erase(it);else++it;
    d["Reactions"].push_back(r);
  }
  return d;
}
std::size_t method(const Json& raw){const std::map<std::string,std::size_t> m={{"iast",0},{"siast",1},{"ei",2},{"sei",3},{"sci",4},{"spi",5},{"mpd",6}};auto s=low(get(raw,"MixturePredictionMethod","IAST").get<std::string>());if(!m.contains(s))throw std::invalid_argument("Main material adapter supports SPI/SCI/EI/SEI/IAST/SIAST/MPD");return m.at(s);}
MixturePrediction mixture(const std::vector<Component>& c,std::size_t prediction,double T){
  std::size_t count=0,index=0;for(std::size_t j=0;j<c.size();++j)if(c[j].isCarrierGas){++count;index=j;}
  return MixturePrediction("Shared main models",c,count,index,T,1.,1e6,10,0,prediction,0);
}
}
struct MainMaterialModels::Impl {
  Json raw;
  bool liquid=false;
  std::size_t prediction=0;
  std::optional<MPDSettings> mpd_settings;
  mutable std::optional<MixturePrediction> cached;
  mutable double cached_temperature=-1.;
};
MainMaterialModels::MainMaterialModels(const Json& raw,const std::string& material):p(std::make_unique<Impl>()){
  p->raw=raw;p->liquid=low(get(raw,"FluidPhase","Gas").get<std::string>())=="liquid";p->prediction=method(raw);
  const InputReader reader("simulation.json",parser_input(raw,material).dump());
  if (p->prediction == 6) {
    // The InputReader constructor validates the public MPD schema and the
    // non-carrier dimension order. Keep one loaded kernel across all T values;
    // its original reweighting routine takes the local temperature on each call.
    p->mpd_settings = reader.mpdSettings;
    p->cached.emplace(reader);
  }
  components=reader.components;reactions=reader.reactions;geometry=reader.geometry;particle_density=reader.particleDensity;
  for(const auto& c:components){chemical_sites=std::max(chemical_sites,c.chemisorption.numberOfSites);for(const auto& s:c.chemisorption.sites)surface_pore|=s.usesSurfacePoreTransport();}
  for(const auto& r:reactions)surface_pore|=r.phase==Reaction::Phase::PoreConcentration;
}
MainMaterialModels::~MainMaterialModels()=default;
MainMaterialModels::MainMaterialModels(const MainMaterialModels& a):components(a.components),reactions(a.reactions),geometry(a.geometry),chemical_sites(a.chemical_sites),surface_pore(a.surface_pore),particle_density(a.particle_density),p(std::make_unique<Impl>(*a.p)){}
MainMaterialModels& MainMaterialModels::operator=(const MainMaterialModels& a){if(this!=&a){components=a.components;reactions=a.reactions;geometry=a.geometry;chemical_sites=a.chemical_sites;surface_pore=a.surface_pore;particle_density=a.particle_density;p=std::make_unique<Impl>(*a.p);}return *this;}
std::optional<Vec> MainMaterialModels::mpd_maximum_loadings() const {
  if (!p->mpd_settings) return std::nullopt;
  constexpr double avogadro = 6.02214076e23;
  const double conversion = 1. / (avogadro * p->mpd_settings->referenceFrameworkMass);
  Vec capacities(components.size(), 0.);
  std::size_t dimension = 0;
  for (std::size_t j = 0; j < components.size(); ++j)
    if (!components[j].isCarrierGas)
      capacities[j] = p->mpd_settings->componentBounds.at(dimension++).nMax * conversion;
  return capacities;
}
double MainMaterialModels::pH(const Vec& concentrations) const {
  const auto mode=low(get(p->raw,"pHMode","Fixed").get<std::string>());
  if(mode=="fixed")return get(p->raw,"pHValue",7.).get<double>();
  const auto index=get(p->raw,"pHComponent",0).get<std::size_t>();
  if(index>=concentrations.size())throw std::invalid_argument("pHComponent is outside Components");
  const double exponent=-std::log10(std::max(concentrations[index]*1e-3,1e-300));
  return mode=="hplus"?exponent:get(p->raw,"pKw",14.).get<double>()-exponent;
}
Vec MainMaterialModels::equilibrium(const Vec& driving,const Vec& temperature,std::size_t nz) const {
  const auto nc=components.size();if(driving.size()!=nc*nz||(temperature.size()!=1&&temperature.size()!=nz))throw std::invalid_argument("Main equilibrium array shape");
  Vec result(driving.size(),0.);
  const auto raw_components=get(p->raw,"Components");
  for(std::size_t i=0;i<nz;++i){
    const double T=temperature[temperature.size()==1?0:i];
    if(p->prediction!=6 && (!p->cached||p->cached_temperature!=T)){
      auto c=components;
      for(std::size_t j=0;j<nc;++j){
        const auto raw=raw_components[j];const auto sites=get(raw,"PhysisorptionSites",Json::array());
        for(std::size_t s=0;s<c[j].isotherm.sites.size();++s){
          // The original parser pads shorter site lists (including carriers)
          // with inactive sites. They have no corresponding raw JSON entry.
          if(s>=sites.size())continue;
          auto& iso=c[j].isotherm.sites[s];const auto r=sites[s];
          const double heat=get(r,"EquilibriumHeat",get(raw,"nonIsothermal",false).get<bool>()?get(raw,"HeatOfAdsorption",0.):Json(0.));
          const auto ref=get(r,"referenceTemperature",get(raw,"referenceTemperature"));
          const double factor=std::exp(std::clamp(heat/8.31446261815324*(1./T-(ref.is_number()?1./ref.get<double>():0.)),-700.,700.));
          if(iso.type==Isotherm::Type::Henry)iso.parameters[0]*=factor;
          else if(iso.parameters.size()>1&&iso.type!=Isotherm::Type::Freundlich&&iso.type!=Isotherm::Type::GAB)iso.parameters[1]*=factor;
          const auto law=get(r,"CapacityTemperature");
          if(law.is_object())iso.parameters[0]*=std::exp(std::clamp(get(law,"Coefficient",0.).get<double>()*(T/get(law,"referenceTemperature",298.15).get<double>()-1.),-700.,700.));
        }
        c[j].nonIsothermal=false;
      }
      p->cached=mixture(c,p->prediction,T);p->cached_temperature=T;
    }
    Vec v(nc),x(nc),q(nc),cp(nc*std::max<std::size_t>(1,p->cached->maxIsothermTerms)),cg(std::max<std::size_t>(1,p->cached->maxIsothermTerms));
    double total=0.;for(std::size_t j=0;j<nc;++j){v[j]=std::max(driving[j*nz+i],0.);total+=v[j];}
    if(total==0.)continue;
    double Ti=T;
    if(p->liquid)p->cached->predictMixture(v,1.,x,q,cp,cg,Ti,pH(v),MixturePrediction::DrivingForceInput::Concentration);
    else {for(double& a:v)a/=total;p->cached->predictMixture(v,total,x,q,cp,cg,Ti,pH(v));}
    for(std::size_t j=0;j<nc;++j){if(!std::isfinite(q[j])||q[j]<-1e-12)throw std::runtime_error("Main equilibrium returned invalid loading");result[j*nz+i]=std::max(q[j],0.);}
  }
  return result;
}
Vec MainMaterialModels::chemical_equilibrium(const Vec& c,double Tg,double Ts) const {
  const auto nc=components.size();Vec result(chemical_sites*nc,0.);
  if(!chemical_sites)return result;
  // As in the original main path, MPD describes physical adsorption; the
  // independent chemical-site prediction uses IAST.
  auto chemical_components=components;
  std::string selected;
  bool physical_source=false;
  for(auto& component:chemical_components)for(auto& site:component.chemisorption.sites) {
    physical_source|=site.physicalEquilibrium;
    if(!site.mixtureMethod.empty()) {
      if(!selected.empty()&&low(selected)!=low(site.mixtureMethod))throw std::invalid_argument("Chemical sites require a consistent MixturePredictionMethod");
      selected=site.mixtureMethod;
    }
    if(site.isotherm && (site.equilibriumReferenceTemperature || site.equilibriumHeat!=0.0)) {
      auto& iso=*site.isotherm;
      const double ref=site.equilibriumReferenceTemperature.value_or(298.15);
      const double factor=std::exp(std::clamp(site.equilibriumHeat/8.31446261815324*(1.0/Ts-1.0/ref),-700.0,700.0));
      if(iso.type==Isotherm::Type::Henry)iso.parameters[0]*=factor;
      else if(iso.parameters.size()>1)iso.parameters[1]*=factor;
      iso.nonIsothermal=false;
    }
  }
  const auto prediction=selected.empty()?(p->prediction==6?0:p->prediction):method(Json{{"MixturePredictionMethod",selected}});
  auto physical=mixture(chemical_components,prediction,Ts);
  auto chem=MixturePrediction::makeChemisorptionPrediction(physical,chemical_components);
  Vec driving=c,x(nc),q(nc),cp(nc*std::max<std::size_t>(1,chem.maxIsothermTerms)),cg(std::max<std::size_t>(1,chem.maxIsothermTerms));
  double total=std::accumulate(c.begin(),c.end(),0.),T=Ts;
  if(p->liquid)chem.predictMixture(driving,1.,x,q,cp,cg,T,pH(c),MixturePrediction::DrivingForceInput::Concentration);
  else if(total>0.){for(auto& v:driving)v=std::max(v,0.)/total;chem.predictMixture(driving,total*8.31446261815324*Tg,x,q,cp,cg,T,pH(c));}
  Vec physical_target;
  if(physical_source){Vec force=c;if(!p->liquid)for(auto& value:force)value*=8.31446261815324*Tg;physical_target=equilibrium(force,Vec{Ts},1);}
  for(std::size_t s=0;s<chemical_sites;++s)for(std::size_t j=0;j<nc;++j){
    if((s*nc+j)<chem.equilibriumSiteLoadings.size())result[s*nc+j]=std::max(chem.equilibriumSiteLoadings[s*nc+j],0.);
    if(s<components[j].chemisorption.sites.size()){
      const auto& k=components[j].chemisorption.sites[s];
      if(k.physicalEquilibrium)result[s*nc+j]=physical_target[j]*k.equilibriumFraction;
      if(!k.usesEquilibriumLoading())result[s*nc+j]=k.maximumLoading;
      if(k.maximumLoading>0.)result[s*nc+j]=std::min(result[s*nc+j],k.maximumLoading);
    }
  }
  return result;
}
}
