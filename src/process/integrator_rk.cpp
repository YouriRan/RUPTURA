#include "cycle/integrator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cycle {
namespace {
using V=OdeState;
constexpr double eps=2.2204460492503131e-16;
bool finite(const V& y) { for(double v:y)if(!std::isfinite(v))return false;return true; }
double atol(const BdfOptions& o,std::size_t i) {return o.absolute_tolerance[o.absolute_tolerance.size()==1?0:i];}
double error_norm(const V& e,const V& a,const V& b,const BdfOptions& o) {
  const auto n=o.physical_dimension?o.physical_dimension:e.size();
  double sum=0.;
  for(std::size_t i=0;i<n;++i) {
    const double q=e[i]/(atol(o,i)+o.relative_tolerance*std::max(std::abs(a[i]),std::abs(b[i])));
    sum+=q*q;
  }
  return std::sqrt(sum/static_cast<double>(n));
}
bool admissible(const V& y,const BdfOptions& o) {
  if(!finite(y))return false;
  if(o.apply_nonnegative_constraints)
    for(auto i:o.nonnegative_indices)if(y.at(i)<-10.*atol(o,i))return false;
  return true;
}
V solve(V a,V b,std::size_t n) {
  for(std::size_t k=0;k<n;++k) {
    auto p=k;
    for(std::size_t i=k+1;i<n;++i)if(std::abs(a[i*n+k])>std::abs(a[p*n+k]))p=i;
    if(!std::isfinite(a[p*n+k])||std::abs(a[p*n+k])<1e-290)throw std::runtime_error("Singular implicit RK stage");
    if(p!=k){for(std::size_t j=k;j<n;++j)std::swap(a[k*n+j],a[p*n+j]);std::swap(b[k],b[p]);}
    for(std::size_t i=k+1;i<n;++i){const double f=a[i*n+k]/a[k*n+k];for(std::size_t j=k+1;j<n;++j)a[i*n+j]-=f*a[k*n+j];b[i]-=f*b[k];}
  }
  for(std::size_t ii=n;ii-->0;){for(std::size_t j=ii+1;j<n;++j)b[ii]-=a[ii*n+j]*b[j];b[ii]/=a[ii*n+ii];}
  return b;
}
struct RK {
  const OdeRhs& full;
  const OdeRhs& explicit_rhs;
  const BdfOptions& o;
  BdfStatistics stats;
  bool semi;
  std::string last_error;
  V cached_matrix;
  double cached_h=-1.;
  int matrix_uses=0;
  V evaluate(const OdeRhs& f,double t,const V& y) {
    if(o.maximum_rhs_evaluations && stats.rhs_evaluations>=o.maximum_rhs_evaluations)
      throw std::runtime_error("Maximum RHS evaluation count exceeded");
    V out(y.size());f(t,y,out);++stats.rhs_evaluations;
    if(out.size()!=y.size()||!finite(out))throw std::runtime_error("Nonfinite or incorrectly sized RK RHS");
    return out;
  }
  V implicit(double t,const V& y) {
    V a=evaluate(full,t,y),b=evaluate(explicit_rhs,t,y);
    for(std::size_t i=0;i<a.size();++i)a[i]-=b[i];return a;
  }
  V stage(const V& base,double t,double h) {
    V y=base;
    const std::size_t n=o.physical_dimension?o.physical_dimension:y.size();
    for(int iteration=0;iteration<12;++iteration) {
      const V g=implicit(t,y);V residual(y.size());
      for(std::size_t i=0;i<y.size();++i)residual[i]=y[i]-base[i]-h*g[i];
      const double norm=error_norm(residual,base,y,o);
      if(norm<.02) {
        // Passive ledger rows do not influence the physical system. Complete
        // their implicit equation after the physical Newton solve.
        for(std::size_t i=n;i<y.size();++i)y[i]=base[i]+h*g[i];
        if(!admissible(y,o))throw std::runtime_error("Implicit RK stage violates nonnegative constraints");
        return y;
      }
      V jac,right(n);
      for(std::size_t i=0;i<n;++i)right[i]=-residual[i];
      const bool rebuild=cached_matrix.size()!=n*n||cached_h!=h||matrix_uses>=12||iteration==6;
      if(rebuild) {
      jac.resize(n*n);
      for(std::size_t j=0;j<n;++j) {
        V trial=y;
        double delta=std::sqrt(eps)*std::max({std::abs(y[j]),atol(o,j)/o.relative_tolerance,1e-12});
        trial[j]+=delta;delta=trial[j]-y[j];const V next=implicit(t,trial);
        for(std::size_t i=0;i<n;++i)jac[i*n+j]=(i==j?1.:0.)-h*(next[i]-g[i])/delta;
      }
      ++stats.jacobian_evaluations;cached_matrix=jac;cached_h=h;matrix_uses=0;
      } else jac=cached_matrix;
      ++matrix_uses;++stats.lu_factorizations;
      const V correction=solve(std::move(jac),std::move(right),n);
      bool improved=false;
      for(double damping=1.;damping>=1./128.;damping*=.5) {
        V trial=y;
        for(std::size_t i=0;i<n;++i)trial[i]+=damping*correction[i];
        try {
          const V next=implicit(t,trial);V r(trial.size());
          for(std::size_t i=0;i<r.size();++i)r[i]=trial[i]-base[i]-h*next[i];
          if(error_norm(r,base,trial,o)<norm){y=std::move(trial);improved=true;break;}
        } catch(const std::exception&) {}
      }
      if(!improved){if(!rebuild){cached_matrix.clear();continue;}throw std::runtime_error("Implicit RK Newton line search failed");}
    }
    throw std::runtime_error("Implicit RK stage Newton iteration failed");
  }
  V step(double t,const V& y,double h) {
    if(!semi) {
      V k=evaluate(full,t,y),a=y;
      for(std::size_t i=0;i<y.size();++i)a[i]+=h*k[i];
      if(!admissible(a,o))throw std::runtime_error("SSPRK3 first stage violates constraints; reduce TimeStep");
      k=evaluate(full,t+h,a);V b=y;
      for(std::size_t i=0;i<y.size();++i)b[i]=.75*y[i]+.25*(a[i]+h*k[i]);
      if(!admissible(b,o))throw std::runtime_error("SSPRK3 second stage violates constraints; reduce TimeStep");
      k=evaluate(full,t+.5*h,b);V out=y;
      for(std::size_t i=0;i<y.size();++i)out[i]=y[i]/3.+(2./3.)*(b[i]+h*k[i]);
      return out;
    }
    // Pareschi & Russo, J. Sci. Comput. 25 (2005), Table 6:
    // SSP3(4,3,3). SSP refers to its explicit/limiting scheme, not an
    // unconditional positivity guarantee for arbitrary nonlinear splitting.
    constexpr double a=.24169426078821,b=.06042356519705,e=.12915286960590;
    const double ai[4][4]={{a,0,0,0},{-a,a,0,0},{0,1-a,a,0},{b,e,.5-b-e-a,a}};
    const double ae[4][4]={{0,0,0,0},{0,0,0,0},{0,1,0,0},{0,.25,.25,0}};
    constexpr double ci[4]={a,0,1,.5},ce[4]={0,0,1,.5},w[4]={0,1./6.,1./6.,2./3.};
    std::array<V,4> f,g;
    for(int s=0;s<4;++s) {
      V base=y;
      for(int j=0;j<s;++j)for(std::size_t k=0;k<y.size();++k)base[k]+=h*(ae[s][j]*f[j][k]+ai[s][j]*g[j][k]);
      const V state=stage(base,t+ci[s]*h,h*ai[s][s]);
      f[s]=evaluate(explicit_rhs,t+ce[s]*h,state);g[s]=implicit(t+ci[s]*h,state);
    }
    V out=y;
    for(int s=0;s<4;++s)for(std::size_t k=0;k<y.size();++k)out[k]+=h*w[s]*(f[s][k]+g[s][k]);
    return out;
  }
  V advance(double t,const V& y,double h) {
    if(!o.adaptive_rk)return step(t,y,h);
    V a=step(t,y,h*.5);return step(t+h*.5,a,h*.5);
  }
};
BdfResult integrate_rk(const OdeRhs& full,const OdeRhs& explicit_rhs,bool semi,double t,const V& initial,
                       double bound,const BdfOptions& o,const OdeAcceptedCallback& accepted,
                       const std::vector<OdeEvent>& events,const OdeSampleCallback& sample) {
  if(!full||(semi&&!explicit_rhs)||!std::isfinite(t)||!std::isfinite(bound)||bound<t||initial.empty()||!finite(initial))
    throw std::invalid_argument("Invalid RK initial-value problem");
  if(!(o.relative_tolerance>0.)||!std::isfinite(o.relative_tolerance)||!(o.maximum_step>0.)||std::isnan(o.maximum_step)||!o.maximum_steps||o.physical_dimension>initial.size())
    throw std::invalid_argument("Invalid RK step/tolerance options");
  if(o.absolute_tolerance.size()!=1&&o.absolute_tolerance.size()!=initial.size())throw std::invalid_argument("RK absolute tolerance shape");
  for(double v:o.absolute_tolerance)if(!(v>0.)||!std::isfinite(v))throw std::invalid_argument("Invalid RK absolute tolerance");
  double previous=t;
  for(double v:o.sample_times){if(v<previous||v>bound||!std::isfinite(v))throw std::invalid_argument("Invalid RK sample time");previous=v;}
  for(const auto& ev:events)if(!ev.value)throw std::invalid_argument("Empty RK event");
  RK rk{full,explicit_rhs,o,{},semi,{}};
  rk.stats.largest_order=3;rk.stats.linear_solver=semi?"native dense Newton SSP3(4,3,3)":"explicit SSPRK3";
  BdfResult result;result.time=t;result.state=initial;
  auto store=[&](double time,const V& value){if(result.times.empty()||time!=result.times.back()){result.times.push_back(time);result.states.push_back(value);}};
  store(t,initial);V y=initial;
  std::size_t si=0;
  while(si<o.sample_times.size()&&o.sample_times[si]==t){if(sample)sample(t,y);++si;}
  double h=o.first_step>0.?o.first_step:std::min(o.maximum_step,std::max((bound-t)*.01,1e-8));
  if(!std::isfinite(h))h=std::max((bound-t)*.01,1e-8);
  V before(events.size());for(std::size_t i=0;i<events.size();++i)before[i]=events[i].value(t,y);
  std::size_t attempts=0;
  try {
    while(t<bound) {
      if(rk.stats.accepted_steps>=o.maximum_steps||++attempts>10*o.maximum_steps+100)throw std::runtime_error("RK maximum step count exceeded");
      double step=std::min({h,o.maximum_step,bound-t});
      if(si<o.sample_times.size())step=std::min(step,o.sample_times[si]-t);
      if(o.initialization_steps&&rk.stats.accepted_steps<o.initialization_steps) {
        const double f=std::pow(1e-4,1.-double(rk.stats.accepted_steps+1)/double(o.initialization_steps));
        step=std::min(step,o.maximum_step*f);
      }
      if(!(step>0.)||t+step==t)throw std::runtime_error("RK step underflow: "+rk.last_error);
      V next;double error=0.;
      try {
        next=rk.advance(t,y,step);
        if(!admissible(next,o))throw std::runtime_error("RK accepted candidate violates constraints");
        if(o.adaptive_rk){const V coarse=rk.step(t,y,step);V d(y.size());for(std::size_t i=0;i<y.size();++i)d[i]=(next[i]-coarse[i])/7.;error=error_norm(d,y,next,o);}
      } catch(const std::exception& ex){rk.last_error=ex.what();error=INFINITY;}
      if(error>1.||!std::isfinite(error)){++rk.stats.rejected_steps;h=step*(std::isfinite(error)?std::clamp(.8*std::pow(error,-.25),.1,.5):.5);continue;}
      double next_t=t+step;
      std::vector<OdeEventHit> hits;
      for(std::size_t i=0;i<events.size();++i) {
        const double after=events[i].value(next_t,next),a=before[i];
        const bool rising=a<0.&&after>=0.,falling=a>0.&&after<=0.;
        if(!((rising&&events[i].direction>=0)||(falling&&events[i].direction<=0)))continue;
        double lo=0.,hi=step;V edge=next;
        for(int it=0;it<80&&hi-lo>std::max(o.event_time_tolerance,32*eps*std::max(1.,std::abs(t)));++it){double mid=.5*(lo+hi);V trial=rk.advance(t,y,mid);double value=events[i].value(t+mid,trial);if((a<0.&&value>=0.)||(a>0.&&value<=0.)){hi=mid;edge=std::move(trial);}else lo=mid;}
        edge=rk.advance(t,y,hi);hits.push_back({i,events[i].name,t+hi,std::move(edge)});
      }
      std::sort(hits.begin(),hits.end(),[](const auto& a,const auto& b){return a.time<b.time;});
      bool terminal=false;
      for(const auto& hit:hits){result.events.push_back(hit);if(events[hit.index].terminal){next_t=hit.time;next=hit.state;terminal=true;break;}}
      t=next_t;y=std::move(next);++rk.stats.accepted_steps;
      result.time=t;result.state=y;
      while(si<o.sample_times.size()&&o.sample_times[si]<=t+16*eps*std::max(1.,std::abs(t))){if(sample)sample(t,y);store(t,y);++si;}
      if(o.store_accepted_steps||terminal)store(t,y);
      if(accepted&&!accepted(t,y,rk.stats)){result.stopped_by_callback=true;break;}
      if(terminal){result.stopped_by_event=true;break;}
      for(std::size_t i=0;i<events.size();++i)before[i]=events[i].value(t,y);
      h=o.adaptive_rk?step*std::clamp(.9*std::pow(std::max(error,1e-12),-.25),.2,3.):o.maximum_step;
    }
    result.success=true;result.message="Integration completed";
  } catch(const std::exception& ex){result.message=ex.what();}
  store(result.time,result.state);result.statistics=rk.stats;return result;
}
} // namespace
BdfResult integrate_bdf1(const OdeRhs& f,double t,const V& y,double end,const BdfOptions& options,
    const OdeAcceptedCallback& a,const std::vector<OdeEvent>& e,const OdeSampleCallback& s) {
  BdfOptions o=options;o.maximum_order=1;o.ndf=false;return integrate_bdf(f,t,y,end,o,a,e,s);
}
BdfResult integrate_ssprk3(const OdeRhs& f,double t,const V& y,double end,const BdfOptions& o,
    const OdeAcceptedCallback& a,const std::vector<OdeEvent>& e,const OdeSampleCallback& s) {
  return integrate_rk(f,{},false,t,y,end,o,a,e,s);
}
BdfResult integrate_semi_implicit_ssprk3(const OdeRhs& f,const OdeRhs& explicit_f,double t,const V& y,double end,
    const BdfOptions& o,const OdeAcceptedCallback& a,const std::vector<OdeEvent>& e,const OdeSampleCallback& s) {
  return integrate_rk(f,explicit_f,true,t,y,end,o,a,e,s);
}
} // namespace cycle
