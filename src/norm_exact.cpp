// Exact streaming enumeration of contiguous normCCS partitions.
// R objects and R API calls are confined to the invoking main thread.
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __FAST_MATH__
#error "normCCS exact scoring must be compiled without fast-math"
#endif
// Preserve legacy trapezoid intermediates even when the target supports FMA.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif
using namespace Rcpp;
namespace ccs_exact {
using U64=std::uint64_t;
struct Model {
  int m,nc;
  std::vector<double> loss,first,cp,cn,base;
  std::vector<int> median,direction,real;
  std::vector<std::string> names,keys;
  static std::string utf8(SEXP s) {
    if(s==NA_STRING)stop("Native exact model strings must not be NA.");
    return std::string(Rf_translateCharUTF8(s));
  }
  inline std::size_t idx(int a,int b,int c) const noexcept {
    return (std::size_t(a-1)+std::size_t(b-1)*m)*nc+c;
  }
  inline std::size_t keyidx(int a,int b) const noexcept {
    return std::size_t(a-1)+std::size_t(b-1)*m;
  }
  explicit Model(List x,bool need_keys=true) {
    m=as<int>(x["m"]);
    if(m<2||m>512)stop("Native exact model m must be between 2 and 512.");
    NumericVector b=x["base_auc"];
    nc=b.size();if(nc<2||nc>512)stop("Native exact model must contain All and 1..511 real cohorts.");
    base.assign(b.begin(),b.end());
    NumericVector l=x["loss"],f=x["first"];
    IntegerVector ld=l.attr("dim"),fd=f.attr("dim");
    if(ld.size()!=3||fd.size()!=3||ld[0]!=nc||ld[1]!=m||ld[2]!=m||fd[0]!=nc||fd[1]!=m||fd[2]!=m)
      stop("Invalid native exact interval array dimensions.");
    loss.assign(l.begin(),l.end());first.assign(f.begin(),f.end());
    NumericMatrix p=x["cp"],n=x["cn"];
    if(p.nrow()!=nc||p.ncol()!=m+1||n.nrow()!=nc||n.ncol()!=m+1)
      stop("Invalid native exact prefix count dimensions.");
    cp.resize(std::size_t(nc)*(m+1));cn.resize(std::size_t(nc)*(m+1));
    for(int c=0;c<nc;c++)for(int j=0;j<=m;j++) {
      double pv=p(c,j),nv=n(c,j);
      if(!std::isfinite(pv)||!std::isfinite(nv)||pv<0||nv<0||pv!=std::floor(pv)||nv!=std::floor(nv)||pv>9007199254740991.0||nv>9007199254740991.0)
        stop("Native exact class counts must be finite nonnegative exactly represented integers.");
      if((j==0&&(pv!=0||nv!=0))||(j>0&&(pv<p(c,j-1)||nv<n(c,j-1))))
        stop("Native exact class prefixes must start at zero and be nondecreasing.");
      cp[c*(m+1)+j]=pv;cn[c*(m+1)+j]=nv;
    }
    IntegerMatrix md=x["median_states"];
    if(md.nrow()!=nc||md.ncol()!=4)stop("Invalid native exact median-state dimensions.");
    median.resize(std::size_t(4)*nc);
    for(int c=0;c<nc;c++)for(int j=0;j<4;j++) {
      int s=md(c,j);if(s<1||s>m)stop("Invalid native exact median state.");median[c*4+j]=s;
    }
    CharacterVector nm=x["cohorts"];
    if(nm.size()!=nc)stop("Invalid native exact cohort name count.");
    std::unordered_set<std::string> seen;
    int pooled=0;
    for(int c=0;c<nc;c++) {
      std::string name=utf8(nm[c]);
      if(name.empty()||!seen.insert(name).second)stop("Native exact cohort names must be nonempty and unique.");
      names.push_back(name);if(name=="All")pooled++;else real.push_back(c);
      double pt=cp[c*(m+1)+m],nt=cn[c*(m+1)+m];
      if(!std::isfinite(base[c])||base[c]<-1e-10||base[c]>1+1e-10||pt<=0||nt<=0)
        stop("Native exact cohorts must have both classes and finite base AUC in [0,1].");
      for(int j=0;j<4;j++) {
        const double* pref=(j<2?&cp[c*(m+1)]:&cn[c*(m+1)]);
        double total=j<2?pt:nt,target=j%2==0?std::ceil(total/2):std::floor(total/2)+1;
        int expected=int(std::lower_bound(pref+1,pref+m+1,target)-pref);
        if(median[c*4+j]!=expected)stop("Native exact median states disagree with prefix counts.");
      }
      for(int a=1;a<=m;a++)for(int z=a;z<=m;z++) {
        double v=loss[idx(a,z,c)],ix=first[idx(a,z,c)];
        double count=cp[c*(m+1)+z]-cp[c*(m+1)+a-1]+cn[c*(m+1)+z]-cn[c*(m+1)+a-1];
        if(!std::isfinite(v)||std::isnan(ix)||ix<1||(count>0&&!std::isfinite(ix))||(count==0&&ix!=std::numeric_limits<double>::infinity())||(std::isfinite(ix)&&(ix!=std::floor(ix)||ix>pt+nt)))
          stop("Invalid native exact interval loss or first-observation index.");
      }
    }
    if(pooled!=1||real.empty())stop("Native exact model must contain exactly one All and at least one real cohort.");
    direction.assign(std::size_t(m+1)*nc,0);
    for(int t=1;t<m;t++)for(int c=0;c<nc;c++)
      direction[t*nc+c]=int(t<median[c*4+2])+int(t<median[c*4+3])-int(t<median[c*4])-int(t<median[c*4+1]);
    if(!need_keys)return;
    keys.resize(std::size_t(m)*m);
    if(x.containsElementNamed("ranked_labels")) {
      CharacterVector labels=x["ranked_labels"];
      if(labels.size()!=m)stop("ranked_labels must have one label per state.");
      std::vector<std::string> label(m);seen.clear();
      for(int i=0;i<m;i++){label[i]=utf8(labels[i]);if(label[i].empty()||!seen.insert(label[i]).second)stop("ranked_labels must be nonempty and unique.");}
      for(int a=1;a<=m;a++){
        std::string joined;
        for(int z=a;z<=m;z++){if(z>a)joined+=",";joined+=label[z-1];keys[keyidx(a,z)]=joined+"; ";}
      }
    } else if(x.containsElementNamed("interval_keys")) {
      CharacterMatrix supplied=x["interval_keys"];
      if(supplied.nrow()!=m||supplied.ncol()!=m)stop("interval_keys must be m by m.");
      for(int a=1;a<=m;a++)for(int z=a;z<=m;z++){
        std::string key=utf8(supplied(a-1,z-1));
        if(key.size()<3||key.compare(key.size()-2,2,"; ")!=0)stop("Native exact interval keys must end in '; '.");
        keys[keyidx(a,z)]=std::move(key);
      }
    } else stop("Native exact model requires ranked_labels or interval_keys.");
  }
};
struct Combinations {
  std::vector<std::array<U64,8>> table;
  explicit Combinations(int n):table(n+1) {
    for(int i=0;i<=n;i++){
      table[i].fill(0);table[i][0]=1;
      for(int j=1;j<=7&&j<=i;j++)table[i][j]=table[i-1][j-1]+table[i-1][j];
    }
  }
  U64 get(int n,int q) const noexcept {return n<0||q<0||q>n||q>7?0:table[n][q];}
  bool unrank(U64 rank,int n,int q,std::array<int,7>& cuts) const noexcept {
    int low=1;
    for(int j=0;j<q;j++){
      bool assigned=false;
      for(int v=low;v<=n-(q-j-1);v++){
        U64 block=get(n-v,q-j-1);
        if(rank<block){cuts[j]=v;low=v+1;assigned=true;break;}rank-=block;
      }
      if(!assigned)return false;
    }
    return true;
  }
};
inline int next(std::array<int,7>& cut,int n,int q) noexcept {
  int j=q-1;
  while(j>=0&&cut[j]==n-(q-1-j))--j;
  if(j<0)return -1;
  ++cut[j];
  for(int a=j+1;a<q;a++)cut[a]=cut[a-1]+1;
  return j;
}
inline double binary_med(double x,double total) noexcept {return double(x>=total/2)/2+double(x>total/2)/2;}
struct Scorer {
  const Model& model; int k; bool rank, optimized;
  std::vector<long double> prefix;
  std::vector<int> direction_prefix;
  std::vector<double> auc,work;
  Scorer(const Model& m,int kk,bool rr,bool op):model(m),k(kk),rank(rr&&kk>2),optimized(op),prefix((kk+1)*m.nc,0),direction_prefix(kk*m.nc,0),auc(m.nc),work(m.nc){}
  void eval(const std::array<int,7>& cut,int changed=0){
    int nc=model.nc,q=k-1,m=model.m;
    if(!rank){
      if(optimized){
        for(int g=changed;g<k;g++){
          int a=g==0?1:cut[g-1]+1,b=g==q?m:cut[g];
          const double* losses=&model.loss[model.idx(a,b,0)];
          for(int c=0;c<nc;c++)prefix[(g+1)*nc+c]=prefix[g*nc+c]+(long double)losses[c];
          if(g<q){
            const int* contribution=&model.direction[cut[g]*nc];
            for(int c=0;c<nc;c++)direction_prefix[(g+1)*nc+c]=direction_prefix[g*nc+c]+contribution[c];
          }
        }
      }
      for(int c=0;c<nc;c++){
        long double sum=0;
        if(optimized)sum=prefix[k*nc+c];
        else for(int g=0;g<k;g++)sum+=(long double)model.loss[model.idx(g==0?1:cut[g-1]+1,g==q?m:cut[g],c)];
        double a=model.base[c]-(double)sum;
        bool reverse;
        if(optimized)reverse=direction_prefix[q*nc+c]<0;
        else {
          int md[4]={0,0,0,0};
          for(int j=0;j<4;j++)for(int g=0;g<q;g++)md[j]+=cut[g]<model.median[c*4+j];
          reverse=md[2]+md[3]<md[0]+md[1];
        }
        auc[c]=reverse?1-a:a;
      }
      return;
    }
    for(int c=0;c<nc;c++){
      double p[8],n[8],ovr[8],first[8];int order[8];
      double pt=model.cp[c*(m+1)+m],nt=model.cn[c*(m+1)+m];
      for(int g=0;g<k;g++){
        int a=g==0?1:cut[g-1]+1,b=g==q?m:cut[g],off=c*(m+1);
        p[g]=model.cp[off+b]-model.cp[off+a-1];n[g]=model.cn[off+b]-model.cn[off+a-1];
        double se=p[g]/pt,sp=(nt-n[g])/nt;
        double score=((1+se)/2)*sp+(se/2)*(1-sp);
        if(binary_med(n[g],nt)>binary_med(p[g],pt)){
          double rs=(pt-p[g])/pt,rsp=n[g]/nt;
          score=((1+rs)/2)*rsp+(rs/2)*(1-rsp);
        }
        ovr[g]=p[g]+n[g]==0?std::numeric_limits<double>::infinity():score;
        first[g]=model.first[model.idx(a,b,c)];order[g]=g;
      }
      // K<=8: insertion sort avoids libstdc++'s speculative 16-item
      // sorting path and preserves legacy AUC/first-observation tie order.
      auto precedes=[&](int a,int b){
        if(ovr[a]!=ovr[b])return ovr[a]<ovr[b];
        if(first[a]!=first[b])return first[a]<first[b];
        return a<b;
      };
      for(int j=1;j<k;j++){
        int item=order[j],position=j;
        while(position>0&&precedes(item,order[position-1])){
          order[position]=order[position-1];--position;
        }
        order[position]=item;
      }
      double cumulative=0,pc=0;long double acc=0;
      double mp=1,mn=1;
      for(int j=0;j<k;j++){
        int g=order[j];cumulative+=n[g];
        double term=p[g]*(cumulative-n[g]/2);acc+=(long double)term;
        // norm_median_bin counts cumulative values strictly below each middle position.
        pc+=p[g];
        mp+=double(pc<std::ceil(pt/2))/2+double(pc<std::floor(pt/2)+1)/2;
        mn+=double(cumulative<std::ceil(nt/2))/2+double(cumulative<std::floor(nt/2)+1)/2;
      }
      double a=(double)acc/(pt*nt);auc[c]=mn>mp?1-a:a;
    }
  }
  static double median_value(std::vector<double>& values,int count) {
    int h=count/2;
    std::nth_element(values.begin(),values.begin()+h,values.begin()+count);
    if(count%2)return values[h];
    double lo=*std::max_element(values.begin(),values.begin()+h);
    return (lo+values[h])/2;
  }
  void aggregate(int objective,double& score,double& mean) {
    long double sum=0;
    for(int c:model.real)sum+=(long double)auc[c];
    int count=(int)model.real.size();
    mean=(double)(sum/(long double)count);
    if(objective==0){std::copy(auc.begin(),auc.end(),work.begin());score=median_value(work,model.nc);return;}
    if(objective==1){score=mean;return;}
    for(int j=0;j<count;j++)work[j]=auc[model.real[j]];
    if(objective==4){score=median_value(work,count);return;}
    std::sort(work.begin(),work.begin()+count);
    double mass=.25*count;int whole=(int)std::floor(mass);double fraction=mass-whole;
    long double lower_sum=0;
    for(int j=0;j<whole;j++)lower_sum+=(long double)work[j];
    double tail=((double)lower_sum+(fraction>0?fraction*work[whole]:0.0))/mass;
    score=objective==2?tail:.5*mean+.5*tail;
  }
};
struct Candidate {double score,mean;std::array<int,7> cuts;};
struct InfoCursor {
  const Model& model;const std::array<int,7>& cuts;int q,grade;std::size_t offset;
  InfoCursor(const Model& m,const std::array<int,7>& c,int qq):model(m),cuts(c),q(qq),grade(0),offset(0){}
  bool next_byte(unsigned char& out) noexcept {
    while(grade<=q){
      if(offset<2){out=offset==0?static_cast<unsigned char>('0'+grade):static_cast<unsigned char>('=');++offset;return true;}
      int g=q-grade,a=g==0?1:cuts[g-1]+1,b=g==q?model.m:cuts[g];
      const std::string& key=model.keys[model.keyidx(a,b)];
      std::size_t size=key.size()-(grade==q?2:0);
      if(offset-2<size){out=static_cast<unsigned char>(key[offset-2]);++offset;return true;}
      ++grade;offset=0;
    }
    return false;
  }
};
struct Better {
  const Model* model;int q;bool secondary;
  bool operator()(const Candidate& a,const Candidate& b) const noexcept {
    if(a.score!=b.score)return a.score>b.score;
    if(secondary&&a.mean!=b.mean)return a.mean>b.mean;
    InfoCursor ac(*model,a.cuts,q),bc(*model,b.cuts,q);
    unsigned char av=0,bv=0;
    while(true){bool has_a=ac.next_byte(av),has_b=bc.next_byte(bv);if(!has_a||!has_b){if(has_a!=has_b)return !has_a;break;}if(av!=bv)return av<bv;}
    return std::lexicographical_compare(a.cuts.begin(),a.cuts.begin()+q,b.cuts.begin(),b.cuts.begin()+q);
  }
};
inline void keep(std::vector<Candidate>& heap,const Candidate& candidate,int limit,const Better& better) {
  // With the better-than comparator, heap.front() is the worst retained item.
  if((int)heap.size()<limit){heap.push_back(candidate);std::push_heap(heap.begin(),heap.end(),better);}
  else if(better(candidate,heap.front())){std::pop_heap(heap.begin(),heap.end(),better);heap.back()=candidate;std::push_heap(heap.begin(),heap.end(),better);}
}
inline int objective_id(const std::string& s){if(s=="median")return 0;if(s=="mean")return 1;if(s=="lower_tail")return 2;if(s=="mean_lower_tail")return 3;if(s=="cohort_median")return 4;stop("Unsupported native exact objective.");return 0;}
inline bool rank_id(const std::string& s){if(s=="Raw")return false;if(s=="OptimizeRank")return true;stop("Native exact method must be Raw or OptimizeRank.");return false;}
inline void validate_cuts(const IntegerMatrix& x,int m){
  if(x.ncol()<1||x.ncol()>7||x.ncol()>=m)stop("Native exact cuts must define K between 2 and min(8,m).");
  for(int i=0;i<x.nrow();i++){int previous=0;for(int j=0;j<x.ncol();j++){int value=x(i,j);if(value<=previous||value>=m)stop("Native exact cuts must increase strictly within 1..m-1.");previous=value;}}
}
struct alignas(64) Worker {
  Scorer scorer;std::vector<Candidate> heap;U64 count;
  Worker(const Model& model,int k,bool rank,int limit):scorer(model,k,rank,true),count(0){heap.reserve(limit);}
};
} // namespace ccs_exact

// [[Rcpp::export]]
NumericMatrix norm_exact_auc_cpp(IntegerMatrix cuts,List model,std::string method="Raw") {
  using namespace ccs_exact;
  Model data(model,false);validate_cuts(cuts,data.m);bool rank=rank_id(method);
  NumericMatrix result(cuts.nrow(),data.nc);Scorer scorer(data,cuts.ncol()+1,rank,false);
  for(int i=0;i<cuts.nrow();i++){
    if(i%4096==0)Rcpp::checkUserInterrupt();
    std::array<int,7> cut{};for(int j=0;j<cuts.ncol();j++)cut[j]=cuts(i,j);
    scorer.eval(cut);for(int c=0;c<data.nc;c++)result(i,c)=scorer.auc[c];
  }
  result.attr("dimnames")=List::create(R_NilValue,wrap(data.names));return result;
}

// [[Rcpp::export]]
List norm_exact_cpp(List model,int k,std::string method,std::string objective,int top_n=2000,int threads=1) {
  using namespace ccs_exact;
  const auto started=std::chrono::steady_clock::now();
  Model data(model);bool rank=rank_id(method);int obj=objective_id(objective);
  if(k<2||k>8||k>data.m)stop("Native exact K must be between 2 and min(8,m).");
  if(threads<1||threads>1024)stop("Native exact thread count must be between 1 and 1024.");
  if(top_n<1)stop("Native exact top_n must be positive.");
  Combinations combinations(data.m-1);int q=k-1;
  U64 total=combinations.get(data.m-1,q);
  if(total==0||total>9007199254740991ULL)stop("Native exact combination count is outside the supported exact range.");
  int limit=int(std::min(U64(top_n),total));
  const U64 batch=4096,wave_size=batch*256;
  U64 all_chunks=(total+batch-1)/batch;
  int worker_count=int(std::min(U64(threads),all_chunks));
  bool openmp=false;
#ifdef _OPENMP
  openmp=true;
#else
  worker_count=1;
#endif
  std::vector<std::unique_ptr<Worker>> workers;
  workers.reserve(worker_count);
  for(int t=0;t<worker_count;t++)workers.emplace_back(new Worker(data,k,rank,limit));
  Better better{&data,q,obj!=0};
  int actual_threads=1;
  for(U64 wave_start=0;wave_start<total;wave_start+=wave_size){
    Rcpp::checkUserInterrupt();
    U64 wave_end=std::min(total,wave_start+wave_size);
    int blocks=int((wave_end-wave_start+batch-1)/batch),wave_actual=1;
    std::atomic<bool> failed(false);
    char worker_error[512]={0};
#ifdef _OPENMP
#pragma omp parallel num_threads(worker_count)
#endif
    {
      int tid=0;
#ifdef _OPENMP
      tid=omp_get_thread_num();
#pragma omp single
      {wave_actual=omp_get_num_threads();}
#endif
      Worker& worker=*workers[tid];Scorer& scorer=worker.scorer;
      std::vector<Candidate> local_heap=std::move(worker.heap);
      U64 local_count=0;std::array<int,7> cuts{};
#ifdef _OPENMP
#pragma omp for schedule(dynamic,1)
#endif
      for(int block=0;block<blocks;block++){
        if(failed.load(std::memory_order_relaxed))continue;
        try {
          U64 lo=wave_start+U64(block)*batch,len=std::min(batch,wave_end-lo);
          if(!combinations.unrank(lo,data.m-1,q,cuts)){failed.store(true,std::memory_order_relaxed);continue;}
          int changed=0;
          for(U64 j=0;j<len;j++){
            scorer.eval(cuts,changed);double score,mean;scorer.aggregate(obj,score,mean);
            if(!std::isfinite(score)||!std::isfinite(mean)){failed.store(true,std::memory_order_relaxed);break;}
            keep(local_heap,Candidate{score,mean,cuts},limit,better);++local_count;
            if(j+1<len){changed=next(cuts,data.m-1,q);if(changed<0){failed.store(true,std::memory_order_relaxed);break;}}
          }
        } catch(const std::exception& error) {
          bool expected=false;
          if(failed.compare_exchange_strong(expected,true,std::memory_order_relaxed))
            std::snprintf(worker_error,sizeof(worker_error),"%s",error.what());
        } catch(...) {failed.store(true,std::memory_order_relaxed);}
      }
      worker.heap=std::move(local_heap);worker.count+=local_count;
    }
    if(failed.load(std::memory_order_relaxed)){
      if(worker_error[0])stop("Native exact worker failed: %s",worker_error);
      stop("Native exact worker failed; no complete search result was returned.");
    }
    actual_threads=std::max(actual_threads,wave_actual);
  }
  Rcpp::checkUserInterrupt();
  std::vector<Candidate> heap;heap.reserve(limit);U64 visited=0;
  for(const auto& worker:workers){visited+=worker->count;for(const Candidate& candidate:worker->heap)keep(heap,candidate,limit,better);}
  if(visited!=total)stop("Native exact search did not cover the full combination domain.");
  std::sort(heap.begin(),heap.end(),better);
  NumericVector scores(heap.size()),means(heap.size());IntegerMatrix cuts(heap.size(),q);
  for(int i=0;i<(int)heap.size();i++){scores[i]=heap[i].score;means[i]=heap[i].mean;for(int j=0;j<q;j++)cuts(i,j)=heap[i].cuts[j];}
  const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  return List::create(_["cuts"]=cuts,_["scores"]=scores,_["means"]=means,_["visited"]=double(visited),_["total_candidates"]=double(total),_["complete"]=true,_["search_elapsed_seconds"]=elapsed,_["actual_threads"]=actual_threads,_["openmp"]=openmp,_["requested_threads"]=threads,_["top_n"]=limit,_["chunk_size"]=double(batch),_["interrupt_wave_size"]=double(wave_size));
}
