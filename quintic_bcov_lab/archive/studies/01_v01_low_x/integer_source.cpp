// Optional GMP accelerator for the unchanged direct BCOV genus convolution.
// Build: g++ -O3 -shared -fPIC integer_source.cpp -o integer_source.so -lgmp
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>
#if __has_include(<gmp.h>)
#include <gmp.h>
#else
// GMP's public ABI, for a runtime that has libgmp but no development header.
extern "C" {
struct __mpz_struct { int _mp_alloc; int _mp_size; unsigned long *_mp_d; };
typedef __mpz_struct mpz_t[1];
void __gmpz_init(mpz_t); void __gmpz_clear(mpz_t);
void __gmpz_set(mpz_t,const mpz_t); int __gmpz_set_str(mpz_t,const char*,int);
void __gmpz_mul(mpz_t,const mpz_t,const mpz_t);
void __gmpz_mul_ui(mpz_t,const mpz_t,unsigned long);
void __gmpz_addmul(mpz_t,const mpz_t,const mpz_t);
void __gmpz_lcm(mpz_t,const mpz_t,const mpz_t);
void __gmpz_gcd(mpz_t,const mpz_t,const mpz_t);
void __gmpz_divexact(mpz_t,const mpz_t,const mpz_t);
void __gmpz_neg(mpz_t,const mpz_t);
size_t __gmpz_sizeinbase(const mpz_t,int);
char *__gmpz_get_str(char*,int,const mpz_t);
}
#define mpz_init __gmpz_init
#define mpz_clear __gmpz_clear
#define mpz_set __gmpz_set
#define mpz_set_str __gmpz_set_str
#define mpz_mul __gmpz_mul
#define mpz_mul_ui __gmpz_mul_ui
#define mpz_addmul __gmpz_addmul
#define mpz_lcm __gmpz_lcm
#define mpz_gcd __gmpz_gcd
#define mpz_divexact __gmpz_divexact
#define mpz_neg __gmpz_neg
#define mpz_sizeinbase __gmpz_sizeinbase
#define mpz_get_str __gmpz_get_str
#endif

struct Z {
 mpz_t z;
 Z(){mpz_init(z);}
 Z(const Z &b){mpz_init(z);mpz_set(z,b.z);}
 Z& operator=(const Z&b){mpz_set(z,b.z);return *this;}
 ~Z(){mpz_clear(z);}
 std::string str() const {
  std::vector<char> buf(mpz_sizeinbase(z,10)+3);
  mpz_get_str(buf.data(),10,z);return std::string(buf.data());
 }
};
struct Polynomial {std::vector<int> codes,degrees;std::vector<Z> values;Z den;};
static std::map<int,Polynomial> vertices;
static int degree(int code){code/=8;int d=0;for(int i=0;i<4;++i){d+=code%8;code/=8;}return d;}
extern "C" void reset_vertices(){vertices.clear();}
extern "C" void add_vertex(int g,int count,const int*codes,const char**numbers,const char*den){
 Polynomial p;p.codes.assign(codes,codes+count);p.values.resize(count);
 for(int i=0;i<count;++i){p.degrees.push_back(degree(codes[i]));mpz_set_str(p.values[i].z,numbers[i],10);}
 mpz_set_str(p.den.z,den,10);vertices[g]=p;
}
extern "C" char* compute_source(int g){
 const auto &previous=vertices.at(-1);Z den,product,factor,temp;
 mpz_set(den.z,previous.den.z);
 for(int h=1;h<=g/2;++h){
  mpz_mul(product.z,vertices.at(h).den.z,vertices.at(g-h).den.z);
  mpz_lcm(den.z,den.z,product.z);
 }
 std::map<int,Z> totals;
 mpz_divexact(factor.z,den.z,previous.den.z);
 for(size_t i=0;i<previous.codes.size();++i)
  mpz_mul(totals[previous.codes[i]].z,previous.values[i].z,factor.z);
 for(int h=1;h<=g/2;++h){
  const auto &a=vertices.at(h),&b=vertices.at(g-h);
  mpz_mul(product.z,a.den.z,b.den.z);mpz_divexact(factor.z,den.z,product.z);
  if(2*h!=g)mpz_mul_ui(factor.z,factor.z,2);
  std::map<int,Z> partial;
  for(size_t i=0;i<a.codes.size();++i)for(size_t j=0;j<b.codes.size();++j){
   if(a.degrees[i]+b.degrees[j]>3)continue;
   mpz_addmul(partial[a.codes[i]+b.codes[j]].z,a.values[i].z,b.values[j].z);
  }
  for(auto &entry:partial)mpz_addmul(totals[entry.first].z,entry.second.z,factor.z);
 }
 mpz_mul_ui(den.z,den.z,2);
 std::string out="REDUCED\n";Z divisor,reduced_den;
 for(auto &entry:totals){
  mpz_neg(temp.z,entry.second.z);
  mpz_gcd(divisor.z,temp.z,den.z);
  mpz_divexact(temp.z,temp.z,divisor.z);
  mpz_divexact(reduced_den.z,den.z,divisor.z);
  out+=std::to_string(entry.first)+" "+temp.str()+" "+reduced_den.str()+"\n";
 }
 char *result=(char*)std::malloc(out.size()+1);std::memcpy(result,out.c_str(),out.size()+1);return result;
}
extern "C" void free_result(char*p){std::free(p);}
