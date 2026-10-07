// Artemis - portable C++11 edition of Artemis3.py
// Original: Jamie J. Alnasir, 2014, Royal Holloway University of London, CSSB.
// Copyright (c) 2014-2026 Jamie J. Alnasir, All Rights Reserved.
// C++ port: October 2026. No third-party libraries are required.
// The only platform-specific call detects whether standard input is a terminal.
// Build: c++ -std=c++11 -O2 -Wall -Wextra -pedantic artemis.cpp -o artemis
// Define ARTEMIS_NO_MAIN to include this file in another C++ program.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cstdio>
#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace artemis {
using std::string;
static const std::size_t npos = std::numeric_limits<std::size_t>::max();

inline string trim(const string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n\f\v");
    return a == string::npos ? "" : s.substr(a, s.find_last_not_of(" \t\r\n\f\v")-a+1);
}
inline string field(const string& s, std::size_t a, std::size_t b) {
    return a < s.size() ? s.substr(a, b-a) : "";
}
inline string upper(string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
inline long integer(const string& value) {
    string s=trim(value); std::size_t i=0;
    if (s.empty()) throw std::invalid_argument("Empty integer field");
    if (s[0]=='+' || s[0]=='-') i=1;
    if (i==s.size()) throw std::invalid_argument("Invalid integer: " + s);
    for (; i<s.size(); ++i)
        if (s[i]<'0' || s[i]>'9') throw std::invalid_argument("Invalid integer: " + s);
    std::istringstream stream(s); stream.imbue(std::locale::classic()); long result;
    if (!(stream >> result)) throw std::invalid_argument("Integer out of range: " + s);
    return result;
}
inline double real(const string& value) {
    std::istringstream stream(trim(value)); stream.imbue(std::locale::classic());
    double result; char extra;
    if (!(stream >> result) || (stream >> extra) || !std::isfinite(result))
        throw std::invalid_argument("Invalid finite numeric field: " + value);
    return result;
}
inline string fixed(double value, int precision) {
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite numeric value");
    std::ostringstream out; out.imbue(std::locale::classic());
    out << std::fixed << std::setprecision(precision) << value; return out.str();
}
inline string canonical(string s) {
    std::replace(s.begin(),s.end(),'*','\'');
    if (s=="O1P") return "OP1";
    if (s=="O2P") return "OP2";
    if (s=="O3P") return "OP3";
    return s;
}
inline string nucleotide(string s) {
    s=upper(trim(s));
    if (s.size()==2 && s[0]=='D' && string("ATGCUI").find(s[1])!=string::npos) s=s.substr(1);
    return s.size()==1 && string("ATGCUI").find(s[0])!=string::npos ? s : "";
}
enum class ResidueType { Unknown, AA, DNA }; // DNA includes recognised RNA, as in Artemis3.
enum class Hybridisation { Unknown, SP1, SP2, SP3 };
inline ResidueType residue_type(const string& name) {
    static const std::set<string> amino={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
    if (!nucleotide(name).empty()) return ResidueType::DNA;
    return amino.count(upper(trim(name))) ? ResidueType::AA : ResidueType::Unknown;
}
using Vec3 = std::array<double,3>;
inline Vec3 subtract(const Vec3& a,const Vec3& b) { return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}}; }
inline double dot(const Vec3& a,const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
inline Vec3 cross(const Vec3& a,const Vec3& b) { return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}}; }
inline double distance(const Vec3& a,const Vec3& b) { Vec3 d=subtract(a,b); return std::sqrt(dot(d,d)); }
inline double degrees(double radians) { return radians*(180.0/std::acos(-1.0)); }
inline double angle(const Vec3& centre,const Vec3& a,const Vec3& b) {
    Vec3 v=subtract(a,centre),w=subtract(b,centre);
    double length=std::sqrt(dot(v,v))*std::sqrt(dot(w,w));
    if (!(length>0.0) || !std::isfinite(length)) throw std::invalid_argument("Undefined angle");
    double cosine=dot(v,w)/length;
    if (!std::isfinite(cosine)) throw std::invalid_argument("Undefined angle");
    return degrees(std::acos(std::max(-1.0,std::min(1.0,cosine))));
}
inline double dihedral(const Vec3& a,const Vec3& b,const Vec3& c,const Vec3& d) {
    Vec3 b0=subtract(a,b),b1=subtract(c,b),b2=subtract(d,c);
    double length=std::sqrt(dot(b1,b1));
    if (!(length>1e-12) || !std::isfinite(length)) throw std::invalid_argument("Undefined central bond");
    Vec3 unit={{b1[0]/length,b1[1]/length,b1[2]/length}},v,w;
    double d0=dot(b0,unit),d2=dot(b2,unit);
    for (std::size_t i=0;i<3;++i) { v[i]=b0[i]-d0*unit[i]; w[i]=b2[i]-d2*unit[i]; }
    if (!(dot(v,v)>1e-24) || !(dot(w,w)>1e-24)) throw std::invalid_argument("Collinear dihedral atoms");
    double x=dot(v,w),y=dot(cross(unit,v),w);
    if (!std::isfinite(x) || !std::isfinite(y)) throw std::invalid_argument("Non-finite dihedral");
    return degrees(std::atan2(y,x));
}

struct Atom {
    string serial,name,alt_loc,chain_id,res_name,res_seq,icode,x,y,z,occ,temp,element,charge;
    Hybridisation hybridisation=Hybridisation::Unknown;
    Vec3 coords() const { return {{real(x),real(y),real(z)}}; }
    static Atom parse(const string& line) {
        if (line.size()<54) throw std::invalid_argument("Short ATOM record");
        Atom a;
        a.serial=trim(field(line,6,11)); a.name=trim(field(line,12,16));
        a.alt_loc=trim(field(line,16,17)); a.res_name=trim(field(line,17,20));
        a.chain_id=trim(field(line,21,22)); a.res_seq=trim(field(line,22,26));
        a.icode=trim(field(line,26,27)); a.x=trim(field(line,30,38));
        a.y=trim(field(line,38,46)); a.z=trim(field(line,46,54));
        a.occ=trim(field(line,54,60)); a.temp=trim(field(line,60,66));
        a.element=trim(field(line,76,78)); a.charge=trim(field(line,78,80));
        integer(a.serial); integer(a.res_seq); a.coords();
        if (!a.occ.empty()) real(a.occ);
        if (!a.temp.empty()) real(a.temp);
        if (a.name.empty() || a.res_name.empty()) throw std::invalid_argument("Missing atom/residue name");
        return a;
    }
    bool same_identity(const Atom& b) const {
        return std::tie(name,alt_loc,chain_id,res_name,res_seq,icode)==
               std::tie(b.name,b.alt_loc,b.chain_id,b.res_name,b.res_seq,b.icode);
    }
};
struct BondSpec { const char* first; const char* second; int order; bool backbone; };
// Bond atom indices refer to the owning residue's atoms, avoiding dangling pointers.
struct Bond { std::size_t atom1,atom2; int order; bool backbone; };
// Residue templates ported from Artemis3.py; no runtime Python dependency.
static const std::map<std::string, std::vector<BondSpec> >& bond_templates() {
    static const std::map<std::string, std::vector<BondSpec> > table = {
        {"A", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"N9", "C8", 1, false},
            {"C8", "N7", 2, false},
            {"N7", "C5", 1, false},
            {"C5", "C6", 1, false},
            {"C6", "N6", 1, false},
            {"C6", "N1", 2, false},
            {"N1", "C2", 1, false},
            {"C2", "N3", 2, false},
            {"C5", "C4", 2, false},
            {"C4", "N9", 1, false},
            {"C4", "N3", 1, false},
            {"N9", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"ALA", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"ARG", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CD", "NE", 1, false},
            {"NE", "CZ", 1, false},
            {"CZ", "NH1", 2, false},
            {"CZ", "NH2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"ASN", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "OD1", 2, false},
            {"CG", "ND2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"ASP", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "OD1", 2, false},
            {"CG", "OD2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"C", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"C6", "C5", 2, false},
            {"C5", "C4", 1, false},
            {"C4", "N3", 2, false},
            {"C4", "N4", 1, false},
            {"N3", "C2", 1, false},
            {"O2", "C2", 2, false},
            {"C2", "N1", 1, false},
            {"N1", "C6", 1, false},
            {"N1", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"CYS", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CB", "SG", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"G", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"N9", "C8", 1, false},
            {"C8", "N7", 2, false},
            {"N7", "C5", 1, false},
            {"C5", "C6", 1, false},
            {"C6", "O6", 2, false},
            {"C6", "N1", 1, false},
            {"N1", "C2", 1, false},
            {"C2", "N3", 2, false},
            {"C5", "C4", 2, false},
            {"C4", "N9", 1, false},
            {"C4", "N3", 1, false},
            {"C2", "N2", 1, false},
            {"N9", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"GLN", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CD", "OE1", 2, false},
            {"CD", "NE2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"GLU", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CD", "OE1", 2, false},
            {"CD", "OE2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"GLY", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"HIS", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "CD2", 2, false},
            {"CG", "ND1", 1, false},
            {"ND1", "CE1", 2, false},
            {"CD2", "NE2", 1, false},
            {"NE2", "CE1", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"I", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"N9", "C8", 1, false},
            {"C8", "N7", 2, false},
            {"N7", "C5", 1, false},
            {"C5", "C6", 1, false},
            {"C6", "O6", 2, false},
            {"C6", "N1", 1, false},
            {"N1", "C2", 1, false},
            {"C2", "N3", 2, false},
            {"C5", "C4", 2, false},
            {"C4", "N9", 1, false},
            {"C4", "N3", 1, false},
            {"N9", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"ILE", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CB", "CG1", 1, false},
            {"CB", "CG2", 1, false},
            {"CG1", "CD1", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"LEU", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "CD1", 1, false},
            {"CG", "CD2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"LYS", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CD", "CE", 1, false},
            {"CE", "NZ", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"MET", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "SD", 1, false},
            {"SD", "CE", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"PHE", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "CD1", 2, false},
            {"CG", "CD2", 1, false},
            {"CD1", "CE1", 1, false},
            {"CD2", "CE2", 2, false},
            {"CE1", "CZ", 2, false},
            {"CE2", "CZ", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"PRO", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CD", "N", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"SER", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CB", "OG", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"T", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"N1", "C2", 1, false},
            {"C2", "O2", 2, false},
            {"C2", "N3", 1, false},
            {"N3", "C4", 1, false},
            {"C4", "C5", 1, false},
            {"C4", "O4", 2, false},
            {"C5", "C6", 2, false},
            {"N1", "C6", 1, false},
            {"N1", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C5", "C7", 1, false},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"THR", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CB", "OG1", 1, false},
            {"CB", "CG2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"TRP", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "CD1", 2, false},
            {"CG", "CD2", 1, false},
            {"CD1", "NE1", 1, false},
            {"NE1", "CE2", 1, false},
            {"CE2", "CD2", 2, false},
            {"CE2", "CZ2", 1, false},
            {"CZ2", "CH2", 2, false},
            {"CH2", "CZ3", 1, false},
            {"CZ3", "CE3", 2, false},
            {"CE3", "CD2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"TYR", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CG", "CD1", 2, false},
            {"CG", "CD2", 1, false},
            {"CD1", "CE1", 1, false},
            {"CD2", "CE2", 2, false},
            {"CE1", "CZ", 2, false},
            {"CE2", "CZ", 1, false},
            {"CZ", "OH", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
        {"U", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"N1", "C2", 1, false},
            {"C2", "O2", 2, false},
            {"C2", "N3", 1, false},
            {"N3", "C4", 1, false},
            {"C4", "C5", 1, false},
            {"C4", "O4", 2, false},
            {"C5", "C6", 2, false},
            {"N1", "C6", 1, false},
            {"N1", "C1'", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C2'", "O2'", 1, false},
            {"P", "OP3", 1, false},
        }},
        {"VAL", {
            {"N", "CA", 1, true},
            {"CA", "C", 1, true},
            {"C", "O", 2, true},
            {"CA", "CB", 1, false},
            {"CB", "CG", 1, false},
            {"CG", "CD", 1, false},
            {"CB", "CG1", 1, false},
            {"CB", "CG2", 1, false},
            {"OP1", "P", 1, true},
            {"OP2", "P", 2, true},
            {"C1'", "C2'", 1, true},
            {"C2'", "C3'", 1, true},
            {"C3'", "C4'", 1, true},
            {"C4'", "O4'", 1, true},
            {"O4'", "C1'", 1, true},
            {"C4'", "C5'", 1, true},
            {"C3'", "O3'", 1, true},
            {"C5'", "O5'", 1, true},
            {"P", "O5'", 1, true},
            {"C", "OXT", 1, false},
        }},
    };
    return table;
}

struct Residue {
    string res_name,res_seq,icode,chain_id,model_id="1";
    int model_index=0,segment_id=0;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    ResidueType type() const { return residue_type(res_name); }
    using Selection=std::map<string,std::size_t>;
    Selection selected_atoms() const {
        std::map<string,std::pair<double,std::size_t> > scores;
        for (const Atom& a:atoms) if (!a.alt_loc.empty()) {
            auto& score=scores[a.alt_loc]; score.first+=a.occ.empty()?0.0:real(a.occ); ++score.second;
        }
        string label; double best=-std::numeric_limits<double>::infinity();
        for (const auto& item:scores) {
            double mean=item.second.first/static_cast<double>(item.second.second);
            if (mean>best || (mean==best && item.first=="A")) { best=mean; label=item.first; }
        }
        Selection selected;
        for (std::size_t i=0;i<atoms.size();++i) {
            const Atom& a=atoms[i];
            if (!a.alt_loc.empty() && a.alt_loc!=label) continue;
            string name=canonical(a.name); auto old=selected.find(name);
            if (old==selected.end() || (!atoms[old->second].alt_loc.empty() && a.alt_loc.empty())) selected[name]=i;
        }
        return selected;
    }
    void compute_bonds() {
        bonds.clear(); for (Atom& a:atoms) a.hybridisation=Hybridisation::Unknown;
        string name=nucleotide(res_name); if (name.empty()) name=upper(trim(res_name));
        auto definition=bond_templates().find(name); if (definition==bond_templates().end()) return;
        Selection selected=selected_atoms(); std::set<std::pair<std::size_t,std::size_t> > seen;
        for (const BondSpec& spec:definition->second) {
            auto first=selected.find(spec.first),second=selected.find(spec.second);
            if (first==selected.end() || second==selected.end()) continue;
            auto key=std::minmax(first->second,second->second);
            if (seen.insert(key).second) bonds.push_back({first->second,second->second,spec.order,spec.backbone});
        }
        static const std::map<string,std::set<string> > planar={
            {"ARG",{"NE","CZ","NH1","NH2"}}, {"ASN",{"CG","OD1","ND2"}},
            {"ASP",{"CG","OD1","OD2"}}, {"GLN",{"CD","OE1","NE2"}},
            {"GLU",{"CD","OE1","OE2"}}, {"PHE",{"CG","CD1","CD2","CE1","CE2","CZ"}},
            {"TYR",{"CG","CD1","CD2","CE1","CE2","CZ","OH"}},
            {"HIS",{"CG","ND1","CD2","CE1","NE2"}},
            {"TRP",{"CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"}}
        };
        static const std::set<string> base_planar={"N1","N2","N3","N4","N6","N7","N9","C2","C4","C5","C6","C8","O2","O4","O6"};
        for (const Bond& bond:bonds) if (bond.order==2) {
            atoms[bond.atom1].hybridisation=Hybridisation::SP2;
            atoms[bond.atom2].hybridisation=Hybridisation::SP2;
        }
        for (const auto& entry:selected) {
            Atom& a=atoms[entry.second]; const string& n=entry.first;
            if (type()==ResidueType::AA) {
                auto group=planar.find(name);
                if (n=="C" || n=="O" || n=="OXT" || (group!=planar.end() && group->second.count(n))) a.hybridisation=Hybridisation::SP2;
                else if (!n.empty() && n[0]=='C') a.hybridisation=Hybridisation::SP3;
            } else {
                if ((!n.empty() && n[0]=='C' && n.find('\'')!=string::npos) || n=="C7") a.hybridisation=Hybridisation::SP3;
                else if (base_planar.count(n)) a.hybridisation=Hybridisation::SP2;
            }
        }
    }
};
inline const Atom* get_atom(const Residue& r,const Residue::Selection& selected,const string& name) {
    auto it=selected.find(name); return it==selected.end()?nullptr:&r.atoms[it->second];
}
inline bool peptide_linked(const Residue& a,const Residue& b,const Residue::Selection& sa,const Residue::Selection& sb) {
    if (a.model_index!=b.model_index || a.segment_id!=b.segment_id || a.chain_id!=b.chain_id ||
        a.type()!=ResidueType::AA || b.type()!=ResidueType::AA) return false;
    const Atom* c=get_atom(a,sa,"C"); const Atom* n=get_atom(b,sb,"N");
    if (!c || !n || (!c->alt_loc.empty() && !n->alt_loc.empty() && c->alt_loc!=n->alt_loc)) return false;
    double length=distance(c->coords(),n->coords()); return length>0.5 && length<=2.0;
}
struct Torsions {
    double phi=std::numeric_limits<double>::quiet_NaN();
    double psi=std::numeric_limits<double>::quiet_NaN();
    double omega=std::numeric_limits<double>::quiet_NaN();
};
inline double atom_dihedral(const Atom* a,const Atom* b,const Atom* c,const Atom* d) {
    if (a && b && c && d) {
        try { return dihedral(a->coords(),b->coords(),c->coords(),d->coords()); }
        catch (const std::invalid_argument&) {}
    }
    return std::numeric_limits<double>::quiet_NaN();
}
using ResiduePtr=std::shared_ptr<Residue>;
inline Torsions residue_dihedrals(const std::vector<ResiduePtr>& residues,std::size_t index) {
    const Residue& r=*residues.at(index); Torsions result;
    if (r.type()!=ResidueType::AA) return result;
    auto selected=r.selected_atoms();
    const Atom* n=get_atom(r,selected,"N"),*ca=get_atom(r,selected,"CA"),*c=get_atom(r,selected,"C");
    if (index>0) {
        const Residue& before=*residues[index-1]; auto s=before.selected_atoms();
        if (peptide_linked(before,r,s,selected)) result.phi=atom_dihedral(get_atom(before,s,"C"),n,ca,c);
    }
    if (index+1<residues.size()) {
        const Residue& after=*residues[index+1]; auto s=after.selected_atoms();
        if (peptide_linked(r,after,selected,s)) {
            result.psi=atom_dihedral(n,ca,c,get_atom(after,s,"N"));
            result.omega=atom_dihedral(ca,c,get_atom(after,s,"N"),get_atom(after,s,"CA"));
        }
    }
    return result;
}
inline void put(string& line,std::size_t start,std::size_t width,const string& value,bool right=false) {
    if (value.size()>width || value.find_first_of("\r\n")!=string::npos)
        throw std::invalid_argument("PDB field exceeds columns " + std::to_string(start+1)+"-"+std::to_string(start+width)+": "+value);
    line.replace(start,width,right?string(width-value.size(),' ')+value:value+string(width-value.size(),' '));
}
inline string format_atom(const Atom& a,long serial,string original="") {
    string line=original; if (line.size()<80) line.resize(80,' ');
    put(line,0,6,"ATOM"); put(line,6,5,std::to_string(serial),true); put(line,11,1,"");
    string name=a.name;
    if (!original.empty() && trim(field(original,12,16))==name) name=field(original,12,16);
    else if (name.size()<4 && !name.empty() && !std::isdigit(static_cast<unsigned char>(name[0])) && trim(a.element).size()!=2) name=" "+name;
    put(line,12,4,name); put(line,16,1,a.alt_loc); put(line,17,3,a.res_name,true);
    put(line,20,1,""); put(line,21,1,a.chain_id); put(line,22,4,std::to_string(integer(a.res_seq)),true);
    put(line,26,1,a.icode); put(line,27,3,"");
    put(line,30,8,fixed(real(a.x),3),true); put(line,38,8,fixed(real(a.y),3),true); put(line,46,8,fixed(real(a.z),3),true);
    put(line,54,6,trim(a.occ).empty()?"":fixed(real(a.occ),2),true);
    put(line,60,6,trim(a.temp).empty()?"":fixed(real(a.temp),2),true);
    put(line,76,2,a.element,true); put(line,78,2,a.charge,true); return line;
}

struct Record {
    string line; int model_index; ResiduePtr residue; Atom original;
    bool is_atom() const { return static_cast<bool>(residue); }
};
class Model {
public:
    std::vector<ResiduePtr> residues;
    std::vector<string> chains;
    // Preserved records are public for callers deliberately editing linked metadata.
    std::vector<Record> records;
    bool explicit_models=false;
    void load(std::istream& input) {
        // Parse into an independent object; a failed reload preserves this model.
        Model next; std::vector<string> lines; string line;
        while (std::getline(input,line)) {
            while (!line.empty() && line.back()=='\r') line.pop_back();
            if (trim(field(line,0,6))=="MODEL") next.explicit_models=true;
            lines.push_back(line);
        }
        if (input.bad() || (input.fail() && !input.eof())) throw std::runtime_error("Error reading PDB stream");
        int model=next.explicit_models?-1:0,segment=0;
        bool inside=!next.explicit_models; string model_id="1"; ResiduePtr previous;
        for (std::size_t i=0;i<lines.size();++i) {
            line=lines[i]; string kind=trim(field(line,0,6));
            try {
                if (kind=="MODEL") {
                    if (inside) throw std::invalid_argument("Nested MODEL");
                    if (model==std::numeric_limits<int>::max()) throw std::invalid_argument("Too many models");
                    ++model; model_id=trim(field(line,10,14)); if (model_id.empty()) model_id=std::to_string(model+1);
                    segment=0; previous.reset(); inside=true;
                } else if (kind=="ENDMDL") {
                    if (!next.explicit_models || !inside) throw std::invalid_argument("Unmatched ENDMDL");
                    inside=false; previous.reset();
                } else if (kind=="TER") {
                    if (segment==std::numeric_limits<int>::max()) throw std::invalid_argument("Too many segments");
                    ++segment; previous.reset();
                }
                if (kind!="ATOM") { next.records.push_back({line,model,ResiduePtr(),Atom()}); continue; }
                if (!inside) throw std::invalid_argument("ATOM outside MODEL");
                Atom a=Atom::parse(line);
                if (!previous || previous->model_index!=model || previous->segment_id!=segment ||
                    std::tie(previous->chain_id,previous->res_seq,previous->icode,previous->res_name)!=std::tie(a.chain_id,a.res_seq,a.icode,a.res_name)) {
                    previous=std::make_shared<Residue>(); previous->chain_id=a.chain_id; previous->res_seq=a.res_seq;
                    previous->icode=a.icode; previous->res_name=a.res_name; previous->model_index=model;
                    previous->model_id=model_id; previous->segment_id=segment; next.residues.push_back(previous);
                }
                previous->atoms.push_back(a);
                if (std::find(next.chains.begin(),next.chains.end(),a.chain_id)==next.chains.end()) next.chains.push_back(a.chain_id);
                next.records.push_back({line,model,previous,a});
            } catch (const std::invalid_argument& e) {
                throw std::invalid_argument("Invalid PDB at line "+std::to_string(i+1)+": "+e.what());
            }
        }
        if (next.explicit_models && inside) throw std::invalid_argument("MODEL block has no ENDMDL");
        *this=std::move(next);
    }
    void load_file(const string& path) {
        std::ifstream input(path.c_str(),std::ios::binary);
        if (!input) throw std::runtime_error("Cannot open input: "+path);
        load(input);
    }
    std::vector<ResiduePtr> residues_by_chain(const string& chain,int model=-1,int segment=-1) const {
        std::vector<ResiduePtr> out;
        for (const auto& r:residues) if (r->chain_id==chain && (model<0 || r->model_index==model) && (segment<0 || r->segment_id==segment)) out.push_back(r);
        return out;
    }
    std::size_t atom_count() const { std::size_t n=0; for (const auto& r:residues) n+=r->atoms.size(); return n; }
    std::vector<string> rebuild() const;
    void save(const string& path) const {
        const std::vector<string> lines=rebuild(); // Validate fully before opening output; input is already in memory.
        std::ofstream output(path.c_str(),std::ios::binary|std::ios::trunc);
        if (!output) throw std::runtime_error("Cannot open output: "+path);
        for (const string& line:lines) output << line << '\n';
        output.close();
        if (!output) throw std::runtime_error("Failed writing output: "+path);
    }
};

std::vector<string> Model::rebuild() const {
    using SerialKey=std::pair<int,long>;
    using AtomKey=std::pair<const Residue*,std::size_t>;
    std::set<const Residue*> old_res,current_res;
    std::set<int> known_models;
    std::map<int,std::set<long> > reserved,seen;
    std::map<SerialKey,Atom> original,current;
    std::map<const Residue*,std::size_t> last_slot,positions;
    std::size_t original_count=0;
    for (std::size_t i=0;i<records.size();++i) {
        const Record& rec=records[i]; string kind=trim(field(rec.line,0,6));
        if (rec.is_atom()) {
            old_res.insert(rec.residue.get()); ++original_count; last_slot[rec.residue.get()]=i;
            long serial=integer(rec.original.serial); SerialKey key(rec.model_index,serial);
            if (!original.emplace(key,rec.original).second) throw std::invalid_argument("Duplicate original atom serial within model");
            reserved[rec.model_index].insert(serial);
        } else {
            if (kind=="MODEL") known_models.insert(rec.model_index);
            if ((kind=="HETATM" || kind=="TER") && !trim(field(rec.line,6,11)).empty()) {
                long serial=integer(field(rec.line,6,11)); reserved[rec.model_index].insert(serial); seen[rec.model_index].insert(serial);
            }
        }
    }
    if (!explicit_models) known_models.insert(0);
    auto automatic=[](const string& s) { return trim(s).empty() || upper(trim(s))=="AUTO"; };
    for (const auto& r:residues) {
        if (!r || !current_res.insert(r.get()).second) throw std::invalid_argument("Null/duplicate residue object");
        if (!known_models.count(r->model_index)) throw std::invalid_argument("Residue refers to unknown model");
        positions[r.get()]=0;
        for (const Atom& a:r->atoms) if (!automatic(a.serial)) {
            long serial=integer(a.serial);
            if (serial<1 || serial>99999) throw std::invalid_argument("Atom serial must be between 1 and 99999");
            if (!seen[r->model_index].insert(serial).second) throw std::invalid_argument("Duplicate serial in model");
            reserved[r->model_index].insert(serial);
        }
    }
    std::map<AtomKey,long> serials;
    for (const auto& r:residues) for (std::size_t i=0;i<r->atoms.size();++i) {
        const Atom& a=r->atoms[i]; long serial;
        if (automatic(a.serial)) {
            auto& used=reserved[r->model_index];
            long maximum=used.empty()?0:*used.rbegin();
            if (maximum>=99999) throw std::invalid_argument("No serial space remaining in PDB model");
            serial=maximum+1; used.insert(serial);
        } else serial=integer(a.serial);
        serials[{r.get(),i}]=serial; current[{r->model_index,serial}]=a;
    }
    for (const Record& rec:records) if (!rec.is_atom()) {
        string kind=trim(field(rec.line,0,6));
        if (kind!="CONECT" && kind!="ANISOU" && kind!="SIGATM" && kind!="SIGUIJ") continue;
        std::size_t end=kind=="CONECT"?rec.line.size():11;
        for (std::size_t pos=6;pos<end;pos+=5) {
            string value=trim(field(rec.line,pos,pos+5)); if (value.empty()) continue;
            long serial=integer(value);
            for (const auto& old:original) if (old.first.second==serial && (kind=="CONECT" || old.first.first==rec.model_index)) {
                auto found=current.find(old.first);
                if (found==current.end() || !old.second.same_identity(found->second))
                    throw std::invalid_argument("Edit invalidates retained "+kind+" record for atom "+value);
            }
        }
    }
    std::vector<ResiduePtr> additions;
    for (const auto& r:residues) if (!old_res.count(r.get())) additions.push_back(r);
    bool changed=original_count!=atom_count() || old_res!=current_res;
    std::set<const Residue*> emitted_new;
    std::vector<string> output;
    auto emit_atom=[&](const ResiduePtr& r,std::size_t index,const string& original_line) {
        output.push_back(format_atom(r->atoms[index],serials.at({r.get(),index}),original_line));
    };
    auto emit_new=[&](int model) {
        std::vector<ResiduePtr> pending;
        for (const auto& r:additions) if (r->model_index==model && !emitted_new.count(r.get())) pending.push_back(r);
        if (pending.empty()) return;
        if (std::any_of(output.begin(),output.end(),[](const string& line){return trim(field(line,0,6))=="ATOM";})) output.push_back("TER");
        ResiduePtr previous;
        for (const auto& r:pending) {
            if (previous && (previous->chain_id!=r->chain_id || previous->segment_id!=r->segment_id)) output.push_back("TER");
            for (std::size_t i=0;i<r->atoms.size();++i) emit_atom(r,i,"");
            previous=r; emitted_new.insert(r.get());
        }
    };
    for (std::size_t i=0;i<records.size();++i) {
        const Record& rec=records[i];
        if (!rec.is_atom()) {
            string line=rec.line,kind=trim(field(line,0,6));
            if (kind=="ENDMDL" || (!explicit_models && (kind=="CONECT" || kind=="MASTER" || kind=="END"))) emit_new(explicit_models?rec.model_index:0);
            if (kind=="MASTER" && changed) continue;
            if (kind=="TER" && line.size()>=27) for (auto it=output.rbegin();it!=output.rend();++it) {
                string k=trim(field(*it,0,6));
                if (k=="ATOM" || k=="MODEL" || k=="ENDMDL" || k=="TER") {
                    if (k=="ATOM") line.replace(17,10,it->substr(17,10));
                    break;
                }
            }
            output.push_back(line); continue;
        }
        if (!current_res.count(rec.residue.get())) continue;
        std::size_t& position=positions[rec.residue.get()];
        if (position<rec.residue->atoms.size()) { emit_atom(rec.residue,position,rec.line); ++position; }
        if (i==last_slot[rec.residue.get()]) while (position<rec.residue->atoms.size()) { emit_atom(rec.residue,position,""); ++position; }
    }
    if (!explicit_models) emit_new(0);
    if (emitted_new.size()!=additions.size()) throw std::invalid_argument("Could not place new residue in model");
    return output;
}

// Familiar aliases for users of the Python version.
using TAtom=Atom;
using TResidue=Residue;
using TBond=Bond;
using TPDBModel=Model;

struct Options {
    string input="-",output,chain; long model=0;
    bool input_provided=false;
    bool chain_set=false,no_residues=false,no_bonds=false,no_dihedrals=false,quiet=false,help=false;
};
struct ArgumentError:std::runtime_error { explicit ArgumentError(const string& message):std::runtime_error(message) {} };
inline void help(std::ostream& out) {
    out << "Usage: artemis [options] [PDB]\n\n"
        "Artemis: PDB residues, intra-residue bonds and backbone dihedral angles.\n"
        "PDB omitted reads piped/redirected input; a bare terminal invocation shows help.\n"
        "Use '-' explicitly to read stdin interactively.\n\n"
        "  --no-residues       Suppress the atom/residue listing; bonds remain independent\n"
        "  --no-bonds          Skip bond calculations and reports\n"
        "  --no-dihedrals      Skip dihedral calculations and reports\n"
        "  --chain ID         Report this case-sensitive chain ('' for blank)\n"
        "  --model N          Report this positive PDB MODEL identifier\n"
        "  -o, --output PDB   Save the full rebuilt structure ('-' for stdout)\n"
        "  -q, --quiet        Suppress all reports and their calculations\n"
        "  -h, --help         Show this help\n"
        "  --                 End options (for filenames beginning with '-')\n\n"
        "Examples:\n"
        "  artemis structure.pdb\n"
        "  artemis structure.pdb --chain A --no-residues --no-bonds\n"
        "  artemis structure.pdb --model 2 --no-dihedrals\n"
        "  cat structure.pdb | artemis - --no-bonds\n"
        "  artemis structure.pdb -q -o rebuilt.pdb\n\n"
        "Chain/model filters affect reports only; saved PDBs retain all models and\n"
        "conformers. Without MODEL records the model identifier is 1. HETATM records\n"
        "are preserved but excluded from analysis. N/A means an undefined angle.\n"
        "With -o -, PDB goes to stdout and reports go to stderr.\n"
        "Exit codes: 0 success, 1 processing/I/O error, 2 argument error.\n";
}
inline Options parse_options(int argc,char** argv) {
    Options o; bool positional=false,end=false;
    for (int i=1;i<argc;++i) {
        string arg=argv[i];
        if (!end && arg=="--") { end=true; continue; }
        if (!end && (arg=="-h" || arg=="--help")) { o.help=true; return o; }
        if (!end && arg=="--no-residues") { o.no_residues=true; continue; }
        if (!end && arg=="--no-bonds") { o.no_bonds=true; continue; }
        if (!end && arg=="--no-dihedrals") { o.no_dihedrals=true; continue; }
        if (!end && (arg=="-q" || arg=="--quiet")) { o.quiet=true; continue; }
        if (!end && arg.size()>1 && arg[0]=='-') {
            string key=arg,value; std::size_t equal=arg.find('='); bool attached=equal!=string::npos;
            if (attached) { key=arg.substr(0,equal); value=arg.substr(equal+1); }
            if (key!="--chain" && key!="--model" && key!="--output" && key!="-o") throw ArgumentError("Unknown option: "+arg);
            if (!attached) {
                if (i+1>=argc || (string(argv[i+1]).size()>1 && argv[i+1][0]=='-')) throw ArgumentError("Missing value for "+key);
                value=argv[++i];
            }
            if (key=="--chain") {
                if (value.size()>1 || (!value.empty() && std::isspace(static_cast<unsigned char>(value[0])))) throw ArgumentError("--chain requires one character or ''");
                o.chain=value; o.chain_set=true;
            } else if (key=="--model") {
                try { o.model=integer(value); } catch (const std::exception&) { throw ArgumentError("--model requires a positive integer"); }
                if (o.model<1) throw ArgumentError("--model requires a positive integer");
            } else { if (value.empty()) throw ArgumentError("Empty output filename"); o.output=value; }
            continue;
        }
        if (positional) throw ArgumentError("Only one input PDB can be supplied");
        o.input=arg; o.input_provided=true; positional=true;
    }
    return o;
}
inline bool matches(const Residue& r,const Options& o) {
    if (o.chain_set && r.chain_id!=o.chain) return false;
    if (o.model) {
        try { return integer(r.model_id)==o.model; } catch (const std::invalid_argument&) { return false; }
    }
    return true;
}
inline string angle_text(double value) { return std::isfinite(value)?fixed(value,3):"N/A"; }
inline void report(Model& model,const Options& o,std::ostream& out) {
    if (o.quiet || (o.no_residues && o.no_bonds && o.no_dihedrals)) return;
    out << "Artemis - C++ molecular library for Brookhaven PDB files\n"
        "Copyright (c) 2014 Jamie J. Al-Nasir, All Rights Reserved\n"
        "Loaded: " << (o.input=="-"?"<stdin>":o.input) << "\n\n";
    for (std::size_t begin=0;begin<model.residues.size();) {
        std::size_t end=begin+1; const Residue& first=*model.residues[begin];
        while (end<model.residues.size()) {
            const Residue& r=*model.residues[end];
            if (std::tie(r.model_index,r.chain_id,r.segment_id)!=std::tie(first.model_index,first.chain_id,first.segment_id)) break;
            ++end;
        }
        if (matches(first,o)) {
            out << "Model: " << first.model_id << "; Chain: " << (first.chain_id.empty()?"(blank)":first.chain_id) << "; Segment: " << first.segment_id+1 << '\n';
            for (std::size_t i=begin;i<end && (!o.no_residues || !o.no_bonds);++i) {
                Residue& r=*model.residues[i]; out << "\nResidue: " << r.res_name << r.res_seq << r.icode << '\n';
                if (!o.no_residues) {
                    out << "chain\tname\tresidue\tx\ty\tz\taltLoc\n";
                    for (const Atom& a:r.atoms) out << (a.chain_id.empty()?"(blank)":a.chain_id) << '\t' << a.name << '\t' << a.res_name << a.res_seq << a.icode << '\t' << a.x << '\t' << a.y << '\t' << a.z << '\t' << a.alt_loc << '\n';
                }
                if (!o.no_bonds) {
                    r.compute_bonds(); out << '\n' << r.bonds.size() << " Covalent bond(s) computed within this residue:\n";
                    for (const Bond& b:r.bonds) out << r.atoms[b.atom1].name << " - " << r.atoms[b.atom2].name << '\n';
                }
            }
            if (!o.no_dihedrals) {
                out << "\nDihedral/Torsional angles (N/A = undefined, missing atoms or chain break):\n";
                std::vector<ResiduePtr> group(model.residues.begin()+begin,model.residues.begin()+end);
                for (std::size_t i=0;i<group.size();++i) if (group[i]->type()==ResidueType::AA) {
                    Torsions t=residue_dihedrals(group,i); const Residue& r=*group[i];
                    out << "Residue " << r.res_name << r.res_seq << r.icode << ": Phi=" << angle_text(t.phi) << ", Psi=" << angle_text(t.psi) << ", Ohmega=" << angle_text(t.omega) << '\n';
                }
            }
            out << '\n';
        }
        begin=end;
    }
}
inline bool stdin_is_terminal() {
#if defined(_WIN32)
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(fileno(stdin)) != 0;
#endif
}
inline int run(int argc,char** argv) {
    try {
        Options options=parse_options(argc,argv);
        if (options.help) { help(std::cout); return 0; }
        if (!options.input_provided && stdin_is_terminal()) {
            if (argc==1) { help(std::cout); return 0; }
            throw ArgumentError("Provide a PDB filename, pipe PDB data, or use '-' explicitly for stdin");
        }
        Model model;
        if (options.input=="-") model.load(std::cin); else model.load_file(options.input);
        if (model.atom_count()==0) throw std::invalid_argument("No ATOM records found in input (HETATM-only files are not analysed)");
        if (!std::any_of(model.residues.begin(),model.residues.end(),[&](const ResiduePtr& r){return matches(*r,options);})) throw ArgumentError("No ATOM residues match the requested chain/model");
        if (options.output=="-") { for (const string& line:model.rebuild()) std::cout << line << '\n'; }
        else if (!options.output.empty()) model.save(options.output);
        report(model,options,options.output=="-"?std::cerr:std::cout);
        std::cout.flush(); std::cerr.flush();
        if (!std::cout || !std::cerr) return 1;
        return 0;
    } catch (const ArgumentError& e) { std::cerr << "Artemis: " << e.what() << "\nUse --help for usage.\n"; return 2; }
      catch (const std::exception& e) { std::cerr << "Artemis: " << e.what() << '\n'; return 1; }
}
} // namespace artemis

#ifndef ARTEMIS_NO_MAIN
int main(int argc,char** argv) { return artemis::run(argc,argv); }
#endif
