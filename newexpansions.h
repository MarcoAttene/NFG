#pragma once

#include "numerics.h"

#pragma intrinsic(fabs)

struct expansion {
	// This is the expansion. If N>4 the last slot stores the actual length
	double *_e;	
	uint64_t _len;

	expansion(uint64_t reserved_size) :
		_e(AllocDoubles((uint32_t)reserved_size)) { }

	~expansion() { FreeDoubles(_e); }

	expansion(double x) :
		_e(AllocDoubles(1)),
		_len(1) {
		*_e = x;
	}

	uint64_t len() const { return _len; }

	void setlen(uint64_t l) { _len = l; }

	double get_d() const { return _e[len() - 1]; }

	operator double() { return *_e; }
	operator const double() const { return *_e; }
	operator double* () { return _e; }
	operator const double* () const { return _e; }
};

inline std::ostream& operator<<(std::ostream& os, const  expansion& e)
{
	for (uint64_t i = 0; i < e.len(); i++) os << e._e[i] << " ";
	return os;
}

void invert(expansion& e) { expansionObject::Gen_Invert((int)e.len(), e); }

expansion operator+(const expansion& a, const expansion& b) {
	const uint64_t N = a.len(), M = b.len();
	expansion result(N+M);

	result.setlen(expansionObject::Gen_Sum((int)N, a, (int)M, b, result));

	return result;
}

expansion operator-(const expansion& a, const expansion& b) {
	const uint64_t N = a.len(), M = b.len();
	expansion result(N + M);

	result.setlen(expansionObject::Gen_Diff((int)N, a, (int)M, b, result));

	return result;
}

expansion operator*(const expansion& a, const expansion& b) {
	const uint64_t N = a.len(), M = b.len();
	expansion result(N*M*2);
	result.setlen(expansionObject::Gen_Product((int)N, a, (int)M, b, result));

	return result;
}

expansion sqr(const expansion& a) {
	const uint64_t N = a.len();
	expansion result(N*4 - 2);
	if (N == 1) expansionObject::Square(a, result);
	else if (N == 2) expansionObject::Two_Square(a, result);
	else ip_error("sqr() for long expansions: unimplemented.\n");

	return result;
}

int sgn(const expansion& e) { 
	const double d = e.get_d();
	return ((d > 0) - (d < 0));
}


class Point3c {
public:
	double x, y, z;

	Point3c() {
		x = rand() / ((double)RAND_MAX);
		y = rand() / ((double)RAND_MAX);
		z = rand() / ((double)RAND_MAX);
	}
};

inline int inSphere_exact(double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
	double aex[2];
	expansionObject::two_Diff(pax, pex, aex);
	double aey[2];
	expansionObject::two_Diff(pay, pey, aey);
	double aez[2];
	expansionObject::two_Diff(paz, pez, aez);
	double bex[2];
	expansionObject::two_Diff(pbx, pex, bex);
	double bey[2];
	expansionObject::two_Diff(pby, pey, bey);
	double bez[2];
	expansionObject::two_Diff(pbz, pez, bez);
	double cex[2];
	expansionObject::two_Diff(pcx, pex, cex);
	double cey[2];
	expansionObject::two_Diff(pcy, pey, cey);
	double cez[2];
	expansionObject::two_Diff(pcz, pez, cez);
	double dex[2];
	expansionObject::two_Diff(pdx, pex, dex);
	double dey[2];
	expansionObject::two_Diff(pdy, pey, dey);
	double dez[2];
	expansionObject::two_Diff(pdz, pez, dez);
	double aexbey[8];
	int aexbey_len = expansionObject::Gen_Product(2, aex, 2, bey, aexbey);
	double bexaey[8];
	int bexaey_len = expansionObject::Gen_Product(2, bex, 2, aey, bexaey);
	double ab[16];
	int ab_len = expansionObject::Gen_Diff(aexbey_len, aexbey, bexaey_len, bexaey, ab);
	double bexcey[8];
	int bexcey_len = expansionObject::Gen_Product(2, bex, 2, cey, bexcey);
	double cexbey[8];
	int cexbey_len = expansionObject::Gen_Product(2, cex, 2, bey, cexbey);
	double bc[16];
	int bc_len = expansionObject::Gen_Diff(bexcey_len, bexcey, cexbey_len, cexbey, bc);
	double cexdey[8];
	int cexdey_len = expansionObject::Gen_Product(2, cex, 2, dey, cexdey);
	double dexcey[8];
	int dexcey_len = expansionObject::Gen_Product(2, dex, 2, cey, dexcey);
	double cd[16];
	int cd_len = expansionObject::Gen_Diff(cexdey_len, cexdey, dexcey_len, dexcey, cd);
	double dexaey[8];
	int dexaey_len = expansionObject::Gen_Product(2, dex, 2, aey, dexaey);
	double aexdey[8];
	int aexdey_len = expansionObject::Gen_Product(2, aex, 2, dey, aexdey);
	double da[16];
	int da_len = expansionObject::Gen_Diff(dexaey_len, dexaey, aexdey_len, aexdey, da);
	double aexcey[8];
	int aexcey_len = expansionObject::Gen_Product(2, aex, 2, cey, aexcey);
	double cexaey[8];
	int cexaey_len = expansionObject::Gen_Product(2, cex, 2, aey, cexaey);
	double ac[16];
	int ac_len = expansionObject::Gen_Diff(aexcey_len, aexcey, cexaey_len, cexaey, ac);
	double bexdey[8];
	int bexdey_len = expansionObject::Gen_Product(2, bex, 2, dey, bexdey);
	double dexbey[8];
	int dexbey_len = expansionObject::Gen_Product(2, dex, 2, bey, dexbey);
	double bd[16];
	int bd_len = expansionObject::Gen_Diff(bexdey_len, bexdey, dexbey_len, dexbey, bd);
	double abc1_p[32], * abc1 = abc1_p;
	int abc1_len = expansionObject::Gen_Product_With_PreAlloc(2, aez, bc_len, bc, &abc1, 32);
	double abc2_p[32], * abc2 = abc2_p;
	int abc2_len = expansionObject::Gen_Product_With_PreAlloc(2, bez, ac_len, ac, &abc2, 32);
	double abc3_p[32], * abc3 = abc3_p;
	int abc3_len = expansionObject::Gen_Product_With_PreAlloc(2, cez, ab_len, ab, &abc3, 32);
	double abc4_p[32], * abc4 = abc4_p;
	int abc4_len = expansionObject::Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 32);
	double abc_p[32], * abc = abc_p;
	int abc_len = expansionObject::Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 32);
	double bcd1_p[32], * bcd1 = bcd1_p;
	int bcd1_len = expansionObject::Gen_Product_With_PreAlloc(2, bez, cd_len, cd, &bcd1, 32);
	double bcd2_p[32], * bcd2 = bcd2_p;
	int bcd2_len = expansionObject::Gen_Product_With_PreAlloc(2, cez, bd_len, bd, &bcd2, 32);
	double bcd3_p[32], * bcd3 = bcd3_p;
	int bcd3_len = expansionObject::Gen_Product_With_PreAlloc(2, dez, bc_len, bc, &bcd3, 32);
	double bcd4_p[32], * bcd4 = bcd4_p;
	int bcd4_len = expansionObject::Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 32);
	double bcd_p[32], * bcd = bcd_p;
	int bcd_len = expansionObject::Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 32);
	double cda1_p[32], * cda1 = cda1_p;
	int cda1_len = expansionObject::Gen_Product_With_PreAlloc(2, cez, da_len, da, &cda1, 32);
	double cda2_p[32], * cda2 = cda2_p;
	int cda2_len = expansionObject::Gen_Product_With_PreAlloc(2, dez, ac_len, ac, &cda2, 32);
	double cda3_p[32], * cda3 = cda3_p;
	int cda3_len = expansionObject::Gen_Product_With_PreAlloc(2, aez, cd_len, cd, &cda3, 32);
	double cda4_p[32], * cda4 = cda4_p;
	int cda4_len = expansionObject::Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 32);
	double cda_p[32], * cda = cda_p;
	int cda_len = expansionObject::Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 32);
	double dab1_p[32], * dab1 = dab1_p;
	int dab1_len = expansionObject::Gen_Product_With_PreAlloc(2, dez, ab_len, ab, &dab1, 32);
	double dab2_p[32], * dab2 = dab2_p;
	int dab2_len = expansionObject::Gen_Product_With_PreAlloc(2, aez, bd_len, bd, &dab2, 32);
	double dab3_p[32], * dab3 = dab3_p;
	int dab3_len = expansionObject::Gen_Product_With_PreAlloc(2, bez, da_len, da, &dab3, 32);
	double dab4_p[32], * dab4 = dab4_p;
	int dab4_len = expansionObject::Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 32);
	double dab_p[32], * dab = dab_p;
	int dab_len = expansionObject::Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 32);
	double al1[8];
	int al1_len = expansionObject::Gen_Product(2, aex, 2, aex, al1);
	double al2[8];
	int al2_len = expansionObject::Gen_Product(2, aey, 2, aey, al2);
	double al3[8];
	int al3_len = expansionObject::Gen_Product(2, aez, 2, aez, al3);
	double al4[16];
	int al4_len = expansionObject::Gen_Sum(al1_len, al1, al2_len, al2, al4);
	double alift[24];
	int alift_len = expansionObject::Gen_Sum(al4_len, al4, al3_len, al3, alift);
	double bl1[8];
	int bl1_len = expansionObject::Gen_Product(2, bex, 2, bex, bl1);
	double bl2[8];
	int bl2_len = expansionObject::Gen_Product(2, bey, 2, bey, bl2);
	double bl3[8];
	int bl3_len = expansionObject::Gen_Product(2, bez, 2, bez, bl3);
	double bl4[16];
	int bl4_len = expansionObject::Gen_Sum(bl1_len, bl1, bl2_len, bl2, bl4);
	double blift[24];
	int blift_len = expansionObject::Gen_Sum(bl4_len, bl4, bl3_len, bl3, blift);
	double cl1[8];
	int cl1_len = expansionObject::Gen_Product(2, cex, 2, cex, cl1);
	double cl2[8];
	int cl2_len = expansionObject::Gen_Product(2, cey, 2, cey, cl2);
	double cl3[8];
	int cl3_len = expansionObject::Gen_Product(2, cez, 2, cez, cl3);
	double cl4[16];
	int cl4_len = expansionObject::Gen_Sum(cl1_len, cl1, cl2_len, cl2, cl4);
	double clift[24];
	int clift_len = expansionObject::Gen_Sum(cl4_len, cl4, cl3_len, cl3, clift);
	double dl1[8];
	int dl1_len = expansionObject::Gen_Product(2, dex, 2, dex, dl1);
	double dl2[8];
	int dl2_len = expansionObject::Gen_Product(2, dey, 2, dey, dl2);
	double dl3[8];
	int dl3_len = expansionObject::Gen_Product(2, dez, 2, dez, dl3);
	double dl4[16];
	int dl4_len = expansionObject::Gen_Sum(dl1_len, dl1, dl2_len, dl2, dl4);
	double dlift[24];
	int dlift_len = expansionObject::Gen_Sum(dl4_len, dl4, dl3_len, dl3, dlift);
	double ds1_p[32], * ds1 = ds1_p;
	int ds1_len = expansionObject::Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 32);
	double ds2_p[32], * ds2 = ds2_p;
	int ds2_len = expansionObject::Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 32);
	double dl_p[32], * dl = dl_p;
	int dl_len = expansionObject::Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1_len, ds1, &dl, 32);
	double dr1_p[32], * dr1 = dr1_p;
	int dr1_len = expansionObject::Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 32);
	double dr2_p[32], * dr2 = dr2_p;
	int dr2_len = expansionObject::Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 32);
	double dr_p[32], * dr = dr_p;
	int dr_len = expansionObject::Gen_Diff_With_PreAlloc(dr2_len, dr2, dr1_len, dr1, &dr, 32);
	double det_p[32], * det = det_p;
	int det_len = expansionObject::Gen_Sum_With_PreAlloc(dl_len, dl, dr_len, dr, &det, 32);

	double return_value = det[det_len - 1];
	if (det_p != det) FreeDoubles(det);
	if (dr_p != dr) FreeDoubles(dr);
	if (dr2_p != dr2) FreeDoubles(dr2);
	if (dr1_p != dr1) FreeDoubles(dr1);
	if (dl_p != dl) FreeDoubles(dl);
	if (ds2_p != ds2) FreeDoubles(ds2);
	if (ds1_p != ds1) FreeDoubles(ds1);
	if (dab_p != dab) FreeDoubles(dab);
	if (dab4_p != dab4) FreeDoubles(dab4);
	if (dab3_p != dab3) FreeDoubles(dab3);
	if (dab2_p != dab2) FreeDoubles(dab2);
	if (dab1_p != dab1) FreeDoubles(dab1);
	if (cda_p != cda) FreeDoubles(cda);
	if (cda4_p != cda4) FreeDoubles(cda4);
	if (cda3_p != cda3) FreeDoubles(cda3);
	if (cda2_p != cda2) FreeDoubles(cda2);
	if (cda1_p != cda1) FreeDoubles(cda1);
	if (bcd_p != bcd) FreeDoubles(bcd);
	if (bcd4_p != bcd4) FreeDoubles(bcd4);
	if (bcd3_p != bcd3) FreeDoubles(bcd3);
	if (bcd2_p != bcd2) FreeDoubles(bcd2);
	if (bcd1_p != bcd1) FreeDoubles(bcd1);
	if (abc_p != abc) FreeDoubles(abc);
	if (abc4_p != abc4) FreeDoubles(abc4);
	if (abc3_p != abc3) FreeDoubles(abc3);
	if (abc2_p != abc2) FreeDoubles(abc2);
	if (abc1_p != abc1) FreeDoubles(abc1);

	if (return_value > 0) return 1;
	if (return_value < 0) return -1;
	return 0;
}

inline int inSphere_bigfloat(bigfloat pax, bigfloat pay, bigfloat paz, bigfloat pbx, bigfloat pby, bigfloat pbz, bigfloat pcx, bigfloat pcy, bigfloat pcz, bigfloat pdx, bigfloat pdy, bigfloat pdz, bigfloat pex, bigfloat pey, bigfloat pez)
{
	const bigfloat aex = pax - pex;
	const bigfloat aey = pay - pey;
	const bigfloat aez = paz - pez;
	const bigfloat bex = pbx - pex;
	const bigfloat bey = pby - pey;
	const bigfloat bez = pbz - pez;
	const bigfloat cex = pcx - pex;
	const bigfloat cey = pcy - pey;
	const bigfloat cez = pcz - pez;
	const bigfloat dex = pdx - pex;
	const bigfloat dey = pdy - pey;
	const bigfloat dez = pdz - pez;
	const bigfloat aexbey = aex * bey;
	const bigfloat bexaey = bex * aey;
	const bigfloat ab = aexbey - bexaey;
	const bigfloat bexcey = bex * cey;
	const bigfloat cexbey = cex * bey;
	const bigfloat bc = bexcey - cexbey;
	const bigfloat cexdey = cex * dey;
	const bigfloat dexcey = dex * cey;
	const bigfloat cd = cexdey - dexcey;
	const bigfloat dexaey = dex * aey;
	const bigfloat aexdey = aex * dey;
	const bigfloat da = dexaey - aexdey;
	const bigfloat aexcey = aex * cey;
	const bigfloat cexaey = cex * aey;
	const bigfloat ac = aexcey - cexaey;
	const bigfloat bexdey = bex * dey;
	const bigfloat dexbey = dex * bey;
	const bigfloat bd = bexdey - dexbey;
	const bigfloat abc1 = aez * bc;
	const bigfloat abc2 = bez * ac;
	const bigfloat abc3 = cez * ab;
	const bigfloat abc4 = abc1 + abc3;
	const bigfloat abc = abc4 - abc2;
	const bigfloat bcd1 = bez * cd;
	const bigfloat bcd2 = cez * bd;
	const bigfloat bcd3 = dez * bc;
	const bigfloat bcd4 = bcd1 + bcd3;
	const bigfloat bcd = bcd4 - bcd2;
	const bigfloat cda1 = cez * da;
	const bigfloat cda2 = dez * ac;
	const bigfloat cda3 = aez * cd;
	const bigfloat cda4 = cda1 + cda3;
	const bigfloat cda = cda4 + cda2;
	const bigfloat dab1 = dez * ab;
	const bigfloat dab2 = aez * bd;
	const bigfloat dab3 = bez * da;
	const bigfloat dab4 = dab1 + dab3;
	const bigfloat dab = dab4 + dab2;
	const bigfloat al1 = aex * aex;
	const bigfloat al2 = aey * aey;
	const bigfloat al3 = aez * aez;
	const bigfloat al4 = al1 + al2;
	const bigfloat alift = al4 + al3;
	const bigfloat bl1 = bex * bex;
	const bigfloat bl2 = bey * bey;
	const bigfloat bl3 = bez * bez;
	const bigfloat bl4 = bl1 + bl2;
	const bigfloat blift = bl4 + bl3;
	const bigfloat cl1 = cex * cex;
	const bigfloat cl2 = cey * cey;
	const bigfloat cl3 = cez * cez;
	const bigfloat cl4 = cl1 + cl2;
	const bigfloat clift = cl4 + cl3;
	const bigfloat dl1 = dex * dex;
	const bigfloat dl2 = dey * dey;
	const bigfloat dl3 = dez * dez;
	const bigfloat dl4 = dl1 + dl2;
	const bigfloat dlift = dl4 + dl3;
	const bigfloat ds1 = dlift * abc;
	const bigfloat ds2 = clift * dab;
	const bigfloat dl = ds2 - ds1;
	const bigfloat dr1 = blift * cda;
	const bigfloat dr2 = alift * bcd;
	const bigfloat dr = dr2 - dr1;
	const bigfloat det = dl + dr;

	return det.sgn();
}

inline int inSphere_expansion(expansion pax, expansion pay, expansion paz, expansion pbx, expansion pby, expansion pbz, expansion pcx, expansion pcy, expansion pcz, expansion pdx, expansion pdy, expansion pdz, expansion pex, expansion pey, expansion pez)
{
	const auto aex = pax - pex;
	const auto aey = pay - pey;
	const auto aez = paz - pez;
	const auto bex = pbx - pex;
	const auto bey = pby - pey;
	const auto bez = pbz - pez;
	const auto cex = pcx - pex;
	const auto cey = pcy - pey;
	const auto cez = pcz - pez;
	const auto dex = pdx - pex;
	const auto dey = pdy - pey;
	const auto dez = pdz - pez;
	const auto aexbey = aex * bey;
	const auto bexaey = bex * aey;
	const auto ab = aexbey - bexaey;
	const auto bexcey = bex * cey;
	const auto cexbey = cex * bey;
	const auto bc = bexcey - cexbey;
	const auto cexdey = cex * dey;
	const auto dexcey = dex * cey;
	const auto cd = cexdey - dexcey;
	const auto dexaey = dex * aey;
	const auto aexdey = aex * dey;
	const auto da = dexaey - aexdey;
	const auto aexcey = aex * cey;
	const auto cexaey = cex * aey;
	const auto ac = aexcey - cexaey;
	const auto bexdey = bex * dey;
	const auto dexbey = dex * bey;
	const auto bd = bexdey - dexbey;
	const auto abc1 = aez * bc;
	const auto abc2 = bez * ac;
	const auto abc3 = cez * ab;
	const auto abc4 = abc1 + abc3;
	const auto abc = abc4 - abc2;
	const auto bcd1 = bez * cd;
	const auto bcd2 = cez * bd;
	const auto bcd3 = dez * bc;
	const auto bcd4 = bcd1 + bcd3;
	const auto bcd = bcd4 - bcd2;
	const auto cda1 = cez * da;
	const auto cda2 = dez * ac;
	const auto cda3 = aez * cd;
	const auto cda4 = cda1 + cda3;
	const auto cda = cda4 + cda2;
	const auto dab1 = dez * ab;
	const auto dab2 = aez * bd;
	const auto dab3 = bez * da;
	const auto dab4 = dab1 + dab3;
	const auto dab = dab4 + dab2;
	const auto al1 = aex * aex;
	const auto al2 = aey * aey;
	const auto al3 = aez * aez;
	const auto al4 = al1 + al2;
	const auto alift = al4 + al3;
	const auto bl1 = bex * bex;
	const auto bl2 = bey * bey;
	const auto bl3 = bez * bez;
	const auto bl4 = bl1 + bl2;
	const auto blift = bl4 + bl3;
	const auto cl1 = cex * cex;
	const auto cl2 = cey * cey;
	const auto cl3 = cez * cez;
	const auto cl4 = cl1 + cl2;
	const auto clift = cl4 + cl3;
	const auto dl1 = dex * dex;
	const auto dl2 = dey * dey;
	const auto dl3 = dez * dez;
	const auto dl4 = dl1 + dl2;
	const auto dlift = dl4 + dl3;
	const auto ds1 = dlift * abc;
	const auto ds2 = clift * dab;
	const auto dl = ds2 - ds1;
	const auto dr1 = blift * cda;
	const auto dr2 = alift * bcd;
	const auto dr = dr2 - dr1;
	const auto det = dl + dr;

	return sgn(det);
}

int inSphere_bigfloat(const Point3c& a, const Point3c& b, const Point3c& c, const Point3c& d, const Point3c& e) {
	return inSphere_bigfloat(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z, e.x, e.y, e.z);
}

int inSphere_expansion(const Point3c& a, const Point3c& b, const Point3c& c, const Point3c& d, const Point3c& e) {
	return inSphere_expansion(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z, e.x, e.y, e.z);
}

int inSphere_exact(const Point3c& a, const Point3c& b, const Point3c& c, const Point3c& d, const Point3c& e) {
	return inSphere_exact(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z, e.x, e.y, e.z);
}



void checkAndProfilePredicates() {
	std::chrono::steady_clock::time_point begin, end;
	long long timed;

	const int vsize = 100000; // Number of predicate calls
	Point3c* arr = new Point3c[vsize];
	int acc = 0; // Dummy accumulator to avoid that the optimizer kills everything

	//// Check
	//for (int i = vsize - 5; --i >= 0; )
	//	if (inSphere_expansion(arr[i], arr[i+1], arr[i+2], arr[i+3], arr[i + 4]) !=
	//		inSphere_bigfloat(arr[i], arr[i + 1], arr[i + 2], arr[i + 3], arr[i + 4])) ip_error("Mismatching Orient3D\n");

	// Profile
	begin = std::chrono::steady_clock::now();
	for (int i = vsize - 5; --i >= 0; ) acc += inSphere_bigfloat(arr[i], arr[i + 1], arr[i + 2], arr[i + 3], arr[i + 4]);
	end = std::chrono::steady_clock::now();
	timed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("%d\rIS>>> bigfloat %lld\n", acc, timed);

	begin = std::chrono::steady_clock::now();
	for (int i = vsize - 5; --i >= 0; ) acc += inSphere_expansion(arr[i], arr[i + 1], arr[i + 2], arr[i + 3], arr[i + 4]);
	end = std::chrono::steady_clock::now();
	timed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("%d\rIS>>> expansion %lld\n", acc, timed);

	begin = std::chrono::steady_clock::now();
	for (int i = vsize - 5; --i >= 0; ) acc += inSphere_exact(arr[i], arr[i + 1], arr[i + 2], arr[i + 3], arr[i + 4]);
	end = std::chrono::steady_clock::now();
	timed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("%d\rIS>>> exact %lld\n", acc, timed);
}
