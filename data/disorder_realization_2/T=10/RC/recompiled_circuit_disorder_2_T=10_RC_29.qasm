OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(5.7603523) q[0];
sx q[0];
rz(8.8011959) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5335124) q[0];
sx q[0];
rz(-2.170616) q[0];
sx q[0];
rz(-1.5754726) q[0];
rz(-pi) q[1];
rz(0.093437336) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(2.0051533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3574672) q[1];
sx q[1];
rz(-2.3298414) q[1];
sx q[1];
rz(-1.1151421) q[1];
rz(-0.15668232) q[3];
sx q[3];
rz(-2.0041668) q[3];
sx q[3];
rz(2.3436848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-0.99213123) q[3];
sx q[3];
rz(2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37110776) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(2.4761377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21214813) q[0];
sx q[0];
rz(-2.0245027) q[0];
sx q[0];
rz(1.2544022) q[0];
rz(-2.6639054) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(1.8600841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7484819) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(1.837681) q[1];
rz(-pi) q[2];
rz(1.6104224) q[3];
sx q[3];
rz(-0.70981662) q[3];
sx q[3];
rz(-2.1982847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-0.75876045) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6436359) q[0];
sx q[0];
rz(-0.47753497) q[0];
sx q[0];
rz(-0.60995539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7064352) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(-0.26091012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6058265) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(-0.95226007) q[1];
rz(-1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(-2.1624485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(0.34805527) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340654) q[0];
sx q[0];
rz(-2.4287927) q[0];
sx q[0];
rz(-0.45760052) q[0];
x q[1];
rz(1.1966755) q[2];
sx q[2];
rz(-1.7110363) q[2];
sx q[2];
rz(1.8082878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.066612331) q[1];
sx q[1];
rz(-0.96251026) q[1];
sx q[1];
rz(3.0327256) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9799558) q[3];
sx q[3];
rz(-1.4378387) q[3];
sx q[3];
rz(-2.0032351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7248914) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87930644) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(0.53043956) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.8431429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2753678) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(0.84262459) q[0];
rz(1.8680044) q[2];
sx q[2];
rz(-2.1564335) q[2];
sx q[2];
rz(2.1538018) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7057719) q[1];
sx q[1];
rz(-0.37322361) q[1];
sx q[1];
rz(-2.865764) q[1];
rz(-pi) q[2];
rz(-0.88921806) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-2.6679664) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(-2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-2.1437058) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-2.6748437) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022910718) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(0.33606152) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3495965) q[2];
sx q[2];
rz(-1.1525407) q[2];
sx q[2];
rz(2.7115371) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0007243) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
x q[2];
rz(-1.407133) q[3];
sx q[3];
rz(-1.6657077) q[3];
sx q[3];
rz(-1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2433462) q[0];
sx q[0];
rz(-1.6582489) q[0];
sx q[0];
rz(1.5792219) q[0];
rz(-pi) q[1];
rz(-0.98165841) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(-0.32064082) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6429813) q[1];
sx q[1];
rz(-1.4652518) q[1];
sx q[1];
rz(1.5555744) q[1];
rz(-pi) q[2];
rz(0.34479721) q[3];
sx q[3];
rz(-2.0119917) q[3];
sx q[3];
rz(-1.7896717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4411366) q[0];
sx q[0];
rz(-1.4608129) q[0];
sx q[0];
rz(2.3814047) q[0];
x q[1];
rz(-0.35347519) q[2];
sx q[2];
rz(-0.77958737) q[2];
sx q[2];
rz(2.0445063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1738759) q[1];
sx q[1];
rz(-1.516725) q[1];
sx q[1];
rz(1.4500344) q[1];
rz(-pi) q[2];
rz(0.30808361) q[3];
sx q[3];
rz(-0.83762533) q[3];
sx q[3];
rz(1.777491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(-2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-0.54661173) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(2.2792363) q[0];
rz(-pi) q[1];
rz(-1.1999646) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(1.3818936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4709028) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(1.1501269) q[1];
rz(-pi) q[2];
rz(0.96484465) q[3];
sx q[3];
rz(-0.75111872) q[3];
sx q[3];
rz(1.5497108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(1.9316797) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(-0.70770121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428064) q[0];
sx q[0];
rz(-1.2278623) q[0];
sx q[0];
rz(2.2862611) q[0];
rz(0.14214469) q[2];
sx q[2];
rz(-1.3488349) q[2];
sx q[2];
rz(2.1665426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0545132) q[1];
sx q[1];
rz(-0.88216773) q[1];
sx q[1];
rz(2.7282532) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(-2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83508867) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(0.46335285) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];