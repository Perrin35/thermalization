OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7492548) q[0];
sx q[0];
rz(5.9049913) q[0];
sx q[0];
rz(9.7752934) q[0];
rz(0.89138436) q[1];
sx q[1];
rz(-2.3494224) q[1];
sx q[1];
rz(-1.9686735) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25305155) q[0];
sx q[0];
rz(-1.3810147) q[0];
sx q[0];
rz(-1.7642598) q[0];
x q[1];
rz(0.15526659) q[2];
sx q[2];
rz(-1.7874996) q[2];
sx q[2];
rz(-2.7086176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80591737) q[1];
sx q[1];
rz(-1.8875673) q[1];
sx q[1];
rz(-1.3992228) q[1];
rz(-pi) q[2];
rz(-2.6698106) q[3];
sx q[3];
rz(-1.4174479) q[3];
sx q[3];
rz(0.23424304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(2.9288536) q[2];
rz(3.0563266) q[3];
sx q[3];
rz(-1.250896) q[3];
sx q[3];
rz(-1.8703478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79412115) q[0];
sx q[0];
rz(-0.81231064) q[0];
sx q[0];
rz(-1.043327) q[0];
rz(-0.13283816) q[1];
sx q[1];
rz(-2.9265407) q[1];
sx q[1];
rz(-1.1588233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83860676) q[0];
sx q[0];
rz(-1.9667454) q[0];
sx q[0];
rz(-1.7274169) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6558596) q[2];
sx q[2];
rz(-2.9797005) q[2];
sx q[2];
rz(-0.4134824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0483425) q[1];
sx q[1];
rz(-2.5777317) q[1];
sx q[1];
rz(-2.6287931) q[1];
rz(2.0015249) q[3];
sx q[3];
rz(-1.8716806) q[3];
sx q[3];
rz(-1.6611851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(0.72466737) q[2];
rz(0.66257462) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(0.29964963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1835943) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(2.2268353) q[0];
rz(-0.58685189) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(2.9920726) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5356333) q[0];
sx q[0];
rz(-2.0894755) q[0];
sx q[0];
rz(2.9796322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9376041) q[2];
sx q[2];
rz(-1.2292394) q[2];
sx q[2];
rz(1.8579409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5352763) q[1];
sx q[1];
rz(-0.80937591) q[1];
sx q[1];
rz(-0.80324976) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4305184) q[3];
sx q[3];
rz(-2.1319763) q[3];
sx q[3];
rz(-0.64519889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84825039) q[2];
sx q[2];
rz(-1.205227) q[2];
sx q[2];
rz(0.80043522) q[2];
rz(-0.46659255) q[3];
sx q[3];
rz(-2.0355909) q[3];
sx q[3];
rz(-0.65565562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7993497) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(-0.85860646) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(-0.98791775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8289017) q[0];
sx q[0];
rz(-0.36417555) q[0];
sx q[0];
rz(-2.0065424) q[0];
rz(1.0756827) q[2];
sx q[2];
rz(-1.6778113) q[2];
sx q[2];
rz(-2.1008976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4639484) q[1];
sx q[1];
rz(-1.8419187) q[1];
sx q[1];
rz(1.2605002) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8186422) q[3];
sx q[3];
rz(-1.3846372) q[3];
sx q[3];
rz(1.3398088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19114384) q[2];
sx q[2];
rz(-1.3942275) q[2];
sx q[2];
rz(-0.44551715) q[2];
rz(0.71349239) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-0.92002404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54818654) q[0];
sx q[0];
rz(-2.0596518) q[0];
sx q[0];
rz(-1.5418381) q[0];
rz(1.2708739) q[1];
sx q[1];
rz(-1.7393232) q[1];
sx q[1];
rz(-0.52938968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177757) q[0];
sx q[0];
rz(-0.43250205) q[0];
sx q[0];
rz(-0.94571094) q[0];
rz(-2.5240493) q[2];
sx q[2];
rz(-1.7360592) q[2];
sx q[2];
rz(-2.3338855) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3965949) q[1];
sx q[1];
rz(-1.0675808) q[1];
sx q[1];
rz(-1.8580336) q[1];
rz(2.9934817) q[3];
sx q[3];
rz(-1.221773) q[3];
sx q[3];
rz(2.7452041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14220898) q[2];
sx q[2];
rz(-1.9500407) q[2];
sx q[2];
rz(-1.3231529) q[2];
rz(0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2926272) q[0];
sx q[0];
rz(-2.3830074) q[0];
sx q[0];
rz(-2.1204156) q[0];
rz(-1.5317597) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-0.80345947) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96413104) q[0];
sx q[0];
rz(-1.877907) q[0];
sx q[0];
rz(1.7384647) q[0];
x q[1];
rz(1.6551465) q[2];
sx q[2];
rz(-1.1599031) q[2];
sx q[2];
rz(1.8253758) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72587126) q[1];
sx q[1];
rz(-0.5827924) q[1];
sx q[1];
rz(-2.9225213) q[1];
x q[2];
rz(-0.62688503) q[3];
sx q[3];
rz(-1.6520368) q[3];
sx q[3];
rz(-0.33259847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44567406) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(2.7109801) q[2];
rz(-0.63498354) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(-1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1465313) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(-1.5978285) q[0];
rz(2.457288) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(2.3421471) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22173026) q[0];
sx q[0];
rz(-1.5654025) q[0];
sx q[0];
rz(0.07660596) q[0];
x q[1];
rz(-3.0657049) q[2];
sx q[2];
rz(-0.85691707) q[2];
sx q[2];
rz(0.80684987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3189074) q[1];
sx q[1];
rz(-0.46055183) q[1];
sx q[1];
rz(0.022629398) q[1];
rz(-0.83752172) q[3];
sx q[3];
rz(-1.9068516) q[3];
sx q[3];
rz(0.94810644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53576175) q[2];
sx q[2];
rz(-0.83767086) q[2];
sx q[2];
rz(2.3568995) q[2];
rz(-3.0253518) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(1.2522662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5410974) q[0];
sx q[0];
rz(-1.9116115) q[0];
sx q[0];
rz(2.1121209) q[0];
rz(3.0211499) q[1];
sx q[1];
rz(-1.30013) q[1];
sx q[1];
rz(2.5708503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1627038) q[0];
sx q[0];
rz(-1.9523639) q[0];
sx q[0];
rz(-1.0779774) q[0];
rz(0.71855259) q[2];
sx q[2];
rz(-1.6552123) q[2];
sx q[2];
rz(-2.8511503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6315122) q[1];
sx q[1];
rz(-1.8951547) q[1];
sx q[1];
rz(1.1034758) q[1];
x q[2];
rz(2.0550957) q[3];
sx q[3];
rz(-0.37621337) q[3];
sx q[3];
rz(-1.0602601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8735147) q[2];
sx q[2];
rz(-0.87468481) q[2];
sx q[2];
rz(0.52379215) q[2];
rz(-1.1526456) q[3];
sx q[3];
rz(-2.5609784) q[3];
sx q[3];
rz(-1.7136278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48083392) q[0];
sx q[0];
rz(-1.2013712) q[0];
sx q[0];
rz(2.7177366) q[0];
rz(-2.2747874) q[1];
sx q[1];
rz(-1.024217) q[1];
sx q[1];
rz(-0.20326916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05860672) q[0];
sx q[0];
rz(-1.5241677) q[0];
sx q[0];
rz(0.65445047) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16746232) q[2];
sx q[2];
rz(-0.31412087) q[2];
sx q[2];
rz(-0.7213074) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8228559) q[1];
sx q[1];
rz(-1.7137495) q[1];
sx q[1];
rz(2.8566384) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8762693) q[3];
sx q[3];
rz(-1.3246228) q[3];
sx q[3];
rz(-2.9873029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0554492) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(-2.8622368) q[2];
rz(1.2934359) q[3];
sx q[3];
rz(-1.8297628) q[3];
sx q[3];
rz(1.2156585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9781037) q[0];
sx q[0];
rz(-2.3390529) q[0];
sx q[0];
rz(-1.4168903) q[0];
rz(2.2143927) q[1];
sx q[1];
rz(-1.5617153) q[1];
sx q[1];
rz(-1.4097811) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7977236) q[0];
sx q[0];
rz(-0.90995212) q[0];
sx q[0];
rz(3.0190574) q[0];
rz(2.0538834) q[2];
sx q[2];
rz(-1.8894686) q[2];
sx q[2];
rz(1.1181732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7076787) q[1];
sx q[1];
rz(-1.1362038) q[1];
sx q[1];
rz(-2.4151925) q[1];
rz(-pi) q[2];
rz(3.0994371) q[3];
sx q[3];
rz(-1.8131527) q[3];
sx q[3];
rz(2.4740296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7659144) q[2];
sx q[2];
rz(-1.3216852) q[2];
sx q[2];
rz(0.970617) q[2];
rz(-2.4961903) q[3];
sx q[3];
rz(-2.161945) q[3];
sx q[3];
rz(1.0055044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29393016) q[0];
sx q[0];
rz(-1.4959338) q[0];
sx q[0];
rz(-1.3069859) q[0];
rz(-2.7597799) q[1];
sx q[1];
rz(-1.7996856) q[1];
sx q[1];
rz(-2.3736384) q[1];
rz(-2.1099595) q[2];
sx q[2];
rz(-0.8349541) q[2];
sx q[2];
rz(-0.63944774) q[2];
rz(-2.0313203) q[3];
sx q[3];
rz(-1.9721748) q[3];
sx q[3];
rz(1.0655793) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
