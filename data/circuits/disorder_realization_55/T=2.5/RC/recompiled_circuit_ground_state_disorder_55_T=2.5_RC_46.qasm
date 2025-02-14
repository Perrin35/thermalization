OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3352873) q[0];
sx q[0];
rz(-1.0399613) q[0];
sx q[0];
rz(2.7916173) q[0];
rz(-1.6236053) q[1];
sx q[1];
rz(-0.94437683) q[1];
sx q[1];
rz(-1.1854393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1040504) q[0];
sx q[0];
rz(-0.54706956) q[0];
sx q[0];
rz(0.68010917) q[0];
rz(-pi) q[1];
rz(-2.7230206) q[2];
sx q[2];
rz(-1.7836193) q[2];
sx q[2];
rz(0.71453071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9361286) q[1];
sx q[1];
rz(-1.3359114) q[1];
sx q[1];
rz(-1.0297736) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9583232) q[3];
sx q[3];
rz(-0.96310593) q[3];
sx q[3];
rz(-2.303108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8650633) q[2];
sx q[2];
rz(-0.41215602) q[2];
sx q[2];
rz(0.40974799) q[2];
rz(0.77541882) q[3];
sx q[3];
rz(-1.5284458) q[3];
sx q[3];
rz(2.8342136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85768098) q[0];
sx q[0];
rz(-0.49864054) q[0];
sx q[0];
rz(2.6890802) q[0];
rz(0.2287989) q[1];
sx q[1];
rz(-2.6401873) q[1];
sx q[1];
rz(0.14678821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2499224) q[0];
sx q[0];
rz(-1.2091214) q[0];
sx q[0];
rz(-2.937243) q[0];
rz(-pi) q[1];
rz(1.3323305) q[2];
sx q[2];
rz(-1.9386282) q[2];
sx q[2];
rz(-0.17300082) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2816946) q[1];
sx q[1];
rz(-1.7051776) q[1];
sx q[1];
rz(-1.6107537) q[1];
rz(2.4815219) q[3];
sx q[3];
rz(-1.1997031) q[3];
sx q[3];
rz(-2.0160272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9734834) q[2];
sx q[2];
rz(-2.8066469) q[2];
sx q[2];
rz(0.87146634) q[2];
rz(-2.787309) q[3];
sx q[3];
rz(-2.2058217) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63534921) q[0];
sx q[0];
rz(-2.5214218) q[0];
sx q[0];
rz(1.1059906) q[0];
rz(2.1982819) q[1];
sx q[1];
rz(-0.7990852) q[1];
sx q[1];
rz(1.4897289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8603926) q[0];
sx q[0];
rz(-1.6985122) q[0];
sx q[0];
rz(-3.1312577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6666127) q[2];
sx q[2];
rz(-1.7783217) q[2];
sx q[2];
rz(0.40643164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.572444) q[1];
sx q[1];
rz(-0.97578543) q[1];
sx q[1];
rz(1.6165363) q[1];
rz(0.23814985) q[3];
sx q[3];
rz(-1.643073) q[3];
sx q[3];
rz(1.0840462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6279383) q[2];
sx q[2];
rz(-1.9460121) q[2];
sx q[2];
rz(0.99620831) q[2];
rz(2.5683688) q[3];
sx q[3];
rz(-0.51155353) q[3];
sx q[3];
rz(0.63428026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298249) q[0];
sx q[0];
rz(-1.8695762) q[0];
sx q[0];
rz(1.5559394) q[0];
rz(-0.32444435) q[1];
sx q[1];
rz(-2.3093846) q[1];
sx q[1];
rz(1.197804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75173855) q[0];
sx q[0];
rz(-2.1445988) q[0];
sx q[0];
rz(1.874254) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.156331) q[2];
sx q[2];
rz(-1.6949495) q[2];
sx q[2];
rz(-0.92701605) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50655925) q[1];
sx q[1];
rz(-0.50093073) q[1];
sx q[1];
rz(1.2975724) q[1];
x q[2];
rz(-1.8327375) q[3];
sx q[3];
rz(-1.7290232) q[3];
sx q[3];
rz(-2.5867274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0962301) q[2];
sx q[2];
rz(-2.5924293) q[2];
sx q[2];
rz(-1.9722923) q[2];
rz(-0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(-2.4128778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28155926) q[0];
sx q[0];
rz(-0.39880729) q[0];
sx q[0];
rz(1.9418203) q[0];
rz(1.3072321) q[1];
sx q[1];
rz(-0.93991005) q[1];
sx q[1];
rz(1.7050381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2685674) q[0];
sx q[0];
rz(-1.5102666) q[0];
sx q[0];
rz(-1.1824754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13135037) q[2];
sx q[2];
rz(-1.7997912) q[2];
sx q[2];
rz(-0.74818767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2695527) q[1];
sx q[1];
rz(-1.8694573) q[1];
sx q[1];
rz(-0.36901094) q[1];
rz(0.41884274) q[3];
sx q[3];
rz(-2.1757033) q[3];
sx q[3];
rz(-2.1640282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3522219) q[2];
sx q[2];
rz(-1.7033966) q[2];
sx q[2];
rz(0.5729202) q[2];
rz(-0.88417435) q[3];
sx q[3];
rz(-2.1638162) q[3];
sx q[3];
rz(2.6463215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5578515) q[0];
sx q[0];
rz(-1.252625) q[0];
sx q[0];
rz(-2.0649233) q[0];
rz(2.63511) q[1];
sx q[1];
rz(-0.61512893) q[1];
sx q[1];
rz(-1.8411676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5426555) q[0];
sx q[0];
rz(-0.050572336) q[0];
sx q[0];
rz(1.2957606) q[0];
rz(-pi) q[1];
rz(0.22407786) q[2];
sx q[2];
rz(-0.59106088) q[2];
sx q[2];
rz(-3.0853809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4757931) q[1];
sx q[1];
rz(-0.25036004) q[1];
sx q[1];
rz(1.0151035) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9415116) q[3];
sx q[3];
rz(-1.9164547) q[3];
sx q[3];
rz(2.8576771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7252245) q[2];
sx q[2];
rz(-0.84364265) q[2];
sx q[2];
rz(-1.1923403) q[2];
rz(-0.013817712) q[3];
sx q[3];
rz(-2.4255224) q[3];
sx q[3];
rz(2.1684087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8927807) q[0];
sx q[0];
rz(-2.2154494) q[0];
sx q[0];
rz(2.7698351) q[0];
rz(0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(-0.54381347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11228526) q[0];
sx q[0];
rz(-0.72192955) q[0];
sx q[0];
rz(0.63755905) q[0];
rz(-2.2901754) q[2];
sx q[2];
rz(-1.2179551) q[2];
sx q[2];
rz(-2.1746496) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2799171) q[1];
sx q[1];
rz(-1.1959198) q[1];
sx q[1];
rz(1.5779297) q[1];
rz(-pi) q[2];
rz(1.4187846) q[3];
sx q[3];
rz(-1.4626039) q[3];
sx q[3];
rz(-0.044807981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.030293327) q[2];
sx q[2];
rz(-1.0617504) q[2];
sx q[2];
rz(-2.566805) q[2];
rz(-2.287367) q[3];
sx q[3];
rz(-1.4762907) q[3];
sx q[3];
rz(-2.0265719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2146571) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(2.8171203) q[0];
rz(0.60943162) q[1];
sx q[1];
rz(-2.29988) q[1];
sx q[1];
rz(-2.6268974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99583944) q[0];
sx q[0];
rz(-2.0059364) q[0];
sx q[0];
rz(-3.1313722) q[0];
x q[1];
rz(-1.9799405) q[2];
sx q[2];
rz(-2.3858641) q[2];
sx q[2];
rz(2.8314396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60694164) q[1];
sx q[1];
rz(-1.3365627) q[1];
sx q[1];
rz(0.57778426) q[1];
rz(-pi) q[2];
rz(1.2013161) q[3];
sx q[3];
rz(-0.29467988) q[3];
sx q[3];
rz(-2.0979543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1631761) q[2];
sx q[2];
rz(-1.6360444) q[2];
sx q[2];
rz(-2.8669538) q[2];
rz(0.93242532) q[3];
sx q[3];
rz(-0.43007389) q[3];
sx q[3];
rz(-0.29621223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4579929) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(0.13548166) q[0];
rz(-0.97700351) q[1];
sx q[1];
rz(-0.75825399) q[1];
sx q[1];
rz(-3.1405084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8196306) q[0];
sx q[0];
rz(-1.0792097) q[0];
sx q[0];
rz(-1.7474773) q[0];
rz(-pi) q[1];
rz(1.7747701) q[2];
sx q[2];
rz(-0.53278953) q[2];
sx q[2];
rz(0.29292695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3848968) q[1];
sx q[1];
rz(-2.2731509) q[1];
sx q[1];
rz(0.76354247) q[1];
rz(0.67768456) q[3];
sx q[3];
rz(-0.88078684) q[3];
sx q[3];
rz(3.1180988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5929426) q[2];
sx q[2];
rz(-2.103001) q[2];
sx q[2];
rz(2.9366117) q[2];
rz(-0.18260469) q[3];
sx q[3];
rz(-0.20668106) q[3];
sx q[3];
rz(-0.73777795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36271998) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(0.47738281) q[0];
rz(-0.1167156) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(2.3027072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8461743) q[0];
sx q[0];
rz(-1.8356165) q[0];
sx q[0];
rz(0.85746641) q[0];
x q[1];
rz(1.9220542) q[2];
sx q[2];
rz(-0.87703401) q[2];
sx q[2];
rz(-1.1379682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6179671) q[1];
sx q[1];
rz(-0.8899506) q[1];
sx q[1];
rz(1.8115854) q[1];
rz(-pi) q[2];
rz(1.4492744) q[3];
sx q[3];
rz(-2.9009933) q[3];
sx q[3];
rz(2.6783063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3597151) q[2];
sx q[2];
rz(-2.8350267) q[2];
sx q[2];
rz(-2.8753885) q[2];
rz(2.6340458) q[3];
sx q[3];
rz(-2.4188953) q[3];
sx q[3];
rz(0.44107309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97892852) q[0];
sx q[0];
rz(-1.4586108) q[0];
sx q[0];
rz(2.3391531) q[0];
rz(-2.2294527) q[1];
sx q[1];
rz(-1.9504539) q[1];
sx q[1];
rz(0.84074195) q[1];
rz(-2.5723614) q[2];
sx q[2];
rz(-0.29217671) q[2];
sx q[2];
rz(3.0859699) q[2];
rz(1.1900564) q[3];
sx q[3];
rz(-1.1968812) q[3];
sx q[3];
rz(1.4637917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
