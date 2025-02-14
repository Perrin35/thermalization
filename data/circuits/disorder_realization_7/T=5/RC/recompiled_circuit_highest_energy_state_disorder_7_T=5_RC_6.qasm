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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(-0.74695865) q[1];
sx q[1];
rz(0.42626122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3328457) q[0];
sx q[0];
rz(-1.7431743) q[0];
sx q[0];
rz(0.7953978) q[0];
x q[1];
rz(2.6381016) q[2];
sx q[2];
rz(-1.3174743) q[2];
sx q[2];
rz(1.5510786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1805686) q[1];
sx q[1];
rz(-0.60446179) q[1];
sx q[1];
rz(3.040041) q[1];
rz(2.0419042) q[3];
sx q[3];
rz(-0.63835164) q[3];
sx q[3];
rz(1.9698576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2446186) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(-2.8986325) q[2];
rz(-0.36863676) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31545562) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(2.499089) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314302) q[0];
sx q[0];
rz(-1.1513984) q[0];
sx q[0];
rz(0.51879306) q[0];
rz(-1.9733866) q[2];
sx q[2];
rz(-1.0559096) q[2];
sx q[2];
rz(1.8414611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31256235) q[1];
sx q[1];
rz(-2.1305269) q[1];
sx q[1];
rz(-2.1055431) q[1];
rz(3.0152937) q[3];
sx q[3];
rz(-2.0492611) q[3];
sx q[3];
rz(-2.0984405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.7712234) q[2];
rz(-0.227452) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(-1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(-1.1981717) q[0];
rz(-1.0802957) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(0.094873039) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5115882) q[0];
sx q[0];
rz(-0.60902022) q[0];
sx q[0];
rz(2.0320351) q[0];
rz(1.0611141) q[2];
sx q[2];
rz(-0.72619263) q[2];
sx q[2];
rz(-0.92483556) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4725634) q[1];
sx q[1];
rz(-1.5298784) q[1];
sx q[1];
rz(-1.0340235) q[1];
rz(-pi) q[2];
rz(-1.9547988) q[3];
sx q[3];
rz(-1.8285164) q[3];
sx q[3];
rz(2.6741762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0692811) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(1.9629924) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(1.1627722) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56226319) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(-2.80559) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(3.0373108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096237) q[0];
sx q[0];
rz(-1.2237566) q[0];
sx q[0];
rz(-0.35023679) q[0];
rz(-2.9851172) q[2];
sx q[2];
rz(-1.556852) q[2];
sx q[2];
rz(0.61323159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1946009) q[1];
sx q[1];
rz(-0.90819383) q[1];
sx q[1];
rz(-0.47269447) q[1];
x q[2];
rz(-0.4533259) q[3];
sx q[3];
rz(-1.6764055) q[3];
sx q[3];
rz(2.3345514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5229554) q[2];
sx q[2];
rz(-2.0433661) q[2];
sx q[2];
rz(-1.9970419) q[2];
rz(-0.54730225) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(-2.0539637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2523786) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(2.2303708) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(0.30141452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2269991) q[0];
sx q[0];
rz(-0.98557878) q[0];
sx q[0];
rz(2.8666158) q[0];
rz(1.7784987) q[2];
sx q[2];
rz(-1.951521) q[2];
sx q[2];
rz(0.31685874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6993048) q[1];
sx q[1];
rz(-1.7831452) q[1];
sx q[1];
rz(-1.9599171) q[1];
rz(-pi) q[2];
rz(0.31378515) q[3];
sx q[3];
rz(-1.2126392) q[3];
sx q[3];
rz(1.4729983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-3.1346698) q[2];
rz(0.93112469) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(2.0324223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865006) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(-1.6078) q[0];
rz(2.4389229) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(-0.30219561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077009) q[0];
sx q[0];
rz(-2.2284628) q[0];
sx q[0];
rz(3.012863) q[0];
x q[1];
rz(1.9230559) q[2];
sx q[2];
rz(-0.88636905) q[2];
sx q[2];
rz(0.5796488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51660358) q[1];
sx q[1];
rz(-1.7440197) q[1];
sx q[1];
rz(2.1982963) q[1];
x q[2];
rz(1.4190361) q[3];
sx q[3];
rz(-1.5192864) q[3];
sx q[3];
rz(2.6035978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89221421) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(-1.9269491) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-0.76699081) q[0];
sx q[0];
rz(1.1085229) q[0];
rz(-0.33379894) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(-1.0858067) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.007581) q[0];
sx q[0];
rz(-1.3659048) q[0];
sx q[0];
rz(-1.6519068) q[0];
x q[1];
rz(-1.9778697) q[2];
sx q[2];
rz(-0.43292871) q[2];
sx q[2];
rz(2.0100398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9232417) q[1];
sx q[1];
rz(-2.3585547) q[1];
sx q[1];
rz(-0.79773517) q[1];
rz(-pi) q[2];
rz(0.13768519) q[3];
sx q[3];
rz(-2.2028649) q[3];
sx q[3];
rz(-1.9683624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51599017) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(-1.9773352) q[2];
rz(0.26816756) q[3];
sx q[3];
rz(-1.8987013) q[3];
sx q[3];
rz(-0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.1489768) q[0];
rz(0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(2.099096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5940018) q[0];
sx q[0];
rz(-2.2896575) q[0];
sx q[0];
rz(-0.44525624) q[0];
rz(1.8015968) q[2];
sx q[2];
rz(-0.79446213) q[2];
sx q[2];
rz(-0.59150254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62952215) q[1];
sx q[1];
rz(-2.7334556) q[1];
sx q[1];
rz(-0.44993181) q[1];
rz(-pi) q[2];
x q[2];
rz(0.028826272) q[3];
sx q[3];
rz(-0.90698637) q[3];
sx q[3];
rz(-3.1026866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59757549) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.3647122) q[2];
rz(0.65258604) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6259916) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(-3.1242477) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(-2.9877072) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8872334) q[0];
sx q[0];
rz(-1.4951997) q[0];
sx q[0];
rz(-0.31799728) q[0];
x q[1];
rz(-2.4272404) q[2];
sx q[2];
rz(-2.0660984) q[2];
sx q[2];
rz(2.60204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9964136) q[1];
sx q[1];
rz(-2.6657678) q[1];
sx q[1];
rz(-2.4563172) q[1];
rz(1.4948972) q[3];
sx q[3];
rz(-1.5468949) q[3];
sx q[3];
rz(2.0233512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1880356) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-0.44450644) q[2];
rz(0.19017531) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1403777) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(1.0706527) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-2.3505223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85002029) q[0];
sx q[0];
rz(-2.794461) q[0];
sx q[0];
rz(-1.8308543) q[0];
rz(-pi) q[1];
rz(-0.58035459) q[2];
sx q[2];
rz(-0.46541801) q[2];
sx q[2];
rz(-0.67829715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.4262096) q[1];
sx q[1];
rz(-1.9265623) q[1];
sx q[1];
rz(0.4224311) q[1];
rz(1.8077352) q[3];
sx q[3];
rz(-0.51722368) q[3];
sx q[3];
rz(0.57843971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.644824) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(-1.4727288) q[2];
rz(1.5276927) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(-1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2033757) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(0.51364246) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(-0.41507872) q[2];
sx q[2];
rz(-0.54176258) q[2];
sx q[2];
rz(-1.1272346) q[2];
rz(-2.0075825) q[3];
sx q[3];
rz(-0.81098771) q[3];
sx q[3];
rz(0.24843957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
