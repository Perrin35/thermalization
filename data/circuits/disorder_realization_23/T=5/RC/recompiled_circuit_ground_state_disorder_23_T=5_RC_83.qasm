OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(3.8138226) q[0];
sx q[0];
rz(11.091118) q[0];
rz(0.9530468) q[1];
sx q[1];
rz(-2.9763728) q[1];
sx q[1];
rz(1.2686977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27086285) q[0];
sx q[0];
rz(-2.6755973) q[0];
sx q[0];
rz(2.6129338) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2242429) q[2];
sx q[2];
rz(-1.1684061) q[2];
sx q[2];
rz(1.0831329) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70763146) q[1];
sx q[1];
rz(-0.58506706) q[1];
sx q[1];
rz(1.332704) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8016863) q[3];
sx q[3];
rz(-2.3143263) q[3];
sx q[3];
rz(-0.32318599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98422009) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(0.79322195) q[2];
rz(2.4685517) q[3];
sx q[3];
rz(-2.0561736) q[3];
sx q[3];
rz(-1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15758812) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(2.5322835) q[0];
rz(0.20978236) q[1];
sx q[1];
rz(-0.79133004) q[1];
sx q[1];
rz(0.9598859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6580556) q[0];
sx q[0];
rz(-1.5262881) q[0];
sx q[0];
rz(-0.07907148) q[0];
x q[1];
rz(-1.7018713) q[2];
sx q[2];
rz(-2.6496151) q[2];
sx q[2];
rz(0.78322369) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9896796) q[1];
sx q[1];
rz(-1.6335347) q[1];
sx q[1];
rz(-0.83338085) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0783362) q[3];
sx q[3];
rz(-1.1843269) q[3];
sx q[3];
rz(1.2852576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3557055) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(-2.6500224) q[2];
rz(0.0056886557) q[3];
sx q[3];
rz(-1.0892884) q[3];
sx q[3];
rz(0.51386851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095057644) q[0];
sx q[0];
rz(-0.53545606) q[0];
sx q[0];
rz(2.3989578) q[0];
rz(0.62478089) q[1];
sx q[1];
rz(-0.94015986) q[1];
sx q[1];
rz(2.0754441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0381706) q[0];
sx q[0];
rz(-0.26525149) q[0];
sx q[0];
rz(-2.1384504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.891648) q[2];
sx q[2];
rz(-1.751747) q[2];
sx q[2];
rz(0.89676434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6934476) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(-0.13685302) q[1];
rz(-pi) q[2];
rz(2.9227299) q[3];
sx q[3];
rz(-1.0893679) q[3];
sx q[3];
rz(-0.66314135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8891958) q[2];
sx q[2];
rz(-2.9434581) q[2];
sx q[2];
rz(0.55257094) q[2];
rz(-0.50734723) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(-2.5721917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6602537) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(-1.8950155) q[0];
rz(-1.4056816) q[1];
sx q[1];
rz(-2.3697) q[1];
sx q[1];
rz(0.71800047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7557186) q[0];
sx q[0];
rz(-1.0654952) q[0];
sx q[0];
rz(-0.85407172) q[0];
rz(-pi) q[1];
rz(1.0077072) q[2];
sx q[2];
rz(-1.7741084) q[2];
sx q[2];
rz(-0.85693923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93054616) q[1];
sx q[1];
rz(-2.3706243) q[1];
sx q[1];
rz(0.20937415) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9830789) q[3];
sx q[3];
rz(-1.6974546) q[3];
sx q[3];
rz(1.4617006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(0.38257515) q[3];
sx q[3];
rz(-1.4380598) q[3];
sx q[3];
rz(-2.2883435) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48069561) q[0];
sx q[0];
rz(-3.1013885) q[0];
sx q[0];
rz(0.30886343) q[0];
rz(-2.1924696) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(-2.5905051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4590722) q[0];
sx q[0];
rz(-2.7210752) q[0];
sx q[0];
rz(2.7608052) q[0];
rz(-pi) q[1];
rz(-1.6597719) q[2];
sx q[2];
rz(-2.4818261) q[2];
sx q[2];
rz(-0.24487409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70610395) q[1];
sx q[1];
rz(-1.7167673) q[1];
sx q[1];
rz(1.7273497) q[1];
x q[2];
rz(0.77648987) q[3];
sx q[3];
rz(-1.6740419) q[3];
sx q[3];
rz(-2.6342027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71517313) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(2.5227762) q[2];
rz(0.62301451) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(0.4272517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039336786) q[0];
sx q[0];
rz(-0.66348851) q[0];
sx q[0];
rz(2.6797507) q[0];
rz(2.6448008) q[1];
sx q[1];
rz(-1.3235612) q[1];
sx q[1];
rz(1.367548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753938) q[0];
sx q[0];
rz(-0.68635041) q[0];
sx q[0];
rz(-1.696287) q[0];
rz(-2.0678287) q[2];
sx q[2];
rz(-1.70514) q[2];
sx q[2];
rz(1.4984992) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16713472) q[1];
sx q[1];
rz(-0.30213812) q[1];
sx q[1];
rz(-1.664851) q[1];
x q[2];
rz(0.9487834) q[3];
sx q[3];
rz(-1.439073) q[3];
sx q[3];
rz(-1.6831116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0355012) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(-2.3512225) q[2];
rz(2.5370989) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-2.7214971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4496434) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(0.41437909) q[0];
rz(-0.62037933) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(-1.5827804) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5904424) q[0];
sx q[0];
rz(-1.4865459) q[0];
sx q[0];
rz(0.1131375) q[0];
rz(-pi) q[1];
rz(0.16838603) q[2];
sx q[2];
rz(-0.10659519) q[2];
sx q[2];
rz(1.3063947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61395105) q[1];
sx q[1];
rz(-1.5802586) q[1];
sx q[1];
rz(1.6285614) q[1];
x q[2];
rz(-1.956432) q[3];
sx q[3];
rz(-1.4298855) q[3];
sx q[3];
rz(0.27045617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9862426) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(-2.9637994) q[2];
rz(0.47884652) q[3];
sx q[3];
rz(-1.1136473) q[3];
sx q[3];
rz(-3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14684045) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(0.11257182) q[0];
rz(-0.62537891) q[1];
sx q[1];
rz(-1.1275147) q[1];
sx q[1];
rz(2.2523527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13876943) q[0];
sx q[0];
rz(-3.0516612) q[0];
sx q[0];
rz(-1.0732717) q[0];
x q[1];
rz(-0.26548141) q[2];
sx q[2];
rz(-1.6471787) q[2];
sx q[2];
rz(0.35778174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1804643) q[1];
sx q[1];
rz(-0.58149177) q[1];
sx q[1];
rz(2.0093469) q[1];
rz(3.0744138) q[3];
sx q[3];
rz(-1.3357196) q[3];
sx q[3];
rz(-0.76418166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2832977) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(-0.27321401) q[2];
rz(0.22935271) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44723085) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(0.91621512) q[0];
rz(0.18237309) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.1590385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0908546) q[0];
sx q[0];
rz(-1.3373475) q[0];
sx q[0];
rz(-1.7660733) q[0];
x q[1];
rz(-1.5196475) q[2];
sx q[2];
rz(-1.4754093) q[2];
sx q[2];
rz(0.21610987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43913821) q[1];
sx q[1];
rz(-1.7299486) q[1];
sx q[1];
rz(1.252878) q[1];
x q[2];
rz(-2.4156225) q[3];
sx q[3];
rz(-1.0683224) q[3];
sx q[3];
rz(0.53086764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75659043) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(2.8263367) q[2];
rz(0.59424019) q[3];
sx q[3];
rz(-2.2599594) q[3];
sx q[3];
rz(-2.5074904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2600128) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(0.15923937) q[0];
rz(-1.0881933) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(0.27096567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9054026) q[0];
sx q[0];
rz(-1.3199249) q[0];
sx q[0];
rz(-1.0948927) q[0];
x q[1];
rz(0.6728716) q[2];
sx q[2];
rz(-1.4904163) q[2];
sx q[2];
rz(-2.5774235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3828363) q[1];
sx q[1];
rz(-0.74302948) q[1];
sx q[1];
rz(-2.6322175) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5529717) q[3];
sx q[3];
rz(-0.56097066) q[3];
sx q[3];
rz(0.64631337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(-0.32269746) q[2];
rz(-3.0835551) q[3];
sx q[3];
rz(-0.80153424) q[3];
sx q[3];
rz(-0.29840741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.4651481) q[0];
sx q[0];
rz(-1.5100751) q[0];
sx q[0];
rz(-0.81612192) q[0];
rz(2.8836518) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(-1.8564687) q[2];
sx q[2];
rz(-1.2930047) q[2];
sx q[2];
rz(1.5563896) q[2];
rz(-1.1318351) q[3];
sx q[3];
rz(-2.1793032) q[3];
sx q[3];
rz(1.3201154) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
