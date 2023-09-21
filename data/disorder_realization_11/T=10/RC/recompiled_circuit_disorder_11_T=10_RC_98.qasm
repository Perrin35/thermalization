OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55460632) q[0];
sx q[0];
rz(-0.89590961) q[0];
sx q[0];
rz(-1.4045665) q[0];
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(-1.6552742) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5777187) q[0];
sx q[0];
rz(-1.3868252) q[0];
sx q[0];
rz(-1.4509393) q[0];
x q[1];
rz(1.4990184) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(2.3274802) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7662498) q[1];
sx q[1];
rz(-1.9613593) q[1];
sx q[1];
rz(1.6745964) q[1];
x q[2];
rz(1.5641883) q[3];
sx q[3];
rz(-1.4918543) q[3];
sx q[3];
rz(0.034949485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(3.063607) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29794824) q[0];
sx q[0];
rz(-2.1878562) q[0];
sx q[0];
rz(-0.88275568) q[0];
rz(-pi) q[1];
rz(-2.4215135) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(-0.0062696487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3188394) q[1];
sx q[1];
rz(-1.0471294) q[1];
sx q[1];
rz(-0.19284064) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1017786) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(-3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(-1.2603166) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-0.46762064) q[0];
rz(0.46472654) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25160015) q[0];
sx q[0];
rz(-1.359299) q[0];
sx q[0];
rz(-0.081758008) q[0];
rz(-2.8094693) q[2];
sx q[2];
rz(-0.87288522) q[2];
sx q[2];
rz(-2.3534564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1130502) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(2.261496) q[1];
rz(-pi) q[2];
rz(1.7160277) q[3];
sx q[3];
rz(-1.8947453) q[3];
sx q[3];
rz(0.63889972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(-2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(-2.4688597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8254291) q[0];
sx q[0];
rz(-2.5220693) q[0];
sx q[0];
rz(2.6114458) q[0];
rz(-pi) q[1];
rz(2.7364199) q[2];
sx q[2];
rz(-2.5537958) q[2];
sx q[2];
rz(1.6853465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.67191254) q[1];
sx q[1];
rz(-1.3393991) q[1];
sx q[1];
rz(-3.0242821) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7476051) q[3];
sx q[3];
rz(-0.62649957) q[3];
sx q[3];
rz(0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(-1.6617552) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(3.1033049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0262888) q[0];
sx q[0];
rz(-2.8235108) q[0];
sx q[0];
rz(-1.173107) q[0];
rz(-pi) q[1];
rz(-0.64735909) q[2];
sx q[2];
rz(-2.416496) q[2];
sx q[2];
rz(-0.65078306) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8748862) q[1];
sx q[1];
rz(-1.453207) q[1];
sx q[1];
rz(-1.7456732) q[1];
rz(-pi) q[2];
rz(-2.0376301) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(1.8339748) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6742498) q[0];
sx q[0];
rz(-1.7102339) q[0];
sx q[0];
rz(1.3583202) q[0];
rz(-pi) q[1];
rz(2.7912748) q[2];
sx q[2];
rz(-0.91674524) q[2];
sx q[2];
rz(2.0803114) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29589614) q[1];
sx q[1];
rz(-2.1138043) q[1];
sx q[1];
rz(3.1177109) q[1];
rz(-1.0530924) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(2.6368486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(0.42022589) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(2.1113254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-2.6269873) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(-2.1234535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87228862) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(1.3809204) q[1];
rz(0.19018634) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(-2.412652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(-2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(1.9497005) q[0];
rz(2.2299956) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(1.0620767) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870283) q[0];
sx q[0];
rz(-1.6977068) q[0];
sx q[0];
rz(-1.0303322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1658792) q[2];
sx q[2];
rz(-1.1761883) q[2];
sx q[2];
rz(-1.6369866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8909104) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(-pi) q[2];
rz(0.90537269) q[3];
sx q[3];
rz(-1.5432127) q[3];
sx q[3];
rz(-1.8679004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97201559) q[0];
sx q[0];
rz(-1.2017842) q[0];
sx q[0];
rz(-0.24032648) q[0];
rz(1.4184985) q[2];
sx q[2];
rz(-1.9405775) q[2];
sx q[2];
rz(2.5068482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4832284) q[1];
sx q[1];
rz(-1.3434957) q[1];
sx q[1];
rz(-1.8140273) q[1];
rz(-2.455515) q[3];
sx q[3];
rz(-0.25732532) q[3];
sx q[3];
rz(-2.0285332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(-0.84154877) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(-0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(2.9291005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1157163) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(-0.91761597) q[0];
rz(-1.882952) q[2];
sx q[2];
rz(-1.5317132) q[2];
sx q[2];
rz(-2.8651819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54086911) q[1];
sx q[1];
rz(-1.0717372) q[1];
sx q[1];
rz(-0.869511) q[1];
x q[2];
rz(0.28189567) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(-2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385948) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(-1.8226345) q[2];
sx q[2];
rz(-3.0694198) q[2];
sx q[2];
rz(1.0138489) q[2];
rz(-2.6247737) q[3];
sx q[3];
rz(-1.9782981) q[3];
sx q[3];
rz(0.78936418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];