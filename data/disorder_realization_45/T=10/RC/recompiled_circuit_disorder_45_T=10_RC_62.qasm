OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8279995) q[0];
sx q[0];
rz(-1.8090973) q[0];
sx q[0];
rz(-0.39512623) q[0];
x q[1];
rz(2.2864561) q[2];
sx q[2];
rz(-0.54358608) q[2];
sx q[2];
rz(0.91993514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3798843) q[1];
sx q[1];
rz(-1.1176846) q[1];
sx q[1];
rz(2.3945827) q[1];
rz(-pi) q[2];
rz(1.9786505) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(-0.70139635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.550094) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(-2.1226728) q[0];
rz(-1.4787176) q[1];
sx q[1];
rz(-0.61518413) q[1];
sx q[1];
rz(-0.63308024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3490152) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(-2.3809459) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1585629) q[2];
sx q[2];
rz(-2.2289742) q[2];
sx q[2];
rz(2.7764729) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9620348) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(2.7707997) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8301864) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(0.11745545) q[2];
rz(0.30101267) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.85787073) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(3.0531847) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-2.450768) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.847825) q[0];
sx q[0];
rz(-2.9710037) q[0];
sx q[0];
rz(-0.45844309) q[0];
x q[1];
rz(1.8059398) q[2];
sx q[2];
rz(-0.81980875) q[2];
sx q[2];
rz(-2.2355134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5882727) q[1];
sx q[1];
rz(-1.5447445) q[1];
sx q[1];
rz(-1.7528898) q[1];
x q[2];
rz(-0.74294528) q[3];
sx q[3];
rz(-0.49693123) q[3];
sx q[3];
rz(-3.1019773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(-2.5879588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.722256) q[0];
sx q[0];
rz(-1.1124848) q[0];
sx q[0];
rz(-0.76461794) q[0];
rz(2.373898) q[2];
sx q[2];
rz(-2.1622217) q[2];
sx q[2];
rz(-0.64085811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8071027) q[1];
sx q[1];
rz(-0.79542167) q[1];
sx q[1];
rz(2.9544178) q[1];
rz(2.3447498) q[3];
sx q[3];
rz(-1.5505655) q[3];
sx q[3];
rz(-0.86864558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(2.6468357) q[2];
rz(0.90302145) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(1.0850798) q[0];
rz(2.1814573) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(2.9575612) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83443123) q[0];
sx q[0];
rz(-2.8330028) q[0];
sx q[0];
rz(2.79074) q[0];
x q[1];
rz(-1.2940302) q[2];
sx q[2];
rz(-1.2185316) q[2];
sx q[2];
rz(2.2500452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80464333) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(-0.6964535) q[1];
rz(-pi) q[2];
rz(-2.4155951) q[3];
sx q[3];
rz(-1.7222002) q[3];
sx q[3];
rz(-2.9263673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(2.0007755) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891156) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(-0.14850798) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214971) q[0];
sx q[0];
rz(-1.7317061) q[0];
sx q[0];
rz(-0.73385977) q[0];
rz(-1.6806904) q[2];
sx q[2];
rz(-1.3434778) q[2];
sx q[2];
rz(-1.1306888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.442246) q[1];
sx q[1];
rz(-1.3117426) q[1];
sx q[1];
rz(3.0287292) q[1];
rz(2.4466483) q[3];
sx q[3];
rz(-0.57999014) q[3];
sx q[3];
rz(0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.274652) q[0];
sx q[0];
rz(-2.7526703) q[0];
sx q[0];
rz(-1.2961943) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.869404) q[2];
sx q[2];
rz(-1.4223137) q[2];
sx q[2];
rz(-2.6872925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69562558) q[1];
sx q[1];
rz(-1.8672767) q[1];
sx q[1];
rz(-2.8545024) q[1];
rz(1.0057698) q[3];
sx q[3];
rz(-0.98494782) q[3];
sx q[3];
rz(0.98423959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0685048) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.9203141) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(1.01064) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0518274) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(0.21659539) q[0];
rz(-2.0881565) q[2];
sx q[2];
rz(-0.55609716) q[2];
sx q[2];
rz(-0.14511395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9606564) q[1];
sx q[1];
rz(-1.3410543) q[1];
sx q[1];
rz(0.22682637) q[1];
x q[2];
rz(-2.741022) q[3];
sx q[3];
rz(-1.6555602) q[3];
sx q[3];
rz(0.29464196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6179787) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(-2.5726035) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0555608) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(2.6468497) q[0];
rz(2.5231979) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2878802) q[0];
sx q[0];
rz(-1.259385) q[0];
sx q[0];
rz(2.1561554) q[0];
rz(-pi) q[1];
rz(-1.2443301) q[2];
sx q[2];
rz(-2.9565161) q[2];
sx q[2];
rz(1.9290123) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6203783) q[1];
sx q[1];
rz(-1.673939) q[1];
sx q[1];
rz(0.68876536) q[1];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(-1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(-2.7005844) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84687418) q[0];
sx q[0];
rz(-2.8537769) q[0];
sx q[0];
rz(0.78886445) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4961692) q[2];
sx q[2];
rz(-1.7241038) q[2];
sx q[2];
rz(0.15664936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8969352) q[1];
sx q[1];
rz(-1.999563) q[1];
sx q[1];
rz(0.37533092) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8701843) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(-2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-0.30187541) q[2];
rz(-0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(3.1148615) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(2.6731861) q[2];
sx q[2];
rz(-0.9006587) q[2];
sx q[2];
rz(-2.547154) q[2];
rz(-2.0704913) q[3];
sx q[3];
rz(-2.069996) q[3];
sx q[3];
rz(1.3531006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
