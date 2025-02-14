OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9031653) q[0];
sx q[0];
rz(2.981346) q[0];
sx q[0];
rz(7.5709406) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(-0.75872672) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37739524) q[0];
sx q[0];
rz(-1.8599417) q[0];
sx q[0];
rz(1.8125066) q[0];
rz(-2.8992462) q[2];
sx q[2];
rz(-2.0449315) q[2];
sx q[2];
rz(0.15210064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.030144) q[1];
sx q[1];
rz(-1.9322825) q[1];
sx q[1];
rz(2.3356781) q[1];
rz(-pi) q[2];
rz(2.4843744) q[3];
sx q[3];
rz(-2.2003678) q[3];
sx q[3];
rz(-2.8384125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(-0.50660261) q[2];
rz(1.5233585) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27187207) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(3.077935) q[0];
rz(1.9438538) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(-2.1228085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2301529) q[0];
sx q[0];
rz(-0.24556118) q[0];
sx q[0];
rz(2.5018238) q[0];
rz(-pi) q[1];
rz(1.8070142) q[2];
sx q[2];
rz(-1.902632) q[2];
sx q[2];
rz(2.163909) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0662811) q[1];
sx q[1];
rz(-2.3763658) q[1];
sx q[1];
rz(2.736453) q[1];
x q[2];
rz(2.3237259) q[3];
sx q[3];
rz(-1.4730771) q[3];
sx q[3];
rz(-2.1678501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7055052) q[2];
sx q[2];
rz(-0.59385308) q[2];
sx q[2];
rz(0.84211055) q[2];
rz(-2.0841133) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24404003) q[0];
sx q[0];
rz(-1.0201447) q[0];
sx q[0];
rz(-1.1919588) q[0];
rz(0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.8398197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90067139) q[0];
sx q[0];
rz(-2.0427454) q[0];
sx q[0];
rz(2.1941585) q[0];
rz(-pi) q[1];
rz(-0.19201943) q[2];
sx q[2];
rz(-0.20279591) q[2];
sx q[2];
rz(-1.6209728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8346602) q[1];
sx q[1];
rz(-1.3921157) q[1];
sx q[1];
rz(0.4731005) q[1];
rz(-0.93940063) q[3];
sx q[3];
rz(-2.7223848) q[3];
sx q[3];
rz(-1.7260176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7630345) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(-3.0540826) q[2];
rz(-1.8267501) q[3];
sx q[3];
rz(-1.0564691) q[3];
sx q[3];
rz(0.99353138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73115504) q[0];
sx q[0];
rz(-1.1719828) q[0];
sx q[0];
rz(0.80899578) q[0];
rz(-2.1058829) q[1];
sx q[1];
rz(-0.46329841) q[1];
sx q[1];
rz(2.1083924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22543487) q[0];
sx q[0];
rz(-1.6098798) q[0];
sx q[0];
rz(0.90644159) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90649267) q[2];
sx q[2];
rz(-0.64721738) q[2];
sx q[2];
rz(-0.83409079) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3174008) q[1];
sx q[1];
rz(-0.5017952) q[1];
sx q[1];
rz(2.7609894) q[1];
rz(-pi) q[2];
rz(-2.4304588) q[3];
sx q[3];
rz(-2.4681849) q[3];
sx q[3];
rz(-2.8720958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0741299) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(-1.6039675) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-2.0310903) q[3];
sx q[3];
rz(-1.0874776) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23987016) q[0];
sx q[0];
rz(-0.19342315) q[0];
sx q[0];
rz(-1.1892009) q[0];
rz(-2.4173648) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(1.6667746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132202) q[0];
sx q[0];
rz(-2.459754) q[0];
sx q[0];
rz(-2.4260055) q[0];
rz(-1.7314265) q[2];
sx q[2];
rz(-2.361627) q[2];
sx q[2];
rz(0.51973625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96383038) q[1];
sx q[1];
rz(-2.6796067) q[1];
sx q[1];
rz(-2.3038008) q[1];
rz(-pi) q[2];
rz(1.0737562) q[3];
sx q[3];
rz(-1.1599564) q[3];
sx q[3];
rz(-0.88676329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38179794) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(2.6749715) q[2];
rz(-1.7032547) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(-0.94021016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4284215) q[0];
sx q[0];
rz(-0.82601014) q[0];
sx q[0];
rz(3.1312422) q[0];
rz(-1.6185282) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(-1.3759605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3288157) q[0];
sx q[0];
rz(-1.8121908) q[0];
sx q[0];
rz(-1.9184979) q[0];
rz(-pi) q[1];
rz(2.9508123) q[2];
sx q[2];
rz(-1.4425124) q[2];
sx q[2];
rz(3.0012263) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37782323) q[1];
sx q[1];
rz(-1.0727912) q[1];
sx q[1];
rz(-1.777202) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3924005) q[3];
sx q[3];
rz(-2.0741047) q[3];
sx q[3];
rz(-2.1319413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7213664) q[2];
sx q[2];
rz(-1.0488291) q[2];
sx q[2];
rz(-2.4143207) q[2];
rz(0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(-1.4179199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(-2.7346101) q[0];
rz(2.175323) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(0.13872096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4086558) q[0];
sx q[0];
rz(-1.6870903) q[0];
sx q[0];
rz(0.035712042) q[0];
rz(-pi) q[1];
rz(-0.27757711) q[2];
sx q[2];
rz(-1.3581561) q[2];
sx q[2];
rz(0.18924282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2695864) q[1];
sx q[1];
rz(-1.2254189) q[1];
sx q[1];
rz(2.7589382) q[1];
rz(-0.14782186) q[3];
sx q[3];
rz(-1.9114219) q[3];
sx q[3];
rz(-0.56440777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8653284) q[2];
sx q[2];
rz(-1.3778957) q[2];
sx q[2];
rz(0.31065568) q[2];
rz(1.3079414) q[3];
sx q[3];
rz(-1.6303948) q[3];
sx q[3];
rz(0.37180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(0.088223591) q[0];
rz(-0.010559646) q[1];
sx q[1];
rz(-1.6190395) q[1];
sx q[1];
rz(2.6710076) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751401) q[0];
sx q[0];
rz(-1.0658403) q[0];
sx q[0];
rz(0.83374896) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3345474) q[2];
sx q[2];
rz(-1.4331487) q[2];
sx q[2];
rz(-1.1040083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6424017) q[1];
sx q[1];
rz(-2.8676053) q[1];
sx q[1];
rz(1.8887331) q[1];
x q[2];
rz(1.9235189) q[3];
sx q[3];
rz(-1.7021057) q[3];
sx q[3];
rz(-2.2985947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9528902) q[2];
sx q[2];
rz(-0.95658335) q[2];
sx q[2];
rz(-1.0799705) q[2];
rz(0.88932577) q[3];
sx q[3];
rz(-1.5417128) q[3];
sx q[3];
rz(-1.5842452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476986) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-2.3241924) q[0];
rz(2.3951702) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(0.93961632) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22720756) q[0];
sx q[0];
rz(-0.017134754) q[0];
sx q[0];
rz(-1.5428154) q[0];
x q[1];
rz(-1.290326) q[2];
sx q[2];
rz(-1.913841) q[2];
sx q[2];
rz(1.5743992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0683953) q[1];
sx q[1];
rz(-1.4540352) q[1];
sx q[1];
rz(-1.319127) q[1];
x q[2];
rz(-1.9838748) q[3];
sx q[3];
rz(-1.373817) q[3];
sx q[3];
rz(-0.22991163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8115936) q[2];
sx q[2];
rz(-0.96572319) q[2];
sx q[2];
rz(-2.7962371) q[2];
rz(-1.5251478) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(0.42720544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3578607) q[0];
sx q[0];
rz(-1.9851728) q[0];
sx q[0];
rz(-0.082948908) q[0];
rz(-2.9792765) q[1];
sx q[1];
rz(-1.4444618) q[1];
sx q[1];
rz(0.69581318) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9385949) q[0];
sx q[0];
rz(-1.2913398) q[0];
sx q[0];
rz(1.5273111) q[0];
x q[1];
rz(-1.5478304) q[2];
sx q[2];
rz(-0.12929475) q[2];
sx q[2];
rz(2.8064319) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79244991) q[1];
sx q[1];
rz(-2.6969243) q[1];
sx q[1];
rz(0.044388958) q[1];
x q[2];
rz(1.2862465) q[3];
sx q[3];
rz(-0.88594809) q[3];
sx q[3];
rz(1.7189792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42056981) q[2];
sx q[2];
rz(-2.9837954) q[2];
sx q[2];
rz(-1.2201307) q[2];
rz(-0.57341352) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(-2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2777916) q[0];
sx q[0];
rz(-1.6391123) q[0];
sx q[0];
rz(3.0085051) q[0];
rz(3.1387023) q[1];
sx q[1];
rz(-2.7255701) q[1];
sx q[1];
rz(2.4400673) q[1];
rz(1.3221424) q[2];
sx q[2];
rz(-1.7380309) q[2];
sx q[2];
rz(-2.0792014) q[2];
rz(0.057249476) q[3];
sx q[3];
rz(-0.80116873) q[3];
sx q[3];
rz(2.7309668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
