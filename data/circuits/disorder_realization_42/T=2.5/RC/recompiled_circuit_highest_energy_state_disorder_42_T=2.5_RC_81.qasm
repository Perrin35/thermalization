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
rz(-1.090156) q[0];
sx q[0];
rz(-0.12026726) q[0];
sx q[0];
rz(-0.076844849) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(0.82077789) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3479697) q[0];
sx q[0];
rz(-1.6467148) q[0];
sx q[0];
rz(0.47077076) q[0];
rz(2.2256741) q[2];
sx q[2];
rz(-1.7015463) q[2];
sx q[2];
rz(-0.98721013) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3617435) q[1];
sx q[1];
rz(-2.0498154) q[1];
sx q[1];
rz(-2.3142525) q[1];
x q[2];
rz(1.7480811) q[3];
sx q[3];
rz(-1.0330135) q[3];
sx q[3];
rz(3.1238525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0277752) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.5296193) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(-2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40365264) q[0];
sx q[0];
rz(-2.0495179) q[0];
sx q[0];
rz(-1.2427484) q[0];
rz(-0.87951648) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.8125199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14318109) q[0];
sx q[0];
rz(-0.77309009) q[0];
sx q[0];
rz(-0.66112681) q[0];
rz(-1.0813367) q[2];
sx q[2];
rz(-1.1266409) q[2];
sx q[2];
rz(1.7760488) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.052292) q[1];
sx q[1];
rz(-2.0366021) q[1];
sx q[1];
rz(-2.0897081) q[1];
rz(2.9475698) q[3];
sx q[3];
rz(-0.7027773) q[3];
sx q[3];
rz(1.3181835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30895147) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(-1.857117) q[2];
rz(0.66713157) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(-0.66481203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(1.941823) q[0];
rz(1.2780227) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(2.3866381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66865488) q[0];
sx q[0];
rz(-1.1140214) q[0];
sx q[0];
rz(0.0044538023) q[0];
rz(-1.4226025) q[2];
sx q[2];
rz(-0.8575927) q[2];
sx q[2];
rz(2.7813965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9082005) q[1];
sx q[1];
rz(-2.3182239) q[1];
sx q[1];
rz(0.047659454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69061227) q[3];
sx q[3];
rz(-2.0194144) q[3];
sx q[3];
rz(-1.2909362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42489147) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.1699886) q[2];
rz(-0.26015002) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(-2.9366176) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3527356) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(2.4810171) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-0.6159997) q[1];
sx q[1];
rz(2.7052243) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0562387) q[0];
sx q[0];
rz(-1.071601) q[0];
sx q[0];
rz(-0.27277314) q[0];
rz(-1.9881093) q[2];
sx q[2];
rz(-1.0247604) q[2];
sx q[2];
rz(-0.87620902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.032724) q[1];
sx q[1];
rz(-1.0671664) q[1];
sx q[1];
rz(2.7313066) q[1];
rz(-1.6693258) q[3];
sx q[3];
rz(-1.2312345) q[3];
sx q[3];
rz(-1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.074162) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(-1.1379498) q[2];
rz(0.84749311) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(-2.9345595) q[0];
rz(-1.3218468) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(2.102899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2406851) q[0];
sx q[0];
rz(-1.0552013) q[0];
sx q[0];
rz(2.2365966) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8530099) q[2];
sx q[2];
rz(-2.5009857) q[2];
sx q[2];
rz(-2.688302) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0814235) q[1];
sx q[1];
rz(-0.93460235) q[1];
sx q[1];
rz(-1.549975) q[1];
rz(-2.8106899) q[3];
sx q[3];
rz(-0.32308137) q[3];
sx q[3];
rz(2.4289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1553104) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(2.5347064) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(1.5978395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5034921) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(0.18216369) q[0];
rz(-2.4496574) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(1.8686132) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385524) q[0];
sx q[0];
rz(-1.6196169) q[0];
sx q[0];
rz(-3.131034) q[0];
rz(-2.9115306) q[2];
sx q[2];
rz(-1.7126691) q[2];
sx q[2];
rz(-2.134568) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0804191) q[1];
sx q[1];
rz(-0.81810942) q[1];
sx q[1];
rz(-1.929275) q[1];
x q[2];
rz(0.16830651) q[3];
sx q[3];
rz(-0.3625488) q[3];
sx q[3];
rz(-2.1358638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(-0.87598652) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(-2.7178154) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9559602) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(0.56450379) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(1.5435262) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9912883) q[0];
sx q[0];
rz(-1.8841916) q[0];
sx q[0];
rz(-1.4652976) q[0];
x q[1];
rz(-2.6993611) q[2];
sx q[2];
rz(-2.6870607) q[2];
sx q[2];
rz(1.0584128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7891414) q[1];
sx q[1];
rz(-1.1256212) q[1];
sx q[1];
rz(-2.6295202) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6552175) q[3];
sx q[3];
rz(-1.0135627) q[3];
sx q[3];
rz(3.1319905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(-1.7730664) q[2];
rz(0.3174828) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484061) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(1.0411221) q[0];
rz(-2.5913473) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(-0.1951898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.081851) q[0];
sx q[0];
rz(-1.7072816) q[0];
sx q[0];
rz(1.5166212) q[0];
rz(-0.11048079) q[2];
sx q[2];
rz(-2.049597) q[2];
sx q[2];
rz(-0.88518054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0620408) q[1];
sx q[1];
rz(-1.5490546) q[1];
sx q[1];
rz(1.2840391) q[1];
rz(-1.0848673) q[3];
sx q[3];
rz(-0.80683364) q[3];
sx q[3];
rz(3.0148447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.8769334) q[2];
rz(1.5876611) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(1.9896265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3569262) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(-2.7333562) q[0];
rz(1.6881855) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(-0.85809392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6397809) q[0];
sx q[0];
rz(-2.1970891) q[0];
sx q[0];
rz(-2.7485952) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0442601) q[2];
sx q[2];
rz(-1.6989577) q[2];
sx q[2];
rz(-0.50099696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8112642) q[1];
sx q[1];
rz(-1.1712046) q[1];
sx q[1];
rz(-2.7113804) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87808064) q[3];
sx q[3];
rz(-1.47444) q[3];
sx q[3];
rz(0.89123301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.475829) q[2];
sx q[2];
rz(-2.2593311) q[2];
sx q[2];
rz(-0.25064722) q[2];
rz(-2.2541239) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(2.7953776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733521) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(-2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(-1.685198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9495372) q[0];
sx q[0];
rz(-2.9350781) q[0];
sx q[0];
rz(1.9028936) q[0];
x q[1];
rz(0.9002879) q[2];
sx q[2];
rz(-1.4354726) q[2];
sx q[2];
rz(1.4462665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8700614) q[1];
sx q[1];
rz(-1.2200331) q[1];
sx q[1];
rz(0.49141541) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5771477) q[3];
sx q[3];
rz(-0.58957499) q[3];
sx q[3];
rz(1.6605476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7437848) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(1.3249116) q[2];
rz(-1.2815453) q[3];
sx q[3];
rz(-0.89373389) q[3];
sx q[3];
rz(-2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21541967) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(-0.46319766) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(-1.2590969) q[2];
sx q[2];
rz(-2.1915956) q[2];
sx q[2];
rz(-1.6804463) q[2];
rz(2.9897965) q[3];
sx q[3];
rz(-0.67852466) q[3];
sx q[3];
rz(-1.8192374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
