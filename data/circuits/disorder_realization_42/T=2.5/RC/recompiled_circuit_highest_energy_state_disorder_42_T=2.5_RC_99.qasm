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
rz(2.0514367) q[0];
sx q[0];
rz(-3.0213254) q[0];
sx q[0];
rz(0.076844849) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(0.82077789) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92511237) q[0];
sx q[0];
rz(-2.6651933) q[0];
sx q[0];
rz(-0.16615482) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9772867) q[2];
sx q[2];
rz(-2.2191357) q[2];
sx q[2];
rz(2.458218) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8153238) q[1];
sx q[1];
rz(-2.2827049) q[1];
sx q[1];
rz(0.91628847) q[1];
x q[2];
rz(-2.5968339) q[3];
sx q[3];
rz(-1.7228456) q[3];
sx q[3];
rz(-1.6800547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0277752) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(1.5296193) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.8703929) q[3];
sx q[3];
rz(2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40365264) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(-1.2427484) q[0];
rz(0.87951648) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.3290728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2220228) q[0];
sx q[0];
rz(-1.1276414) q[0];
sx q[0];
rz(-2.485347) q[0];
rz(-1.0813367) q[2];
sx q[2];
rz(-1.1266409) q[2];
sx q[2];
rz(1.7760488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9936867) q[1];
sx q[1];
rz(-2.4588486) q[1];
sx q[1];
rz(2.3629689) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19402282) q[3];
sx q[3];
rz(-0.7027773) q[3];
sx q[3];
rz(1.8234092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(-1.2844757) q[2];
rz(2.4744611) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(0.66481203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2892147) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(1.1997696) q[0];
rz(-1.86357) q[1];
sx q[1];
rz(-0.27690241) q[1];
sx q[1];
rz(-2.3866381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2374868) q[0];
sx q[0];
rz(-1.5747935) q[0];
sx q[0];
rz(-1.1140175) q[0];
rz(-pi) q[1];
rz(-1.4226025) q[2];
sx q[2];
rz(-0.8575927) q[2];
sx q[2];
rz(-0.36019614) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36981407) q[1];
sx q[1];
rz(-1.5358471) q[1];
sx q[1];
rz(0.82280226) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1289914) q[3];
sx q[3];
rz(-0.95944384) q[3];
sx q[3];
rz(0.064289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7167012) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.971604) q[2];
rz(-2.8814426) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(-2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78885704) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(0.66057551) q[0];
rz(0.13064101) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(-2.7052243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35231464) q[0];
sx q[0];
rz(-1.3320062) q[0];
sx q[0];
rz(2.0859857) q[0];
rz(2.5549468) q[2];
sx q[2];
rz(-1.9244951) q[2];
sx q[2];
rz(0.92087899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1088686) q[1];
sx q[1];
rz(-1.0671664) q[1];
sx q[1];
rz(-2.7313066) q[1];
rz(-pi) q[2];
x q[2];
rz(2.800501) q[3];
sx q[3];
rz(-1.6636831) q[3];
sx q[3];
rz(-2.9117381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.067430647) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(1.1379498) q[2];
rz(0.84749311) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.4165037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013019) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(-2.9345595) q[0];
rz(1.3218468) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(-2.102899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90090758) q[0];
sx q[0];
rz(-2.0863914) q[0];
sx q[0];
rz(-0.90499608) q[0];
rz(-pi) q[1];
rz(-2.8530099) q[2];
sx q[2];
rz(-2.5009857) q[2];
sx q[2];
rz(0.45329061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0463883) q[1];
sx q[1];
rz(-2.5051053) q[1];
sx q[1];
rz(0.028179006) q[1];
rz(2.8106899) q[3];
sx q[3];
rz(-2.8185113) q[3];
sx q[3];
rz(2.4289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1553104) q[2];
sx q[2];
rz(-1.0460514) q[2];
sx q[2];
rz(0.60688621) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(1.5978395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63810054) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(2.959429) q[0];
rz(-0.69193524) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(-1.8686132) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48995724) q[0];
sx q[0];
rz(-3.0916442) q[0];
sx q[0];
rz(-1.7836216) q[0];
x q[1];
rz(-1.4251377) q[2];
sx q[2];
rz(-1.7985059) q[2];
sx q[2];
rz(-2.6109254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0611736) q[1];
sx q[1];
rz(-2.3234832) q[1];
sx q[1];
rz(1.2123176) q[1];
rz(-pi) q[2];
rz(1.6342513) q[3];
sx q[3];
rz(-1.9279908) q[3];
sx q[3];
rz(0.82596367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(2.2656061) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(2.7178154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1856325) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(-0.53674269) q[0];
rz(-0.56450379) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(1.5980665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7537345) q[0];
sx q[0];
rz(-1.6711387) q[0];
sx q[0];
rz(-0.31503408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.364643) q[2];
sx q[2];
rz(-1.9788303) q[2];
sx q[2];
rz(-1.598151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98087039) q[1];
sx q[1];
rz(-2.0288784) q[1];
sx q[1];
rz(1.0700109) q[1];
x q[2];
rz(-2.1847637) q[3];
sx q[3];
rz(-1.162863) q[3];
sx q[3];
rz(1.3077426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(-1.3685263) q[2];
rz(2.8241099) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(-0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19318652) q[0];
sx q[0];
rz(-0.53006154) q[0];
sx q[0];
rz(-2.1004706) q[0];
rz(-2.5913473) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(-0.1951898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379162) q[0];
sx q[0];
rz(-1.5171255) q[0];
sx q[0];
rz(-0.1366833) q[0];
rz(-1.3615029) q[2];
sx q[2];
rz(-0.49041623) q[2];
sx q[2];
rz(1.121466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7239388) q[1];
sx q[1];
rz(-0.28755763) q[1];
sx q[1];
rz(-1.4940667) q[1];
x q[2];
rz(0.45359277) q[3];
sx q[3];
rz(-2.2633584) q[3];
sx q[3];
rz(2.3627687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0027085) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(-1.2646593) q[2];
rz(1.5876611) q[3];
sx q[3];
rz(-1.2041661) q[3];
sx q[3];
rz(1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3569262) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(-1.4534072) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(2.2834987) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88622296) q[0];
sx q[0];
rz(-2.4165034) q[0];
sx q[0];
rz(-1.0839456) q[0];
x q[1];
rz(-3.0442601) q[2];
sx q[2];
rz(-1.6989577) q[2];
sx q[2];
rz(-2.6405957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3303285) q[1];
sx q[1];
rz(-1.9703881) q[1];
sx q[1];
rz(0.43021221) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4205877) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66576362) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(-2.2541239) q[3];
sx q[3];
rz(-1.9991425) q[3];
sx q[3];
rz(-2.7953776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733521) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(-0.52262753) q[1];
sx q[1];
rz(-2.3326383) q[1];
sx q[1];
rz(-1.685198) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53082836) q[0];
sx q[0];
rz(-1.7658773) q[0];
sx q[0];
rz(-0.06819703) q[0];
rz(-1.3550884) q[2];
sx q[2];
rz(-2.459639) q[2];
sx q[2];
rz(-0.29302675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27148068) q[1];
sx q[1];
rz(-0.59529724) q[1];
sx q[1];
rz(2.4820296) q[1];
rz(-pi) q[2];
rz(-1.9144451) q[3];
sx q[3];
rz(-1.081774) q[3];
sx q[3];
rz(-2.1320313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3978079) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(1.816681) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21541967) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(2.678395) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(2.4972476) q[2];
sx q[2];
rz(-1.8229137) q[2];
sx q[2];
rz(0.075621123) q[2];
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
