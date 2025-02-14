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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(-0.92010486) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(-0.97715598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29229627) q[0];
sx q[0];
rz(-1.1998645) q[0];
sx q[0];
rz(-1.7975259) q[0];
x q[1];
rz(0.10213587) q[2];
sx q[2];
rz(-3.0086824) q[2];
sx q[2];
rz(-1.5166211) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10763532) q[1];
sx q[1];
rz(-0.9743685) q[1];
sx q[1];
rz(-1.4075085) q[1];
x q[2];
rz(1.8771873) q[3];
sx q[3];
rz(-2.236111) q[3];
sx q[3];
rz(0.09395919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(2.5207632) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8949378) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(-2.5625693) q[0];
rz(-2.31965) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(-0.45375219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762208) q[0];
sx q[0];
rz(-2.1848618) q[0];
sx q[0];
rz(-1.8472415) q[0];
rz(1.4582107) q[2];
sx q[2];
rz(-1.1069008) q[2];
sx q[2];
rz(-2.2245906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0099185506) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(2.591955) q[1];
rz(-2.1165127) q[3];
sx q[3];
rz(-1.1835872) q[3];
sx q[3];
rz(-0.51147616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73432505) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(-0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(-1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9722209) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(0.55150223) q[0];
rz(0.36477271) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(0.50484467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802949) q[0];
sx q[0];
rz(-0.43879959) q[0];
sx q[0];
rz(-0.3831692) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7288923) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(-2.7049261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4703731) q[1];
sx q[1];
rz(-0.80802901) q[1];
sx q[1];
rz(0.68617448) q[1];
rz(2.4931413) q[3];
sx q[3];
rz(-2.5621427) q[3];
sx q[3];
rz(3.0045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(-2.1164472) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.717201) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(-1.4855509) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43760532) q[0];
sx q[0];
rz(-0.018089596) q[0];
sx q[0];
rz(-1.4996858) q[0];
rz(-0.35648326) q[2];
sx q[2];
rz(-1.0986137) q[2];
sx q[2];
rz(-3.0038578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8532456) q[1];
sx q[1];
rz(-1.6430055) q[1];
sx q[1];
rz(1.9068524) q[1];
x q[2];
rz(2.5136389) q[3];
sx q[3];
rz(-2.0139004) q[3];
sx q[3];
rz(-1.7114044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.7690313) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-0.6438846) q[3];
sx q[3];
rz(-0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32886252) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(-0.10511705) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-0.9143908) q[1];
sx q[1];
rz(2.0786659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30468291) q[0];
sx q[0];
rz(-1.5848918) q[0];
sx q[0];
rz(-2.430116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.379871) q[2];
sx q[2];
rz(-1.3889379) q[2];
sx q[2];
rz(-1.8363425) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.059710709) q[1];
sx q[1];
rz(-0.060259911) q[1];
sx q[1];
rz(-2.9032272) q[1];
x q[2];
rz(-1.5751198) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(0.81567848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(-1.3366535) q[2];
rz(-0.054232728) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9688251) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-0.85375839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19148286) q[0];
sx q[0];
rz(-1.6304509) q[0];
sx q[0];
rz(-1.7050939) q[0];
rz(-0.75550236) q[2];
sx q[2];
rz(-1.1248338) q[2];
sx q[2];
rz(-1.8309636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(2.3385919) q[1];
rz(-2.742203) q[3];
sx q[3];
rz(-0.37617427) q[3];
sx q[3];
rz(1.1272421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13009109) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(0.74756527) q[2];
rz(-2.2085564) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(-2.172442) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(2.0279121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232491) q[0];
sx q[0];
rz(-1.9750496) q[0];
sx q[0];
rz(0.1057616) q[0];
rz(-pi) q[1];
rz(2.6081354) q[2];
sx q[2];
rz(-2.1367624) q[2];
sx q[2];
rz(3.0135558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.021999849) q[1];
sx q[1];
rz(-0.55460677) q[1];
sx q[1];
rz(1.9153992) q[1];
rz(-pi) q[2];
rz(-0.05984743) q[3];
sx q[3];
rz(-0.52111292) q[3];
sx q[3];
rz(-1.9684362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(0.92195177) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7477357) q[0];
sx q[0];
rz(-1.4841635) q[0];
sx q[0];
rz(1.9282233) q[0];
x q[1];
rz(2.5002648) q[2];
sx q[2];
rz(-1.3654764) q[2];
sx q[2];
rz(-2.5032798) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53271097) q[1];
sx q[1];
rz(-1.5517762) q[1];
sx q[1];
rz(-0.089781922) q[1];
rz(1.1128759) q[3];
sx q[3];
rz(-2.1919214) q[3];
sx q[3];
rz(-0.21645138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(1.1456683) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.3409748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670973) q[0];
sx q[0];
rz(-1.4112368) q[0];
sx q[0];
rz(-1.480353) q[0];
rz(-pi) q[1];
rz(-3.0485504) q[2];
sx q[2];
rz(-0.75936717) q[2];
sx q[2];
rz(-1.134128) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0611524) q[1];
sx q[1];
rz(-2.3576479) q[1];
sx q[1];
rz(-2.817201) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.789647) q[3];
sx q[3];
rz(-2.7776981) q[3];
sx q[3];
rz(2.4879251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0535023) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(-0.11605334) q[2];
rz(-1.8807489) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(-1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5905404) q[0];
sx q[0];
rz(-0.11369471) q[0];
sx q[0];
rz(-2.9569448) q[0];
rz(-0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.437423) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4686615) q[0];
sx q[0];
rz(-1.3672921) q[0];
sx q[0];
rz(-2.0528021) q[0];
rz(-3.060861) q[2];
sx q[2];
rz(-2.417832) q[2];
sx q[2];
rz(2.2249976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.559222) q[1];
sx q[1];
rz(-0.6383903) q[1];
sx q[1];
rz(1.8965782) q[1];
rz(-pi) q[2];
rz(-0.67622306) q[3];
sx q[3];
rz(-1.3075324) q[3];
sx q[3];
rz(-0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47716466) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(2.3264558) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(-2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.7246134) q[0];
sx q[0];
rz(-1.8920349) q[0];
sx q[0];
rz(-2.98988) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(-1.2909605) q[2];
sx q[2];
rz(-1.6879514) q[2];
sx q[2];
rz(0.25456706) q[2];
rz(-1.5315957) q[3];
sx q[3];
rz(-1.3269674) q[3];
sx q[3];
rz(1.5404601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
