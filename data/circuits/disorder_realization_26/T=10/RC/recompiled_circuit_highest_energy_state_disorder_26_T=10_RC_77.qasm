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
rz(-0.24902046) q[0];
sx q[0];
rz(-0.95872107) q[0];
sx q[0];
rz(-0.22554654) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(3.5993242) q[1];
sx q[1];
rz(9.2158894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1516722) q[0];
sx q[0];
rz(-2.6197126) q[0];
sx q[0];
rz(-0.56445091) q[0];
x q[1];
rz(1.8525847) q[2];
sx q[2];
rz(-2.1827841) q[2];
sx q[2];
rz(0.91152292) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3882313) q[1];
sx q[1];
rz(-1.2665788) q[1];
sx q[1];
rz(-2.6016584) q[1];
rz(-0.31583438) q[3];
sx q[3];
rz(-1.2441934) q[3];
sx q[3];
rz(2.0706319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1170342) q[2];
sx q[2];
rz(-2.890675) q[2];
sx q[2];
rz(0.92570242) q[2];
rz(-2.3671345) q[3];
sx q[3];
rz(-1.2239933) q[3];
sx q[3];
rz(-1.2831877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16933146) q[0];
sx q[0];
rz(-0.65003482) q[0];
sx q[0];
rz(-1.4647123) q[0];
rz(-2.2526534) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(-1.3533786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.618093) q[0];
sx q[0];
rz(-0.67636469) q[0];
sx q[0];
rz(0.47250611) q[0];
rz(3.0841295) q[2];
sx q[2];
rz(-2.7749535) q[2];
sx q[2];
rz(-0.55546782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.9026437) q[1];
sx q[1];
rz(-1.3629991) q[1];
sx q[1];
rz(-2.5420339) q[1];
x q[2];
rz(-0.7781182) q[3];
sx q[3];
rz(-1.7675085) q[3];
sx q[3];
rz(2.9784378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9951524) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(0.97274441) q[2];
rz(1.815833) q[3];
sx q[3];
rz(-1.1498068) q[3];
sx q[3];
rz(0.61746922) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165118) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(-0.02956477) q[0];
rz(2.120453) q[1];
sx q[1];
rz(-2.857326) q[1];
sx q[1];
rz(-0.31612843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9122304) q[0];
sx q[0];
rz(-1.5925613) q[0];
sx q[0];
rz(0.020177186) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7720376) q[2];
sx q[2];
rz(-1.6627208) q[2];
sx q[2];
rz(-2.6215009) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21507922) q[1];
sx q[1];
rz(-1.2383019) q[1];
sx q[1];
rz(-0.18400561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.088313266) q[3];
sx q[3];
rz(-2.4188015) q[3];
sx q[3];
rz(2.4105154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11008392) q[2];
sx q[2];
rz(-1.2981334) q[2];
sx q[2];
rz(2.617344) q[2];
rz(2.5414069) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(-1.0711099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6780739) q[0];
sx q[0];
rz(-1.9215895) q[0];
sx q[0];
rz(-0.36233166) q[0];
rz(-2.433297) q[1];
sx q[1];
rz(-2.7999122) q[1];
sx q[1];
rz(-1.9963616) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6714685) q[0];
sx q[0];
rz(-0.13397476) q[0];
sx q[0];
rz(2.2132316) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6537204) q[2];
sx q[2];
rz(-0.93624253) q[2];
sx q[2];
rz(1.3417336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64001361) q[1];
sx q[1];
rz(-0.7487491) q[1];
sx q[1];
rz(-1.1925936) q[1];
x q[2];
rz(-0.16412603) q[3];
sx q[3];
rz(-1.6200674) q[3];
sx q[3];
rz(-2.6409145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5273253) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(3.1329727) q[2];
rz(-2.4873867) q[3];
sx q[3];
rz(-2.4253186) q[3];
sx q[3];
rz(-2.8692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006186) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(2.3531438) q[0];
rz(-1.0139326) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(-3.0533155) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006314) q[0];
sx q[0];
rz(-0.92871237) q[0];
sx q[0];
rz(0.93812801) q[0];
x q[1];
rz(-1.6730673) q[2];
sx q[2];
rz(-2.132093) q[2];
sx q[2];
rz(3.1401538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4406529) q[1];
sx q[1];
rz(-1.0723317) q[1];
sx q[1];
rz(0.32157125) q[1];
rz(1.621219) q[3];
sx q[3];
rz(-2.6459624) q[3];
sx q[3];
rz(-2.4448066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2225515) q[2];
sx q[2];
rz(-1.4107979) q[2];
sx q[2];
rz(-0.50773531) q[2];
rz(0.12568411) q[3];
sx q[3];
rz(-1.1832184) q[3];
sx q[3];
rz(-0.79508933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8102201) q[0];
sx q[0];
rz(-2.3851244) q[0];
sx q[0];
rz(-3.0561225) q[0];
rz(2.5147009) q[1];
sx q[1];
rz(-1.1566999) q[1];
sx q[1];
rz(1.8524106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4816718) q[0];
sx q[0];
rz(-2.916475) q[0];
sx q[0];
rz(-1.0876924) q[0];
rz(-pi) q[1];
rz(2.3979509) q[2];
sx q[2];
rz(-1.6847451) q[2];
sx q[2];
rz(-2.7211962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8823008) q[1];
sx q[1];
rz(-1.3332102) q[1];
sx q[1];
rz(2.9670694) q[1];
rz(-pi) q[2];
rz(3.1269642) q[3];
sx q[3];
rz(-1.808015) q[3];
sx q[3];
rz(0.62495172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36279303) q[2];
sx q[2];
rz(-2.1976566) q[2];
sx q[2];
rz(-1.2550521) q[2];
rz(-3.0826027) q[3];
sx q[3];
rz(-1.2041644) q[3];
sx q[3];
rz(-2.1571933) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7207709) q[0];
sx q[0];
rz(-1.2728007) q[0];
sx q[0];
rz(0.76501784) q[0];
rz(-1.7401241) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(-0.23745647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119806) q[0];
sx q[0];
rz(-1.6061826) q[0];
sx q[0];
rz(2.1982212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.075421926) q[2];
sx q[2];
rz(-0.89176501) q[2];
sx q[2];
rz(2.0184269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4478893) q[1];
sx q[1];
rz(-0.2940601) q[1];
sx q[1];
rz(0.15090461) q[1];
rz(-3.1116099) q[3];
sx q[3];
rz(-0.87453547) q[3];
sx q[3];
rz(3.1122623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6852297) q[2];
sx q[2];
rz(-0.80358973) q[2];
sx q[2];
rz(2.2881499) q[2];
rz(-1.338909) q[3];
sx q[3];
rz(-1.572255) q[3];
sx q[3];
rz(-2.9760823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.45561871) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(2.9562505) q[0];
rz(-0.39400426) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(1.0967163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029848969) q[0];
sx q[0];
rz(-1.9557448) q[0];
sx q[0];
rz(-2.0864331) q[0];
rz(-pi) q[1];
rz(1.4273066) q[2];
sx q[2];
rz(-2.4713785) q[2];
sx q[2];
rz(1.0619947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1450069) q[1];
sx q[1];
rz(-1.2357986) q[1];
sx q[1];
rz(0.62535357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1878566) q[3];
sx q[3];
rz(-2.1309867) q[3];
sx q[3];
rz(-2.0479353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3700221) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(0.30878511) q[2];
rz(-1.8874946) q[3];
sx q[3];
rz(-1.2718688) q[3];
sx q[3];
rz(-1.1312283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.18987385) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(-0.21981123) q[0];
rz(-1.9728569) q[1];
sx q[1];
rz(-1.7857779) q[1];
sx q[1];
rz(-1.8692325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74368661) q[0];
sx q[0];
rz(-0.96677654) q[0];
sx q[0];
rz(0.029725909) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6556622) q[2];
sx q[2];
rz(-1.0615292) q[2];
sx q[2];
rz(2.9743663) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1049526) q[1];
sx q[1];
rz(-2.7014241) q[1];
sx q[1];
rz(2.7326581) q[1];
rz(-pi) q[2];
rz(-1.06274) q[3];
sx q[3];
rz(-2.0174167) q[3];
sx q[3];
rz(-1.5039218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9610338) q[2];
sx q[2];
rz(-0.84332931) q[2];
sx q[2];
rz(-0.45832222) q[2];
rz(2.787163) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(-1.1752769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71427041) q[0];
sx q[0];
rz(-1.2027807) q[0];
sx q[0];
rz(-2.400193) q[0];
rz(-1.4330014) q[1];
sx q[1];
rz(-0.42867908) q[1];
sx q[1];
rz(-0.0892078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2682104) q[0];
sx q[0];
rz(-2.6307627) q[0];
sx q[0];
rz(1.9477773) q[0];
rz(-pi) q[1];
rz(1.1871502) q[2];
sx q[2];
rz(-1.0133146) q[2];
sx q[2];
rz(-1.5526183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7877823) q[1];
sx q[1];
rz(-1.5292166) q[1];
sx q[1];
rz(1.2316684) q[1];
rz(-1.1600288) q[3];
sx q[3];
rz(-2.2645519) q[3];
sx q[3];
rz(-0.18405562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1409113) q[2];
sx q[2];
rz(-0.050364308) q[2];
sx q[2];
rz(0.64765206) q[2];
rz(-0.40376136) q[3];
sx q[3];
rz(-1.8881366) q[3];
sx q[3];
rz(0.13127479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972926) q[0];
sx q[0];
rz(-2.1680752) q[0];
sx q[0];
rz(-2.0808921) q[0];
rz(-0.79509673) q[1];
sx q[1];
rz(-2.2207694) q[1];
sx q[1];
rz(2.3650852) q[1];
rz(-2.2262103) q[2];
sx q[2];
rz(-1.6554828) q[2];
sx q[2];
rz(0.2632904) q[2];
rz(1.6609876) q[3];
sx q[3];
rz(-1.492751) q[3];
sx q[3];
rz(0.94730151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
