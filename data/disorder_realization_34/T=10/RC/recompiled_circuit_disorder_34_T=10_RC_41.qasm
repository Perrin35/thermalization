OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62277943) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(-1.326484) q[0];
x q[1];
rz(-2.9847203) q[2];
sx q[2];
rz(-1.9777159) q[2];
sx q[2];
rz(-3.0770609) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4408735) q[1];
sx q[1];
rz(-0.45342017) q[1];
sx q[1];
rz(1.2749626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8257915) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(-2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(1.5995021) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(0.19533531) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(0.23981747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0175184) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(0.56970861) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2113091) q[2];
sx q[2];
rz(-1.8654056) q[2];
sx q[2];
rz(0.78795563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51057306) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(-1.6454979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5544942) q[3];
sx q[3];
rz(-1.5849515) q[3];
sx q[3];
rz(2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(1.0659165) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33879963) q[0];
sx q[0];
rz(-0.25164139) q[0];
sx q[0];
rz(1.7358857) q[0];
x q[1];
rz(2.9124444) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(2.0958401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8733858) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-1.4170253) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-0.032547396) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(-1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33071163) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(0.34513721) q[0];
rz(1.2232259) q[2];
sx q[2];
rz(-2.1516557) q[2];
sx q[2];
rz(2.6207404) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.023991) q[1];
sx q[1];
rz(-1.4211434) q[1];
sx q[1];
rz(-2.8531122) q[1];
rz(-1.5781457) q[3];
sx q[3];
rz(-0.74439936) q[3];
sx q[3];
rz(-3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-0.016074093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.21503) q[0];
sx q[0];
rz(-0.83850551) q[0];
sx q[0];
rz(2.3123884) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1400914) q[2];
sx q[2];
rz(-1.5632196) q[2];
sx q[2];
rz(-1.3318782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88735089) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(-2.2345047) q[1];
rz(-pi) q[2];
rz(0.36230476) q[3];
sx q[3];
rz(-1.5665896) q[3];
sx q[3];
rz(1.6844974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(-2.9212852) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82299267) q[0];
sx q[0];
rz(-0.66245125) q[0];
sx q[0];
rz(-1.2160765) q[0];
rz(-2.0133063) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(-2.92885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.280937) q[1];
sx q[1];
rz(-0.60815629) q[1];
sx q[1];
rz(1.4742736) q[1];
rz(-pi) q[2];
rz(2.740432) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(1.7743558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.1619103) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(0.70077983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535671) q[0];
sx q[0];
rz(-2.0206847) q[0];
sx q[0];
rz(-0.1247503) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8220196) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(0.061352913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3157981) q[1];
sx q[1];
rz(-2.1343532) q[1];
sx q[1];
rz(0.10679306) q[1];
x q[2];
rz(-0.84007646) q[3];
sx q[3];
rz(-1.6335765) q[3];
sx q[3];
rz(1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-2.2498806) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448924) q[0];
sx q[0];
rz(-2.3377044) q[0];
sx q[0];
rz(-0.21960396) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0217998) q[2];
sx q[2];
rz(-1.0378671) q[2];
sx q[2];
rz(1.9110796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8430431) q[1];
sx q[1];
rz(-1.9798568) q[1];
sx q[1];
rz(0.72960735) q[1];
rz(-pi) q[2];
rz(-2.8250152) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(-0.77565912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(2.419557) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.7766215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6745973) q[0];
sx q[0];
rz(-0.82775263) q[0];
sx q[0];
rz(1.6270431) q[0];
rz(-pi) q[1];
rz(-1.7609673) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(1.6585569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0627928) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(-2.417056) q[1];
rz(-pi) q[2];
rz(-1.2737208) q[3];
sx q[3];
rz(-1.7617412) q[3];
sx q[3];
rz(1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.9101248) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(2.058303) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(2.1733984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4817081) q[0];
sx q[0];
rz(-1.9511001) q[0];
sx q[0];
rz(2.7754521) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41843157) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(3.067311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72659513) q[1];
sx q[1];
rz(-1.8907049) q[1];
sx q[1];
rz(1.725561) q[1];
rz(0.00023437436) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-2.8075519) q[2];
sx q[2];
rz(-2.3636849) q[2];
sx q[2];
rz(-1.452527) q[2];
rz(2.7868174) q[3];
sx q[3];
rz(-1.0412024) q[3];
sx q[3];
rz(2.1333221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
