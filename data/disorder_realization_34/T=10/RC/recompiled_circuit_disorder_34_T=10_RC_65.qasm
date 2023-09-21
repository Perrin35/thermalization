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
rz(-0.43963471) q[0];
sx q[0];
rz(3.0602732) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(0.78279701) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62277943) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.326484) q[0];
x q[1];
rz(1.2230258) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(-0.44473106) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3738457) q[1];
sx q[1];
rz(-1.1384374) q[1];
sx q[1];
rz(-3.0004764) q[1];
rz(-1.6601895) q[3];
sx q[3];
rz(-1.8854183) q[3];
sx q[3];
rz(2.6320626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-2.9017752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12407423) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(2.571884) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9302836) q[2];
sx q[2];
rz(-1.8654056) q[2];
sx q[2];
rz(-0.78795563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0071734) q[1];
sx q[1];
rz(-1.6234142) q[1];
sx q[1];
rz(0.79018553) q[1];
x q[2];
rz(1.5870985) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-0.53888884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.3298539) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
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
rz(-pi) q[3];
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
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(-1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(1.0659165) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695755) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(1.3224365) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4270646) q[2];
sx q[2];
rz(-1.020069) q[2];
sx q[2];
rz(1.8260173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8733858) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-pi) q[2];
rz(-1.7245674) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(0.94661498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-0.032547396) q[2];
rz(2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.6548086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958172) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(2.4367024) q[0];
x q[1];
rz(0.60947946) q[2];
sx q[2];
rz(-1.2820499) q[2];
sx q[2];
rz(-0.8537054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.023991) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(-2.8531122) q[1];
rz(-pi) q[2];
rz(-1.5634469) q[3];
sx q[3];
rz(-0.74439936) q[3];
sx q[3];
rz(3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(0.016074093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094729312) q[0];
sx q[0];
rz(-2.097058) q[0];
sx q[0];
rz(-0.88386436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(0.24239937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6962471) q[1];
sx q[1];
rz(-2.1228585) q[1];
sx q[1];
rz(0.50260431) q[1];
rz(-2.7792879) q[3];
sx q[3];
rz(-1.575003) q[3];
sx q[3];
rz(-1.6844974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428403) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(1.2247359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82299267) q[0];
sx q[0];
rz(-2.4791414) q[0];
sx q[0];
rz(1.9255161) q[0];
rz(-pi) q[1];
rz(2.0133063) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(2.92885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.280937) q[1];
sx q[1];
rz(-0.60815629) q[1];
sx q[1];
rz(-1.6673191) q[1];
x q[2];
rz(-2.740432) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(-1.7743558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99888745) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66877767) q[0];
sx q[0];
rz(-0.46572177) q[0];
sx q[0];
rz(-1.8229683) q[0];
x q[1];
rz(2.8220196) q[2];
sx q[2];
rz(-1.6817037) q[2];
sx q[2];
rz(-3.0802397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5138445) q[1];
sx q[1];
rz(-2.5690837) q[1];
sx q[1];
rz(1.7379012) q[1];
rz(1.4768836) q[3];
sx q[3];
rz(-0.73291534) q[3];
sx q[3];
rz(-2.7914999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(3.11943) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(2.2498806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1336846) q[0];
sx q[0];
rz(-2.3500751) q[0];
sx q[0];
rz(1.7931116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1197929) q[2];
sx q[2];
rz(-2.1037256) q[2];
sx q[2];
rz(-1.9110796) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2878694) q[1];
sx q[1];
rz(-2.3239377) q[1];
sx q[1];
rz(-2.5649648) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1710839) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(1.8613929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(0.23323664) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.7766215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6745973) q[0];
sx q[0];
rz(-2.31384) q[0];
sx q[0];
rz(-1.6270431) q[0];
x q[1];
rz(0.083655595) q[2];
sx q[2];
rz(-1.1612648) q[2];
sx q[2];
rz(-1.4510029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2024467) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(-2.0161611) q[1];
x q[2];
rz(-2.1544103) q[3];
sx q[3];
rz(-0.35161388) q[3];
sx q[3];
rz(-0.83948638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.4779133) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4817081) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(0.36614059) q[0];
x q[1];
rz(2.0696938) q[2];
sx q[2];
rz(-1.1985156) q[2];
sx q[2];
rz(1.6921713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72659513) q[1];
sx q[1];
rz(-1.8907049) q[1];
sx q[1];
rz(1.725561) q[1];
x q[2];
rz(1.5700941) q[3];
sx q[3];
rz(-0.32214221) q[3];
sx q[3];
rz(0.57442564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(1.1408268) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(0.36623476) q[1];
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
