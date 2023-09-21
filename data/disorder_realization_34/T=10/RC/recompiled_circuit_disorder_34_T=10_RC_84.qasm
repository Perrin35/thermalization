OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(4.4250017) q[1];
sx q[1];
rz(11.783574) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047526) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(0.063637861) q[0];
x q[1];
rz(1.9185669) q[2];
sx q[2];
rz(-2.7070621) q[2];
sx q[2];
rz(-0.44473106) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.767747) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(0.14111622) q[1];
x q[2];
rz(-0.31580117) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5995021) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(-2.9462573) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(-2.9017752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12407423) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(2.571884) q[0];
rz(-pi) q[1];
rz(2.1755978) q[2];
sx q[2];
rz(-0.36075337) q[2];
sx q[2];
rz(-0.15168562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61566831) q[1];
sx q[1];
rz(-2.3500372) q[1];
sx q[1];
rz(-0.073992373) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1274356) q[3];
sx q[3];
rz(-1.5544958) q[3];
sx q[3];
rz(-1.0316767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.3247066) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.4734316) q[0];
rz(-1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(1.0659165) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695755) q[0];
sx q[0];
rz(-1.6117275) q[0];
sx q[0];
rz(1.3224365) q[0];
rz(2.9124444) q[2];
sx q[2];
rz(-0.56729588) q[2];
sx q[2];
rz(1.0457525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2324893) q[1];
sx q[1];
rz(-2.4965582) q[1];
sx q[1];
rz(2.049936) q[1];
rz(-pi) q[2];
rz(1.4170253) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(0.94661498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(2.7815946) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.4867841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33071163) q[0];
sx q[0];
rz(-1.2905621) q[0];
sx q[0];
rz(-0.34513721) q[0];
rz(-0.60947946) q[2];
sx q[2];
rz(-1.2820499) q[2];
sx q[2];
rz(-2.2878873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1176016) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(0.28848044) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82641043) q[3];
sx q[3];
rz(-1.5757757) q[3];
sx q[3];
rz(1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(-1.460357) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(0.016074093) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094729312) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(-0.88386436) q[0];
rz(-pi) q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(-0.24239937) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2542418) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(-0.90708797) q[1];
x q[2];
rz(3.1297242) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2621433) q[0];
sx q[0];
rz(-0.95603846) q[0];
sx q[0];
rz(-0.26457796) q[0];
x q[1];
rz(0.20829006) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(1.2693894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.280937) q[1];
sx q[1];
rz(-2.5334364) q[1];
sx q[1];
rz(1.6673191) q[1];
x q[2];
rz(-2.740432) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(1.3672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(-1.7512158) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.128292) q[0];
sx q[0];
rz(-1.6830781) q[0];
sx q[0];
rz(-1.1178455) q[0];
rz(-1.6875661) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(1.5955398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3157981) q[1];
sx q[1];
rz(-2.1343532) q[1];
sx q[1];
rz(3.0347996) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4768836) q[3];
sx q[3];
rz(-2.4086773) q[3];
sx q[3];
rz(-2.7914999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-1.4027275) q[1];
sx q[1];
rz(-0.8917121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69670024) q[0];
sx q[0];
rz(-2.3377044) q[0];
sx q[0];
rz(-2.9219887) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0217998) q[2];
sx q[2];
rz(-1.0378671) q[2];
sx q[2];
rz(-1.230513) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29854952) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(-2.4119853) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9705087) q[3];
sx q[3];
rz(-0.69972316) q[3];
sx q[3];
rz(1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4760251) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.3649712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9997172) q[0];
sx q[0];
rz(-1.6122072) q[0];
sx q[0];
rz(2.3977604) q[0];
rz(-1.1599837) q[2];
sx q[2];
rz(-1.4940726) q[2];
sx q[2];
rz(3.0551747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.939146) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(2.0161611) q[1];
rz(-2.1544103) q[3];
sx q[3];
rz(-2.7899788) q[3];
sx q[3];
rz(-2.3021063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(0.96819425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4614842) q[0];
sx q[0];
rz(-2.619954) q[0];
sx q[0];
rz(-2.3011544) q[0];
rz(-2.7231611) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(-0.07428169) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.89322461) q[1];
sx q[1];
rz(-1.7176504) q[1];
sx q[1];
rz(-0.3235154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5714985) q[3];
sx q[3];
rz(-2.8194504) q[3];
sx q[3];
rz(0.57442564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(2.8075519) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
rz(1.0352186) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
