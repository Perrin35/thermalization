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
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62277943) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.8151087) q[0];
rz(-pi) q[1];
rz(-2.9847203) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(3.0770609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0040972) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(2.0069684) q[1];
x q[2];
rz(2.8257915) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(-2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-2.9462573) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(-0.23981747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0175184) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(-0.56970861) q[0];
rz(1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(2.2965455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.61566831) q[1];
sx q[1];
rz(-2.3500372) q[1];
sx q[1];
rz(3.0676003) q[1];
rz(-pi) q[2];
rz(-1.5870985) q[3];
sx q[3];
rz(-1.5849515) q[3];
sx q[3];
rz(2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720172) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(1.3224365) q[0];
x q[1];
rz(-2.9124444) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(-2.0958401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4855811) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-pi) q[2];
rz(0.12001868) q[3];
sx q[3];
rz(-0.91306049) q[3];
sx q[3];
rz(-2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(3.1090453) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(-1.4867841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958172) q[0];
sx q[0];
rz(-2.7005807) q[0];
sx q[0];
rz(0.70489024) q[0];
rz(-pi) q[1];
rz(0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(2.037231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6390037) q[1];
sx q[1];
rz(-1.8559615) q[1];
sx q[1];
rz(1.7267978) q[1];
rz(-pi) q[2];
rz(-0.82641043) q[3];
sx q[3];
rz(-1.5757757) q[3];
sx q[3];
rz(-1.5116364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.0356045) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.460357) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(0.016074093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0468633) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(2.2577283) q[0];
x q[1];
rz(1.5526505) q[2];
sx q[2];
rz(-0.43076736) q[2];
sx q[2];
rz(0.22242966) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9863723) q[1];
sx q[1];
rz(-1.1481789) q[1];
sx q[1];
rz(-0.95814725) q[1];
x q[2];
rz(-0.011868422) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2621433) q[0];
sx q[0];
rz(-2.1855542) q[0];
sx q[0];
rz(2.8770147) q[0];
rz(-0.20829006) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(-1.2693894) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21056255) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(-2.1767666) q[1];
x q[2];
rz(1.6812427) q[3];
sx q[3];
rz(-1.1718281) q[3];
sx q[3];
rz(-0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(0.21329221) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(2.4408128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.128292) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-2.0237472) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8009311) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(-1.9549191) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.339387) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(-2.1369364) q[1];
x q[2];
rz(0.84007646) q[3];
sx q[3];
rz(-1.5080161) q[3];
sx q[3];
rz(-1.1508133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(2.7409592) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-2.2498806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.007908) q[0];
sx q[0];
rz(-2.3500751) q[0];
sx q[0];
rz(-1.7931116) q[0];
rz(-2.5614369) q[2];
sx q[2];
rz(-1.1859425) q[2];
sx q[2];
rz(-2.5600524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29854952) q[1];
sx q[1];
rz(-1.9798568) q[1];
sx q[1];
rz(2.4119853) q[1];
x q[2];
rz(-2.8250152) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(-0.77565912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6531758) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(2.908356) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418755) q[0];
sx q[0];
rz(-1.6122072) q[0];
sx q[0];
rz(2.3977604) q[0];
rz(-pi) q[1];
rz(1.3806254) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(-1.4830358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0627928) q[1];
sx q[1];
rz(-1.9140869) q[1];
sx q[1];
rz(-0.72453665) q[1];
rz(-pi) q[2];
rz(2.1544103) q[3];
sx q[3];
rz(-2.7899788) q[3];
sx q[3];
rz(-0.83948638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23046872) q[0];
sx q[0];
rz(-1.231912) q[0];
sx q[0];
rz(-1.1662657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41843157) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(3.067311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4149975) q[1];
sx q[1];
rz(-1.2508878) q[1];
sx q[1];
rz(1.4160316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.00023437436) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(0.57516592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0782464) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.39682) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-1.2583854) q[2];
sx q[2];
rz(-0.84597833) q[2];
sx q[2];
rz(-0.99920263) q[2];
rz(-1.0352186) q[3];
sx q[3];
rz(-0.62789161) q[3];
sx q[3];
rz(-0.37554489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
