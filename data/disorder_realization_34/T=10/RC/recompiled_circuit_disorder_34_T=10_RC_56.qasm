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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570219) q[0];
sx q[0];
rz(-1.5091358) q[0];
sx q[0];
rz(1.3205359) q[0];
rz(1.9185669) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(0.44473106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.767747) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(-3.0004764) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31580117) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-2.9017752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0175184) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(-0.56970861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2698783) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(-0.84504715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51057306) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(-1.6454979) q[1];
x q[2];
rz(-1.5544942) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-0.53888884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(-1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-2.0756762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.802793) q[0];
sx q[0];
rz(-2.8899513) q[0];
sx q[0];
rz(1.405707) q[0];
rz(-1.4270646) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(-1.3155754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90910339) q[1];
sx q[1];
rz(-0.64503446) q[1];
sx q[1];
rz(-1.0916566) q[1];
rz(-pi) q[2];
rz(-3.021574) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(-0.79157296) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2457755) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(-0.70489024) q[0];
rz(-1.9183667) q[2];
sx q[2];
rz(-2.1516557) q[2];
sx q[2];
rz(2.6207404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0123803) q[1];
sx q[1];
rz(-0.32402363) q[1];
sx q[1];
rz(-2.654241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82641043) q[3];
sx q[3];
rz(-1.5757757) q[3];
sx q[3];
rz(-1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(0.016074093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094729312) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(2.2577283) q[0];
rz(-3.1332544) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(2.8991933) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(0.95814725) q[1];
x q[2];
rz(-0.36230476) q[3];
sx q[3];
rz(-1.5665896) q[3];
sx q[3];
rz(1.4570953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428403) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.2247359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3186) q[0];
sx q[0];
rz(-0.66245125) q[0];
sx q[0];
rz(1.2160765) q[0];
rz(1.9903509) q[2];
sx q[2];
rz(-0.47861368) q[2];
sx q[2];
rz(-1.4065557) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9310301) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(-0.964826) q[1];
rz(2.740432) q[3];
sx q[3];
rz(-1.6725371) q[3];
sx q[3];
rz(1.7743558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619103) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535671) q[0];
sx q[0];
rz(-2.0206847) q[0];
sx q[0];
rz(3.0168424) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4540265) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(-1.5955398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80220561) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(1.0046563) q[1];
rz(-pi) q[2];
rz(-1.4768836) q[3];
sx q[3];
rz(-2.4086773) q[3];
sx q[3];
rz(-2.7914999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(-1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-0.8917121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72043334) q[0];
sx q[0];
rz(-1.7283068) q[0];
sx q[0];
rz(-2.349855) q[0];
rz(-pi) q[1];
rz(2.5614369) q[2];
sx q[2];
rz(-1.1859425) q[2];
sx q[2];
rz(2.5600524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5276287) q[1];
sx q[1];
rz(-2.2288537) q[1];
sx q[1];
rz(-1.0440473) q[1];
rz(1.9705087) q[3];
sx q[3];
rz(-0.69972316) q[3];
sx q[3];
rz(1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(1.7766215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997172) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(2.3977604) q[0];
x q[1];
rz(1.7609673) q[2];
sx q[2];
rz(-0.41751465) q[2];
sx q[2];
rz(-1.4830358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.078799876) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(2.417056) q[1];
rz(-2.1544103) q[3];
sx q[3];
rz(-2.7899788) q[3];
sx q[3];
rz(0.83948638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(-0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.0868527) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614842) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(-2.3011544) q[0];
x q[1];
rz(-0.41843157) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(-3.067311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4149975) q[1];
sx q[1];
rz(-1.8907049) q[1];
sx q[1];
rz(-1.725561) q[1];
x q[2];
rz(3.1413583) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(-2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(1.1408268) q[3];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-0.36623476) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(0.74950779) q[2];
sx q[2];
rz(-1.8029677) q[2];
sx q[2];
rz(0.36063902) q[2];
rz(-1.0127388) q[3];
sx q[3];
rz(-1.2663208) q[3];
sx q[3];
rz(-2.3940621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
