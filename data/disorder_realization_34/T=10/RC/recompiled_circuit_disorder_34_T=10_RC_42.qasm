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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711174) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(0.063637861) q[0];
x q[1];
rz(-1.159367) q[2];
sx q[2];
rz(-1.7147658) q[2];
sx q[2];
rz(1.4437445) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.6988519) q[1];
sx q[1];
rz(2.0069684) q[1];
rz(-1.6601895) q[3];
sx q[3];
rz(-1.2561744) q[3];
sx q[3];
rz(-2.6320626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(0.19533531) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-0.23981747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12407423) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(0.56970861) q[0];
rz(-pi) q[1];
rz(1.8717143) q[2];
sx q[2];
rz(-1.3687203) q[2];
sx q[2];
rz(0.84504715) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1344192) q[1];
sx q[1];
rz(-1.6234142) q[1];
sx q[1];
rz(0.79018553) q[1];
x q[2];
rz(-1.5870985) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-1.714528) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(-1.8260173) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7245674) q[3];
sx q[3];
rz(-2.4745998) q[3];
sx q[3];
rz(0.94661498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6548086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958172) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(-0.70489024) q[0];
rz(2.5321132) q[2];
sx q[2];
rz(-1.8595427) q[2];
sx q[2];
rz(-0.8537054) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1292124) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(-0.48735168) q[1];
rz(0.0067700245) q[3];
sx q[3];
rz(-0.82642185) q[3];
sx q[3];
rz(0.063746728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(-2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-3.1255186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.21503) q[0];
sx q[0];
rz(-2.3030871) q[0];
sx q[0];
rz(-0.82920427) q[0];
x q[1];
rz(-1.5889421) q[2];
sx q[2];
rz(-2.7108253) q[2];
sx q[2];
rz(-0.22242966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(2.1834454) q[1];
rz(-pi) q[2];
rz(0.011868422) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(3.0167937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82299267) q[0];
sx q[0];
rz(-0.66245125) q[0];
sx q[0];
rz(1.2160765) q[0];
x q[1];
rz(-0.20829006) q[2];
sx q[2];
rz(-1.1366833) q[2];
sx q[2];
rz(-1.8722033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9310301) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(2.1767666) q[1];
rz(2.8858658) q[3];
sx q[3];
rz(-0.41318196) q[3];
sx q[3];
rz(2.7030088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.7512158) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-1.1178455) q[0];
rz(-pi) q[1];
rz(1.4540265) q[2];
sx q[2];
rz(-1.8883369) q[2];
sx q[2];
rz(1.5460528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62774819) q[1];
sx q[1];
rz(-0.57250896) q[1];
sx q[1];
rz(-1.7379012) q[1];
x q[2];
rz(-0.84007646) q[3];
sx q[3];
rz(-1.6335765) q[3];
sx q[3];
rz(-1.1508133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(-0.8917121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211593) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-2.349855) q[0];
rz(-pi) q[1];
rz(-0.58015577) q[2];
sx q[2];
rz(-1.9556502) q[2];
sx q[2];
rz(-2.5600524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(2.4119853) q[1];
x q[2];
rz(2.8250152) q[3];
sx q[3];
rz(-2.205875) q[3];
sx q[3];
rz(2.3659335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-0.78424224) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.7766215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9997172) q[0];
sx q[0];
rz(-1.6122072) q[0];
sx q[0];
rz(2.3977604) q[0];
rz(-pi) q[1];
rz(1.7609673) q[2];
sx q[2];
rz(-0.41751465) q[2];
sx q[2];
rz(-1.4830358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2024467) q[1];
sx q[1];
rz(-0.89679634) q[1];
sx q[1];
rz(1.1254315) q[1];
rz(-0.19946675) q[3];
sx q[3];
rz(-1.8623127) q[3];
sx q[3];
rz(2.9150972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6598845) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(2.7754521) q[0];
rz(-pi) q[1];
rz(2.7231611) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(0.07428169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89322461) q[1];
sx q[1];
rz(-1.4239422) q[1];
sx q[1];
rz(-2.8180772) q[1];
rz(-pi) q[2];
rz(-1.5714985) q[3];
sx q[3];
rz(-0.32214221) q[3];
sx q[3];
rz(0.57442564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(-0.36623476) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(-0.33404074) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
rz(-0.35477521) q[3];
sx q[3];
rz(-1.0412024) q[3];
sx q[3];
rz(2.1333221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
