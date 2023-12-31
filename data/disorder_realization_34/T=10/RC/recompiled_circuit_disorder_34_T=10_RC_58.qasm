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
rz(-2.4913139) q[1];
sx q[1];
rz(4.4250017) q[1];
sx q[1];
rz(11.783574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5188132) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(-1.326484) q[0];
x q[1];
rz(-0.15687234) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(0.064531782) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3738457) q[1];
sx q[1];
rz(-1.1384374) q[1];
sx q[1];
rz(0.14111622) q[1];
rz(-1.4814032) q[3];
sx q[3];
rz(-1.8854183) q[3];
sx q[3];
rz(0.50953007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-0.23981747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(-0.56970861) q[0];
rz(-pi) q[1];
rz(1.2698783) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(0.84504715) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0071734) q[1];
sx q[1];
rz(-1.5181784) q[1];
sx q[1];
rz(-0.79018553) q[1];
rz(2.2858743) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(1.74687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(-2.0756762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6324368) q[0];
sx q[0];
rz(-1.8189438) q[0];
sx q[0];
rz(3.0993673) q[0];
x q[1];
rz(1.4270646) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(1.3155754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2324893) q[1];
sx q[1];
rz(-0.64503446) q[1];
sx q[1];
rz(2.049936) q[1];
rz(-pi) q[2];
rz(-0.12001868) q[3];
sx q[3];
rz(-0.91306049) q[3];
sx q[3];
rz(-1.1413871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.4867841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810881) q[0];
sx q[0];
rz(-1.2905621) q[0];
sx q[0];
rz(2.7964554) q[0];
x q[1];
rz(-0.60947946) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.023991) q[1];
sx q[1];
rz(-1.4211434) q[1];
sx q[1];
rz(-0.28848044) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(3.0778459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(2.0452943) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(-3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094729312) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(-0.88386436) q[0];
x q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(-0.24239937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2542418) q[1];
sx q[1];
rz(-2.4130531) q[1];
sx q[1];
rz(2.2345047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(3.0294861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(1.9168568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8794494) q[0];
sx q[0];
rz(-0.95603846) q[0];
sx q[0];
rz(-2.8770147) q[0];
rz(-pi) q[1];
rz(-2.0133063) q[2];
sx q[2];
rz(-1.759521) q[2];
sx q[2];
rz(-0.21274266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8606557) q[1];
sx q[1];
rz(-0.60815629) q[1];
sx q[1];
rz(1.6673191) q[1];
rz(-1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-2.9810867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619103) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38802559) q[0];
sx q[0];
rz(-2.0206847) q[0];
sx q[0];
rz(-0.1247503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4540265) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(-1.5460528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82579457) q[1];
sx q[1];
rz(-1.0072395) q[1];
sx q[1];
rz(0.10679306) q[1];
rz(-3.0573781) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(-2.6654411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(-0.8917121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69670024) q[0];
sx q[0];
rz(-2.3377044) q[0];
sx q[0];
rz(-2.9219887) q[0];
rz(-pi) q[1];
rz(-2.5052091) q[2];
sx q[2];
rz(-2.4578265) q[2];
sx q[2];
rz(-2.6725339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5276287) q[1];
sx q[1];
rz(-2.2288537) q[1];
sx q[1];
rz(-1.0440473) q[1];
rz(-pi) q[2];
rz(0.31657747) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(-0.77565912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-0.72203565) q[0];
rz(-2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-2.3968292) q[0];
sx q[0];
rz(0.061116771) q[0];
x q[1];
rz(-1.3806254) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(1.4830358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2024467) q[1];
sx q[1];
rz(-0.89679634) q[1];
sx q[1];
rz(-2.0161611) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1544103) q[3];
sx q[3];
rz(-2.7899788) q[3];
sx q[3];
rz(-0.83948638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1086796) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(-3.1075297) q[3];
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
rz(2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(2.058303) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9111239) q[0];
sx q[0];
rz(-1.9096806) q[0];
sx q[0];
rz(1.1662657) q[0];
rz(-2.7231611) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(-3.067311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26607516) q[1];
sx q[1];
rz(-0.35421687) q[1];
sx q[1];
rz(0.43550272) q[1];
rz(-1.5714985) q[3];
sx q[3];
rz(-2.8194504) q[3];
sx q[3];
rz(2.567167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-0.74950779) q[2];
sx q[2];
rz(-1.3386249) q[2];
sx q[2];
rz(-2.7809536) q[2];
rz(1.0127388) q[3];
sx q[3];
rz(-1.8752718) q[3];
sx q[3];
rz(0.74753052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
