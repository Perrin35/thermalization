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
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5188132) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.8151087) q[0];
x q[1];
rz(0.15687234) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(-0.064531782) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1374955) q[1];
sx q[1];
rz(-1.6988519) q[1];
sx q[1];
rz(1.1346243) q[1];
rz(-pi) q[2];
rz(-2.8257915) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(-2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.9462573) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-2.9017752) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12407423) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(-2.571884) q[0];
rz(-pi) q[1];
rz(0.2113091) q[2];
sx q[2];
rz(-1.276187) q[2];
sx q[2];
rz(0.78795563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61566831) q[1];
sx q[1];
rz(-0.79155542) q[1];
sx q[1];
rz(-0.073992373) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.014157045) q[3];
sx q[3];
rz(-1.5544958) q[3];
sx q[3];
rz(1.0316767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6372765) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.2770237) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(-1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(1.0659165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695755) q[0];
sx q[0];
rz(-1.6117275) q[0];
sx q[0];
rz(1.8191562) q[0];
rz(-pi) q[1];
rz(2.5862323) q[2];
sx q[2];
rz(-1.4484324) q[2];
sx q[2];
rz(-0.33081474) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4855811) q[1];
sx q[1];
rz(-1.0080358) q[1];
sx q[1];
rz(0.33388168) q[1];
rz(0.90956456) q[3];
sx q[3];
rz(-1.665691) q[3];
sx q[3];
rz(0.50300099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(0.70670635) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.810881) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(0.34513721) q[0];
x q[1];
rz(-1.9183667) q[2];
sx q[2];
rz(-0.98993694) q[2];
sx q[2];
rz(0.52085224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.023991) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(2.8531122) q[1];
rz(2.3151822) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(-1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(3.1255186) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.21503) q[0];
sx q[0];
rz(-2.3030871) q[0];
sx q[0];
rz(-2.3123884) q[0];
rz(-1.5889421) q[2];
sx q[2];
rz(-2.7108253) q[2];
sx q[2];
rz(2.919163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(2.1834454) q[1];
rz(-pi) q[2];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.9330977) q[3];
sx q[3];
rz(3.0294861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(-0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(1.2247359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46366102) q[0];
sx q[0];
rz(-1.7860798) q[0];
sx q[0];
rz(2.2023375) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0133063) q[2];
sx q[2];
rz(-1.759521) q[2];
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
rz(2.9310301) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(2.1767666) q[1];
rz(-pi) q[2];
rz(-2.8858658) q[3];
sx q[3];
rz(-0.41318196) q[3];
sx q[3];
rz(-2.7030088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(-2.2018946) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66877767) q[0];
sx q[0];
rz(-2.6758709) q[0];
sx q[0];
rz(-1.8229683) q[0];
x q[1];
rz(0.31957303) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(-3.0802397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82579457) q[1];
sx q[1];
rz(-2.1343532) q[1];
sx q[1];
rz(-0.10679306) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0573781) q[3];
sx q[3];
rz(-2.299752) q[3];
sx q[3];
rz(-0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(0.022162612) q[2];
rz(2.3968905) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(-2.1180034) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(-2.2498806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72043334) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-0.79173761) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5052091) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(0.46905876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5276287) q[1];
sx q[1];
rz(-0.91273897) q[1];
sx q[1];
rz(1.0440473) q[1];
rz(-pi) q[2];
rz(0.31657747) q[3];
sx q[3];
rz(-2.205875) q[3];
sx q[3];
rz(-2.3659335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(2.908356) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(2.419557) q[0];
rz(-2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.3649712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(3.0804759) q[0];
x q[1];
rz(1.3806254) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(1.6585569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.939146) q[1];
sx q[1];
rz(-0.89679634) q[1];
sx q[1];
rz(-1.1254315) q[1];
rz(-2.9421259) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(2.9150972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.6636794) q[0];
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
rz(-1.4817081) q[0];
sx q[0];
rz(-1.9511001) q[0];
sx q[0];
rz(-2.7754521) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(2.248368) q[1];
sx q[1];
rz(-1.4239422) q[1];
sx q[1];
rz(-0.3235154) q[1];
rz(-3.1413583) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.8075519) q[2];
sx q[2];
rz(-2.3636849) q[2];
sx q[2];
rz(-1.452527) q[2];
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
