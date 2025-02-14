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
rz(3.0577793) q[0];
sx q[0];
rz(2.7901791) q[0];
sx q[0];
rz(8.3793381) q[0];
rz(-5.4391556) q[1];
sx q[1];
rz(0.56517711) q[1];
sx q[1];
rz(11.204389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3484074) q[0];
sx q[0];
rz(-2.4688666) q[0];
sx q[0];
rz(0.5714709) q[0];
rz(-0.092081618) q[2];
sx q[2];
rz(-1.9972587) q[2];
sx q[2];
rz(-1.895176) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5246084) q[1];
sx q[1];
rz(-1.9739264) q[1];
sx q[1];
rz(1.3289846) q[1];
x q[2];
rz(2.8027439) q[3];
sx q[3];
rz(-1.0940427) q[3];
sx q[3];
rz(0.92425215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7094946) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(-2.438365) q[2];
rz(0.75913366) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(2.4659992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264939) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(0.79459992) q[0];
rz(-2.3130747) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(-0.4812831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3414087) q[0];
sx q[0];
rz(-1.5427378) q[0];
sx q[0];
rz(-3.0223439) q[0];
rz(-pi) q[1];
rz(1.0727302) q[2];
sx q[2];
rz(-2.2553372) q[2];
sx q[2];
rz(1.7430151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0161288) q[1];
sx q[1];
rz(-1.4939346) q[1];
sx q[1];
rz(1.8015339) q[1];
x q[2];
rz(-1.0167502) q[3];
sx q[3];
rz(-2.1475128) q[3];
sx q[3];
rz(-0.63652906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98722297) q[2];
sx q[2];
rz(-2.1798446) q[2];
sx q[2];
rz(-0.50755429) q[2];
rz(-0.65732035) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.9067732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364612) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(0.24612799) q[0];
rz(2.4283465) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(2.2770142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9893892) q[0];
sx q[0];
rz(-1.2626022) q[0];
sx q[0];
rz(-2.9257923) q[0];
x q[1];
rz(2.6106493) q[2];
sx q[2];
rz(-1.0304255) q[2];
sx q[2];
rz(-0.8671538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9731701) q[1];
sx q[1];
rz(-1.4900786) q[1];
sx q[1];
rz(1.1087085) q[1];
rz(-2.354218) q[3];
sx q[3];
rz(-1.678576) q[3];
sx q[3];
rz(-0.4103578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(1.9351287) q[2];
rz(0.11780277) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(-0.12731586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8449611) q[0];
sx q[0];
rz(-1.6806108) q[0];
sx q[0];
rz(0.46517459) q[0];
rz(2.2370715) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(1.7101589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44712191) q[0];
sx q[0];
rz(-1.9794802) q[0];
sx q[0];
rz(-2.7848836) q[0];
rz(-1.0776005) q[2];
sx q[2];
rz(-1.7494798) q[2];
sx q[2];
rz(2.5692232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5618973) q[1];
sx q[1];
rz(-0.919678) q[1];
sx q[1];
rz(-1.1267594) q[1];
rz(-0.72700084) q[3];
sx q[3];
rz(-1.3073587) q[3];
sx q[3];
rz(0.7924197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5005834) q[2];
sx q[2];
rz(-0.21037978) q[2];
sx q[2];
rz(1.6453936) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(-0.7588318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043623) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(-1.5414365) q[0];
rz(-1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(1.5024332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1609057) q[0];
sx q[0];
rz(-1.1146) q[0];
sx q[0];
rz(-2.6260757) q[0];
rz(-0.2404332) q[2];
sx q[2];
rz(-2.5747262) q[2];
sx q[2];
rz(-2.4971003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17325832) q[1];
sx q[1];
rz(-1.9957665) q[1];
sx q[1];
rz(2.6476184) q[1];
x q[2];
rz(1.6385025) q[3];
sx q[3];
rz(-1.4622524) q[3];
sx q[3];
rz(0.53311611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0043103546) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(0.38055554) q[2];
rz(-2.6015094) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(-1.165747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.6441017) q[0];
sx q[0];
rz(-3.0971425) q[0];
sx q[0];
rz(-2.7569726) q[0];
rz(-3.0907471) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(-1.6772038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58005262) q[0];
sx q[0];
rz(-1.1275575) q[0];
sx q[0];
rz(1.1053941) q[0];
rz(-pi) q[1];
rz(-0.88846859) q[2];
sx q[2];
rz(-2.0231477) q[2];
sx q[2];
rz(0.41146989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44838833) q[1];
sx q[1];
rz(-0.12388661) q[1];
sx q[1];
rz(2.3883567) q[1];
x q[2];
rz(0.25050779) q[3];
sx q[3];
rz(-2.4412812) q[3];
sx q[3];
rz(-0.61628646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1579608) q[2];
sx q[2];
rz(-2.3775358) q[2];
sx q[2];
rz(0.91510406) q[2];
rz(-2.8369331) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604928) q[0];
sx q[0];
rz(-0.89291328) q[0];
sx q[0];
rz(-2.0258946) q[0];
rz(-2.5005493) q[1];
sx q[1];
rz(-1.8843001) q[1];
sx q[1];
rz(1.3263652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68841877) q[0];
sx q[0];
rz(-0.50473112) q[0];
sx q[0];
rz(-2.8674528) q[0];
x q[1];
rz(-0.4037598) q[2];
sx q[2];
rz(-0.40382622) q[2];
sx q[2];
rz(0.60690875) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5976482) q[1];
sx q[1];
rz(-2.1955288) q[1];
sx q[1];
rz(-0.90170699) q[1];
x q[2];
rz(-0.78115233) q[3];
sx q[3];
rz(-2.255126) q[3];
sx q[3];
rz(-0.39881166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.556584) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(-1.8741685) q[2];
rz(3.0804539) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(-2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1214445) q[0];
sx q[0];
rz(-0.76524884) q[0];
sx q[0];
rz(-2.6159317) q[0];
rz(-1.6540182) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(-2.9590327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8480685) q[0];
sx q[0];
rz(-1.7901359) q[0];
sx q[0];
rz(2.7128986) q[0];
x q[1];
rz(-2.8307249) q[2];
sx q[2];
rz(-2.0070878) q[2];
sx q[2];
rz(2.7300961) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1614496) q[1];
sx q[1];
rz(-1.399446) q[1];
sx q[1];
rz(2.2049987) q[1];
rz(2.0917039) q[3];
sx q[3];
rz(-2.6080797) q[3];
sx q[3];
rz(0.81797879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(-3.0328499) q[2];
rz(3.0336618) q[3];
sx q[3];
rz(-0.35522541) q[3];
sx q[3];
rz(-1.7598553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906319) q[0];
sx q[0];
rz(-1.2010295) q[0];
sx q[0];
rz(-2.4884124) q[0];
rz(-1.8494891) q[1];
sx q[1];
rz(-2.8958246) q[1];
sx q[1];
rz(3.0652769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.171577) q[0];
sx q[0];
rz(-1.8948012) q[0];
sx q[0];
rz(-2.8879763) q[0];
x q[1];
rz(3.1082319) q[2];
sx q[2];
rz(-1.7387058) q[2];
sx q[2];
rz(1.0238992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9093949) q[1];
sx q[1];
rz(-1.4320717) q[1];
sx q[1];
rz(-0.75896427) q[1];
x q[2];
rz(2.2003284) q[3];
sx q[3];
rz(-1.8615033) q[3];
sx q[3];
rz(2.6846026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2448347) q[2];
sx q[2];
rz(-1.1847757) q[2];
sx q[2];
rz(2.1484788) q[2];
rz(1.7874329) q[3];
sx q[3];
rz(-2.6024151) q[3];
sx q[3];
rz(-1.4208992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83574522) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(-2.848023) q[0];
rz(-2.1096032) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(-0.31732496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5164233) q[0];
sx q[0];
rz(-0.95808555) q[0];
sx q[0];
rz(-1.4132383) q[0];
rz(-pi) q[1];
rz(1.7909088) q[2];
sx q[2];
rz(-1.0270455) q[2];
sx q[2];
rz(-1.7644644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60734487) q[1];
sx q[1];
rz(-0.56612724) q[1];
sx q[1];
rz(0.67791636) q[1];
rz(-pi) q[2];
rz(-2.844048) q[3];
sx q[3];
rz(-1.625522) q[3];
sx q[3];
rz(-1.7136337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5130634) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(-0.31488669) q[2];
rz(3.1079187) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(-1.4382188) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92937975) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(0.20044151) q[1];
sx q[1];
rz(-2.0574175) q[1];
sx q[1];
rz(2.5331694) q[1];
rz(-0.54666109) q[2];
sx q[2];
rz(-0.056686747) q[2];
sx q[2];
rz(-2.0591339) q[2];
rz(2.759886) q[3];
sx q[3];
rz(-2.5504938) q[3];
sx q[3];
rz(2.9989178) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
