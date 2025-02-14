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
rz(-0.35141355) q[0];
sx q[0];
rz(-2.0961528) q[0];
rz(-2.2975629) q[1];
sx q[1];
rz(-0.56517711) q[1];
sx q[1];
rz(-1.3619818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1050674) q[0];
sx q[0];
rz(-2.1224667) q[0];
sx q[0];
rz(1.9776634) q[0];
x q[1];
rz(1.3711113) q[2];
sx q[2];
rz(-0.43569316) q[2];
sx q[2];
rz(-2.1148122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2842362) q[1];
sx q[1];
rz(-1.7928837) q[1];
sx q[1];
rz(-0.41389334) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0697738) q[3];
sx q[3];
rz(-1.27099) q[3];
sx q[3];
rz(-2.6553947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7094946) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(-0.70322767) q[2];
rz(-0.75913366) q[3];
sx q[3];
rz(-2.1742564) q[3];
sx q[3];
rz(2.4659992) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21509875) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(-2.3469927) q[0];
rz(-0.82851797) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(0.4812831) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3743418) q[0];
sx q[0];
rz(-1.6899979) q[0];
sx q[0];
rz(1.5425372) q[0];
rz(-pi) q[1];
rz(-0.74864804) q[2];
sx q[2];
rz(-1.949913) q[2];
sx q[2];
rz(0.15896713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0121721) q[1];
sx q[1];
rz(-2.8986065) q[1];
sx q[1];
rz(-1.8956196) q[1];
x q[2];
rz(-2.4886143) q[3];
sx q[3];
rz(-2.0275473) q[3];
sx q[3];
rz(1.8819609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98722297) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(0.65732035) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(-1.9067732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2051314) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(-2.8954647) q[0];
rz(-0.71324619) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(-0.86457843) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3636091) q[0];
sx q[0];
rz(-0.37425259) q[0];
sx q[0];
rz(0.97866367) q[0];
rz(-2.1786392) q[2];
sx q[2];
rz(-1.1216444) q[2];
sx q[2];
rz(2.1445865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8997174) q[1];
sx q[1];
rz(-2.6730099) q[1];
sx q[1];
rz(1.750293) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78737463) q[3];
sx q[3];
rz(-1.678576) q[3];
sx q[3];
rz(-2.7312349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(1.2064639) q[2];
rz(3.0237899) q[3];
sx q[3];
rz(-0.67540568) q[3];
sx q[3];
rz(-0.12731586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.6806108) q[0];
sx q[0];
rz(-0.46517459) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-1.0135087) q[1];
sx q[1];
rz(-1.7101589) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9412044) q[0];
sx q[0];
rz(-0.53576195) q[0];
sx q[0];
rz(-0.8922116) q[0];
rz(-2.0639922) q[2];
sx q[2];
rz(-1.3921129) q[2];
sx q[2];
rz(2.5692232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2273754) q[1];
sx q[1];
rz(-2.3721116) q[1];
sx q[1];
rz(2.6282267) q[1];
rz(-pi) q[2];
rz(-1.2243829) q[3];
sx q[3];
rz(-2.2675121) q[3];
sx q[3];
rz(0.55075607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5005834) q[2];
sx q[2];
rz(-0.21037978) q[2];
sx q[2];
rz(1.496199) q[2];
rz(0.03838852) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(0.7588318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043623) q[0];
sx q[0];
rz(-2.4431603) q[0];
sx q[0];
rz(1.5414365) q[0];
rz(1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(1.6391594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1609057) q[0];
sx q[0];
rz(-1.1146) q[0];
sx q[0];
rz(-2.6260757) q[0];
rz(-pi) q[1];
rz(2.9011594) q[2];
sx q[2];
rz(-2.5747262) q[2];
sx q[2];
rz(-2.4971003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9683343) q[1];
sx q[1];
rz(-1.1458261) q[1];
sx q[1];
rz(2.6476184) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5030902) q[3];
sx q[3];
rz(-1.4622524) q[3];
sx q[3];
rz(-2.6084765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1372823) q[2];
sx q[2];
rz(-2.8305125) q[2];
sx q[2];
rz(2.7610371) q[2];
rz(0.54008326) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(-1.165747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.49749097) q[0];
sx q[0];
rz(-3.0971425) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(3.0907471) q[1];
sx q[1];
rz(-2.6652002) q[1];
sx q[1];
rz(1.4643889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58005262) q[0];
sx q[0];
rz(-2.0140352) q[0];
sx q[0];
rz(1.1053941) q[0];
rz(2.5821677) q[2];
sx q[2];
rz(-0.96772268) q[2];
sx q[2];
rz(0.81808264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7685827) q[1];
sx q[1];
rz(-1.4861732) q[1];
sx q[1];
rz(-3.0510034) q[1];
rz(-pi) q[2];
rz(0.6847295) q[3];
sx q[3];
rz(-1.7312418) q[3];
sx q[3];
rz(1.1477276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1579608) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(-0.91510406) q[2];
rz(0.30465952) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(-1.4753531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604928) q[0];
sx q[0];
rz(-2.2486794) q[0];
sx q[0];
rz(-1.1156981) q[0];
rz(2.5005493) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(-1.8152274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4531739) q[0];
sx q[0];
rz(-2.6368615) q[0];
sx q[0];
rz(-2.8674528) q[0];
rz(-0.4037598) q[2];
sx q[2];
rz(-0.40382622) q[2];
sx q[2];
rz(0.60690875) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6815344) q[1];
sx q[1];
rz(-2.0979954) q[1];
sx q[1];
rz(0.74337007) q[1];
x q[2];
rz(0.71621098) q[3];
sx q[3];
rz(-2.1478601) q[3];
sx q[3];
rz(-0.61208597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5850087) q[2];
sx q[2];
rz(-2.1592087) q[2];
sx q[2];
rz(1.2674241) q[2];
rz(-3.0804539) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214445) q[0];
sx q[0];
rz(-2.3763438) q[0];
sx q[0];
rz(0.52566093) q[0];
rz(1.4875745) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(0.18255998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.96344) q[0];
sx q[0];
rz(-1.1530253) q[0];
sx q[0];
rz(1.3304292) q[0];
rz(0.31086773) q[2];
sx q[2];
rz(-2.0070878) q[2];
sx q[2];
rz(2.7300961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9801431) q[1];
sx q[1];
rz(-1.7421466) q[1];
sx q[1];
rz(-0.93659393) q[1];
x q[2];
rz(2.0917039) q[3];
sx q[3];
rz(-2.6080797) q[3];
sx q[3];
rz(-2.3236139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(-0.1087428) q[2];
rz(-0.10793081) q[3];
sx q[3];
rz(-2.7863672) q[3];
sx q[3];
rz(1.7598553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906319) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(2.4884124) q[0];
rz(1.8494891) q[1];
sx q[1];
rz(-0.2457681) q[1];
sx q[1];
rz(-0.076315708) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51173333) q[0];
sx q[0];
rz(-2.7329067) q[0];
sx q[0];
rz(-0.9291533) q[0];
x q[1];
rz(-3.1082319) q[2];
sx q[2];
rz(-1.4028869) q[2];
sx q[2];
rz(-2.1176934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46901033) q[1];
sx q[1];
rz(-2.3206876) q[1];
sx q[1];
rz(-1.7609079) q[1];
rz(-pi) q[2];
rz(2.2003284) q[3];
sx q[3];
rz(-1.2800893) q[3];
sx q[3];
rz(0.45699003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2448347) q[2];
sx q[2];
rz(-1.9568169) q[2];
sx q[2];
rz(0.99311382) q[2];
rz(-1.7874329) q[3];
sx q[3];
rz(-0.5391776) q[3];
sx q[3];
rz(-1.4208992) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83574522) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(-0.29356965) q[0];
rz(-2.1096032) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(2.8242677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5164233) q[0];
sx q[0];
rz(-2.1835071) q[0];
sx q[0];
rz(-1.4132383) q[0];
rz(-1.7909088) q[2];
sx q[2];
rz(-1.0270455) q[2];
sx q[2];
rz(-1.3771283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7723497) q[1];
sx q[1];
rz(-2.0017821) q[1];
sx q[1];
rz(-1.950077) q[1];
x q[2];
rz(-1.6280319) q[3];
sx q[3];
rz(-1.2737107) q[3];
sx q[3];
rz(0.12606584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5130634) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(-2.826706) q[2];
rz(0.033673938) q[3];
sx q[3];
rz(-1.0957054) q[3];
sx q[3];
rz(-1.4382188) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2122129) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(-0.20044151) q[1];
sx q[1];
rz(-1.0841752) q[1];
sx q[1];
rz(-0.60842327) q[1];
rz(3.0931531) q[2];
sx q[2];
rz(-1.5413399) q[2];
sx q[2];
rz(0.057609859) q[2];
rz(1.8157806) q[3];
sx q[3];
rz(-2.1143338) q[3];
sx q[3];
rz(-0.59296617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
