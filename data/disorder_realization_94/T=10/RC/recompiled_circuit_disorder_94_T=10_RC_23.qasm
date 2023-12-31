OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(-1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(0.61520666) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078755137) q[0];
sx q[0];
rz(-0.7987928) q[0];
sx q[0];
rz(0.90057217) q[0];
rz(1.3846606) q[2];
sx q[2];
rz(-1.2190483) q[2];
sx q[2];
rz(-0.56439161) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8021009) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(-0.93356737) q[1];
x q[2];
rz(-1.7513566) q[3];
sx q[3];
rz(-1.872634) q[3];
sx q[3];
rz(-0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-2.8033076) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.5240086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.108792) q[0];
sx q[0];
rz(-1.5724702) q[0];
sx q[0];
rz(-1.3711506) q[0];
rz(1.5288058) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(-3.0595879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6807032) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(-2.1417888) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0475572) q[3];
sx q[3];
rz(-1.3841076) q[3];
sx q[3];
rz(-1.7969014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(2.8895203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789327) q[0];
sx q[0];
rz(-1.1061215) q[0];
sx q[0];
rz(-2.4436823) q[0];
x q[1];
rz(0.88148586) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(-1.880868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1356126) q[1];
sx q[1];
rz(-0.95893919) q[1];
sx q[1];
rz(2.9115885) q[1];
x q[2];
rz(-2.4113703) q[3];
sx q[3];
rz(-2.632004) q[3];
sx q[3];
rz(-2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0079460572) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(-0.72511073) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44666501) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(2.5926673) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9634339) q[1];
sx q[1];
rz(-1.6903094) q[1];
sx q[1];
rz(-0.91576373) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2624192) q[3];
sx q[3];
rz(-1.3652507) q[3];
sx q[3];
rz(1.8431078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2924071) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(0.25462338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.03527) q[0];
sx q[0];
rz(-1.5570886) q[0];
sx q[0];
rz(1.4916219) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5931555) q[2];
sx q[2];
rz(-2.5639113) q[2];
sx q[2];
rz(-0.81800848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29171195) q[1];
sx q[1];
rz(-0.92392081) q[1];
sx q[1];
rz(-1.2667659) q[1];
rz(-pi) q[2];
rz(-2.7084064) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(-1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(0.20656955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26039133) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(0.3221237) q[0];
rz(2.0986404) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(-2.8258459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088591136) q[1];
sx q[1];
rz(-2.2676761) q[1];
sx q[1];
rz(0.23115302) q[1];
rz(-pi) q[2];
rz(-2.8387186) q[3];
sx q[3];
rz(-0.93749638) q[3];
sx q[3];
rz(-2.5686222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(-2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(-0.38152951) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.3279703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0921811) q[0];
sx q[0];
rz(-1.6850527) q[0];
sx q[0];
rz(1.1294424) q[0];
x q[1];
rz(-2.1778657) q[2];
sx q[2];
rz(-0.71787314) q[2];
sx q[2];
rz(1.8075862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5495587) q[1];
sx q[1];
rz(-0.70211071) q[1];
sx q[1];
rz(-1.3105427) q[1];
rz(-0.6789356) q[3];
sx q[3];
rz(-2.5151765) q[3];
sx q[3];
rz(-0.39144799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(-0.95091933) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61688214) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(-pi) q[1];
rz(1.7958926) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(0.30579145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9329405) q[1];
sx q[1];
rz(-1.3622074) q[1];
sx q[1];
rz(1.5044466) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37787921) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9178847) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10044554) q[0];
sx q[0];
rz(-2.1491096) q[0];
sx q[0];
rz(2.5006177) q[0];
rz(1.3525891) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(-1.6831236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1253423) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(-0.36228212) q[1];
x q[2];
rz(1.6660059) q[3];
sx q[3];
rz(-2.6009437) q[3];
sx q[3];
rz(1.8298061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10782345) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(-2.4643722) q[0];
x q[1];
rz(-2.5752441) q[2];
sx q[2];
rz(-2.1841335) q[2];
sx q[2];
rz(1.4373506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1199011) q[1];
sx q[1];
rz(-2.7530115) q[1];
sx q[1];
rz(2.4619224) q[1];
rz(-pi) q[2];
rz(2.8519248) q[3];
sx q[3];
rz(-2.3710459) q[3];
sx q[3];
rz(2.88248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62853652) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.097461854) q[2];
sx q[2];
rz(-0.41257358) q[2];
sx q[2];
rz(-1.67795) q[2];
rz(1.5619754) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
