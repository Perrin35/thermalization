OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(4.6580553) q[0];
sx q[0];
rz(9.1604995) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(-0.73524737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7967448) q[0];
sx q[0];
rz(-1.7186527) q[0];
sx q[0];
rz(2.4792838) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60230435) q[2];
sx q[2];
rz(-0.77568433) q[2];
sx q[2];
rz(1.3210981) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3056065) q[1];
sx q[1];
rz(-1.4164093) q[1];
sx q[1];
rz(-1.2812213) q[1];
rz(-1.5427038) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(2.7147646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-1.1038587) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(0.43637481) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(-2.8754821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5824495) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(2.0931431) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1462609) q[2];
sx q[2];
rz(-1.4333712) q[2];
sx q[2];
rz(-0.95057887) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74782158) q[1];
sx q[1];
rz(-2.1854679) q[1];
sx q[1];
rz(-0.6786896) q[1];
rz(-pi) q[2];
rz(-1.3105884) q[3];
sx q[3];
rz(-0.62285138) q[3];
sx q[3];
rz(-2.7130896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(-1.3471289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12362326) q[0];
sx q[0];
rz(-0.53640134) q[0];
sx q[0];
rz(0.53039741) q[0];
x q[1];
rz(-0.15813078) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(0.49027298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.320142) q[1];
sx q[1];
rz(-0.67678932) q[1];
sx q[1];
rz(1.0980781) q[1];
rz(1.7473162) q[3];
sx q[3];
rz(-2.4268097) q[3];
sx q[3];
rz(2.178758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1469664) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(-1.5220801) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20630079) q[0];
sx q[0];
rz(-1.1456523) q[0];
sx q[0];
rz(1.4817119) q[0];
rz(-pi) q[1];
rz(-0.70978769) q[2];
sx q[2];
rz(-0.94883942) q[2];
sx q[2];
rz(-0.70993916) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2019129) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(1.4918785) q[1];
rz(-pi) q[2];
rz(-1.3341321) q[3];
sx q[3];
rz(-0.79573123) q[3];
sx q[3];
rz(-0.79470293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(0.4310472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-0.34367925) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491964) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(-2.2618494) q[0];
rz(-2.5630066) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(1.7631284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5695565) q[1];
sx q[1];
rz(-1.0735638) q[1];
sx q[1];
rz(0.083942457) q[1];
rz(-pi) q[2];
x q[2];
rz(2.214659) q[3];
sx q[3];
rz(-1.5537795) q[3];
sx q[3];
rz(3.0683558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(0.25973928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7344959) q[0];
sx q[0];
rz(-1.2571113) q[0];
sx q[0];
rz(-1.9514027) q[0];
rz(-pi) q[1];
rz(1.514228) q[2];
sx q[2];
rz(-0.33677378) q[2];
sx q[2];
rz(-2.3677504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.427634) q[1];
sx q[1];
rz(-0.40363388) q[1];
sx q[1];
rz(2.1778818) q[1];
rz(-pi) q[2];
rz(1.1214244) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(-1.8024973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-0.10061131) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(1.7939059) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(-2.535634) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46442407) q[0];
sx q[0];
rz(-0.66093984) q[0];
sx q[0];
rz(-2.7695157) q[0];
rz(-pi) q[1];
rz(1.7714959) q[2];
sx q[2];
rz(-1.8569274) q[2];
sx q[2];
rz(-1.2933033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0026605) q[1];
sx q[1];
rz(-1.2298889) q[1];
sx q[1];
rz(2.1041811) q[1];
x q[2];
rz(-0.80375399) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(-0.758981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24017748) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(-2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-3.0775552) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.3834073) q[0];
sx q[0];
rz(-2.1894356) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0520347) q[2];
sx q[2];
rz(-2.5737737) q[2];
sx q[2];
rz(1.0713112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1053972) q[1];
sx q[1];
rz(-0.35696778) q[1];
sx q[1];
rz(3.0371975) q[1];
rz(3.1021032) q[3];
sx q[3];
rz(-1.2590623) q[3];
sx q[3];
rz(-1.5453651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(2.2122038) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14426194) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147946) q[0];
sx q[0];
rz(-2.9920122) q[0];
sx q[0];
rz(3.0339255) q[0];
x q[1];
rz(-0.53290368) q[2];
sx q[2];
rz(-2.1295296) q[2];
sx q[2];
rz(-1.2677873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3582663) q[1];
sx q[1];
rz(-1.3291385) q[1];
sx q[1];
rz(0.011947167) q[1];
x q[2];
rz(-2.1298218) q[3];
sx q[3];
rz(-1.8852919) q[3];
sx q[3];
rz(2.3712456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(-2.9679427) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-0.39150795) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(0.82410518) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3819645) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(2.1419924) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27406759) q[2];
sx q[2];
rz(-2.0147418) q[2];
sx q[2];
rz(-2.519671) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13092566) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(-0.24433498) q[1];
rz(-pi) q[2];
rz(-0.70865788) q[3];
sx q[3];
rz(-1.7950247) q[3];
sx q[3];
rz(1.538016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-0.097975227) q[2];
sx q[2];
rz(-0.91839472) q[2];
sx q[2];
rz(0.67994946) q[2];
rz(-2.0588856) q[3];
sx q[3];
rz(-1.6255717) q[3];
sx q[3];
rz(-2.0902904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
