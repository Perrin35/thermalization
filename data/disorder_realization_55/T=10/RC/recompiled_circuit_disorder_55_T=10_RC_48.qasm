OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(0.2642785) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(-0.73524737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-2.465415) q[0];
sx q[0];
rz(-2.9039608) q[0];
rz(0.60230435) q[2];
sx q[2];
rz(-0.77568433) q[2];
sx q[2];
rz(1.8204945) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9303317) q[1];
sx q[1];
rz(-0.32713612) q[1];
sx q[1];
rz(-1.0717908) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5427038) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(-2.7147646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-2.0377339) q[2];
rz(1.8707229) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(2.7052178) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(0.26611051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3811831) q[0];
sx q[0];
rz(-1.1927483) q[0];
sx q[0];
rz(0.80947431) q[0];
rz(-pi) q[1];
rz(1.8196462) q[2];
sx q[2];
rz(-0.58983931) q[2];
sx q[2];
rz(-2.313254) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38763603) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(-2.3073879) q[1];
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
rz(-0.46488547) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(0.28765837) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1613306) q[0];
sx q[0];
rz(-1.3093003) q[0];
sx q[0];
rz(-0.47388347) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2855661) q[2];
sx q[2];
rz(-1.7226379) q[2];
sx q[2];
rz(1.0361995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82145065) q[1];
sx q[1];
rz(-0.67678932) q[1];
sx q[1];
rz(-2.0435145) q[1];
x q[2];
rz(-2.9903528) q[3];
sx q[3];
rz(-2.2721707) q[3];
sx q[3];
rz(-2.410694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(2.6456397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(0.55975634) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9352919) q[0];
sx q[0];
rz(-1.9959404) q[0];
sx q[0];
rz(-1.4817119) q[0];
rz(0.81360929) q[2];
sx q[2];
rz(-1.0126197) q[2];
sx q[2];
rz(1.8166325) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2019129) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(-1.6497142) q[1];
rz(-0.23493725) q[3];
sx q[3];
rz(-2.3384691) q[3];
sx q[3];
rz(1.1266176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(-1.7170061) q[2];
rz(-2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-0.4310472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(2.7979134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8415547) q[0];
sx q[0];
rz(-2.2404788) q[0];
sx q[0];
rz(0.71398736) q[0];
rz(-0.57858606) q[2];
sx q[2];
rz(-2.1861976) q[2];
sx q[2];
rz(1.7631284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3949563) q[1];
sx q[1];
rz(-0.50368217) q[1];
sx q[1];
rz(1.4175182) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92693365) q[3];
sx q[3];
rz(-1.5878131) q[3];
sx q[3];
rz(0.073236853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(2.8500309) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(3.1164363) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(0.25973928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.320967) q[0];
sx q[0];
rz(-2.6532986) q[0];
sx q[0];
rz(0.85296209) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6273646) q[2];
sx q[2];
rz(-0.33677378) q[2];
sx q[2];
rz(0.77384225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7139587) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(-0.9637109) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0201683) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(1.3390954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-3.0409813) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-2.8314262) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(0.60595864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180442) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(-1.8463085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3700968) q[2];
sx q[2];
rz(-1.2846652) q[2];
sx q[2];
rz(-1.8482894) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3785432) q[1];
sx q[1];
rz(-1.0710643) q[1];
sx q[1];
rz(2.7507957) q[1];
x q[2];
rz(-2.6113854) q[3];
sx q[3];
rz(-1.0876417) q[3];
sx q[3];
rz(0.14164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(0.85285464) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-0.064037474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1916699) q[0];
sx q[0];
rz(-2.4987698) q[0];
sx q[0];
rz(-1.886801) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30631752) q[2];
sx q[2];
rz(-1.0848572) q[2];
sx q[2];
rz(1.475032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.774051) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(-0.35518412) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4488892) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14426194) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(1.7096747) q[0];
rz(-2.572708) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(-1.127839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8059212) q[0];
sx q[0];
rz(-1.7195042) q[0];
sx q[0];
rz(1.5869897) q[0];
x q[1];
rz(0.88840975) q[2];
sx q[2];
rz(-0.75196224) q[2];
sx q[2];
rz(1.0351406) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4081501) q[1];
sx q[1];
rz(-2.8996455) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0117709) q[3];
sx q[3];
rz(-1.2563007) q[3];
sx q[3];
rz(-0.77034706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.613712) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(2.7500847) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7596282) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(2.1419924) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27406759) q[2];
sx q[2];
rz(-1.1268508) q[2];
sx q[2];
rz(2.519671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13092566) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(-2.8972577) q[1];
x q[2];
rz(-2.8045373) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(0.22102236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(0.53722107) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(0.46335012) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(-2.2255185) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(-1.6871917) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
