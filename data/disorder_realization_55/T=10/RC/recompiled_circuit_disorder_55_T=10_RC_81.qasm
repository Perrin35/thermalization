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
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(-1.2101313) q[1];
sx q[1];
rz(0.73524737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7967448) q[0];
sx q[0];
rz(-1.7186527) q[0];
sx q[0];
rz(0.66230886) q[0];
x q[1];
rz(-0.60230435) q[2];
sx q[2];
rz(-0.77568433) q[2];
sx q[2];
rz(1.3210981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.211261) q[1];
sx q[1];
rz(-0.32713612) q[1];
sx q[1];
rz(1.0717908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5611476) q[3];
sx q[3];
rz(-1.5942897) q[3];
sx q[3];
rz(1.1593727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(-0.27403533) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
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
rz(0.46288681) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-0.26611051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55914315) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(1.0484496) q[0];
rz(0.16337784) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(0.53158224) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7539566) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(2.3073879) q[1];
rz(-pi) q[2];
rz(-1.8310043) q[3];
sx q[3];
rz(-0.62285138) q[3];
sx q[3];
rz(-0.42850307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-0.28765837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4192393) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(1.2786352) q[0];
x q[1];
rz(-0.15813078) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(-2.6513197) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9007064) q[1];
sx q[1];
rz(-2.1623003) q[1];
sx q[1];
rz(-2.790931) q[1];
rz(-1.3942765) q[3];
sx q[3];
rz(-0.71478292) q[3];
sx q[3];
rz(-2.178758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(-2.6456397) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(-2.4194338) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(0.55975634) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4013195) q[0];
sx q[0];
rz(-1.6519321) q[0];
sx q[0];
rz(-0.42663891) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83300029) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(-2.8766362) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9396797) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(-1.6497142) q[1];
rz(-1.3341321) q[3];
sx q[3];
rz(-2.3458614) q[3];
sx q[3];
rz(-2.3468897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42797783) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-0.4310472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(2.7979134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30003798) q[0];
sx q[0];
rz(-2.2404788) q[0];
sx q[0];
rz(-0.71398736) q[0];
rz(-0.86954388) q[2];
sx q[2];
rz(-1.1079259) q[2];
sx q[2];
rz(-2.9733544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1829454) q[1];
sx q[1];
rz(-1.6445541) q[1];
sx q[1];
rz(-2.069509) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.120317) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(1.4847886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-0.053744944) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(2.8818534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.040859) q[0];
sx q[0];
rz(-1.9319527) q[0];
sx q[0];
rz(0.33613899) q[0];
rz(0.01979205) q[2];
sx q[2];
rz(-1.2345825) q[2];
sx q[2];
rz(-2.3078231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.427634) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(-0.9637109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5335347) q[3];
sx q[3];
rz(-0.44964368) q[3];
sx q[3];
rz(-0.26526181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(0.10061131) q[2];
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
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1473734) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(-0.60595864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92354846) q[0];
sx q[0];
rz(-0.96203066) q[0];
sx q[0];
rz(-1.2952842) q[0];
x q[1];
rz(-2.8499243) q[2];
sx q[2];
rz(-1.3783611) q[2];
sx q[2];
rz(0.33484909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1389321) q[1];
sx q[1];
rz(-1.2298889) q[1];
sx q[1];
rz(2.1041811) q[1];
rz(2.3378387) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(-0.758981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(1.3700221) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-3.0775552) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768351) q[0];
sx q[0];
rz(-1.7581853) q[0];
sx q[0];
rz(2.1894356) q[0];
x q[1];
rz(-1.0649101) q[2];
sx q[2];
rz(-1.8407028) q[2];
sx q[2];
rz(3.0907061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9248326) q[1];
sx q[1];
rz(-1.9257345) q[1];
sx q[1];
rz(-1.6096398) q[1];
rz(-1.4488892) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(-0.88820109) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-2.2122038) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(1.4319179) q[0];
rz(-2.572708) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(1.127839) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9147946) q[0];
sx q[0];
rz(-2.9920122) q[0];
sx q[0];
rz(3.0339255) q[0];
rz(-pi) q[1];
rz(0.88840975) q[2];
sx q[2];
rz(-2.3896304) q[2];
sx q[2];
rz(-1.0351406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4081501) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-0.3663775) q[3];
sx q[3];
rz(-2.0994224) q[3];
sx q[3];
rz(-0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-2.7500847) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.4655112) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(2.5691659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51047072) q[0];
sx q[0];
rz(-1.0902176) q[0];
sx q[0];
rz(-2.5170588) q[0];
rz(-pi) q[1];
rz(2.8675251) q[2];
sx q[2];
rz(-2.0147418) q[2];
sx q[2];
rz(-2.519671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4116284) q[1];
sx q[1];
rz(-0.3246383) q[1];
sx q[1];
rz(-2.4050729) q[1];
rz(-pi) q[2];
rz(-2.4329348) q[3];
sx q[3];
rz(-1.7950247) q[3];
sx q[3];
rz(-1.538016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(0.53722107) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(-2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.854241) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(0.46335012) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(-3.0436174) q[2];
sx q[2];
rz(-2.2231979) q[2];
sx q[2];
rz(-2.4616432) q[2];
rz(1.6871917) q[3];
sx q[3];
rz(-0.4909066) q[3];
sx q[3];
rz(2.7248513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
