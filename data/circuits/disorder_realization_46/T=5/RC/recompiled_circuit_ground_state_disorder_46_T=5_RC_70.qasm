OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.159654) q[0];
sx q[0];
rz(-1.5647793) q[0];
sx q[0];
rz(1.2793581) q[0];
rz(-2.1518985) q[1];
sx q[1];
rz(-1.9446179) q[1];
sx q[1];
rz(-2.9762414) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5996313) q[0];
sx q[0];
rz(-1.8867869) q[0];
sx q[0];
rz(0.91008474) q[0];
x q[1];
rz(-2.3897091) q[2];
sx q[2];
rz(-0.57752973) q[2];
sx q[2];
rz(-2.5381388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27439427) q[1];
sx q[1];
rz(-1.1667248) q[1];
sx q[1];
rz(-2.1073125) q[1];
x q[2];
rz(1.8824485) q[3];
sx q[3];
rz(-1.8168279) q[3];
sx q[3];
rz(0.450799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2692261) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(1.1671789) q[2];
rz(2.732318) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(-0.95612139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-2.9874711) q[0];
sx q[0];
rz(1.7896205) q[0];
rz(1.4714454) q[1];
sx q[1];
rz(-2.2189326) q[1];
sx q[1];
rz(-0.89554375) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.937718) q[0];
sx q[0];
rz(-0.26709891) q[0];
sx q[0];
rz(-0.0077203053) q[0];
rz(-pi) q[1];
rz(-0.15094916) q[2];
sx q[2];
rz(-0.47892919) q[2];
sx q[2];
rz(-2.4189309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91949192) q[1];
sx q[1];
rz(-0.89029145) q[1];
sx q[1];
rz(0.51503985) q[1];
rz(-1.9905521) q[3];
sx q[3];
rz(-1.354387) q[3];
sx q[3];
rz(-1.2866502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46513778) q[2];
sx q[2];
rz(-2.319591) q[2];
sx q[2];
rz(-1.6949534) q[2];
rz(-2.1220574) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(-2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119706) q[0];
sx q[0];
rz(-1.1638887) q[0];
sx q[0];
rz(3.1373366) q[0];
rz(-2.440522) q[1];
sx q[1];
rz(-0.509976) q[1];
sx q[1];
rz(3.1140936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.741318) q[0];
sx q[0];
rz(-1.7257462) q[0];
sx q[0];
rz(1.1688031) q[0];
rz(-pi) q[1];
rz(-1.4239935) q[2];
sx q[2];
rz(-0.73063421) q[2];
sx q[2];
rz(-0.50792009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3361275) q[1];
sx q[1];
rz(-1.256811) q[1];
sx q[1];
rz(-2.9538395) q[1];
rz(-1.9070088) q[3];
sx q[3];
rz(-1.5839521) q[3];
sx q[3];
rz(-2.1845078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6444401) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(-2.7437317) q[2];
rz(-1.0128939) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(-2.7329172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48552805) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(-0.84554607) q[0];
rz(-0.0088508765) q[1];
sx q[1];
rz(-0.84819853) q[1];
sx q[1];
rz(1.1525851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51868472) q[0];
sx q[0];
rz(-0.9810514) q[0];
sx q[0];
rz(2.6481762) q[0];
rz(-1.1197309) q[2];
sx q[2];
rz(-1.3374723) q[2];
sx q[2];
rz(-0.87750669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2611672) q[1];
sx q[1];
rz(-2.5486055) q[1];
sx q[1];
rz(-0.050131948) q[1];
x q[2];
rz(1.5043855) q[3];
sx q[3];
rz(-1.8838657) q[3];
sx q[3];
rz(-2.294189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74959603) q[2];
sx q[2];
rz(-1.0276724) q[2];
sx q[2];
rz(-1.2006302) q[2];
rz(2.6642753) q[3];
sx q[3];
rz(-0.47003191) q[3];
sx q[3];
rz(-0.70394713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93477997) q[0];
sx q[0];
rz(-0.88345695) q[0];
sx q[0];
rz(2.2373037) q[0];
rz(-2.1705056) q[1];
sx q[1];
rz(-1.9458408) q[1];
sx q[1];
rz(-0.28360525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0440031) q[0];
sx q[0];
rz(-1.6416802) q[0];
sx q[0];
rz(-1.5547836) q[0];
rz(-pi) q[1];
rz(0.798337) q[2];
sx q[2];
rz(-0.58272782) q[2];
sx q[2];
rz(-0.59205627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3963822) q[1];
sx q[1];
rz(-1.77715) q[1];
sx q[1];
rz(-2.2554805) q[1];
x q[2];
rz(-1.0002076) q[3];
sx q[3];
rz(-2.3688487) q[3];
sx q[3];
rz(-2.7955713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35563719) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(-2.9034485) q[2];
rz(2.1613878) q[3];
sx q[3];
rz(-1.1583867) q[3];
sx q[3];
rz(-0.75077209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.39764431) q[0];
sx q[0];
rz(-1.605796) q[0];
sx q[0];
rz(0.407298) q[0];
rz(0.40335718) q[1];
sx q[1];
rz(-0.74059161) q[1];
sx q[1];
rz(2.2743646) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9089389) q[0];
sx q[0];
rz(-2.5728658) q[0];
sx q[0];
rz(1.4599811) q[0];
x q[1];
rz(1.7795785) q[2];
sx q[2];
rz(-1.0241514) q[2];
sx q[2];
rz(-0.41266135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7379511) q[1];
sx q[1];
rz(-2.3372136) q[1];
sx q[1];
rz(-0.40108313) q[1];
rz(-pi) q[2];
rz(-0.12680291) q[3];
sx q[3];
rz(-1.276166) q[3];
sx q[3];
rz(-2.08088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0200218) q[2];
sx q[2];
rz(-2.3131804) q[2];
sx q[2];
rz(1.4309481) q[2];
rz(2.6089148) q[3];
sx q[3];
rz(-1.0360274) q[3];
sx q[3];
rz(1.7238341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.1581887) q[0];
sx q[0];
rz(-2.991365) q[0];
sx q[0];
rz(-2.0576553) q[0];
rz(3.0896507) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(0.025024978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0662224) q[0];
sx q[0];
rz(-2.0522723) q[0];
sx q[0];
rz(-1.5388778) q[0];
rz(-2.689471) q[2];
sx q[2];
rz(-2.2703301) q[2];
sx q[2];
rz(-1.6665065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1781749) q[1];
sx q[1];
rz(-1.1421848) q[1];
sx q[1];
rz(0.18842285) q[1];
rz(-pi) q[2];
rz(0.42633485) q[3];
sx q[3];
rz(-1.4407743) q[3];
sx q[3];
rz(0.015817964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0673125) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(2.32453) q[2];
rz(0.95064154) q[3];
sx q[3];
rz(-2.1248415) q[3];
sx q[3];
rz(-1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.10512146) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(2.5407319) q[0];
rz(3.1071682) q[1];
sx q[1];
rz(-1.8079146) q[1];
sx q[1];
rz(1.5404125) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87266028) q[0];
sx q[0];
rz(-1.0300228) q[0];
sx q[0];
rz(-2.2898104) q[0];
rz(-pi) q[1];
x q[1];
rz(2.41018) q[2];
sx q[2];
rz(-0.81899511) q[2];
sx q[2];
rz(0.42289823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29112962) q[1];
sx q[1];
rz(-2.7013868) q[1];
sx q[1];
rz(-3.090211) q[1];
rz(-0.65798379) q[3];
sx q[3];
rz(-0.65449981) q[3];
sx q[3];
rz(2.2006048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8704661) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(-1.6477443) q[2];
rz(-0.77914023) q[3];
sx q[3];
rz(-1.4804163) q[3];
sx q[3];
rz(-2.486865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7036024) q[0];
sx q[0];
rz(-2.0605189) q[0];
sx q[0];
rz(-0.78824747) q[0];
rz(-1.8986656) q[1];
sx q[1];
rz(-1.582076) q[1];
sx q[1];
rz(-2.8020614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57836048) q[0];
sx q[0];
rz(-1.4548911) q[0];
sx q[0];
rz(-1.7645592) q[0];
rz(-pi) q[1];
rz(-2.5794284) q[2];
sx q[2];
rz(-1.8149281) q[2];
sx q[2];
rz(-0.31329271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7701246) q[1];
sx q[1];
rz(-1.1684019) q[1];
sx q[1];
rz(-1.9103019) q[1];
x q[2];
rz(2.4838264) q[3];
sx q[3];
rz(-2.2493746) q[3];
sx q[3];
rz(-0.15746169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.01123151) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.8846903) q[2];
rz(-1.0907178) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(-0.91712657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85482368) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(1.0877892) q[0];
rz(-2.4913359) q[1];
sx q[1];
rz(-1.6842664) q[1];
sx q[1];
rz(0.021171721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41474671) q[0];
sx q[0];
rz(-1.124589) q[0];
sx q[0];
rz(-1.3394974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4707321) q[2];
sx q[2];
rz(-1.3952878) q[2];
sx q[2];
rz(-2.144475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0291527) q[1];
sx q[1];
rz(-0.8275367) q[1];
sx q[1];
rz(0.2874916) q[1];
rz(-pi) q[2];
rz(-0.043280799) q[3];
sx q[3];
rz(-2.1028404) q[3];
sx q[3];
rz(-0.30902853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1374958) q[2];
sx q[2];
rz(-1.7986412) q[2];
sx q[2];
rz(2.1583978) q[2];
rz(-0.99299562) q[3];
sx q[3];
rz(-1.1119305) q[3];
sx q[3];
rz(-2.907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642333) q[0];
sx q[0];
rz(-0.79240427) q[0];
sx q[0];
rz(-2.0808676) q[0];
rz(-2.0296774) q[1];
sx q[1];
rz(-0.30525515) q[1];
sx q[1];
rz(0.12259604) q[1];
rz(-0.23774679) q[2];
sx q[2];
rz(-1.100111) q[2];
sx q[2];
rz(-3.0636655) q[2];
rz(1.5817011) q[3];
sx q[3];
rz(-2.0937216) q[3];
sx q[3];
rz(1.0918319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
