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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(-0.48409387) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(-3.0866403) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8380175) q[0];
sx q[0];
rz(-2.2406881) q[0];
sx q[0];
rz(3.1017257) q[0];
x q[1];
rz(-1.3447433) q[2];
sx q[2];
rz(-2.885427) q[2];
sx q[2];
rz(-1.0423755) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72803264) q[1];
sx q[1];
rz(-1.7031341) q[1];
sx q[1];
rz(-1.669426) q[1];
rz(2.0255768) q[3];
sx q[3];
rz(-0.6729799) q[3];
sx q[3];
rz(2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.310828) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(0.27466276) q[2];
rz(-0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(2.8680475) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65983588) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(2.0358548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53351346) q[0];
sx q[0];
rz(-1.7987804) q[0];
sx q[0];
rz(2.6883459) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91715889) q[2];
sx q[2];
rz(-2.0431402) q[2];
sx q[2];
rz(1.2231959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0827786) q[1];
sx q[1];
rz(-1.5681439) q[1];
sx q[1];
rz(3.1383731) q[1];
x q[2];
rz(-0.12064528) q[3];
sx q[3];
rz(-1.835726) q[3];
sx q[3];
rz(0.93075965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(-1.0750809) q[2];
rz(0.69617802) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(0.16304326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27580801) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19720896) q[0];
sx q[0];
rz(-1.5329307) q[0];
sx q[0];
rz(0.029649563) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0946629) q[2];
sx q[2];
rz(-1.4478193) q[2];
sx q[2];
rz(2.8356981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8113831) q[1];
sx q[1];
rz(-0.37322497) q[1];
sx q[1];
rz(-0.95287816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3893395) q[3];
sx q[3];
rz(-2.3779388) q[3];
sx q[3];
rz(1.6490761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8850024) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(2.1021252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19379481) q[0];
sx q[0];
rz(-2.373824) q[0];
sx q[0];
rz(-0.34996535) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8245071) q[2];
sx q[2];
rz(-1.3519962) q[2];
sx q[2];
rz(-1.7148866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1828945) q[1];
sx q[1];
rz(-2.9258203) q[1];
sx q[1];
rz(-1.853626) q[1];
rz(-pi) q[2];
rz(1.3576034) q[3];
sx q[3];
rz(-1.8490095) q[3];
sx q[3];
rz(1.9813862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(0.49631611) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(0.12174363) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(1.3294539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861391) q[0];
sx q[0];
rz(-0.3124961) q[0];
sx q[0];
rz(-1.7715095) q[0];
rz(0.41556032) q[2];
sx q[2];
rz(-1.1141127) q[2];
sx q[2];
rz(-1.6650852) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7763627) q[1];
sx q[1];
rz(-0.91971469) q[1];
sx q[1];
rz(1.0427703) q[1];
rz(1.4546118) q[3];
sx q[3];
rz(-3.0619762) q[3];
sx q[3];
rz(-1.5594359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(-0.15138781) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695892) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(0.53120652) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(2.5154617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75979489) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(1.8438898) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3254204) q[2];
sx q[2];
rz(-0.98773709) q[2];
sx q[2];
rz(-2.8193605) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6408566) q[1];
sx q[1];
rz(-1.8759512) q[1];
sx q[1];
rz(-0.63668294) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(-2.5753491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-1.0142856) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(-1.0534508) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(-0.92051202) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(-1.1309518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3217472) q[0];
sx q[0];
rz(-1.5764569) q[0];
sx q[0];
rz(-2.5438755) q[0];
rz(0.80318309) q[2];
sx q[2];
rz(-0.72942643) q[2];
sx q[2];
rz(0.27952295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56173813) q[1];
sx q[1];
rz(-0.67719141) q[1];
sx q[1];
rz(1.3586292) q[1];
x q[2];
rz(2.6252296) q[3];
sx q[3];
rz(-1.6404021) q[3];
sx q[3];
rz(-0.35610896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7941234) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(-0.40880173) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(-1.1387775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278397) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(2.193702) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(-1.9416169) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5644218) q[0];
sx q[0];
rz(-2.0977712) q[0];
sx q[0];
rz(-1.1061263) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4418996) q[2];
sx q[2];
rz(-2.202575) q[2];
sx q[2];
rz(0.80365411) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34741286) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.7296687) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7539976) q[3];
sx q[3];
rz(-1.6699413) q[3];
sx q[3];
rz(0.20618901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.348749) q[2];
rz(-1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.1143484) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(0.46357402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6025507) q[0];
sx q[0];
rz(-2.9304404) q[0];
sx q[0];
rz(2.0402914) q[0];
rz(-pi) q[1];
rz(1.3691291) q[2];
sx q[2];
rz(-1.3140162) q[2];
sx q[2];
rz(-0.078291206) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.025109865) q[1];
sx q[1];
rz(-1.9466725) q[1];
sx q[1];
rz(2.2459338) q[1];
rz(-1.660668) q[3];
sx q[3];
rz(-1.0139795) q[3];
sx q[3];
rz(1.1695605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8672436) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(-0.80844936) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(-2.8448232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1574157) q[0];
sx q[0];
rz(-2.5317051) q[0];
sx q[0];
rz(-1.0737674) q[0];
rz(2.0149258) q[2];
sx q[2];
rz(-1.9812366) q[2];
sx q[2];
rz(0.12516147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0368288) q[1];
sx q[1];
rz(-1.6460168) q[1];
sx q[1];
rz(-0.22237088) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39408306) q[3];
sx q[3];
rz(-2.218045) q[3];
sx q[3];
rz(2.955472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9639637) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(-0.56619823) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(1.1894777) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2531256) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(0.26168564) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(1.5952806) q[2];
sx q[2];
rz(-1.6527805) q[2];
sx q[2];
rz(2.3775227) q[2];
rz(-0.85121831) q[3];
sx q[3];
rz(-0.85243445) q[3];
sx q[3];
rz(1.7029521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
