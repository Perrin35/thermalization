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
rz(2.6574988) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(-3.0866403) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.84961) q[0];
sx q[0];
rz(-1.5395482) q[0];
sx q[0];
rz(-2.241075) q[0];
rz(-1.3208766) q[2];
sx q[2];
rz(-1.5139765) q[2];
sx q[2];
rz(-0.74733464) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7701227) q[1];
sx q[1];
rz(-0.16487637) q[1];
sx q[1];
rz(2.5046964) q[1];
rz(-pi) q[2];
rz(-0.3368042) q[3];
sx q[3];
rz(-2.1651377) q[3];
sx q[3];
rz(-1.6860733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83076465) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(2.8669299) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(-2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.2363385) q[1];
sx q[1];
rz(-2.0358548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92764054) q[0];
sx q[0];
rz(-1.1301148) q[0];
sx q[0];
rz(-1.3182314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2244338) q[2];
sx q[2];
rz(-2.0431402) q[2];
sx q[2];
rz(-1.9183967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2011021) q[1];
sx q[1];
rz(-3.1374212) q[1];
sx q[1];
rz(-0.68912403) q[1];
x q[2];
rz(-1.9883184) q[3];
sx q[3];
rz(-0.29051775) q[3];
sx q[3];
rz(0.49714303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4411321) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(1.0750809) q[2];
rz(0.69617802) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(-2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769128) q[0];
sx q[0];
rz(-1.541168) q[0];
sx q[0];
rz(1.6086786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9335494) q[2];
sx q[2];
rz(-3.0100076) q[2];
sx q[2];
rz(3.0818444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8113831) q[1];
sx q[1];
rz(-0.37322497) q[1];
sx q[1];
rz(-2.1887145) q[1];
x q[2];
rz(-0.9914753) q[3];
sx q[3];
rz(-1.0414755) q[3];
sx q[3];
rz(2.4059699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(-3.0510862) q[2];
rz(2.6835486) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(2.0335782) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(-2.2099387) q[0];
rz(0.16904198) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(-1.0394675) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19379481) q[0];
sx q[0];
rz(-2.373824) q[0];
sx q[0];
rz(0.34996535) q[0];
x q[1];
rz(-2.2957689) q[2];
sx q[2];
rz(-2.8081144) q[2];
sx q[2];
rz(0.84104702) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9586981) q[1];
sx q[1];
rz(-2.9258203) q[1];
sx q[1];
rz(-1.2879667) q[1];
rz(1.7839892) q[3];
sx q[3];
rz(-1.2925832) q[3];
sx q[3];
rz(-1.1602064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6197551) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(-1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8747044) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(0.12174363) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(-1.8121388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9864101) q[0];
sx q[0];
rz(-1.8768132) q[0];
sx q[0];
rz(-0.06432342) q[0];
x q[1];
rz(0.8829863) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(-0.6907874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5426104) q[1];
sx q[1];
rz(-0.81331454) q[1];
sx q[1];
rz(-2.5573362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0092486898) q[3];
sx q[3];
rz(-1.4917177) q[3];
sx q[3];
rz(1.6759863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(2.9902048) q[2];
rz(-0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720035) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(-2.0701011) q[0];
rz(0.53120652) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(0.62613097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1086639) q[0];
sx q[0];
rz(-1.5538006) q[0];
sx q[0];
rz(3.1368318) q[0];
x q[1];
rz(0.59711908) q[2];
sx q[2];
rz(-1.7749987) q[2];
sx q[2];
rz(-2.0300421) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5007361) q[1];
sx q[1];
rz(-1.2656414) q[1];
sx q[1];
rz(0.63668294) q[1];
rz(-2.6426396) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(-2.5753491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-1.0142856) q[2];
rz(1.2554393) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96100539) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(0.92051202) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-2.0106409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3823272) q[0];
sx q[0];
rz(-0.59774071) q[0];
sx q[0];
rz(-0.010058479) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55552395) q[2];
sx q[2];
rz(-2.0709403) q[2];
sx q[2];
rz(-1.9490567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5798545) q[1];
sx q[1];
rz(-0.67719141) q[1];
sx q[1];
rz(1.3586292) q[1];
rz(2.6252296) q[3];
sx q[3];
rz(-1.5011906) q[3];
sx q[3];
rz(-2.7854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(2.9583904) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-1.1387775) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278397) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(-0.94789061) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(1.9416169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5771709) q[0];
sx q[0];
rz(-1.0438215) q[0];
sx q[0];
rz(-2.0354664) q[0];
x q[1];
rz(-2.3340204) q[2];
sx q[2];
rz(-1.0243197) q[2];
sx q[2];
rz(-1.228491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34741286) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(-1.7296687) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.348749) q[2];
rz(-1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0272442) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(2.0674904) q[0];
rz(2.7117924) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(2.6780186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6493452) q[0];
sx q[0];
rz(-1.6657636) q[0];
sx q[0];
rz(-1.3819206) q[0];
rz(-0.65166574) q[2];
sx q[2];
rz(-2.8164688) q[2];
sx q[2];
rz(0.59949694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.025109865) q[1];
sx q[1];
rz(-1.9466725) q[1];
sx q[1];
rz(0.89565887) q[1];
rz(-pi) q[2];
rz(0.55863278) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(0.35364756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(2.3331433) q[2];
rz(2.6570053) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(-0.038381902) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(0.29676944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0050186) q[0];
sx q[0];
rz(-1.2941735) q[0];
sx q[0];
rz(1.0199976) q[0];
rz(-pi) q[1];
rz(-0.44906868) q[2];
sx q[2];
rz(-1.9757604) q[2];
sx q[2];
rz(1.8835889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44904199) q[1];
sx q[1];
rz(-1.7925279) q[1];
sx q[1];
rz(-1.4936844) q[1];
x q[2];
rz(0.88480437) q[3];
sx q[3];
rz(-1.8821041) q[3];
sx q[3];
rz(1.6303568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9639637) q[2];
sx q[2];
rz(-2.0202426) q[2];
sx q[2];
rz(0.56619823) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(1.1894777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(-0.26168564) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(0.08200866) q[2];
sx q[2];
rz(-1.5463943) q[2];
sx q[2];
rz(0.80873185) q[2];
rz(0.86033173) q[3];
sx q[3];
rz(-2.0900149) q[3];
sx q[3];
rz(2.7505977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
