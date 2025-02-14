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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(1.1991062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48079607) q[0];
sx q[0];
rz(-1.4624339) q[0];
sx q[0];
rz(-1.5078791) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4680696) q[2];
sx q[2];
rz(-1.4269514) q[2];
sx q[2];
rz(-2.0604482) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71237129) q[1];
sx q[1];
rz(-1.0069798) q[1];
sx q[1];
rz(2.9122206) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5254123) q[3];
sx q[3];
rz(-2.1839704) q[3];
sx q[3];
rz(0.67894906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5045515) q[2];
sx q[2];
rz(-1.4013638) q[2];
sx q[2];
rz(1.6999647) q[2];
rz(-2.1291034) q[3];
sx q[3];
rz(-1.1463405) q[3];
sx q[3];
rz(-2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3016894) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(-1.0294234) q[0];
rz(1.3599716) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(0.50672466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.473712) q[0];
sx q[0];
rz(-3.0356044) q[0];
sx q[0];
rz(1.7539658) q[0];
rz(-pi) q[1];
rz(-2.7756491) q[2];
sx q[2];
rz(-0.90922132) q[2];
sx q[2];
rz(0.55398527) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9190127) q[1];
sx q[1];
rz(-1.6363385) q[1];
sx q[1];
rz(-0.99187851) q[1];
rz(0.84552879) q[3];
sx q[3];
rz(-1.7983899) q[3];
sx q[3];
rz(1.1448154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3158675) q[2];
sx q[2];
rz(-0.59810144) q[2];
sx q[2];
rz(0.21710795) q[2];
rz(2.121296) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8751136) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-0.76817051) q[0];
rz(-2.8504596) q[1];
sx q[1];
rz(-1.365064) q[1];
sx q[1];
rz(-1.1870144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6757557) q[0];
sx q[0];
rz(-1.4872941) q[0];
sx q[0];
rz(2.9299987) q[0];
x q[1];
rz(0.044069604) q[2];
sx q[2];
rz(-1.4975093) q[2];
sx q[2];
rz(1.2512887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95988454) q[1];
sx q[1];
rz(-1.925029) q[1];
sx q[1];
rz(-1.6003057) q[1];
x q[2];
rz(0.95245016) q[3];
sx q[3];
rz(-2.420331) q[3];
sx q[3];
rz(0.14474584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8698296) q[2];
sx q[2];
rz(-0.84438476) q[2];
sx q[2];
rz(2.1738906) q[2];
rz(2.9076231) q[3];
sx q[3];
rz(-0.7015737) q[3];
sx q[3];
rz(2.7847737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135603) q[0];
sx q[0];
rz(-1.5910633) q[0];
sx q[0];
rz(1.0796984) q[0];
rz(0.30403852) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(-1.3160204) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2000421) q[0];
sx q[0];
rz(-2.594875) q[0];
sx q[0];
rz(2.9018875) q[0];
rz(-pi) q[1];
x q[1];
rz(0.082985254) q[2];
sx q[2];
rz(-2.1647762) q[2];
sx q[2];
rz(1.9165438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41982728) q[1];
sx q[1];
rz(-0.98079762) q[1];
sx q[1];
rz(0.71728398) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7075482) q[3];
sx q[3];
rz(-2.7117808) q[3];
sx q[3];
rz(0.97625247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5217846) q[2];
sx q[2];
rz(-1.9079756) q[2];
sx q[2];
rz(-3.0002777) q[2];
rz(2.5610793) q[3];
sx q[3];
rz(-1.6655191) q[3];
sx q[3];
rz(-0.22835246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2572131) q[0];
sx q[0];
rz(-0.79638201) q[0];
sx q[0];
rz(1.6656531) q[0];
rz(2.2085704) q[1];
sx q[1];
rz(-0.7111744) q[1];
sx q[1];
rz(1.7500386) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45804241) q[0];
sx q[0];
rz(-1.6876966) q[0];
sx q[0];
rz(2.3085672) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33220187) q[2];
sx q[2];
rz(-0.4449639) q[2];
sx q[2];
rz(1.5277758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36223026) q[1];
sx q[1];
rz(-2.0915151) q[1];
sx q[1];
rz(-3.1073276) q[1];
rz(-1.0443893) q[3];
sx q[3];
rz(-2.5431051) q[3];
sx q[3];
rz(0.55803669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89620227) q[2];
sx q[2];
rz(-0.37084493) q[2];
sx q[2];
rz(2.89213) q[2];
rz(-1.6438515) q[3];
sx q[3];
rz(-1.4062873) q[3];
sx q[3];
rz(-2.7264285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45667085) q[0];
sx q[0];
rz(-0.49600729) q[0];
sx q[0];
rz(-1.3859092) q[0];
rz(-1.2617525) q[1];
sx q[1];
rz(-1.2077786) q[1];
sx q[1];
rz(-0.52880803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.38976) q[0];
sx q[0];
rz(-2.6963137) q[0];
sx q[0];
rz(-1.5543429) q[0];
rz(-pi) q[1];
rz(-2.4086921) q[2];
sx q[2];
rz(-0.7846047) q[2];
sx q[2];
rz(0.32151383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13695194) q[1];
sx q[1];
rz(-1.5037422) q[1];
sx q[1];
rz(-2.2100558) q[1];
x q[2];
rz(-0.68499039) q[3];
sx q[3];
rz(-1.6719975) q[3];
sx q[3];
rz(1.8678968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95278198) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(1.7702276) q[2];
rz(0.2025226) q[3];
sx q[3];
rz(-1.6731508) q[3];
sx q[3];
rz(-1.1892345) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33223575) q[0];
sx q[0];
rz(-1.3487331) q[0];
sx q[0];
rz(-0.15383823) q[0];
rz(-2.4065252) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(-0.14911266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64407883) q[0];
sx q[0];
rz(-0.91877684) q[0];
sx q[0];
rz(-0.10494167) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0015099) q[2];
sx q[2];
rz(-0.39719279) q[2];
sx q[2];
rz(0.45969648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2231517) q[1];
sx q[1];
rz(-2.0479124) q[1];
sx q[1];
rz(2.3972307) q[1];
rz(-pi) q[2];
rz(0.090773067) q[3];
sx q[3];
rz(-2.0730264) q[3];
sx q[3];
rz(1.9050364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4601124) q[2];
sx q[2];
rz(-1.511938) q[2];
sx q[2];
rz(0.45664772) q[2];
rz(-1.418142) q[3];
sx q[3];
rz(-2.2599615) q[3];
sx q[3];
rz(-0.68896967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6316471) q[0];
sx q[0];
rz(-0.25196415) q[0];
sx q[0];
rz(0.052635996) q[0];
rz(-1.3098199) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(-1.2059258) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9034957) q[0];
sx q[0];
rz(-2.4795929) q[0];
sx q[0];
rz(2.2308626) q[0];
x q[1];
rz(1.8032372) q[2];
sx q[2];
rz(-2.1191868) q[2];
sx q[2];
rz(-0.91969925) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0830906) q[1];
sx q[1];
rz(-1.41377) q[1];
sx q[1];
rz(0.67879063) q[1];
x q[2];
rz(2.1931529) q[3];
sx q[3];
rz(-1.3788169) q[3];
sx q[3];
rz(-2.6485659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2974818) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(0.55997509) q[2];
rz(-1.0567793) q[3];
sx q[3];
rz(-2.4323075) q[3];
sx q[3];
rz(-0.81282508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.63447222) q[0];
sx q[0];
rz(-1.3993323) q[0];
sx q[0];
rz(-1.5768453) q[0];
rz(1.5693846) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(2.3326468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2574111) q[0];
sx q[0];
rz(-2.0657259) q[0];
sx q[0];
rz(-1.6135501) q[0];
x q[1];
rz(-2.3303512) q[2];
sx q[2];
rz(-1.8643856) q[2];
sx q[2];
rz(-0.21101235) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.28918) q[1];
sx q[1];
rz(-2.1816263) q[1];
sx q[1];
rz(-3.1256366) q[1];
rz(-pi) q[2];
rz(-1.393717) q[3];
sx q[3];
rz(-1.5054323) q[3];
sx q[3];
rz(-3.1013427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6605777) q[2];
sx q[2];
rz(-0.95260859) q[2];
sx q[2];
rz(-3.0926404) q[2];
rz(-1.281721) q[3];
sx q[3];
rz(-0.47544026) q[3];
sx q[3];
rz(-1.4648645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.319741) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(1.9191746) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.9404575) q[1];
sx q[1];
rz(0.45752057) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37402651) q[0];
sx q[0];
rz(-1.4806642) q[0];
sx q[0];
rz(-0.73390168) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2234688) q[2];
sx q[2];
rz(-2.5312573) q[2];
sx q[2];
rz(-0.68896919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19472518) q[1];
sx q[1];
rz(-1.3683142) q[1];
sx q[1];
rz(2.7848886) q[1];
rz(-pi) q[2];
rz(-0.49623747) q[3];
sx q[3];
rz(-0.75025815) q[3];
sx q[3];
rz(-1.539278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3931302) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(2.1012696) q[2];
rz(2.6535502) q[3];
sx q[3];
rz(-1.7010242) q[3];
sx q[3];
rz(0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671065) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(-1.8472916) q[1];
sx q[1];
rz(-2.558567) q[1];
sx q[1];
rz(-0.8898215) q[1];
rz(-2.0001844) q[2];
sx q[2];
rz(-2.5038237) q[2];
sx q[2];
rz(0.43546168) q[2];
rz(-1.9593976) q[3];
sx q[3];
rz(-0.91872707) q[3];
sx q[3];
rz(1.6891458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
