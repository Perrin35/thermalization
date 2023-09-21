OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(-2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18626285) q[0];
sx q[0];
rz(-1.7832527) q[0];
sx q[0];
rz(-2.0583378) q[0];
rz(-0.21284717) q[2];
sx q[2];
rz(-0.93570645) q[2];
sx q[2];
rz(-2.0095306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72370428) q[1];
sx q[1];
rz(-2.1065518) q[1];
sx q[1];
rz(-2.8292311) q[1];
x q[2];
rz(0.25986259) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5713455) q[0];
sx q[0];
rz(-2.3803108) q[0];
sx q[0];
rz(1.135457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5580514) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(1.5577424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9559905) q[1];
sx q[1];
rz(-2.5163109) q[1];
sx q[1];
rz(0.76693265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63668164) q[3];
sx q[3];
rz(-0.59906206) q[3];
sx q[3];
rz(2.3707795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3765592) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(-0.88699938) q[0];
rz(-pi) q[1];
rz(1.4372196) q[2];
sx q[2];
rz(-1.8250416) q[2];
sx q[2];
rz(-0.1453407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1287071) q[1];
sx q[1];
rz(-1.7696847) q[1];
sx q[1];
rz(-0.30602869) q[1];
x q[2];
rz(2.4967381) q[3];
sx q[3];
rz(-1.2554902) q[3];
sx q[3];
rz(0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47953654) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(-2.7583073) q[0];
rz(-pi) q[1];
rz(-1.6323339) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(2.5663944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5220118) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(0.79230688) q[1];
rz(-pi) q[2];
rz(-1.0158402) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(-0.89158981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(-3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-2.1972426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7992295) q[0];
sx q[0];
rz(-0.21393299) q[0];
sx q[0];
rz(-1.3571204) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3659533) q[2];
sx q[2];
rz(-1.2795942) q[2];
sx q[2];
rz(-1.8679384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22630616) q[1];
sx q[1];
rz(-1.5729135) q[1];
sx q[1];
rz(-1.6318984) q[1];
x q[2];
rz(-2.7086908) q[3];
sx q[3];
rz(-1.3372984) q[3];
sx q[3];
rz(1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(3.1392858) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(3.086673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391101) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(-1.8399747) q[0];
rz(1.247585) q[2];
sx q[2];
rz(-1.5885457) q[2];
sx q[2];
rz(-0.36349597) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37327787) q[1];
sx q[1];
rz(-2.1149966) q[1];
sx q[1];
rz(-2.0139704) q[1];
rz(-pi) q[2];
rz(-0.98723282) q[3];
sx q[3];
rz(-1.1757441) q[3];
sx q[3];
rz(-0.46056718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.6761557) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1032573) q[0];
sx q[0];
rz(-3.0568125) q[0];
sx q[0];
rz(0.8707365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1548642) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(1.8679801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42156223) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(-2.6471495) q[1];
rz(-pi) q[2];
rz(2.1479285) q[3];
sx q[3];
rz(-0.93615195) q[3];
sx q[3];
rz(-1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(1.0459082) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6408765) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(1.4027632) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(0.47972958) q[1];
rz(-pi) q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(2.382747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(0.018652648) q[0];
rz(-0.27300948) q[2];
sx q[2];
rz(-0.62676478) q[2];
sx q[2];
rz(0.11944709) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5938877) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(-1.3681075) q[1];
x q[2];
rz(0.014563668) q[3];
sx q[3];
rz(-0.90523883) q[3];
sx q[3];
rz(1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41994914) q[0];
sx q[0];
rz(-0.83166612) q[0];
sx q[0];
rz(1.182611) q[0];
rz(-pi) q[1];
rz(1.2953193) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(1.7737349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1249485) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(-1.8925341) q[1];
rz(-pi) q[2];
rz(-0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(-1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(1.7667608) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(1.0999023) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
