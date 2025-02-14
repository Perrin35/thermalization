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
rz(-2.1343159) q[0];
sx q[0];
rz(-2.4957823) q[0];
sx q[0];
rz(2.7752152) q[0];
rz(0.74908787) q[1];
sx q[1];
rz(4.7685342) q[1];
sx q[1];
rz(9.2821477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34829545) q[0];
sx q[0];
rz(-2.1014338) q[0];
sx q[0];
rz(-2.2133618) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22444522) q[2];
sx q[2];
rz(-1.7068752) q[2];
sx q[2];
rz(-2.4405757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9318568) q[1];
sx q[1];
rz(-1.0967347) q[1];
sx q[1];
rz(-0.13412383) q[1];
x q[2];
rz(-2.6680115) q[3];
sx q[3];
rz(-1.2788805) q[3];
sx q[3];
rz(2.0278668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7819405) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(-1.5521607) q[2];
rz(1.3945256) q[3];
sx q[3];
rz(-0.37140578) q[3];
sx q[3];
rz(0.68215218) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8438016) q[0];
sx q[0];
rz(-2.5318635) q[0];
sx q[0];
rz(2.8363256) q[0];
rz(1.3097395) q[1];
sx q[1];
rz(-2.7204308) q[1];
sx q[1];
rz(0.052068204) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0355988) q[0];
sx q[0];
rz(-1.4609504) q[0];
sx q[0];
rz(-2.99361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8709677) q[2];
sx q[2];
rz(-1.9294293) q[2];
sx q[2];
rz(1.4305565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4365719) q[1];
sx q[1];
rz(-1.3051045) q[1];
sx q[1];
rz(0.17036856) q[1];
x q[2];
rz(-2.2388715) q[3];
sx q[3];
rz(-1.4666712) q[3];
sx q[3];
rz(1.6834761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0569968) q[2];
sx q[2];
rz(-1.8085542) q[2];
sx q[2];
rz(1.4196654) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-0.82908583) q[3];
sx q[3];
rz(1.4066633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.1664355) q[0];
sx q[0];
rz(-1.9571914) q[0];
sx q[0];
rz(-1.3982406) q[0];
rz(-0.94683975) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(1.85087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7311659) q[0];
sx q[0];
rz(-1.7058347) q[0];
sx q[0];
rz(-3.0655087) q[0];
rz(-pi) q[1];
rz(-0.78611908) q[2];
sx q[2];
rz(-1.9232456) q[2];
sx q[2];
rz(-2.4967683) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5980893) q[1];
sx q[1];
rz(-0.072797983) q[1];
sx q[1];
rz(-1.1612215) q[1];
x q[2];
rz(-2.6638899) q[3];
sx q[3];
rz(-1.3906428) q[3];
sx q[3];
rz(-2.2733062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1281841) q[2];
sx q[2];
rz(-0.42542294) q[2];
sx q[2];
rz(-0.95995963) q[2];
rz(2.8175682) q[3];
sx q[3];
rz(-2.6257381) q[3];
sx q[3];
rz(-2.7014151) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0078761) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(-0.20009759) q[0];
rz(0.61087418) q[1];
sx q[1];
rz(-2.4433544) q[1];
sx q[1];
rz(-2.3470338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88761974) q[0];
sx q[0];
rz(-1.4591503) q[0];
sx q[0];
rz(-0.15038862) q[0];
rz(-pi) q[1];
rz(-2.2438917) q[2];
sx q[2];
rz(-1.4294927) q[2];
sx q[2];
rz(0.92800498) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3172052) q[1];
sx q[1];
rz(-0.35046589) q[1];
sx q[1];
rz(2.0262655) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95501542) q[3];
sx q[3];
rz(-0.94388591) q[3];
sx q[3];
rz(-2.0364498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8841298) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(-1.5696611) q[2];
rz(-1.3563159) q[3];
sx q[3];
rz(-1.9808199) q[3];
sx q[3];
rz(3.0089488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86303478) q[0];
sx q[0];
rz(-0.25535169) q[0];
sx q[0];
rz(2.417946) q[0];
rz(3.0408995) q[1];
sx q[1];
rz(-1.0483402) q[1];
sx q[1];
rz(3.0752799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6479939) q[0];
sx q[0];
rz(-2.1382522) q[0];
sx q[0];
rz(0.074851224) q[0];
rz(-2.2411614) q[2];
sx q[2];
rz(-0.77610972) q[2];
sx q[2];
rz(2.4845207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7620341) q[1];
sx q[1];
rz(-2.4707816) q[1];
sx q[1];
rz(1.8422153) q[1];
x q[2];
rz(0.67541404) q[3];
sx q[3];
rz(-2.6468122) q[3];
sx q[3];
rz(-2.8515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.940332) q[2];
sx q[2];
rz(-1.9278434) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(2.7471733) q[3];
sx q[3];
rz(-2.4557178) q[3];
sx q[3];
rz(-0.39943892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.047693096) q[0];
sx q[0];
rz(-0.91966367) q[0];
sx q[0];
rz(-3.1374875) q[0];
rz(-0.63402367) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(-0.9563458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078236899) q[0];
sx q[0];
rz(-0.18785297) q[0];
sx q[0];
rz(0.22042033) q[0];
x q[1];
rz(2.6317503) q[2];
sx q[2];
rz(-1.9467453) q[2];
sx q[2];
rz(-2.816566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2507657) q[1];
sx q[1];
rz(-0.27784744) q[1];
sx q[1];
rz(-0.020093159) q[1];
rz(-pi) q[2];
rz(-1.6878674) q[3];
sx q[3];
rz(-1.2244976) q[3];
sx q[3];
rz(-2.9786144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.52028209) q[2];
sx q[2];
rz(-1.4352398) q[2];
sx q[2];
rz(-2.3069646) q[2];
rz(0.37072119) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(2.4247775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653559) q[0];
sx q[0];
rz(-3.1270202) q[0];
sx q[0];
rz(-2.6915349) q[0];
rz(0.4484446) q[1];
sx q[1];
rz(-0.83201718) q[1];
sx q[1];
rz(-0.73743302) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56313228) q[0];
sx q[0];
rz(-1.9942584) q[0];
sx q[0];
rz(0.20276258) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4294523) q[2];
sx q[2];
rz(-1.5307883) q[2];
sx q[2];
rz(-1.227898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0998675) q[1];
sx q[1];
rz(-1.2507572) q[1];
sx q[1];
rz(-2.7766224) q[1];
x q[2];
rz(0.78011765) q[3];
sx q[3];
rz(-1.3778049) q[3];
sx q[3];
rz(1.5548201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2927148) q[2];
sx q[2];
rz(-0.60056168) q[2];
sx q[2];
rz(-0.71504492) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.4628937) q[3];
sx q[3];
rz(1.7432632) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364863) q[0];
sx q[0];
rz(-0.54623258) q[0];
sx q[0];
rz(2.7832094) q[0];
rz(-0.37250039) q[1];
sx q[1];
rz(-1.7820216) q[1];
sx q[1];
rz(0.43982664) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43482237) q[0];
sx q[0];
rz(-2.3447038) q[0];
sx q[0];
rz(3.0838941) q[0];
x q[1];
rz(0.40435515) q[2];
sx q[2];
rz(-1.99843) q[2];
sx q[2];
rz(-0.15844395) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29319841) q[1];
sx q[1];
rz(-1.3064791) q[1];
sx q[1];
rz(-1.3879723) q[1];
x q[2];
rz(-2.8183455) q[3];
sx q[3];
rz(-1.6222937) q[3];
sx q[3];
rz(-0.53158376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0899352) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(0.70866054) q[2];
rz(2.8900201) q[3];
sx q[3];
rz(-2.1422062) q[3];
sx q[3];
rz(0.034041762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8468903) q[0];
sx q[0];
rz(-2.1069694) q[0];
sx q[0];
rz(2.5501472) q[0];
rz(3.1239608) q[1];
sx q[1];
rz(-2.089274) q[1];
sx q[1];
rz(-3.0406521) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288947) q[0];
sx q[0];
rz(-0.34358812) q[0];
sx q[0];
rz(1.1593007) q[0];
rz(-pi) q[1];
rz(0.63669651) q[2];
sx q[2];
rz(-1.5392188) q[2];
sx q[2];
rz(-2.4898138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7314607) q[1];
sx q[1];
rz(-1.6900151) q[1];
sx q[1];
rz(-1.8006667) q[1];
rz(-pi) q[2];
rz(0.89102192) q[3];
sx q[3];
rz(-1.4263065) q[3];
sx q[3];
rz(1.4308892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65323222) q[2];
sx q[2];
rz(-2.4182352) q[2];
sx q[2];
rz(-1.4867268) q[2];
rz(-0.15445736) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(-0.36623335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3514997) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(2.067814) q[0];
rz(-2.0692661) q[1];
sx q[1];
rz(-0.25389478) q[1];
sx q[1];
rz(0.059018746) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8813777) q[0];
sx q[0];
rz(-0.20199595) q[0];
sx q[0];
rz(2.9750573) q[0];
rz(2.8779936) q[2];
sx q[2];
rz(-1.0754183) q[2];
sx q[2];
rz(-1.4857298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7714951) q[1];
sx q[1];
rz(-2.1019249) q[1];
sx q[1];
rz(2.4350776) q[1];
x q[2];
rz(2.2667472) q[3];
sx q[3];
rz(-0.84934649) q[3];
sx q[3];
rz(-1.9375999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8129639) q[2];
sx q[2];
rz(-1.2920516) q[2];
sx q[2];
rz(-0.22895075) q[2];
rz(-1.1936584) q[3];
sx q[3];
rz(-2.8173859) q[3];
sx q[3];
rz(-1.6875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345632) q[0];
sx q[0];
rz(-2.3352191) q[0];
sx q[0];
rz(2.4051608) q[0];
rz(1.3855343) q[1];
sx q[1];
rz(-1.7693188) q[1];
sx q[1];
rz(1.7973695) q[1];
rz(1.6177542) q[2];
sx q[2];
rz(-2.6018819) q[2];
sx q[2];
rz(0.56503208) q[2];
rz(1.7731316) q[3];
sx q[3];
rz(-2.1740395) q[3];
sx q[3];
rz(-2.5185829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
