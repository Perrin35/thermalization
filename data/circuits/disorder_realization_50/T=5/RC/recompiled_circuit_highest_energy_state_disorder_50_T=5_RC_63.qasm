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
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(1.1573428) q[0];
rz(-2.3565489) q[1];
sx q[1];
rz(-2.4554689) q[1];
sx q[1];
rz(0.52678984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136381) q[0];
sx q[0];
rz(-2.1305973) q[0];
sx q[0];
rz(0.35388057) q[0];
rz(-2.3700222) q[2];
sx q[2];
rz(-0.89961806) q[2];
sx q[2];
rz(3.0085466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1810468) q[1];
sx q[1];
rz(-1.6072453) q[1];
sx q[1];
rz(-0.44238449) q[1];
rz(-pi) q[2];
rz(-0.029954682) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-2.4966979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36246768) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(0.15164068) q[2];
rz(2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4769984) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(2.3240996) q[0];
rz(3.0290161) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(-0.89961019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(-1.3221571) q[0];
rz(-pi) q[1];
rz(-0.95006534) q[2];
sx q[2];
rz(-2.143444) q[2];
sx q[2];
rz(0.060163035) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6625329) q[1];
sx q[1];
rz(-2.3182437) q[1];
sx q[1];
rz(1.4940628) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95590653) q[3];
sx q[3];
rz(-0.77707505) q[3];
sx q[3];
rz(1.2082548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67109913) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(2.1136843) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.4772011) q[3];
sx q[3];
rz(3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(-1.8434803) q[1];
sx q[1];
rz(-1.4793414) q[1];
sx q[1];
rz(1.6568291) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910716) q[0];
sx q[0];
rz(-1.5245887) q[0];
sx q[0];
rz(0.025006983) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12642352) q[2];
sx q[2];
rz(-1.7806541) q[2];
sx q[2];
rz(-1.017638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72276593) q[1];
sx q[1];
rz(-1.9120941) q[1];
sx q[1];
rz(1.4272293) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9336352) q[3];
sx q[3];
rz(-1.4936183) q[3];
sx q[3];
rz(3.0953323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16126157) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(-2.4397591) q[2];
rz(0.069132239) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.055701) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(2.5162146) q[0];
rz(2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(-1.8249003) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33915658) q[0];
sx q[0];
rz(-0.56471741) q[0];
sx q[0];
rz(-1.3012278) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8181192) q[2];
sx q[2];
rz(-1.6626366) q[2];
sx q[2];
rz(-0.30486456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55344501) q[1];
sx q[1];
rz(-0.7281174) q[1];
sx q[1];
rz(2.8690228) q[1];
rz(-pi) q[2];
rz(-2.1573477) q[3];
sx q[3];
rz(-0.86185019) q[3];
sx q[3];
rz(1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40654287) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(3.1114846) q[2];
rz(-0.41664577) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(-3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3629214) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.2177421) q[0];
rz(0.22449224) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(-0.40103689) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56479708) q[0];
sx q[0];
rz(-2.5677997) q[0];
sx q[0];
rz(2.3032715) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6729092) q[2];
sx q[2];
rz(-0.71063738) q[2];
sx q[2];
rz(0.24220322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2485473) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(-2.9045312) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9791361) q[3];
sx q[3];
rz(-1.6704644) q[3];
sx q[3];
rz(-1.398744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(0.92042813) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4949263) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(1.0194417) q[0];
rz(-1.1722209) q[1];
sx q[1];
rz(-1.9688789) q[1];
sx q[1];
rz(-1.2219465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1555699) q[0];
sx q[0];
rz(-1.9180074) q[0];
sx q[0];
rz(-2.272764) q[0];
rz(-pi) q[1];
rz(3.1062791) q[2];
sx q[2];
rz(-1.0562293) q[2];
sx q[2];
rz(3.0701612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1091369) q[1];
sx q[1];
rz(-1.0278406) q[1];
sx q[1];
rz(1.2486267) q[1];
x q[2];
rz(-0.29767146) q[3];
sx q[3];
rz(-1.0165914) q[3];
sx q[3];
rz(-2.8870074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-1.0931724) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(-2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.2375803) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(-1.4599266) q[0];
rz(1.6319252) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(-3.0622838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722141) q[0];
sx q[0];
rz(-1.3339808) q[0];
sx q[0];
rz(0.083061465) q[0];
rz(-pi) q[1];
rz(-2.5517188) q[2];
sx q[2];
rz(-0.24352077) q[2];
sx q[2];
rz(-0.68984725) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0876336) q[1];
sx q[1];
rz(-1.1529536) q[1];
sx q[1];
rz(1.8709976) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71770151) q[3];
sx q[3];
rz(-1.4677748) q[3];
sx q[3];
rz(0.75005248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1377533) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(0.4129146) q[2];
rz(-1.2952992) q[3];
sx q[3];
rz(-0.92419878) q[3];
sx q[3];
rz(-2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(-0.28537634) q[0];
rz(3.0889619) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(2.9579128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3899468) q[0];
sx q[0];
rz(-0.11706287) q[0];
sx q[0];
rz(-1.879891) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48641522) q[2];
sx q[2];
rz(-0.32054311) q[2];
sx q[2];
rz(-0.17580168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7990098) q[1];
sx q[1];
rz(-1.4945806) q[1];
sx q[1];
rz(0.59554312) q[1];
rz(-pi) q[2];
rz(1.2792789) q[3];
sx q[3];
rz(-2.1983578) q[3];
sx q[3];
rz(-1.732839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8949184) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(1.0651917) q[2];
rz(-1.4554321) q[3];
sx q[3];
rz(-1.8543517) q[3];
sx q[3];
rz(-0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0378549) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.5392342) q[0];
rz(0.94929758) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(-1.4303713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8447782) q[0];
sx q[0];
rz(-1.0460633) q[0];
sx q[0];
rz(1.0206166) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6550421) q[2];
sx q[2];
rz(-0.73053321) q[2];
sx q[2];
rz(-0.73405311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1518703) q[1];
sx q[1];
rz(-0.29726004) q[1];
sx q[1];
rz(3.0402335) q[1];
rz(-pi) q[2];
rz(-1.999302) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(-2.9203897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-0.12651786) q[2];
rz(0.33454076) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(-2.6244923) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(-0.089769207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9341016) q[0];
sx q[0];
rz(-2.4330407) q[0];
sx q[0];
rz(-2.2183378) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0788306) q[2];
sx q[2];
rz(-2.1307935) q[2];
sx q[2];
rz(-2.7095209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43347142) q[1];
sx q[1];
rz(-2.8420487) q[1];
sx q[1];
rz(2.5764861) q[1];
rz(-pi) q[2];
rz(-1.0221972) q[3];
sx q[3];
rz(-0.41893235) q[3];
sx q[3];
rz(-1.054712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53139293) q[2];
sx q[2];
rz(-2.0072082) q[2];
sx q[2];
rz(-0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(-2.9288647) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1995734) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(-1.749281) q[1];
sx q[1];
rz(-1.6284457) q[1];
sx q[1];
rz(1.5096691) q[1];
rz(-0.86633273) q[2];
sx q[2];
rz(-1.489218) q[2];
sx q[2];
rz(1.0897286) q[2];
rz(-2.5531664) q[3];
sx q[3];
rz(-1.0259494) q[3];
sx q[3];
rz(-1.1001724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
