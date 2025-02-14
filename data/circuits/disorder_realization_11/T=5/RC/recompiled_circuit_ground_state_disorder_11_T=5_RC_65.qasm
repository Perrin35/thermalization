OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1240368) q[0];
sx q[0];
rz(-0.15260829) q[0];
sx q[0];
rz(2.9583162) q[0];
rz(0.89219379) q[1];
sx q[1];
rz(5.1434864) q[1];
sx q[1];
rz(8.8465717) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5472422) q[0];
sx q[0];
rz(-0.092012398) q[0];
sx q[0];
rz(-0.1081744) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(0.13730857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.86204927) q[1];
sx q[1];
rz(-2.0506564) q[1];
sx q[1];
rz(2.1016438) q[1];
rz(-pi) q[2];
x q[2];
rz(2.036506) q[3];
sx q[3];
rz(-1.6034725) q[3];
sx q[3];
rz(2.1316043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9649428) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(0.098026015) q[2];
rz(-2.3227504) q[3];
sx q[3];
rz(-0.10418532) q[3];
sx q[3];
rz(-0.63481832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1408511) q[0];
sx q[0];
rz(-0.31743693) q[0];
sx q[0];
rz(0.66824085) q[0];
rz(-2.5375598) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(-0.86413962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174422) q[0];
sx q[0];
rz(-1.6675341) q[0];
sx q[0];
rz(0.58944943) q[0];
rz(-pi) q[1];
rz(1.0572724) q[2];
sx q[2];
rz(-0.74716669) q[2];
sx q[2];
rz(-2.7208089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.493638) q[1];
sx q[1];
rz(-1.3820547) q[1];
sx q[1];
rz(-0.24540875) q[1];
rz(-pi) q[2];
rz(-2.2142124) q[3];
sx q[3];
rz(-2.815554) q[3];
sx q[3];
rz(0.044635208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16022564) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(2.3571864) q[2];
rz(-1.9085599) q[3];
sx q[3];
rz(-2.8630856) q[3];
sx q[3];
rz(-1.2109141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311555) q[0];
sx q[0];
rz(-1.5657319) q[0];
sx q[0];
rz(2.1203777) q[0];
rz(1.637623) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(0.57600299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29991462) q[0];
sx q[0];
rz(-1.155071) q[0];
sx q[0];
rz(-1.6982416) q[0];
rz(-pi) q[1];
rz(0.25635135) q[2];
sx q[2];
rz(-1.1927989) q[2];
sx q[2];
rz(-0.5452273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29167029) q[1];
sx q[1];
rz(-1.909316) q[1];
sx q[1];
rz(-0.68847076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27184527) q[3];
sx q[3];
rz(-1.3619641) q[3];
sx q[3];
rz(-2.6110716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7716498) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(0.096262781) q[2];
rz(-2.7104968) q[3];
sx q[3];
rz(-0.96612203) q[3];
sx q[3];
rz(2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-2.3189321) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(2.8462963) q[0];
rz(2.2604306) q[1];
sx q[1];
rz(-0.22717871) q[1];
sx q[1];
rz(0.49240246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8618362) q[0];
sx q[0];
rz(-1.4040213) q[0];
sx q[0];
rz(1.357973) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0855851) q[2];
sx q[2];
rz(-3.036068) q[2];
sx q[2];
rz(-0.72805041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4223692) q[1];
sx q[1];
rz(-1.8998659) q[1];
sx q[1];
rz(-3.0039139) q[1];
rz(-pi) q[2];
rz(1.2232699) q[3];
sx q[3];
rz(-1.1545968) q[3];
sx q[3];
rz(-2.2496441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1262576) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(-0.34273657) q[2];
rz(-0.98894173) q[3];
sx q[3];
rz(-1.1611232) q[3];
sx q[3];
rz(0.69800085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88687503) q[0];
sx q[0];
rz(-2.8509792) q[0];
sx q[0];
rz(1.2683723) q[0];
rz(2.7791924) q[1];
sx q[1];
rz(-0.86530322) q[1];
sx q[1];
rz(-1.5836345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11626205) q[0];
sx q[0];
rz(-1.1209348) q[0];
sx q[0];
rz(2.5309358) q[0];
rz(1.8475632) q[2];
sx q[2];
rz(-1.74969) q[2];
sx q[2];
rz(1.3053633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0456699) q[1];
sx q[1];
rz(-1.5557003) q[1];
sx q[1];
rz(-2.5445055) q[1];
x q[2];
rz(-2.6222428) q[3];
sx q[3];
rz(-1.8867114) q[3];
sx q[3];
rz(0.43866065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4871939) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-0.85713345) q[2];
rz(2.1060139) q[3];
sx q[3];
rz(-0.77719694) q[3];
sx q[3];
rz(-0.63151675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71984464) q[0];
sx q[0];
rz(-2.6578465) q[0];
sx q[0];
rz(2.8068722) q[0];
rz(0.98182976) q[1];
sx q[1];
rz(-1.5900759) q[1];
sx q[1];
rz(0.88012153) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2379358) q[0];
sx q[0];
rz(-1.4159044) q[0];
sx q[0];
rz(1.2747033) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99079972) q[2];
sx q[2];
rz(-0.31337619) q[2];
sx q[2];
rz(2.6277781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8093258) q[1];
sx q[1];
rz(-0.49254815) q[1];
sx q[1];
rz(-1.2713856) q[1];
rz(-pi) q[2];
rz(-2.2273034) q[3];
sx q[3];
rz(-0.56826545) q[3];
sx q[3];
rz(1.9075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52241391) q[2];
sx q[2];
rz(-0.48888561) q[2];
sx q[2];
rz(-1.4948814) q[2];
rz(-1.8374247) q[3];
sx q[3];
rz(-2.2924278) q[3];
sx q[3];
rz(0.18613923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3592767) q[0];
sx q[0];
rz(-0.19141153) q[0];
sx q[0];
rz(3.0707448) q[0];
rz(0.62295667) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(0.43359044) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1529652) q[0];
sx q[0];
rz(-2.835243) q[0];
sx q[0];
rz(-0.10220893) q[0];
rz(-pi) q[1];
rz(2.8500227) q[2];
sx q[2];
rz(-1.5159515) q[2];
sx q[2];
rz(0.26569732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2411464) q[1];
sx q[1];
rz(-1.4747319) q[1];
sx q[1];
rz(2.8166031) q[1];
rz(-pi) q[2];
rz(1.6852156) q[3];
sx q[3];
rz(-1.1297207) q[3];
sx q[3];
rz(2.499388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6314466) q[2];
sx q[2];
rz(-0.53312174) q[2];
sx q[2];
rz(2.5041637) q[2];
rz(2.0006477) q[3];
sx q[3];
rz(-0.44064042) q[3];
sx q[3];
rz(-0.91658896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7618074) q[0];
sx q[0];
rz(-0.26476911) q[0];
sx q[0];
rz(-0.26527143) q[0];
rz(-0.028566407) q[1];
sx q[1];
rz(-2.9539689) q[1];
sx q[1];
rz(2.2144337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68944478) q[0];
sx q[0];
rz(-1.5981705) q[0];
sx q[0];
rz(-0.22449799) q[0];
x q[1];
rz(-0.36998387) q[2];
sx q[2];
rz(-1.4346037) q[2];
sx q[2];
rz(-0.16384928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78955626) q[1];
sx q[1];
rz(-1.6733783) q[1];
sx q[1];
rz(0.7898777) q[1];
rz(0.24683471) q[3];
sx q[3];
rz(-2.7066019) q[3];
sx q[3];
rz(0.37295476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.022543) q[2];
sx q[2];
rz(-2.0930585) q[2];
sx q[2];
rz(1.0724462) q[2];
rz(-2.8065245) q[3];
sx q[3];
rz(-3.0266893) q[3];
sx q[3];
rz(2.9160299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91113126) q[0];
sx q[0];
rz(-1.9879531) q[0];
sx q[0];
rz(-2.352584) q[0];
rz(-2.543653) q[1];
sx q[1];
rz(-1.1006678) q[1];
sx q[1];
rz(2.423563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2418023) q[0];
sx q[0];
rz(-0.16092096) q[0];
sx q[0];
rz(1.3113381) q[0];
rz(1.687019) q[2];
sx q[2];
rz(-1.5801443) q[2];
sx q[2];
rz(-2.7264894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78290547) q[1];
sx q[1];
rz(-0.99596874) q[1];
sx q[1];
rz(-0.7723104) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47649033) q[3];
sx q[3];
rz(-2.9417178) q[3];
sx q[3];
rz(1.0651922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52704063) q[2];
sx q[2];
rz(-2.8895832) q[2];
sx q[2];
rz(-1.4177812) q[2];
rz(-2.0333911) q[3];
sx q[3];
rz(-2.7751444) q[3];
sx q[3];
rz(-0.38780701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.336816) q[0];
sx q[0];
rz(-0.62938654) q[0];
sx q[0];
rz(-2.8861073) q[0];
rz(-0.77087036) q[1];
sx q[1];
rz(-1.5522542) q[1];
sx q[1];
rz(-1.593387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51963888) q[0];
sx q[0];
rz(-2.0072091) q[0];
sx q[0];
rz(3.0163517) q[0];
rz(-pi) q[1];
rz(2.911627) q[2];
sx q[2];
rz(-0.88092025) q[2];
sx q[2];
rz(-1.3892737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0615427) q[1];
sx q[1];
rz(-1.1876186) q[1];
sx q[1];
rz(-1.8661168) q[1];
rz(-pi) q[2];
rz(1.9491626) q[3];
sx q[3];
rz(-2.1244266) q[3];
sx q[3];
rz(-2.7146482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8678681) q[2];
sx q[2];
rz(-2.3994583) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(-2.6015688) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(-2.36256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9191606) q[0];
sx q[0];
rz(-1.5800911) q[0];
sx q[0];
rz(2.0240361) q[0];
rz(-0.22668214) q[1];
sx q[1];
rz(-0.10348987) q[1];
sx q[1];
rz(-1.6685974) q[1];
rz(1.031735) q[2];
sx q[2];
rz(-0.65209263) q[2];
sx q[2];
rz(1.9917464) q[2];
rz(2.4189714) q[3];
sx q[3];
rz(-1.8137877) q[3];
sx q[3];
rz(2.0210101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
