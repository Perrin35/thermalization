OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.190783) q[0];
sx q[0];
rz(-1.3124183) q[0];
sx q[0];
rz(-1.4825312) q[0];
rz(-pi) q[1];
rz(-0.30774967) q[2];
sx q[2];
rz(-2.0648742) q[2];
sx q[2];
rz(-2.5623164) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6482918) q[1];
sx q[1];
rz(-2.9420605) q[1];
sx q[1];
rz(-2.3471911) q[1];
rz(-pi) q[2];
rz(2.9631859) q[3];
sx q[3];
rz(-2.3462786) q[3];
sx q[3];
rz(0.91245302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67291659) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(-2.1762302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4345877) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(-2.7804651) q[0];
rz(1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(-0.8180058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084413962) q[0];
sx q[0];
rz(-1.955843) q[0];
sx q[0];
rz(1.2731228) q[0];
x q[1];
rz(0.12273522) q[2];
sx q[2];
rz(-1.9397748) q[2];
sx q[2];
rz(-2.4462647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2172912) q[1];
sx q[1];
rz(-2.0089307) q[1];
sx q[1];
rz(-2.2469254) q[1];
rz(1.7198256) q[3];
sx q[3];
rz(-0.36747284) q[3];
sx q[3];
rz(3.0847103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(-1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4017568) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(-0.31618205) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(-0.2972163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1145828) q[0];
sx q[0];
rz(-1.6342589) q[0];
sx q[0];
rz(-2.7364536) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4386114) q[2];
sx q[2];
rz(-2.2361122) q[2];
sx q[2];
rz(-2.2192628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6857022) q[1];
sx q[1];
rz(-2.9487902) q[1];
sx q[1];
rz(1.3193858) q[1];
rz(-2.475297) q[3];
sx q[3];
rz(-2.1994281) q[3];
sx q[3];
rz(2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(-2.4598222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101689) q[0];
sx q[0];
rz(-1.7792601) q[0];
sx q[0];
rz(1.8375977) q[0];
x q[1];
rz(1.8800456) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(2.4583465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7766429) q[1];
sx q[1];
rz(-2.4185191) q[1];
sx q[1];
rz(0.88100453) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9346312) q[3];
sx q[3];
rz(-1.2483276) q[3];
sx q[3];
rz(-0.22842562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(0.99984461) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(1.6220185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1006267) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(0.53951453) q[0];
x q[1];
rz(-2.4155596) q[2];
sx q[2];
rz(-2.3880929) q[2];
sx q[2];
rz(-1.7475278) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2625761) q[1];
sx q[1];
rz(-0.82619709) q[1];
sx q[1];
rz(1.3000342) q[1];
rz(-pi) q[2];
x q[2];
rz(1.394746) q[3];
sx q[3];
rz(-2.5213443) q[3];
sx q[3];
rz(-0.32905096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(0.2229283) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.3870753) q[3];
sx q[3];
rz(-0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.8575352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4053455) q[0];
sx q[0];
rz(-0.9760455) q[0];
sx q[0];
rz(-0.88881641) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.861704) q[2];
sx q[2];
rz(-1.9473238) q[2];
sx q[2];
rz(-1.0371475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2388873) q[1];
sx q[1];
rz(-1.1088015) q[1];
sx q[1];
rz(1.2483031) q[1];
rz(-2.1070126) q[3];
sx q[3];
rz(-1.8367193) q[3];
sx q[3];
rz(1.769161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(-2.2033851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8021916) q[0];
sx q[0];
rz(-1.8759449) q[0];
sx q[0];
rz(1.9992273) q[0];
rz(-pi) q[1];
rz(-2.5222048) q[2];
sx q[2];
rz(-2.3084547) q[2];
sx q[2];
rz(0.39820652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92749121) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(2.6919041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35258099) q[3];
sx q[3];
rz(-2.4028824) q[3];
sx q[3];
rz(-2.0032361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(-1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(2.3103255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76038218) q[0];
sx q[0];
rz(-1.7179278) q[0];
sx q[0];
rz(2.7358664) q[0];
rz(-1.3966884) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(-2.1070534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7716678) q[1];
sx q[1];
rz(-2.2450788) q[1];
sx q[1];
rz(0.64529459) q[1];
rz(-pi) q[2];
rz(2.3723699) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(0.10570082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.9780805) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(2.1642165) q[0];
rz(-1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(2.0358553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0375263) q[0];
sx q[0];
rz(-2.6797047) q[0];
sx q[0];
rz(-2.877263) q[0];
rz(3.0366304) q[2];
sx q[2];
rz(-2.2347921) q[2];
sx q[2];
rz(0.86519372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2178206) q[1];
sx q[1];
rz(-2.8433228) q[1];
sx q[1];
rz(-0.89426269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7902778) q[3];
sx q[3];
rz(-1.8422541) q[3];
sx q[3];
rz(-0.16004496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5332807) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-0.14979714) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0906592) q[0];
sx q[0];
rz(-1.5516084) q[0];
sx q[0];
rz(-3.1113383) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7255515) q[2];
sx q[2];
rz(-0.96826474) q[2];
sx q[2];
rz(-0.11520152) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64393109) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(0.58845206) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4526669) q[3];
sx q[3];
rz(-1.6932634) q[3];
sx q[3];
rz(2.1446705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(0.76114571) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54820838) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(2.7535915) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(-0.64994561) q[2];
sx q[2];
rz(-2.2192498) q[2];
sx q[2];
rz(0.50512209) q[2];
rz(1.0317867) q[3];
sx q[3];
rz(-1.4236593) q[3];
sx q[3];
rz(0.66766213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
