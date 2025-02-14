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
rz(-0.52738458) q[0];
sx q[0];
rz(-0.57317797) q[0];
sx q[0];
rz(0.61668026) q[0];
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(0.33831212) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341921) q[0];
sx q[0];
rz(-1.3284042) q[0];
sx q[0];
rz(1.6079748) q[0];
rz(-pi) q[1];
rz(-0.42023797) q[2];
sx q[2];
rz(-2.3048688) q[2];
sx q[2];
rz(-2.9014996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4614726) q[1];
sx q[1];
rz(-1.0306038) q[1];
sx q[1];
rz(-1.7278746) q[1];
rz(-pi) q[2];
rz(-0.81387384) q[3];
sx q[3];
rz(-1.4424837) q[3];
sx q[3];
rz(1.3303043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3931291) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(2.4836922) q[2];
rz(-2.6864478) q[3];
sx q[3];
rz(-0.15076605) q[3];
sx q[3];
rz(-1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8667792) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(0.28208062) q[0];
rz(0.2433978) q[1];
sx q[1];
rz(-1.8744105) q[1];
sx q[1];
rz(-0.74877053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80303148) q[0];
sx q[0];
rz(-0.73006064) q[0];
sx q[0];
rz(-2.8499313) q[0];
rz(1.6924627) q[2];
sx q[2];
rz(-1.5632544) q[2];
sx q[2];
rz(0.53754025) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0091331) q[1];
sx q[1];
rz(-2.0722425) q[1];
sx q[1];
rz(3.0363068) q[1];
rz(-pi) q[2];
rz(1.5010819) q[3];
sx q[3];
rz(-1.7697213) q[3];
sx q[3];
rz(0.45123842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8613209) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(-1.3410428) q[2];
rz(1.9848112) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(-0.7635428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(-0.65735835) q[0];
rz(1.9961458) q[1];
sx q[1];
rz(-2.1871388) q[1];
sx q[1];
rz(-0.81833902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9377146) q[0];
sx q[0];
rz(-2.7174207) q[0];
sx q[0];
rz(2.5293777) q[0];
x q[1];
rz(1.4253699) q[2];
sx q[2];
rz(-2.8231695) q[2];
sx q[2];
rz(-2.2216036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2215449) q[1];
sx q[1];
rz(-0.86600477) q[1];
sx q[1];
rz(0.66001604) q[1];
x q[2];
rz(1.42055) q[3];
sx q[3];
rz(-0.45805061) q[3];
sx q[3];
rz(0.86298215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1546617) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(-1.7884802) q[2];
rz(0.96495572) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(-2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.9406142) q[0];
rz(-0.0032084223) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(-0.23689717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1334522) q[0];
sx q[0];
rz(-0.72785512) q[0];
sx q[0];
rz(2.2766964) q[0];
x q[1];
rz(-0.63726823) q[2];
sx q[2];
rz(-1.6760525) q[2];
sx q[2];
rz(-0.95544514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8767406) q[1];
sx q[1];
rz(-0.63794604) q[1];
sx q[1];
rz(-1.9015584) q[1];
x q[2];
rz(0.8221738) q[3];
sx q[3];
rz(-2.195925) q[3];
sx q[3];
rz(1.7879947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0420456) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(0.42638865) q[2];
rz(1.03553) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(2.5199913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.34430382) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(-1.9328851) q[1];
sx q[1];
rz(-1.8588926) q[1];
sx q[1];
rz(3.1399609) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8953319) q[0];
sx q[0];
rz(-0.29731942) q[0];
sx q[0];
rz(-2.5410482) q[0];
rz(0.29745488) q[2];
sx q[2];
rz(-1.5337197) q[2];
sx q[2];
rz(-0.86264474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.092068191) q[1];
sx q[1];
rz(-0.18316575) q[1];
sx q[1];
rz(1.8335672) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91301133) q[3];
sx q[3];
rz(-2.8051441) q[3];
sx q[3];
rz(-2.2483765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95167595) q[2];
sx q[2];
rz(-0.92477551) q[2];
sx q[2];
rz(0.16119371) q[2];
rz(-2.7023884) q[3];
sx q[3];
rz(-0.74857155) q[3];
sx q[3];
rz(-2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8301903) q[0];
sx q[0];
rz(-1.4051733) q[0];
sx q[0];
rz(-0.79469529) q[0];
rz(-1.7757802) q[1];
sx q[1];
rz(-2.3410485) q[1];
sx q[1];
rz(-2.6672003) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8891539) q[0];
sx q[0];
rz(-1.948101) q[0];
sx q[0];
rz(2.120156) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55665908) q[2];
sx q[2];
rz(-1.2534598) q[2];
sx q[2];
rz(2.5450626) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51481828) q[1];
sx q[1];
rz(-1.1650135) q[1];
sx q[1];
rz(-0.6231413) q[1];
x q[2];
rz(0.9871363) q[3];
sx q[3];
rz(-1.2452092) q[3];
sx q[3];
rz(-1.1914355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5279493) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(-0.5274241) q[2];
rz(-2.534965) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(1.0724148) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793176) q[0];
sx q[0];
rz(-2.9476705) q[0];
sx q[0];
rz(-0.79757565) q[0];
rz(-2.2480615) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(2.8614047) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93793106) q[0];
sx q[0];
rz(-1.1755921) q[0];
sx q[0];
rz(-0.18919887) q[0];
x q[1];
rz(0.94455384) q[2];
sx q[2];
rz(-1.0172067) q[2];
sx q[2];
rz(-1.1819201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5671605) q[1];
sx q[1];
rz(-0.48266294) q[1];
sx q[1];
rz(-1.5437276) q[1];
rz(-pi) q[2];
rz(-0.35983054) q[3];
sx q[3];
rz(-2.1323279) q[3];
sx q[3];
rz(0.4175182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7402652) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(1.9291482) q[2];
rz(-0.23981833) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(-0.56813204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78385335) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(1.8200112) q[0];
rz(-0.18121885) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(-3.0785353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1940799) q[0];
sx q[0];
rz(-0.79160684) q[0];
sx q[0];
rz(-0.20157728) q[0];
rz(-pi) q[1];
rz(-0.087788344) q[2];
sx q[2];
rz(-1.286418) q[2];
sx q[2];
rz(-1.2149917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.588187) q[1];
sx q[1];
rz(-1.3253731) q[1];
sx q[1];
rz(2.7977944) q[1];
rz(-pi) q[2];
rz(-3.1386915) q[3];
sx q[3];
rz(-1.408268) q[3];
sx q[3];
rz(-1.3292461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(1.6723527) q[2];
rz(-1.3937048) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(-2.1446153) q[0];
rz(-0.793055) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(1.1096032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5376268) q[0];
sx q[0];
rz(-1.0111799) q[0];
sx q[0];
rz(1.9163314) q[0];
rz(1.1446321) q[2];
sx q[2];
rz(-0.89489102) q[2];
sx q[2];
rz(-3.0578095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21002467) q[1];
sx q[1];
rz(-2.3801675) q[1];
sx q[1];
rz(0.85827338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4023191) q[3];
sx q[3];
rz(-2.2003897) q[3];
sx q[3];
rz(-1.633267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41254607) q[2];
sx q[2];
rz(-1.46526) q[2];
sx q[2];
rz(-1.2726146) q[2];
rz(2.0160969) q[3];
sx q[3];
rz(-0.19386217) q[3];
sx q[3];
rz(0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37128714) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(-3.0433997) q[0];
rz(-1.498361) q[1];
sx q[1];
rz(-1.9941092) q[1];
sx q[1];
rz(-2.3032545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107008) q[0];
sx q[0];
rz(-2.3641788) q[0];
sx q[0];
rz(-2.8882508) q[0];
rz(-0.21282332) q[2];
sx q[2];
rz(-0.25770074) q[2];
sx q[2];
rz(0.759173) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0467207) q[1];
sx q[1];
rz(-1.1390242) q[1];
sx q[1];
rz(2.8304965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7915947) q[3];
sx q[3];
rz(-1.9933369) q[3];
sx q[3];
rz(-2.411946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3457634) q[2];
sx q[2];
rz(-2.2553307) q[2];
sx q[2];
rz(-0.32540992) q[2];
rz(0.77643967) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(2.5344892) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008739) q[0];
sx q[0];
rz(-1.2043395) q[0];
sx q[0];
rz(2.1111063) q[0];
rz(0.28618947) q[1];
sx q[1];
rz(-1.20594) q[1];
sx q[1];
rz(1.1084569) q[1];
rz(2.7306225) q[2];
sx q[2];
rz(-0.80425089) q[2];
sx q[2];
rz(-1.2488386) q[2];
rz(-1.2318883) q[3];
sx q[3];
rz(-2.1327426) q[3];
sx q[3];
rz(-0.54973092) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
