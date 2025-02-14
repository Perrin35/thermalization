OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(-1.2644816) q[1];
sx q[1];
rz(-2.2883132) q[1];
sx q[1];
rz(2.5929911) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596244) q[0];
sx q[0];
rz(-2.7197356) q[0];
sx q[0];
rz(2.3530988) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.110496) q[2];
sx q[2];
rz(-0.59760909) q[2];
sx q[2];
rz(-2.5903167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59034172) q[1];
sx q[1];
rz(-1.1635374) q[1];
sx q[1];
rz(1.007651) q[1];
x q[2];
rz(1.886142) q[3];
sx q[3];
rz(-1.3419125) q[3];
sx q[3];
rz(1.5034408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41142622) q[2];
sx q[2];
rz(-0.41894087) q[2];
sx q[2];
rz(-1.0116928) q[2];
rz(2.2569979) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(-1.8959034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678134) q[0];
sx q[0];
rz(-2.1429006) q[0];
sx q[0];
rz(-0.78773898) q[0];
rz(2.1630321) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(1.2501134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7594371) q[0];
sx q[0];
rz(-0.50001745) q[0];
sx q[0];
rz(-1.3445205) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8626094) q[2];
sx q[2];
rz(-1.1820542) q[2];
sx q[2];
rz(-1.1065567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3652894) q[1];
sx q[1];
rz(-0.76500101) q[1];
sx q[1];
rz(1.6056152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3504418) q[3];
sx q[3];
rz(-1.7856132) q[3];
sx q[3];
rz(-2.1493916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1193739) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(3.019943) q[2];
rz(0.71074784) q[3];
sx q[3];
rz(-2.9080279) q[3];
sx q[3];
rz(0.60230437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1043333) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(-0.65993586) q[0];
rz(-2.624699) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(-2.0491811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683974) q[0];
sx q[0];
rz(-1.930113) q[0];
sx q[0];
rz(-0.27984377) q[0];
x q[1];
rz(2.8404929) q[2];
sx q[2];
rz(-0.89576713) q[2];
sx q[2];
rz(1.1464034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2017434) q[1];
sx q[1];
rz(-0.99786109) q[1];
sx q[1];
rz(2.4948289) q[1];
rz(-pi) q[2];
rz(0.33340065) q[3];
sx q[3];
rz(-0.67615792) q[3];
sx q[3];
rz(2.859883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-2.7984012) q[2];
sx q[2];
rz(-1.421831) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.657594) q[3];
sx q[3];
rz(1.4181731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5487109) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(-1.8362057) q[0];
rz(3.1343754) q[1];
sx q[1];
rz(-0.22166285) q[1];
sx q[1];
rz(0.73297393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78295499) q[0];
sx q[0];
rz(-1.452139) q[0];
sx q[0];
rz(1.4374742) q[0];
x q[1];
rz(-0.29455955) q[2];
sx q[2];
rz(-1.6798875) q[2];
sx q[2];
rz(-0.2738758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0187015) q[1];
sx q[1];
rz(-2.3614493) q[1];
sx q[1];
rz(2.200526) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4184138) q[3];
sx q[3];
rz(-1.4878325) q[3];
sx q[3];
rz(-1.7570329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8492154) q[2];
sx q[2];
rz(-1.4039618) q[2];
sx q[2];
rz(0.89548573) q[2];
rz(-0.53564566) q[3];
sx q[3];
rz(-1.5786542) q[3];
sx q[3];
rz(3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2484922) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(2.6336811) q[0];
rz(-1.4503362) q[1];
sx q[1];
rz(-1.0117057) q[1];
sx q[1];
rz(1.3571665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0471418) q[0];
sx q[0];
rz(-3.019382) q[0];
sx q[0];
rz(-2.727174) q[0];
x q[1];
rz(-1.1512027) q[2];
sx q[2];
rz(-1.9997678) q[2];
sx q[2];
rz(1.2871413) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30631599) q[1];
sx q[1];
rz(-2.7359254) q[1];
sx q[1];
rz(1.8464645) q[1];
rz(2.1693863) q[3];
sx q[3];
rz(-1.2523041) q[3];
sx q[3];
rz(-1.9690258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7148529) q[2];
sx q[2];
rz(-1.4740976) q[2];
sx q[2];
rz(2.0392141) q[2];
rz(1.1357931) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(2.3701325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8847467) q[0];
sx q[0];
rz(-0.9698292) q[0];
sx q[0];
rz(-2.2127175) q[0];
rz(-2.7595787) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(1.6928203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9858157) q[0];
sx q[0];
rz(-0.57812968) q[0];
sx q[0];
rz(-0.57684071) q[0];
x q[1];
rz(-2.6567078) q[2];
sx q[2];
rz(-1.6248091) q[2];
sx q[2];
rz(1.206813) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6544899) q[1];
sx q[1];
rz(-1.1740604) q[1];
sx q[1];
rz(0.11782067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9881416) q[3];
sx q[3];
rz(-2.7355237) q[3];
sx q[3];
rz(-1.4514597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53943071) q[2];
sx q[2];
rz(-2.7950037) q[2];
sx q[2];
rz(-2.3740785) q[2];
rz(-1.8481988) q[3];
sx q[3];
rz(-2.4496205) q[3];
sx q[3];
rz(-2.9366711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(-0.25099227) q[0];
rz(0.17768606) q[1];
sx q[1];
rz(-1.4439986) q[1];
sx q[1];
rz(2.4846855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4068864) q[0];
sx q[0];
rz(-1.7771134) q[0];
sx q[0];
rz(1.330349) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4127633) q[2];
sx q[2];
rz(-0.92301805) q[2];
sx q[2];
rz(2.4142746) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4106573) q[1];
sx q[1];
rz(-1.1852131) q[1];
sx q[1];
rz(-1.8944086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64539306) q[3];
sx q[3];
rz(-1.925996) q[3];
sx q[3];
rz(-1.3076289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80829197) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(2.2283238) q[2];
rz(2.463786) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(-1.0031797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.102757) q[0];
sx q[0];
rz(-1.9982194) q[0];
sx q[0];
rz(2.5586149) q[0];
rz(1.4684756) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(0.34117064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.720168) q[0];
sx q[0];
rz(-2.9377794) q[0];
sx q[0];
rz(1.782062) q[0];
rz(-pi) q[1];
rz(-2.8994045) q[2];
sx q[2];
rz(-1.2953087) q[2];
sx q[2];
rz(-1.7603859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64168186) q[1];
sx q[1];
rz(-1.9371095) q[1];
sx q[1];
rz(-1.9670427) q[1];
x q[2];
rz(1.0811133) q[3];
sx q[3];
rz(-2.4016909) q[3];
sx q[3];
rz(0.56671732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94656452) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(-0.80671802) q[2];
rz(0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(0.88361067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.44613999) q[0];
sx q[0];
rz(-1.0574295) q[0];
sx q[0];
rz(2.9651508) q[0];
rz(-0.31013075) q[1];
sx q[1];
rz(-1.4896432) q[1];
sx q[1];
rz(-0.93607059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98956489) q[0];
sx q[0];
rz(-1.7195722) q[0];
sx q[0];
rz(0.29971931) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1175198) q[2];
sx q[2];
rz(-1.0220811) q[2];
sx q[2];
rz(-2.0224151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.34039482) q[1];
sx q[1];
rz(-0.38005334) q[1];
sx q[1];
rz(1.317522) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0447254) q[3];
sx q[3];
rz(-1.6781665) q[3];
sx q[3];
rz(-2.3176258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.04756847) q[2];
sx q[2];
rz(-1.8393643) q[2];
sx q[2];
rz(0.0085208323) q[2];
rz(0.73355567) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(0.12695299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200639) q[0];
sx q[0];
rz(-2.7286752) q[0];
sx q[0];
rz(0.76989663) q[0];
rz(-3.0072615) q[1];
sx q[1];
rz(-2.3513992) q[1];
sx q[1];
rz(-1.2199527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111356) q[0];
sx q[0];
rz(-0.38549462) q[0];
sx q[0];
rz(-0.66118447) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5346181) q[2];
sx q[2];
rz(-1.3262981) q[2];
sx q[2];
rz(0.90842694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1618859) q[1];
sx q[1];
rz(-2.5570013) q[1];
sx q[1];
rz(-2.4120055) q[1];
rz(-pi) q[2];
rz(-0.071000428) q[3];
sx q[3];
rz(-1.419908) q[3];
sx q[3];
rz(1.3726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28025815) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(-2.3560143) q[2];
rz(-0.33513364) q[3];
sx q[3];
rz(-2.1527055) q[3];
sx q[3];
rz(0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817374) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(-0.68589504) q[1];
sx q[1];
rz(-2.5527725) q[1];
sx q[1];
rz(1.7935161) q[1];
rz(-1.5790719) q[2];
sx q[2];
rz(-2.5604421) q[2];
sx q[2];
rz(1.5887518) q[2];
rz(0.036964702) q[3];
sx q[3];
rz(-1.2782469) q[3];
sx q[3];
rz(-0.20991355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
