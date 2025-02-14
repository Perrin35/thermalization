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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34615883) q[0];
sx q[0];
rz(-1.8654658) q[0];
sx q[0];
rz(2.8351889) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.031096641) q[2];
sx q[2];
rz(-0.59760909) q[2];
sx q[2];
rz(-0.55127599) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42014112) q[1];
sx q[1];
rz(-0.68183091) q[1];
sx q[1];
rz(-0.89116606) q[1];
rz(0.92655801) q[3];
sx q[3];
rz(-0.38739714) q[3];
sx q[3];
rz(-0.54033632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7301664) q[2];
sx q[2];
rz(-0.41894087) q[2];
sx q[2];
rz(1.0116928) q[2];
rz(0.88459477) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(-1.2456892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678134) q[0];
sx q[0];
rz(-0.99869204) q[0];
sx q[0];
rz(0.78773898) q[0];
rz(-0.9785606) q[1];
sx q[1];
rz(-1.9906882) q[1];
sx q[1];
rz(-1.2501134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.387991) q[0];
sx q[0];
rz(-1.6785673) q[0];
sx q[0];
rz(-2.0600256) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97888005) q[2];
sx q[2];
rz(-2.667281) q[2];
sx q[2];
rz(-1.7537376) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82455222) q[1];
sx q[1];
rz(-2.3352156) q[1];
sx q[1];
rz(-0.033407465) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9216209) q[3];
sx q[3];
rz(-1.7860054) q[3];
sx q[3];
rz(-0.53088354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1043333) q[0];
sx q[0];
rz(-0.72838825) q[0];
sx q[0];
rz(2.4816568) q[0];
rz(-0.51689369) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(-1.0924115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3583864) q[0];
sx q[0];
rz(-0.45166812) q[0];
sx q[0];
rz(2.2048401) q[0];
x q[1];
rz(0.87320341) q[2];
sx q[2];
rz(-1.8044458) q[2];
sx q[2];
rz(-0.2327118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1498949) q[1];
sx q[1];
rz(-0.83577613) q[1];
sx q[1];
rz(0.81945547) q[1];
x q[2];
rz(2.808192) q[3];
sx q[3];
rz(-2.4654347) q[3];
sx q[3];
rz(-0.28170965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92521459) q[2];
sx q[2];
rz(-0.34319147) q[2];
sx q[2];
rz(-1.421831) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.4839987) q[3];
sx q[3];
rz(1.7234195) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5487109) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(-1.8362057) q[0];
rz(-0.0072172324) q[1];
sx q[1];
rz(-0.22166285) q[1];
sx q[1];
rz(-2.4086187) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.510988) q[0];
sx q[0];
rz(-0.17824358) q[0];
sx q[0];
rz(-2.3018738) q[0];
rz(-pi) q[1];
rz(-0.36075212) q[2];
sx q[2];
rz(-0.3135598) q[2];
sx q[2];
rz(1.5001198) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8165255) q[1];
sx q[1];
rz(-0.96615929) q[1];
sx q[1];
rz(-2.6139392) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4603314) q[3];
sx q[3];
rz(-0.85064954) q[3];
sx q[3];
rz(3.0283749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2923773) q[2];
sx q[2];
rz(-1.7376309) q[2];
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
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2484922) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(2.6336811) q[0];
rz(1.6912564) q[1];
sx q[1];
rz(-1.0117057) q[1];
sx q[1];
rz(-1.7844261) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535808) q[0];
sx q[0];
rz(-1.6199027) q[0];
sx q[0];
rz(-3.0296369) q[0];
x q[1];
rz(2.4139348) q[2];
sx q[2];
rz(-2.5508159) q[2];
sx q[2];
rz(-1.0340921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5366282) q[1];
sx q[1];
rz(-1.1812897) q[1];
sx q[1];
rz(-3.0252181) q[1];
x q[2];
rz(2.761854) q[3];
sx q[3];
rz(-1.0061227) q[3];
sx q[3];
rz(2.5329451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42673972) q[2];
sx q[2];
rz(-1.667495) q[2];
sx q[2];
rz(2.0392141) q[2];
rz(-2.0057996) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(2.3701325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8847467) q[0];
sx q[0];
rz(-2.1717635) q[0];
sx q[0];
rz(-0.92887512) q[0];
rz(0.38201395) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(-1.4487723) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083915868) q[0];
sx q[0];
rz(-1.2681715) q[0];
sx q[0];
rz(0.50047366) q[0];
rz(-1.5097627) q[2];
sx q[2];
rz(-1.0866797) q[2];
sx q[2];
rz(2.8060437) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1898489) q[1];
sx q[1];
rz(-2.7286224) q[1];
sx q[1];
rz(-1.8443405) q[1];
rz(-pi) q[2];
rz(1.1962863) q[3];
sx q[3];
rz(-1.7315961) q[3];
sx q[3];
rz(2.874115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6021619) q[2];
sx q[2];
rz(-0.346589) q[2];
sx q[2];
rz(-2.3740785) q[2];
rz(-1.8481988) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(2.9366711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9354189) q[0];
sx q[0];
rz(-0.17429166) q[0];
sx q[0];
rz(-0.25099227) q[0];
rz(-2.9639066) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(-2.4846855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7347062) q[0];
sx q[0];
rz(-1.3644793) q[0];
sx q[0];
rz(-1.8112436) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4877704) q[2];
sx q[2];
rz(-1.6966239) q[2];
sx q[2];
rz(0.74761151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73093534) q[1];
sx q[1];
rz(-1.9563795) q[1];
sx q[1];
rz(1.247184) q[1];
x q[2];
rz(2.0054988) q[3];
sx q[3];
rz(-2.1699749) q[3];
sx q[3];
rz(2.6223132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80829197) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(0.91326886) q[2];
rz(0.67780668) q[3];
sx q[3];
rz(-1.5014239) q[3];
sx q[3];
rz(2.138413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0388357) q[0];
sx q[0];
rz(-1.9982194) q[0];
sx q[0];
rz(2.5586149) q[0];
rz(-1.4684756) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(-0.34117064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35638973) q[0];
sx q[0];
rz(-1.613253) q[0];
sx q[0];
rz(1.3713942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.274373) q[2];
sx q[2];
rz(-2.7768306) q[2];
sx q[2];
rz(-2.4976969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64168186) q[1];
sx q[1];
rz(-1.2044831) q[1];
sx q[1];
rz(1.9670427) q[1];
rz(-0.40557762) q[3];
sx q[3];
rz(-2.2080407) q[3];
sx q[3];
rz(0.058407053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94656452) q[2];
sx q[2];
rz(-2.3174758) q[2];
sx q[2];
rz(-2.3348746) q[2];
rz(-0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(-0.88361067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6954527) q[0];
sx q[0];
rz(-1.0574295) q[0];
sx q[0];
rz(-0.17644185) q[0];
rz(-0.31013075) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(-2.2055221) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061358) q[0];
sx q[0];
rz(-1.8671037) q[0];
sx q[0];
rz(1.7264051) q[0];
rz(-pi) q[1];
rz(1.1175198) q[2];
sx q[2];
rz(-2.1195115) q[2];
sx q[2];
rz(1.1191776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1470799) q[1];
sx q[1];
rz(-1.4777061) q[1];
sx q[1];
rz(-1.9398111) q[1];
rz(-pi) q[2];
rz(2.0447254) q[3];
sx q[3];
rz(-1.6781665) q[3];
sx q[3];
rz(-2.3176258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0940242) q[2];
sx q[2];
rz(-1.8393643) q[2];
sx q[2];
rz(3.1330718) q[2];
rz(0.73355567) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(-3.0146397) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200639) q[0];
sx q[0];
rz(-0.41291741) q[0];
sx q[0];
rz(2.371696) q[0];
rz(-3.0072615) q[1];
sx q[1];
rz(-2.3513992) q[1];
sx q[1];
rz(1.92164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1320187) q[0];
sx q[0];
rz(-1.2694799) q[0];
sx q[0];
rz(-1.8150041) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8656857) q[2];
sx q[2];
rz(-2.1572626) q[2];
sx q[2];
rz(-0.49582729) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7999477) q[1];
sx q[1];
rz(-1.994767) q[1];
sx q[1];
rz(1.1553702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.071000428) q[3];
sx q[3];
rz(-1.419908) q[3];
sx q[3];
rz(-1.7689266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28025815) q[2];
sx q[2];
rz(-2.4836149) q[2];
sx q[2];
rz(0.78557837) q[2];
rz(-0.33513364) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(-0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32421865) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(-2.4556976) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(3.1361572) q[2];
sx q[2];
rz(-0.98966826) q[2];
sx q[2];
rz(1.5788509) q[2];
rz(-1.2780581) q[3];
sx q[3];
rz(-1.6061898) q[3];
sx q[3];
rz(1.3502179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
