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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(-1.3728859) q[0];
rz(-1.1276487) q[1];
sx q[1];
rz(-1.375066) q[1];
sx q[1];
rz(-2.3553203) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776913) q[0];
sx q[0];
rz(-2.7161975) q[0];
sx q[0];
rz(-1.0546272) q[0];
x q[1];
rz(2.588618) q[2];
sx q[2];
rz(-2.1510501) q[2];
sx q[2];
rz(-1.0667104) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9278501) q[1];
sx q[1];
rz(-1.8600271) q[1];
sx q[1];
rz(-3.1104282) q[1];
rz(-0.97162515) q[3];
sx q[3];
rz(-1.9135981) q[3];
sx q[3];
rz(2.9620217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11975153) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(1.8915668) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64330548) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.9885063) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(0.52871314) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77069717) q[0];
sx q[0];
rz(-0.57790725) q[0];
sx q[0];
rz(0.96630567) q[0];
rz(-pi) q[1];
rz(-1.3350851) q[2];
sx q[2];
rz(-1.3870174) q[2];
sx q[2];
rz(-1.2313796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0549677) q[1];
sx q[1];
rz(-1.0180915) q[1];
sx q[1];
rz(2.5235559) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1558481) q[3];
sx q[3];
rz(-1.364721) q[3];
sx q[3];
rz(-1.2690488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(3.1129692) q[2];
rz(0.35401595) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(0.89135528) q[0];
rz(0.50093961) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0100493) q[0];
sx q[0];
rz(-1.6185404) q[0];
sx q[0];
rz(-3.0555807) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0436613) q[2];
sx q[2];
rz(-1.8190776) q[2];
sx q[2];
rz(2.6206526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0173024) q[1];
sx q[1];
rz(-0.92355403) q[1];
sx q[1];
rz(-2.1976297) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15283382) q[3];
sx q[3];
rz(-0.18251576) q[3];
sx q[3];
rz(1.0618126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(-2.9050262) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(-1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22223602) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(0.87523571) q[0];
rz(-1.596176) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-0.045305591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6117258) q[0];
sx q[0];
rz(-1.5054724) q[0];
sx q[0];
rz(-1.6732814) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90317921) q[2];
sx q[2];
rz(-1.1588237) q[2];
sx q[2];
rz(0.93288405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44573024) q[1];
sx q[1];
rz(-2.4668184) q[1];
sx q[1];
rz(-1.4314907) q[1];
rz(-pi) q[2];
rz(-1.8604467) q[3];
sx q[3];
rz(-0.60263205) q[3];
sx q[3];
rz(-1.8528191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28216013) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(1.947594) q[2];
rz(-2.2512839) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10876656) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(-1.3312181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4109681) q[0];
sx q[0];
rz(-1.5544116) q[0];
sx q[0];
rz(-2.2975305) q[0];
rz(-pi) q[1];
rz(0.26010988) q[2];
sx q[2];
rz(-1.1887728) q[2];
sx q[2];
rz(2.0224689) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9115548) q[1];
sx q[1];
rz(-1.4304407) q[1];
sx q[1];
rz(0.36169238) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2636637) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(3.0832269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0266626) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(2.5866348) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(-0.70372787) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.5237944) q[1];
sx q[1];
rz(-2.9827548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1554398) q[0];
sx q[0];
rz(-1.6912529) q[0];
sx q[0];
rz(-1.7140165) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67262291) q[2];
sx q[2];
rz(-0.19572283) q[2];
sx q[2];
rz(-1.4162404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.006268) q[1];
sx q[1];
rz(-2.6416831) q[1];
sx q[1];
rz(2.8630524) q[1];
x q[2];
rz(-1.7991285) q[3];
sx q[3];
rz(-2.1960746) q[3];
sx q[3];
rz(-2.5061553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77330294) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(2.2840195) q[2];
rz(-1.1693906) q[3];
sx q[3];
rz(-1.2811067) q[3];
sx q[3];
rz(-0.20839553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9661949) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(2.4229557) q[0];
rz(0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-2.4729572) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0575585) q[0];
sx q[0];
rz(-0.72958045) q[0];
sx q[0];
rz(-1.0490824) q[0];
x q[1];
rz(-2.1988772) q[2];
sx q[2];
rz(-1.8657902) q[2];
sx q[2];
rz(1.6571972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1437415) q[1];
sx q[1];
rz(-0.5023379) q[1];
sx q[1];
rz(2.3010002) q[1];
rz(-pi) q[2];
rz(-0.63124048) q[3];
sx q[3];
rz(-1.5641912) q[3];
sx q[3];
rz(2.1932878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4947027) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(2.4155937) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(1.9451709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57121074) q[0];
sx q[0];
rz(-1.397126) q[0];
sx q[0];
rz(0.028868361) q[0];
x q[1];
rz(3.093946) q[2];
sx q[2];
rz(-1.2638014) q[2];
sx q[2];
rz(2.7757211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4774306) q[1];
sx q[1];
rz(-0.71738418) q[1];
sx q[1];
rz(0.75023164) q[1];
x q[2];
rz(2.7128588) q[3];
sx q[3];
rz(-0.7168684) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2339345) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(2.2080803) q[2];
rz(-0.61698169) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(3.0883375) q[0];
rz(0.28757295) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(-0.25921777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0236146) q[0];
sx q[0];
rz(-3.0833457) q[0];
sx q[0];
rz(-1.9969166) q[0];
rz(-pi) q[1];
rz(-1.338663) q[2];
sx q[2];
rz(-1.7001503) q[2];
sx q[2];
rz(-0.30226135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51356835) q[1];
sx q[1];
rz(-1.6322989) q[1];
sx q[1];
rz(2.3223367) q[1];
x q[2];
rz(1.7227371) q[3];
sx q[3];
rz(-2.3309618) q[3];
sx q[3];
rz(-1.6068271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58664924) q[2];
sx q[2];
rz(-1.8879075) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(-0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(0.980353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(2.7594866) q[0];
rz(0.29640472) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.2688676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60703245) q[0];
sx q[0];
rz(-1.1447313) q[0];
sx q[0];
rz(1.0491648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13258719) q[2];
sx q[2];
rz(-2.2397723) q[2];
sx q[2];
rz(2.7821409) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.813432) q[1];
sx q[1];
rz(-2.8012271) q[1];
sx q[1];
rz(1.5068547) q[1];
x q[2];
rz(1.6927229) q[3];
sx q[3];
rz(-2.0215748) q[3];
sx q[3];
rz(-1.6850231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2962013) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(2.7517547) q[2];
rz(1.8440394) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525477) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(2.2139868) q[2];
sx q[2];
rz(-1.19258) q[2];
sx q[2];
rz(1.6473991) q[2];
rz(1.9080347) q[3];
sx q[3];
rz(-0.81546406) q[3];
sx q[3];
rz(-2.245695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
