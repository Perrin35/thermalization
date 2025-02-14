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
rz(0.78627237) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776913) q[0];
sx q[0];
rz(-2.7161975) q[0];
sx q[0];
rz(2.0869654) q[0];
rz(-pi) q[1];
rz(-0.55297466) q[2];
sx q[2];
rz(-0.99054256) q[2];
sx q[2];
rz(1.0667104) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.21374258) q[1];
sx q[1];
rz(-1.2815655) q[1];
sx q[1];
rz(3.1104282) q[1];
rz(-pi) q[2];
rz(0.97162515) q[3];
sx q[3];
rz(-1.2279945) q[3];
sx q[3];
rz(-0.179571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0218411) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(-1.2500259) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.64330548) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.9885063) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(-0.52871314) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3246177) q[0];
sx q[0];
rz(-1.2551089) q[0];
sx q[0];
rz(-1.0782918) q[0];
rz(-pi) q[1];
rz(0.89854449) q[2];
sx q[2];
rz(-2.8437584) q[2];
sx q[2];
rz(0.98984776) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0866249) q[1];
sx q[1];
rz(-2.1235012) q[1];
sx q[1];
rz(-0.61803671) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0924322) q[3];
sx q[3];
rz(-2.6809664) q[3];
sx q[3];
rz(3.0085029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(-0.028623494) q[2];
rz(0.35401595) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(0.55907512) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(-0.89135528) q[0];
rz(-2.640653) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0552154) q[0];
sx q[0];
rz(-0.098345938) q[0];
sx q[0];
rz(0.50758657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0436613) q[2];
sx q[2];
rz(-1.8190776) q[2];
sx q[2];
rz(-0.52094007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14172983) q[1];
sx q[1];
rz(-1.0838306) q[1];
sx q[1];
rz(0.75090639) q[1];
rz(-pi) q[2];
rz(-1.5988878) q[3];
sx q[3];
rz(-1.7511611) q[3];
sx q[3];
rz(-1.9244058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7848876) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(-2.9050262) q[2];
rz(0.92075721) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9193566) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(0.87523571) q[0];
rz(-1.596176) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(3.0962871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6117258) q[0];
sx q[0];
rz(-1.5054724) q[0];
sx q[0];
rz(-1.6732814) q[0];
rz(-pi) q[1];
rz(0.90317921) q[2];
sx q[2];
rz(-1.1588237) q[2];
sx q[2];
rz(-2.2087086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.873535) q[1];
sx q[1];
rz(-2.2378439) q[1];
sx q[1];
rz(-3.0309543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8604467) q[3];
sx q[3];
rz(-0.60263205) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8594325) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.947594) q[2];
rz(-2.2512839) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(-1.8531063) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85859834) q[0];
sx q[0];
rz(-0.72688519) q[0];
sx q[0];
rz(-1.5461393) q[0];
x q[1];
rz(-1.17679) q[2];
sx q[2];
rz(-1.3298243) q[2];
sx q[2];
rz(-0.35277982) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9115548) q[1];
sx q[1];
rz(-1.4304407) q[1];
sx q[1];
rz(-0.36169238) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87792895) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(-0.058365783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1149301) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(2.4746573) q[2];
rz(-0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-0.2989029) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126467) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(-0.70372787) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(2.9827548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1554398) q[0];
sx q[0];
rz(-1.4503398) q[0];
sx q[0];
rz(1.4275761) q[0];
rz(0.67262291) q[2];
sx q[2];
rz(-2.9458698) q[2];
sx q[2];
rz(-1.7253523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68137121) q[1];
sx q[1];
rz(-1.4386144) q[1];
sx q[1];
rz(2.6580407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63780906) q[3];
sx q[3];
rz(-1.3862228) q[3];
sx q[3];
rz(-1.0705494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3682897) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(-0.85757315) q[2];
rz(-1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(0.20839553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9661949) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(-2.4229557) q[0];
rz(0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-2.4729572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0840341) q[0];
sx q[0];
rz(-2.4120122) q[0];
sx q[0];
rz(2.0925103) q[0];
x q[1];
rz(-2.7823607) q[2];
sx q[2];
rz(-0.97373913) q[2];
sx q[2];
rz(0.12166858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.23827502) q[1];
sx q[1];
rz(-1.897745) q[1];
sx q[1];
rz(-1.1823149) q[1];
rz(1.578978) q[3];
sx q[3];
rz(-2.2020209) q[3];
sx q[3];
rz(-0.62731987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4947027) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(-1.0364214) q[2];
rz(2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(-2.4155937) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(-1.1964218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5703819) q[0];
sx q[0];
rz(-1.7444667) q[0];
sx q[0];
rz(-3.1127243) q[0];
x q[1];
rz(1.8781186) q[2];
sx q[2];
rz(-1.5253789) q[2];
sx q[2];
rz(-1.9222593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4359522) q[1];
sx q[1];
rz(-2.0355823) q[1];
sx q[1];
rz(2.57354) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7128588) q[3];
sx q[3];
rz(-0.7168684) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90765816) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-0.93351239) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(-0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1652949) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(0.053255178) q[0];
rz(0.28757295) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(0.25921777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1179781) q[0];
sx q[0];
rz(-0.058246944) q[0];
sx q[0];
rz(-1.144676) q[0];
rz(1.8029297) q[2];
sx q[2];
rz(-1.7001503) q[2];
sx q[2];
rz(2.8393313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51356835) q[1];
sx q[1];
rz(-1.6322989) q[1];
sx q[1];
rz(0.81925591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3756392) q[3];
sx q[3];
rz(-1.6807091) q[3];
sx q[3];
rz(-3.0725057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58664924) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-0.1813691) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3996537) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(0.3821061) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.872725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5871966) q[0];
sx q[0];
rz(-0.66080873) q[0];
sx q[0];
rz(-2.3095678) q[0];
x q[1];
rz(-0.13258719) q[2];
sx q[2];
rz(-0.90182038) q[2];
sx q[2];
rz(-0.35945177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8386744) q[1];
sx q[1];
rz(-1.5921291) q[1];
sx q[1];
rz(1.9105186) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8954072) q[3];
sx q[3];
rz(-2.6757112) q[3];
sx q[3];
rz(1.7307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84539139) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(0.38983795) q[2];
rz(1.8440394) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98904499) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(-0.46089725) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(2.3572902) q[3];
sx q[3];
rz(-1.3275066) q[3];
sx q[3];
rz(2.2307997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
