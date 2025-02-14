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
rz(-0.99875206) q[0];
sx q[0];
rz(-0.771703) q[0];
sx q[0];
rz(-1.2326711) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(2.2751459) q[1];
sx q[1];
rz(12.187764) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9199149) q[0];
sx q[0];
rz(-0.23338813) q[0];
sx q[0];
rz(1.4063157) q[0];
rz(2.1577055) q[2];
sx q[2];
rz(-2.4574033) q[2];
sx q[2];
rz(-0.049287576) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12593658) q[1];
sx q[1];
rz(-1.7864461) q[1];
sx q[1];
rz(0.68839781) q[1];
x q[2];
rz(-0.572851) q[3];
sx q[3];
rz(-1.8924973) q[3];
sx q[3];
rz(-2.6953146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1401225) q[2];
sx q[2];
rz(-1.9688316) q[2];
sx q[2];
rz(2.6918461) q[2];
rz(-0.55284119) q[3];
sx q[3];
rz(-2.5832085) q[3];
sx q[3];
rz(-2.9302178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387892) q[0];
sx q[0];
rz(-1.9961822) q[0];
sx q[0];
rz(0.75622028) q[0];
rz(0.65762949) q[1];
sx q[1];
rz(-0.88784528) q[1];
sx q[1];
rz(-0.58849803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.156608) q[0];
sx q[0];
rz(-0.78240186) q[0];
sx q[0];
rz(2.6020312) q[0];
x q[1];
rz(2.9220044) q[2];
sx q[2];
rz(-2.3052518) q[2];
sx q[2];
rz(-1.7569923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1900764) q[1];
sx q[1];
rz(-1.5786247) q[1];
sx q[1];
rz(-1.2782793) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28266763) q[3];
sx q[3];
rz(-2.1506607) q[3];
sx q[3];
rz(-1.5527703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53380352) q[2];
sx q[2];
rz(-1.4709996) q[2];
sx q[2];
rz(3.0652453) q[2];
rz(2.7401183) q[3];
sx q[3];
rz(-0.52291003) q[3];
sx q[3];
rz(1.128986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56979316) q[0];
sx q[0];
rz(-0.58572584) q[0];
sx q[0];
rz(-2.8049923) q[0];
rz(-1.0353237) q[1];
sx q[1];
rz(-0.88408771) q[1];
sx q[1];
rz(-0.050315637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49281989) q[0];
sx q[0];
rz(-1.1002212) q[0];
sx q[0];
rz(2.5008766) q[0];
rz(-pi) q[1];
rz(0.59479721) q[2];
sx q[2];
rz(-1.7414879) q[2];
sx q[2];
rz(1.6015944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5118647) q[1];
sx q[1];
rz(-1.9105163) q[1];
sx q[1];
rz(1.7930255) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99934484) q[3];
sx q[3];
rz(-1.3908104) q[3];
sx q[3];
rz(-2.6213329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3745554) q[2];
sx q[2];
rz(-0.87696806) q[2];
sx q[2];
rz(-0.98131895) q[2];
rz(2.0902925) q[3];
sx q[3];
rz(-1.6853354) q[3];
sx q[3];
rz(3.1276957) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0392847) q[0];
sx q[0];
rz(-1.4816875) q[0];
sx q[0];
rz(-0.977595) q[0];
rz(-2.5915937) q[1];
sx q[1];
rz(-0.98098749) q[1];
sx q[1];
rz(2.1994793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4956168) q[0];
sx q[0];
rz(-1.8919404) q[0];
sx q[0];
rz(0.34713388) q[0];
rz(-pi) q[1];
rz(0.35772985) q[2];
sx q[2];
rz(-1.6211281) q[2];
sx q[2];
rz(-2.3711287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4244183) q[1];
sx q[1];
rz(-0.62655973) q[1];
sx q[1];
rz(1.8612629) q[1];
rz(-pi) q[2];
rz(2.5994456) q[3];
sx q[3];
rz(-2.039969) q[3];
sx q[3];
rz(2.5960603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.00086870988) q[2];
sx q[2];
rz(-0.372118) q[2];
sx q[2];
rz(1.5127399) q[2];
rz(-2.3673529) q[3];
sx q[3];
rz(-0.95572487) q[3];
sx q[3];
rz(2.9775508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0999488) q[0];
sx q[0];
rz(-2.994717) q[0];
sx q[0];
rz(0.80648333) q[0];
rz(2.3355314) q[1];
sx q[1];
rz(-1.9652003) q[1];
sx q[1];
rz(0.75256601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6735288) q[0];
sx q[0];
rz(-1.6861334) q[0];
sx q[0];
rz(2.2542984) q[0];
rz(-pi) q[1];
rz(0.089293496) q[2];
sx q[2];
rz(-2.7460423) q[2];
sx q[2];
rz(1.4380921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6347152) q[1];
sx q[1];
rz(-1.992883) q[1];
sx q[1];
rz(-1.4603126) q[1];
rz(-pi) q[2];
rz(2.7769776) q[3];
sx q[3];
rz(-0.43496736) q[3];
sx q[3];
rz(1.2193958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8807184) q[2];
sx q[2];
rz(-1.4148834) q[2];
sx q[2];
rz(-0.073337642) q[2];
rz(-0.65555769) q[3];
sx q[3];
rz(-2.1857502) q[3];
sx q[3];
rz(2.5478794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23621479) q[0];
sx q[0];
rz(-2.0423934) q[0];
sx q[0];
rz(-0.69860506) q[0];
rz(-1.3911432) q[1];
sx q[1];
rz(-0.51934424) q[1];
sx q[1];
rz(-3.0379675) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9148283) q[0];
sx q[0];
rz(-1.5323167) q[0];
sx q[0];
rz(-1.5917042) q[0];
rz(1.6334055) q[2];
sx q[2];
rz(-2.6117628) q[2];
sx q[2];
rz(-1.7018715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0857398) q[1];
sx q[1];
rz(-0.21634783) q[1];
sx q[1];
rz(-2.0592709) q[1];
rz(-pi) q[2];
rz(-1.7781939) q[3];
sx q[3];
rz(-1.816664) q[3];
sx q[3];
rz(1.2896001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.708272) q[2];
sx q[2];
rz(-2.2430113) q[2];
sx q[2];
rz(1.1517322) q[2];
rz(-3.057462) q[3];
sx q[3];
rz(-2.3470272) q[3];
sx q[3];
rz(-0.35552037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56650913) q[0];
sx q[0];
rz(-0.36324781) q[0];
sx q[0];
rz(-1.890924) q[0];
rz(1.4051215) q[1];
sx q[1];
rz(-2.5650918) q[1];
sx q[1];
rz(1.0692495) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0435283) q[0];
sx q[0];
rz(-2.4777924) q[0];
sx q[0];
rz(-2.7190156) q[0];
rz(-pi) q[1];
rz(-1.5965514) q[2];
sx q[2];
rz(-2.6633778) q[2];
sx q[2];
rz(1.1104294) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1508559) q[1];
sx q[1];
rz(-1.6526319) q[1];
sx q[1];
rz(-2.2800755) q[1];
rz(-pi) q[2];
rz(1.8054784) q[3];
sx q[3];
rz(-1.6460287) q[3];
sx q[3];
rz(-0.13401517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8413267) q[2];
sx q[2];
rz(-2.1372676) q[2];
sx q[2];
rz(2.0036073) q[2];
rz(-0.56002069) q[3];
sx q[3];
rz(-1.0617826) q[3];
sx q[3];
rz(1.2573857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1244125) q[0];
sx q[0];
rz(-2.7668598) q[0];
sx q[0];
rz(2.4770233) q[0];
rz(-0.38732227) q[1];
sx q[1];
rz(-0.97159425) q[1];
sx q[1];
rz(-1.9688781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5779424) q[0];
sx q[0];
rz(-1.1563753) q[0];
sx q[0];
rz(2.2060288) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46229273) q[2];
sx q[2];
rz(-1.7708808) q[2];
sx q[2];
rz(0.26319661) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93843469) q[1];
sx q[1];
rz(-0.76701346) q[1];
sx q[1];
rz(-2.2545283) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4807289) q[3];
sx q[3];
rz(-0.58739118) q[3];
sx q[3];
rz(1.7174073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6465801) q[2];
sx q[2];
rz(-2.6436372) q[2];
sx q[2];
rz(0.23769561) q[2];
rz(-2.2869535) q[3];
sx q[3];
rz(-1.5187289) q[3];
sx q[3];
rz(2.2304227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2482727) q[0];
sx q[0];
rz(-1.3398291) q[0];
sx q[0];
rz(0.43333897) q[0];
rz(-1.3029107) q[1];
sx q[1];
rz(-1.7830667) q[1];
sx q[1];
rz(-3.0824331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24640326) q[0];
sx q[0];
rz(-1.329591) q[0];
sx q[0];
rz(3.1144322) q[0];
rz(-1.3174345) q[2];
sx q[2];
rz(-2.4640016) q[2];
sx q[2];
rz(-0.80582011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98461565) q[1];
sx q[1];
rz(-1.2885151) q[1];
sx q[1];
rz(2.8033189) q[1];
rz(-pi) q[2];
rz(0.27579422) q[3];
sx q[3];
rz(-1.0101476) q[3];
sx q[3];
rz(-1.8704853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2032418) q[2];
sx q[2];
rz(-1.4138736) q[2];
sx q[2];
rz(0.69796491) q[2];
rz(-2.5159154) q[3];
sx q[3];
rz(-0.2514078) q[3];
sx q[3];
rz(1.3231369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980001) q[0];
sx q[0];
rz(-2.5830144) q[0];
sx q[0];
rz(2.887605) q[0];
rz(0.99098539) q[1];
sx q[1];
rz(-1.5373983) q[1];
sx q[1];
rz(-0.00016798642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5701262) q[0];
sx q[0];
rz(-1.1227221) q[0];
sx q[0];
rz(-0.49093935) q[0];
x q[1];
rz(-0.047880574) q[2];
sx q[2];
rz(-1.6872395) q[2];
sx q[2];
rz(-0.5393942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40305576) q[1];
sx q[1];
rz(-1.3118032) q[1];
sx q[1];
rz(1.8243755) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9659986) q[3];
sx q[3];
rz(-0.43365824) q[3];
sx q[3];
rz(-0.31037229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.99074546) q[2];
sx q[2];
rz(-1.3517697) q[2];
sx q[2];
rz(-3.0899437) q[2];
rz(1.026574) q[3];
sx q[3];
rz(-0.54205042) q[3];
sx q[3];
rz(-0.76977229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50292618) q[0];
sx q[0];
rz(-2.2421056) q[0];
sx q[0];
rz(-2.5594287) q[0];
rz(2.7405986) q[1];
sx q[1];
rz(-0.6656701) q[1];
sx q[1];
rz(2.3013339) q[1];
rz(-1.0377961) q[2];
sx q[2];
rz(-1.1134951) q[2];
sx q[2];
rz(-2.122428) q[2];
rz(2.544315) q[3];
sx q[3];
rz(-0.18571449) q[3];
sx q[3];
rz(-1.0635536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
