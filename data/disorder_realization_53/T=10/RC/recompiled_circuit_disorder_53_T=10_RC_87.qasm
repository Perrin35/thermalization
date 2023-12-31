OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9525113) q[0];
sx q[0];
rz(-3.014325) q[0];
sx q[0];
rz(-0.98841086) q[0];
rz(-0.53108162) q[1];
sx q[1];
rz(-2.792882) q[1];
sx q[1];
rz(1.7928064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3492337) q[0];
sx q[0];
rz(-1.4548737) q[0];
sx q[0];
rz(-2.9754292) q[0];
rz(-pi) q[1];
rz(2.8785107) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(0.50815998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(-0.82378806) q[1];
x q[2];
rz(0.80356055) q[3];
sx q[3];
rz(-0.27419146) q[3];
sx q[3];
rz(-2.1823332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3794043) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(2.7963426) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(0.96864831) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17125601) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-0.41980699) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(2.125724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43243877) q[0];
sx q[0];
rz(-1.0854939) q[0];
sx q[0];
rz(2.0251459) q[0];
rz(0.20611368) q[2];
sx q[2];
rz(-1.7942675) q[2];
sx q[2];
rz(-0.4414562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3381172) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(2.7944399) q[1];
x q[2];
rz(2.2855371) q[3];
sx q[3];
rz(-0.62969172) q[3];
sx q[3];
rz(-0.39470181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(-3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9839086) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(1.8799211) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(0.18149158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6430969) q[0];
sx q[0];
rz(-2.4884014) q[0];
sx q[0];
rz(-2.3343587) q[0];
x q[1];
rz(2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0814221) q[1];
sx q[1];
rz(-0.27407384) q[1];
sx q[1];
rz(1.4826135) q[1];
rz(2.7441032) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(-1.8163619) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(2.0480115) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(0.68853199) q[0];
rz(0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-0.15904388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7796411) q[0];
sx q[0];
rz(-0.96966302) q[0];
sx q[0];
rz(2.4848293) q[0];
rz(-2.7646388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(-0.24573791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58069431) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-0.025440865) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52491297) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(-1.5167985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(-0.68583268) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-2.593186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5655366) q[0];
sx q[0];
rz(-2.5649374) q[0];
sx q[0];
rz(1.6422436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7565281) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(-3.082049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7171399) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(-1.1007924) q[1];
x q[2];
rz(-1.8091651) q[3];
sx q[3];
rz(-0.41394627) q[3];
sx q[3];
rz(-0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(-0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4431865) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(1.959527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5744393) q[2];
sx q[2];
rz(-1.0672896) q[2];
sx q[2];
rz(-2.7267981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1514046) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(-0.0068454725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1676222) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(-0.50658222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(2.2454967) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98260005) q[0];
sx q[0];
rz(-1.9293961) q[0];
sx q[0];
rz(0.067111777) q[0];
rz(-pi) q[1];
rz(0.54589097) q[2];
sx q[2];
rz(-0.87807579) q[2];
sx q[2];
rz(0.32165124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9662465) q[1];
sx q[1];
rz(-0.65755492) q[1];
sx q[1];
rz(-1.6772126) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6610664) q[3];
sx q[3];
rz(-0.98283813) q[3];
sx q[3];
rz(2.5939536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(-3.098439) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1535004) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(2.7283227) q[0];
rz(0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(2.4024898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844334) q[0];
sx q[0];
rz(-1.5214349) q[0];
sx q[0];
rz(-2.5779448) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1355255) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(-2.9059682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6369789) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(1.8085338) q[1];
rz(-1.8352175) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(1.1018264) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-2.326899) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(0.6189515) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(-2.6143262) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0256955) q[0];
sx q[0];
rz(-0.39931116) q[0];
sx q[0];
rz(0.77948178) q[0];
rz(-0.69006069) q[2];
sx q[2];
rz(-2.0482716) q[2];
sx q[2];
rz(2.8183187) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4317324) q[1];
sx q[1];
rz(-1.1255956) q[1];
sx q[1];
rz(0.47275895) q[1];
x q[2];
rz(-0.91068565) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(1.1235352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(-2.2100892) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(2.9634109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44706599) q[0];
sx q[0];
rz(-2.1994626) q[0];
sx q[0];
rz(2.4987614) q[0];
rz(-pi) q[1];
rz(-1.222625) q[2];
sx q[2];
rz(-2.2286219) q[2];
sx q[2];
rz(-0.94217506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1106995) q[1];
sx q[1];
rz(-0.57037121) q[1];
sx q[1];
rz(-2.0026026) q[1];
rz(-pi) q[2];
rz(1.448477) q[3];
sx q[3];
rz(-1.9983851) q[3];
sx q[3];
rz(-1.3631671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2330033) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(-1.9227149) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(0.030011054) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
rz(2.0291438) q[3];
sx q[3];
rz(-2.4670798) q[3];
sx q[3];
rz(-2.3453875) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
