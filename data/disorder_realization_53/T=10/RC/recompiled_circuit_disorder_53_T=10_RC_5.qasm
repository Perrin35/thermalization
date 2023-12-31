OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(0.98841086) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(1.3487863) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006346) q[0];
sx q[0];
rz(-1.7358343) q[0];
sx q[0];
rz(-1.4532695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8336352) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(-0.76560417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29404902) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(-0.82378806) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19282135) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(-1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-0.34525004) q[2];
rz(0.18940645) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-0.96864831) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(0.41980699) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-2.125724) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3623747) q[0];
sx q[0];
rz(-1.9694766) q[0];
sx q[0];
rz(-2.6108512) q[0];
x q[1];
rz(-1.3426571) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(1.1756431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8034755) q[1];
sx q[1];
rz(-0.40121597) q[1];
sx q[1];
rz(-2.7944399) q[1];
x q[2];
rz(0.44554168) q[3];
sx q[3];
rz(-1.1098776) q[3];
sx q[3];
rz(1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(-2.9601011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719291) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(-2.0757872) q[0];
rz(-pi) q[1];
rz(-2.491465) q[2];
sx q[2];
rz(-1.0368772) q[2];
sx q[2];
rz(-0.18661737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.03141244) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(-0.024755342) q[1];
rz(-pi) q[2];
rz(2.7441032) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.3428358) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(-2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-0.15904388) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.521579) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(0.85711993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7775457) q[2];
sx q[2];
rz(-2.6345207) q[2];
sx q[2];
rz(2.5271202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.3590707) q[1];
sx q[1];
rz(-3.1161518) q[1];
x q[2];
rz(-1.791648) q[3];
sx q[3];
rz(-1.0564809) q[3];
sx q[3];
rz(3.0855892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(-2.7159178) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(0.69452906) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(3.0131969) q[0];
rz(0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(0.54840666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054664139) q[0];
sx q[0];
rz(-1.5318649) q[0];
sx q[0];
rz(0.995308) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7565281) q[2];
sx q[2];
rz(-2.3302632) q[2];
sx q[2];
rz(-0.059543691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17715958) q[1];
sx q[1];
rz(-2.0398643) q[1];
sx q[1];
rz(-3.0735077) q[1];
rz(-1.9742825) q[3];
sx q[3];
rz(-1.4756804) q[3];
sx q[3];
rz(-2.000589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8473062) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(2.2646358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(-1.0992682) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5925245) q[0];
sx q[0];
rz(-1.7905092) q[0];
sx q[0];
rz(-2.1474804) q[0];
rz(-1.5671533) q[2];
sx q[2];
rz(-1.0672896) q[2];
sx q[2];
rz(0.41479455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99018807) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(-0.0068454725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9739705) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(0.50658222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2991335) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-2.6966406) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(0.075008579) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61178111) q[0];
sx q[0];
rz(-1.5079594) q[0];
sx q[0];
rz(-1.2114552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54589097) q[2];
sx q[2];
rz(-2.2635169) q[2];
sx q[2];
rz(-2.8199414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.311104) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(2.225609) q[1];
rz(-0.48052629) q[3];
sx q[3];
rz(-2.1587545) q[3];
sx q[3];
rz(0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1535004) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-0.73910284) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9773213) q[0];
sx q[0];
rz(-2.5760206) q[0];
sx q[0];
rz(-0.092202734) q[0];
rz(-pi) q[1];
rz(-2.0268388) q[2];
sx q[2];
rz(-1.5653492) q[2];
sx q[2];
rz(-1.3378439) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0308471) q[1];
sx q[1];
rz(-1.4111641) q[1];
sx q[1];
rz(2.2975132) q[1];
rz(-3.0855582) q[3];
sx q[3];
rz(-1.8348215) q[3];
sx q[3];
rz(0.92428401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(1.1018264) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(-1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(0.5272665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70452481) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.8591451) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69006069) q[2];
sx q[2];
rz(-1.0933211) q[2];
sx q[2];
rz(-2.8183187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83890115) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(-2.332815) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23993581) q[3];
sx q[3];
rz(-2.2168808) q[3];
sx q[3];
rz(-2.5480888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(0.17818174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3522211) q[0];
sx q[0];
rz(-2.2749315) q[0];
sx q[0];
rz(2.2602918) q[0];
rz(0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(-2.7351565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0900967) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(1.0432613) q[1];
rz(-pi) q[2];
rz(0.26161216) q[3];
sx q[3];
rz(-2.6978921) q[3];
sx q[3];
rz(-1.651368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-2.1774489) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(-1.2116829) q[2];
sx q[2];
rz(-1.5426987) q[2];
sx q[2];
rz(0.60088746) q[2];
rz(-2.0291438) q[3];
sx q[3];
rz(-0.67451285) q[3];
sx q[3];
rz(0.79620517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
