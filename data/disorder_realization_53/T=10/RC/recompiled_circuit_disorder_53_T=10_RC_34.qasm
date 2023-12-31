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
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
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
rz(0.16616343) q[0];
rz(-pi) q[1];
rz(2.8785107) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(-2.6334327) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4781293) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(-0.66901916) q[1];
rz(-pi) q[2];
rz(-0.19282135) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(-1.3960658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(2.7963426) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17125601) q[0];
sx q[0];
rz(-2.4617564) q[0];
sx q[0];
rz(0.41980699) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3623747) q[0];
sx q[0];
rz(-1.9694766) q[0];
sx q[0];
rz(2.6108512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83769669) q[2];
sx q[2];
rz(-0.3028377) q[2];
sx q[2];
rz(-2.8267415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8034755) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(0.34715279) q[1];
rz(-pi) q[2];
rz(1.0677098) q[3];
sx q[3];
rz(-1.1745319) q[3];
sx q[3];
rz(-0.5644507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9839086) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-0.18149158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719291) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(-2.0757872) q[0];
rz(2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1101802) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(0.024755342) q[1];
rz(-0.39748945) q[3];
sx q[3];
rz(-0.51525138) q[3];
sx q[3];
rz(-1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(-2.9825488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200136) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(2.2844727) q[0];
x q[1];
rz(-2.7646388) q[2];
sx q[2];
rz(-1.2231584) q[2];
sx q[2];
rz(-2.8958547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98475458) q[1];
sx q[1];
rz(-1.595669) q[1];
sx q[1];
rz(1.3590042) q[1];
x q[2];
rz(0.52491297) q[3];
sx q[3];
rz(-1.3789) q[3];
sx q[3];
rz(1.6247941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(-0.42567483) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869285) q[0];
sx q[0];
rz(-1.5318649) q[0];
sx q[0];
rz(2.1462847) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3734599) q[2];
sx q[2];
rz(-1.7051201) q[2];
sx q[2];
rz(1.7589993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4244528) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(-2.0408003) q[1];
rz(-pi) q[2];
rz(1.3324276) q[3];
sx q[3];
rz(-0.41394627) q[3];
sx q[3];
rz(2.4928989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(-0.87695688) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0393031) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-1.0992682) q[0];
rz(-0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431865) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(1.1820656) q[0];
rz(-1.5744393) q[2];
sx q[2];
rz(-1.0672896) q[2];
sx q[2];
rz(2.7267981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1514046) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(0.0068454725) q[1];
rz(-pi) q[2];
rz(-1.1676222) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(0.50658222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2991335) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(0.075008579) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61178111) q[0];
sx q[0];
rz(-1.5079594) q[0];
sx q[0];
rz(1.9301374) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54589097) q[2];
sx q[2];
rz(-2.2635169) q[2];
sx q[2];
rz(2.8199414) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9662465) q[1];
sx q[1];
rz(-2.4840377) q[1];
sx q[1];
rz(-1.6772126) q[1];
x q[2];
rz(2.6610664) q[3];
sx q[3];
rz(-2.1587545) q[3];
sx q[3];
rz(0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-2.7283227) q[0];
rz(0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(2.4024898) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055187125) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(1.6291717) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0268388) q[2];
sx q[2];
rz(-1.5762435) q[2];
sx q[2];
rz(-1.3378439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3196657) q[1];
sx q[1];
rz(-0.85532665) q[1];
sx q[1];
rz(0.2121851) q[1];
rz(-3.0855582) q[3];
sx q[3];
rz(-1.3067712) q[3];
sx q[3];
rz(2.2173086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-1.1018264) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-2.326899) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11589719) q[0];
sx q[0];
rz(-2.7422815) q[0];
sx q[0];
rz(2.3621109) q[0];
rz(-0.69006069) q[2];
sx q[2];
rz(-2.0482716) q[2];
sx q[2];
rz(-0.32327393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.077721715) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(2.0627756) q[1];
rz(-pi) q[2];
rz(-1.2654952) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(-0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(-2.2100892) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(-2.2562064) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-2.9634109) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3522211) q[0];
sx q[0];
rz(-2.2749315) q[0];
sx q[0];
rz(2.2602918) q[0];
rz(-2.7257596) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(2.7351565) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0514959) q[1];
sx q[1];
rz(-1.7987383) q[1];
sx q[1];
rz(1.0432613) q[1];
rz(-pi) q[2];
rz(1.6931157) q[3];
sx q[3];
rz(-1.1432075) q[3];
sx q[3];
rz(1.7784255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085893) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(1.2188777) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-1.2116829) q[2];
sx q[2];
rz(-1.5426987) q[2];
sx q[2];
rz(0.60088746) q[2];
rz(-1.1124489) q[3];
sx q[3];
rz(-2.4670798) q[3];
sx q[3];
rz(-2.3453875) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
