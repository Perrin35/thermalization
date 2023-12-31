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
rz(3.014325) q[0];
sx q[0];
rz(11.57796) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(1.3487863) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.759401) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(2.5281744) q[0];
x q[1];
rz(-1.6562535) q[2];
sx q[2];
rz(-1.8329617) q[2];
sx q[2];
rz(1.0848349) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(-2.3178046) q[1];
x q[2];
rz(-1.3710255) q[3];
sx q[3];
rz(-1.3817182) q[3];
sx q[3];
rz(-0.13669361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(-2.7833815) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(2.125724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37594721) q[0];
sx q[0];
rz(-0.65212661) q[0];
sx q[0];
rz(2.447522) q[0];
x q[1];
rz(2.303896) q[2];
sx q[2];
rz(-0.3028377) q[2];
sx q[2];
rz(0.31485117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3381172) q[1];
sx q[1];
rz(-0.40121597) q[1];
sx q[1];
rz(2.7944399) q[1];
x q[2];
rz(-0.85605551) q[3];
sx q[3];
rz(-2.5119009) q[3];
sx q[3];
rz(-2.7468908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7887855) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(1.1091728) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(0.060401827) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719291) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(1.0658054) q[0];
rz(-2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(-1.1676163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1101802) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(0.024755342) q[1];
x q[2];
rz(-0.39748945) q[3];
sx q[3];
rz(-0.51525138) q[3];
sx q[3];
rz(-1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.3428358) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(-2.0480115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(-0.15904388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656887) q[0];
sx q[0];
rz(-0.85907798) q[0];
sx q[0];
rz(-2.2982161) q[0];
rz(-0.7775457) q[2];
sx q[2];
rz(-2.6345207) q[2];
sx q[2];
rz(2.5271202) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1568381) q[1];
sx q[1];
rz(-1.595669) q[1];
sx q[1];
rz(-1.7825885) q[1];
rz(-pi) q[2];
rz(2.7717416) q[3];
sx q[3];
rz(-0.55579138) q[3];
sx q[3];
rz(-0.37214798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(-2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.46538019) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(-0.54840666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054664139) q[0];
sx q[0];
rz(-1.6097277) q[0];
sx q[0];
rz(-2.1462847) q[0];
x q[1];
rz(-1.3850645) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(3.082049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4244528) q[1];
sx q[1];
rz(-1.5100749) q[1];
sx q[1];
rz(-2.0408003) q[1];
rz(-1.9742825) q[3];
sx q[3];
rz(-1.6659123) q[3];
sx q[3];
rz(-1.1410037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(1.4667286) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(2.2646358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1625527) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(0.26035505) q[0];
rz(-pi) q[1];
rz(1.5744393) q[2];
sx q[2];
rz(-1.0672896) q[2];
sx q[2];
rz(0.41479455) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0108346) q[1];
sx q[1];
rz(-0.33797435) q[1];
sx q[1];
rz(-1.5902751) q[1];
x q[2];
rz(-2.4410938) q[3];
sx q[3];
rz(-1.2555712) q[3];
sx q[3];
rz(1.3198927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(0.80727243) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(-2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(0.075008579) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.5079594) q[0];
sx q[0];
rz(1.2114552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0117815) q[2];
sx q[2];
rz(-2.2885239) q[2];
sx q[2];
rz(2.7028524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8304886) q[1];
sx q[1];
rz(-1.5058335) q[1];
sx q[1];
rz(0.91598367) q[1];
rz(-2.6610664) q[3];
sx q[3];
rz(-2.1587545) q[3];
sx q[3];
rz(-0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16671495) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(-3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-2.7283227) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(0.73910284) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6571592) q[0];
sx q[0];
rz(-1.6201577) q[0];
sx q[0];
rz(2.5779448) q[0];
x q[1];
rz(1.1147538) q[2];
sx q[2];
rz(-1.5653492) q[2];
sx q[2];
rz(-1.3378439) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5046138) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(1.8085338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8352175) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535764) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(2.1221819) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(0.5272665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4370678) q[0];
sx q[0];
rz(-1.8509812) q[0];
sx q[0];
rz(1.8591451) q[0];
rz(-pi) q[1];
x q[1];
rz(2.451532) q[2];
sx q[2];
rz(-2.0482716) q[2];
sx q[2];
rz(-0.32327393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.077721715) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(1.0788171) q[1];
rz(-pi) q[2];
rz(1.8760975) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(-0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(-0.93150345) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-2.6623181) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-0.17818174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(-2.2602918) q[0];
rz(1.222625) q[2];
sx q[2];
rz(-0.91297075) q[2];
sx q[2];
rz(2.1994176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0514959) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(1.0432613) q[1];
rz(-pi) q[2];
rz(-0.26161216) q[3];
sx q[3];
rz(-0.44370053) q[3];
sx q[3];
rz(1.4902247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(2.318312) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(1.2188777) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-3.1115816) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
rz(0.34006313) q[3];
sx q[3];
rz(-0.9763413) q[3];
sx q[3];
rz(-1.7819596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
