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
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(-1.7928064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006346) q[0];
sx q[0];
rz(-1.4057584) q[0];
sx q[0];
rz(-1.4532695) q[0];
x q[1];
rz(1.4853391) q[2];
sx q[2];
rz(-1.8329617) q[2];
sx q[2];
rz(-2.0567577) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66346332) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(0.66901916) q[1];
x q[2];
rz(-0.19282135) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(-1.3960658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-0.34525004) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(2.125724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37594721) q[0];
sx q[0];
rz(-0.65212661) q[0];
sx q[0];
rz(0.69407065) q[0];
rz(0.20611368) q[2];
sx q[2];
rz(-1.3473251) q[2];
sx q[2];
rz(0.4414562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0527555) q[1];
sx q[1];
rz(-1.4375327) q[1];
sx q[1];
rz(-0.37957508) q[1];
x q[2];
rz(-2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(-1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9839086) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-2.9601011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42230168) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(1.0658054) q[0];
rz(2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5955334) q[1];
sx q[1];
rz(-1.5469578) q[1];
sx q[1];
rz(-1.2977352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4811901) q[3];
sx q[3];
rz(-1.7627197) q[3];
sx q[3];
rz(2.6709675) q[3];
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
rz(1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(2.9825488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36195155) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(2.4848293) q[0];
x q[1];
rz(0.7775457) q[2];
sx q[2];
rz(-2.6345207) q[2];
sx q[2];
rz(-2.5271202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4403968) q[1];
sx q[1];
rz(-2.9283668) q[1];
sx q[1];
rz(1.4529983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6166797) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(-1.6247941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29545414) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(-0.019429026) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(-2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
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
rz(-2.3734599) q[2];
sx q[2];
rz(-1.7051201) q[2];
sx q[2];
rz(1.3825934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17715958) q[1];
sx q[1];
rz(-1.1017283) q[1];
sx q[1];
rz(3.0735077) q[1];
x q[2];
rz(1.1673101) q[3];
sx q[3];
rz(-1.4756804) q[3];
sx q[3];
rz(-2.000589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0393031) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-1.0992682) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1625527) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(0.26035505) q[0];
x q[1];
rz(-0.0066130916) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(-2.7343482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99018807) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(3.1347472) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70049882) q[3];
sx q[3];
rz(-1.8860215) q[3];
sx q[3];
rz(-1.3198927) q[3];
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
rz(-0.10891309) q[3];
sx q[3];
rz(0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(-1.2114552) q[0];
x q[1];
rz(2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.6192186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.311104) q[1];
sx q[1];
rz(-1.5058335) q[1];
sx q[1];
rz(-0.91598367) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6610664) q[3];
sx q[3];
rz(-2.1587545) q[3];
sx q[3];
rz(0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(2.7283227) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(-0.73910284) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16427134) q[0];
sx q[0];
rz(-2.5760206) q[0];
sx q[0];
rz(3.0493899) q[0];
rz(-pi) q[1];
rz(2.0268388) q[2];
sx q[2];
rz(-1.5653492) q[2];
sx q[2];
rz(1.3378439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0308471) q[1];
sx q[1];
rz(-1.7304286) q[1];
sx q[1];
rz(-0.84407945) q[1];
rz(1.3063752) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(-2.0397662) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(2.6143262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4370678) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.8591451) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4591044) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.3860821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7098602) q[1];
sx q[1];
rz(-1.1255956) q[1];
sx q[1];
rz(-2.6688337) q[1];
rz(2.230907) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(-2.0180574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(-0.93150345) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(-2.9634109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885926) q[0];
sx q[0];
rz(-1.0646001) q[0];
sx q[0];
rz(-0.83336713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9189677) q[2];
sx q[2];
rz(-2.2286219) q[2];
sx q[2];
rz(-0.94217506) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.0983314) q[1];
rz(-pi) q[2];
rz(-0.43042572) q[3];
sx q[3];
rz(-1.6820551) q[3];
sx q[3];
rz(-0.15669565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(2.318312) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.9085893) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(-1.2188777) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(1.9299097) q[2];
sx q[2];
rz(-1.5426987) q[2];
sx q[2];
rz(0.60088746) q[2];
rz(-2.1929019) q[3];
sx q[3];
rz(-1.2908251) q[3];
sx q[3];
rz(2.7348107) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
