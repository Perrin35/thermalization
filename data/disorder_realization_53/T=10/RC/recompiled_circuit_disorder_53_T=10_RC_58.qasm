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
rz(0.38219163) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(-2.5281744) q[0];
rz(-pi) q[1];
rz(2.8785107) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(-2.6334327) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66346332) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(2.4725735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9487713) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(1.3960658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(2.7963426) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17125601) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(1.0158687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7656454) q[0];
sx q[0];
rz(-0.65212661) q[0];
sx q[0];
rz(0.69407065) q[0];
rz(1.3426571) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(1.9659496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0888372) q[1];
sx q[1];
rz(-1.4375327) q[1];
sx q[1];
rz(0.37957508) q[1];
rz(-pi) q[2];
rz(0.44554168) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(-1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(2.9601011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42230168) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(-2.0757872) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0814221) q[1];
sx q[1];
rz(-0.27407384) q[1];
sx q[1];
rz(1.6589792) q[1];
x q[2];
rz(1.7865903) q[3];
sx q[3];
rz(-2.0424235) q[3];
sx q[3];
rz(1.0009047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(-2.0480115) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(0.68853199) q[0];
rz(-0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(0.15904388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656887) q[0];
sx q[0];
rz(-0.85907798) q[0];
sx q[0];
rz(0.84337658) q[0];
rz(-pi) q[1];
rz(2.7646388) q[2];
sx q[2];
rz(-1.2231584) q[2];
sx q[2];
rz(2.8958547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4403968) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(-1.4529983) q[1];
x q[2];
rz(-1.3499447) q[3];
sx q[3];
rz(-2.0851118) q[3];
sx q[3];
rz(3.0855892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(2.7159178) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(-2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6762125) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(2.593186) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054664139) q[0];
sx q[0];
rz(-1.5318649) q[0];
sx q[0];
rz(2.1462847) q[0];
rz(-pi) q[1];
rz(1.7565281) q[2];
sx q[2];
rz(-2.3302632) q[2];
sx q[2];
rz(0.059543691) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4244528) q[1];
sx q[1];
rz(-1.5100749) q[1];
sx q[1];
rz(-2.0408003) q[1];
x q[2];
rz(-1.9742825) q[3];
sx q[3];
rz(-1.6659123) q[3];
sx q[3];
rz(-1.1410037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(1.0855801) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69840616) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(-1.1820656) q[0];
x q[1];
rz(1.5744393) q[2];
sx q[2];
rz(-2.074303) q[2];
sx q[2];
rz(-0.41479455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99018807) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(-3.1347472) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46835132) q[3];
sx q[3];
rz(-0.75707179) q[3];
sx q[3];
rz(-0.10145951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2991335) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(2.6966406) q[2];
rz(-0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(-2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79338193) q[0];
sx q[0];
rz(-0.36455867) q[0];
sx q[0];
rz(1.7478463) q[0];
x q[1];
rz(-2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(1.6192186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8320964) q[1];
sx q[1];
rz(-0.91760228) q[1];
sx q[1];
rz(0.081835882) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1771031) q[3];
sx q[3];
rz(-0.74092591) q[3];
sx q[3];
rz(-0.20674202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(2.4024898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0864055) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(1.5124209) q[0];
rz(-pi) q[1];
rz(0.0060672005) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(-2.9059682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5046138) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(1.8085338) q[1];
rz(-3.0855582) q[3];
sx q[3];
rz(-1.8348215) q[3];
sx q[3];
rz(0.92428401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3070613) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-2.0397662) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(0.5272665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0256955) q[0];
sx q[0];
rz(-2.7422815) q[0];
sx q[0];
rz(2.3621109) q[0];
rz(-pi) q[1];
rz(-2.4591044) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.3860821) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0638709) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(-1.0788171) q[1];
x q[2];
rz(-1.8760975) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(-2.9340569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(-0.93150345) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-2.6623181) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(2.9634109) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44706599) q[0];
sx q[0];
rz(-0.94213001) q[0];
sx q[0];
rz(0.64283128) q[0];
rz(-2.4536132) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(-2.2945987) q[2];
sx q[3];
rz(-pi/2) q[3];
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
x q[2];
rz(-0.43042572) q[3];
sx q[3];
rz(-1.4595375) q[3];
sx q[3];
rz(-2.984897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2330033) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(-1.2188777) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(1.490996) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
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
