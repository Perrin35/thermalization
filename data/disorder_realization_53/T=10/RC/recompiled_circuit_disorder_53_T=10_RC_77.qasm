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
rz(-0.53108162) q[1];
sx q[1];
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24095806) q[0];
sx q[0];
rz(-1.7358343) q[0];
sx q[0];
rz(1.4532695) q[0];
x q[1];
rz(-0.30795745) q[2];
sx q[2];
rz(-0.27543682) q[2];
sx q[2];
rz(0.76560417) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66346332) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(2.4725735) q[1];
rz(0.19282135) q[3];
sx q[3];
rz(-1.3746327) q[3];
sx q[3];
rz(1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-0.72221243) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(0.41980699) q[0];
rz(-2.7833815) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(2.125724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7656454) q[0];
sx q[0];
rz(-0.65212661) q[0];
sx q[0];
rz(-2.447522) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7989356) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(1.9659496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8034755) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(2.7944399) q[1];
rz(-pi) q[2];
rz(2.0738828) q[3];
sx q[3];
rz(-1.9670608) q[3];
sx q[3];
rz(2.577142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35280716) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719291) q[0];
sx q[0];
rz(-2.0045067) q[0];
sx q[0];
rz(-2.0757872) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7736189) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(1.1676163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5460593) q[1];
sx q[1];
rz(-1.5946348) q[1];
sx q[1];
rz(1.8438575) q[1];
rz(-pi) q[2];
rz(-0.4811901) q[3];
sx q[3];
rz(-1.378873) q[3];
sx q[3];
rz(-2.6709675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.3428358) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(2.0480115) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(0.21903285) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-0.15904388) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.521579) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(0.85711993) q[0];
rz(-1.9424058) q[2];
sx q[2];
rz(-1.2174165) q[2];
sx q[2];
rz(1.6824739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98475458) q[1];
sx q[1];
rz(-1.5459237) q[1];
sx q[1];
rz(1.7825885) q[1];
x q[2];
rz(0.36985107) q[3];
sx q[3];
rz(-2.5858013) q[3];
sx q[3];
rz(-0.37214798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(-2.7159178) q[2];
rz(3.1221636) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5760561) q[0];
sx q[0];
rz(-0.57665529) q[0];
sx q[0];
rz(-1.499349) q[0];
rz(1.3850645) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(0.059543691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1141495) q[1];
sx q[1];
rz(-0.47361923) q[1];
sx q[1];
rz(-1.7042392) q[1];
rz(1.3324276) q[3];
sx q[3];
rz(-0.41394627) q[3];
sx q[3];
rz(-0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(1.0992682) q[0];
rz(-0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(2.1469595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4431865) q[0];
sx q[0];
rz(-0.61265677) q[0];
sx q[0];
rz(-1.959527) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1349796) q[2];
sx q[2];
rz(-2.6380739) q[2];
sx q[2];
rz(2.7343482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1514046) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(0.0068454725) q[1];
rz(2.4410938) q[3];
sx q[3];
rz(-1.2555712) q[3];
sx q[3];
rz(1.8217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(-0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(-2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-2.2454967) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(1.9301374) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.5223741) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30949621) q[1];
sx q[1];
rz(-0.91760228) q[1];
sx q[1];
rz(-3.0597568) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.215407) q[3];
sx q[3];
rz(-1.965596) q[3];
sx q[3];
rz(1.3046164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(1.9141076) q[2];
rz(-3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(-2.7283227) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(2.4024898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0864055) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(-1.5124209) q[0];
x q[1];
rz(1.5831645) q[2];
sx q[2];
rz(-0.45607273) q[2];
sx q[2];
rz(2.9197444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.11074556) q[1];
sx q[1];
rz(-1.7304286) q[1];
sx q[1];
rz(0.84407945) q[1];
rz(1.8352175) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(2.5097178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-1.1018264) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(2.326899) q[3];
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
rz(0.48801625) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-0.5272665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0256955) q[0];
sx q[0];
rz(-0.39931116) q[0];
sx q[0];
rz(-0.77948178) q[0];
rz(0.97986603) q[2];
sx q[2];
rz(-0.96989378) q[2];
sx q[2];
rz(-2.2566233) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3026915) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(-0.80877766) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2654952) q[3];
sx q[3];
rz(-2.4584208) q[3];
sx q[3];
rz(2.9340569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
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
rz(0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-2.9634109) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(-0.88130086) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(2.7351565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5299915) q[1];
sx q[1];
rz(-1.0582663) q[1];
sx q[1];
rz(0.26228735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8799805) q[3];
sx q[3];
rz(-0.44370053) q[3];
sx q[3];
rz(1.651368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085893) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(3.1115816) q[2];
sx q[2];
rz(-1.2118311) q[2];
sx q[2];
rz(-0.98045469) q[2];
rz(2.1929019) q[3];
sx q[3];
rz(-1.8507675) q[3];
sx q[3];
rz(-0.40678195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
