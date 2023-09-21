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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792359) q[0];
sx q[0];
rz(-1.4548737) q[0];
sx q[0];
rz(0.16616343) q[0];
rz(1.6562535) q[2];
sx q[2];
rz(-1.308631) q[2];
sx q[2];
rz(1.0848349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(0.82378806) q[1];
rz(-pi) q[2];
rz(0.80356055) q[3];
sx q[3];
rz(-0.27419146) q[3];
sx q[3];
rz(-2.1823332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(-0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-2.125724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37594721) q[0];
sx q[0];
rz(-2.489466) q[0];
sx q[0];
rz(-2.447522) q[0];
x q[1];
rz(0.83769669) q[2];
sx q[2];
rz(-0.3028377) q[2];
sx q[2];
rz(2.8267415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8034755) q[1];
sx q[1];
rz(-0.40121597) q[1];
sx q[1];
rz(2.7944399) q[1];
rz(-pi) q[2];
x q[2];
rz(2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(-1.9259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.9839086) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(0.88009673) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(-0.18149158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719291) q[0];
sx q[0];
rz(-2.0045067) q[0];
sx q[0];
rz(1.0658054) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9319839) q[2];
sx q[2];
rz(-1.0227232) q[2];
sx q[2];
rz(2.1266448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.03141244) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(3.1168373) q[1];
rz(-1.3550024) q[3];
sx q[3];
rz(-2.0424235) q[3];
sx q[3];
rz(-2.140688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(0.15904388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656887) q[0];
sx q[0];
rz(-2.2825147) q[0];
sx q[0];
rz(0.84337658) q[0];
x q[1];
rz(2.364047) q[2];
sx q[2];
rz(-0.50707196) q[2];
sx q[2];
rz(-2.5271202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-0.025440865) q[1];
rz(-0.36985107) q[3];
sx q[3];
rz(-0.55579138) q[3];
sx q[3];
rz(2.7694447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(2.7159178) q[2];
rz(-0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(-2.593186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6507051) q[0];
sx q[0];
rz(-0.99579949) q[0];
sx q[0];
rz(3.0951963) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9494909) q[2];
sx q[2];
rz(-2.3641799) q[2];
sx q[2];
rz(-2.8156413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7171399) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(-1.1007924) q[1];
rz(0.10336419) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(-0.38927024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(-1.0855801) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(1.0992682) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.97904) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(-0.26035505) q[0];
x q[1];
rz(-2.6380831) q[2];
sx q[2];
rz(-1.5676055) q[2];
sx q[2];
rz(1.1577595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1514046) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(0.0068454725) q[1];
rz(-1.1676222) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-0.50658222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2991335) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(-2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(-0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(0.89609599) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61178111) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(1.9301374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5957017) q[2];
sx q[2];
rz(-0.87807579) q[2];
sx q[2];
rz(2.8199414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.311104) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(-0.91598367) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1771031) q[3];
sx q[3];
rz(-2.4006667) q[3];
sx q[3];
rz(-0.20674202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.9141076) q[2];
rz(3.098439) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9880923) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(2.7283227) q[0];
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
rz(0.16427134) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(-3.0493899) q[0];
x q[1];
rz(-0.0060672005) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(2.9059682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3196657) q[1];
sx q[1];
rz(-0.85532665) q[1];
sx q[1];
rz(-0.2121851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3063752) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(-1.0194107) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(2.6143262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94811234) q[0];
sx q[0];
rz(-1.8475979) q[0];
sx q[0];
rz(-0.29159082) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97986603) q[2];
sx q[2];
rz(-2.1716989) q[2];
sx q[2];
rz(0.88496937) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.077721715) q[1];
sx q[1];
rz(-1.1472902) q[1];
sx q[1];
rz(1.0788171) q[1];
rz(2.230907) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(-2.0180574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(-2.2100892) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60944027) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(-2.6623181) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-0.17818174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44706599) q[0];
sx q[0];
rz(-2.1994626) q[0];
sx q[0];
rz(2.4987614) q[0];
rz(0.68797942) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(-2.2945987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0514959) q[1];
sx q[1];
rz(-1.7987383) q[1];
sx q[1];
rz(-1.0432613) q[1];
rz(-pi) q[2];
rz(-0.26161216) q[3];
sx q[3];
rz(-2.6978921) q[3];
sx q[3];
rz(1.651368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(1.490996) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
rz(2.8015295) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
