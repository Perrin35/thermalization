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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3492337) q[0];
sx q[0];
rz(-1.6867189) q[0];
sx q[0];
rz(0.16616343) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26308194) q[2];
sx q[2];
rz(-1.4882659) q[2];
sx q[2];
rz(-0.50815998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3218282) q[1];
sx q[1];
rz(-2.1121703) q[1];
sx q[1];
rz(-2.2775047) q[1];
x q[2];
rz(-2.3380321) q[3];
sx q[3];
rz(-0.27419146) q[3];
sx q[3];
rz(0.95925946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3794043) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
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
rz(0.17125601) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(-2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-2.125724) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7091539) q[0];
sx q[0];
rz(-2.0560987) q[0];
sx q[0];
rz(1.1164467) q[0];
rz(-pi) q[1];
rz(1.3426571) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(1.9659496) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8034755) q[1];
sx q[1];
rz(-0.40121597) q[1];
sx q[1];
rz(0.34715279) q[1];
rz(-pi) q[2];
rz(2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-3.0811908) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(0.18149158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6430969) q[0];
sx q[0];
rz(-2.4884014) q[0];
sx q[0];
rz(0.80723395) q[0];
x q[1];
rz(2.3679738) q[2];
sx q[2];
rz(-2.3256989) q[2];
sx q[2];
rz(1.9739763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.060170598) q[1];
sx q[1];
rz(-2.8675188) q[1];
sx q[1];
rz(1.4826135) q[1];
rz(1.7865903) q[3];
sx q[3];
rz(-2.0424235) q[3];
sx q[3];
rz(-2.140688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.024753831) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(-1.8163619) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(-1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670068) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(-2.4530607) q[0];
rz(0.21903285) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(0.15904388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57590398) q[0];
sx q[0];
rz(-2.2825147) q[0];
sx q[0];
rz(-0.84337658) q[0];
x q[1];
rz(0.7775457) q[2];
sx q[2];
rz(-0.50707196) q[2];
sx q[2];
rz(-0.61447243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4403968) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(-1.4529983) q[1];
x q[2];
rz(1.791648) q[3];
sx q[3];
rz(-2.0851118) q[3];
sx q[3];
rz(3.0855892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8461385) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(-3.1221636) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(0.69452906) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(0.12839578) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-0.54840666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4908876) q[0];
sx q[0];
rz(-0.99579949) q[0];
sx q[0];
rz(-3.0951963) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3850645) q[2];
sx q[2];
rz(-2.3302632) q[2];
sx q[2];
rz(-0.059543691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4244528) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(-1.1007924) q[1];
rz(-1.9742825) q[3];
sx q[3];
rz(-1.6659123) q[3];
sx q[3];
rz(2.000589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5925245) q[0];
sx q[0];
rz(-1.3510834) q[0];
sx q[0];
rz(-0.99411221) q[0];
rz(2.6380831) q[2];
sx q[2];
rz(-1.5676055) q[2];
sx q[2];
rz(1.9838331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57833886) q[1];
sx q[1];
rz(-1.5772547) q[1];
sx q[1];
rz(1.9087113) q[1];
rz(-pi) q[2];
rz(2.6732413) q[3];
sx q[3];
rz(-0.75707179) q[3];
sx q[3];
rz(3.0401331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(-0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826236) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(-3.0665841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98260005) q[0];
sx q[0];
rz(-1.2121965) q[0];
sx q[0];
rz(-3.0744809) q[0];
rz(-0.80008614) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.6192186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.311104) q[1];
sx q[1];
rz(-1.5058335) q[1];
sx q[1];
rz(2.225609) q[1];
x q[2];
rz(-0.9644896) q[3];
sx q[3];
rz(-2.4006667) q[3];
sx q[3];
rz(2.9348506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16671495) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1535004) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(2.7283227) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(-0.73910284) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9773213) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(-3.0493899) q[0];
x q[1];
rz(1.1147538) q[2];
sx q[2];
rz(-1.5762435) q[2];
sx q[2];
rz(1.3378439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6369789) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(1.8085338) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3063752) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-1.1018264) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(2.5226412) q[0];
rz(-1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-2.6143262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1934803) q[0];
sx q[0];
rz(-1.2939948) q[0];
sx q[0];
rz(2.8500018) q[0];
rz(-pi) q[1];
rz(0.68248827) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(1.7555106) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83890115) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(-2.332815) q[1];
rz(-pi) q[2];
rz(-2.9016568) q[3];
sx q[3];
rz(-2.2168808) q[3];
sx q[3];
rz(-0.59350384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5321524) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-2.9634109) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4327334) q[0];
sx q[0];
rz(-2.0769925) q[0];
sx q[0];
rz(-0.83336713) q[0];
rz(-0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(2.7351565) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0308932) q[1];
sx q[1];
rz(-2.5712214) q[1];
sx q[1];
rz(1.13899) q[1];
x q[2];
rz(0.26161216) q[3];
sx q[3];
rz(-2.6978921) q[3];
sx q[3];
rz(1.4902247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(-2.318312) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(-1.9227149) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-3.1115816) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
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