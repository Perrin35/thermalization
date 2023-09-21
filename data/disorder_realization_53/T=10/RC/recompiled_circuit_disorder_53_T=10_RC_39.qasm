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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38219163) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(2.5281744) q[0];
x q[1];
rz(-1.6562535) q[2];
sx q[2];
rz(-1.308631) q[2];
sx q[2];
rz(2.0567577) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66346332) q[1];
sx q[1];
rz(-2.1609398) q[1];
sx q[1];
rz(-0.66901916) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80356055) q[3];
sx q[3];
rz(-2.8674012) q[3];
sx q[3];
rz(-0.95925946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(-2.7963426) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.63136357) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.779218) q[0];
sx q[0];
rz(-1.172116) q[0];
sx q[0];
rz(-0.53074145) q[0];
x q[1];
rz(0.20611368) q[2];
sx q[2];
rz(-1.3473251) q[2];
sx q[2];
rz(0.4414562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42900436) q[1];
sx q[1];
rz(-1.9468369) q[1];
sx q[1];
rz(1.7141378) q[1];
rz(2.2855371) q[3];
sx q[3];
rz(-2.5119009) q[3];
sx q[3];
rz(0.39470181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(1.8799211) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(0.18149158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719291) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(-1.0658054) q[0];
rz(-pi) q[1];
rz(-0.9319839) q[2];
sx q[2];
rz(-1.0227232) q[2];
sx q[2];
rz(2.1266448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1101802) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(-3.1168373) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6604026) q[3];
sx q[3];
rz(-1.7627197) q[3];
sx q[3];
rz(0.47062518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.3428358) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57590398) q[0];
sx q[0];
rz(-0.85907798) q[0];
sx q[0];
rz(-2.2982161) q[0];
rz(-pi) q[1];
rz(-2.7646388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(2.8958547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.3590707) q[1];
sx q[1];
rz(0.025440865) q[1];
rz(2.6166797) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(-1.5167985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(0.12839578) q[0];
rz(-0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(-0.54840666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5760561) q[0];
sx q[0];
rz(-0.57665529) q[0];
sx q[0];
rz(-1.6422436) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76813278) q[2];
sx q[2];
rz(-1.4364725) q[2];
sx q[2];
rz(-1.7589993) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17715958) q[1];
sx q[1];
rz(-2.0398643) q[1];
sx q[1];
rz(-0.068084929) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0382285) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(0.38927024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(-2.2646358) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10228957) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(1.0992682) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54906819) q[0];
sx q[0];
rz(-1.7905092) q[0];
sx q[0];
rz(2.1474804) q[0];
x q[1];
rz(-1.5744393) q[2];
sx q[2];
rz(-2.074303) q[2];
sx q[2];
rz(-2.7267981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57833886) q[1];
sx q[1];
rz(-1.5772547) q[1];
sx q[1];
rz(-1.9087113) q[1];
x q[2];
rz(2.6732413) q[3];
sx q[3];
rz(-2.3845209) q[3];
sx q[3];
rz(-3.0401331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2991335) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(2.6966406) q[2];
rz(0.80727243) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826236) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(2.2454967) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(-3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3482107) q[0];
sx q[0];
rz(-2.777034) q[0];
sx q[0];
rz(1.7478463) q[0];
rz(-0.54589097) q[2];
sx q[2];
rz(-0.87807579) q[2];
sx q[2];
rz(-0.32165124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8320964) q[1];
sx q[1];
rz(-0.91760228) q[1];
sx q[1];
rz(0.081835882) q[1];
rz(0.92618561) q[3];
sx q[3];
rz(-1.1759967) q[3];
sx q[3];
rz(1.8369762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.2274851) q[2];
rz(-0.04315367) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1535004) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(2.7283227) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055187125) q[0];
sx q[0];
rz(-1.0079181) q[0];
sx q[0];
rz(-1.5124209) q[0];
rz(-pi) q[1];
rz(1.5831645) q[2];
sx q[2];
rz(-0.45607273) q[2];
sx q[2];
rz(2.9197444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11074556) q[1];
sx q[1];
rz(-1.7304286) q[1];
sx q[1];
rz(-2.2975132) q[1];
rz(-pi) q[2];
rz(-1.3063752) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(-0.63187481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(1.1018264) q[2];
rz(-2.7448591) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(2.5226412) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(0.5272665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0256955) q[0];
sx q[0];
rz(-0.39931116) q[0];
sx q[0];
rz(2.3621109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68248827) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.7555106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4317324) q[1];
sx q[1];
rz(-2.0159971) q[1];
sx q[1];
rz(-0.47275895) q[1];
x q[2];
rz(-0.91068565) q[3];
sx q[3];
rz(-1.761697) q[3];
sx q[3];
rz(2.0180574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(2.6623181) q[0];
rz(0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(2.9634109) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6945267) q[0];
sx q[0];
rz(-2.1994626) q[0];
sx q[0];
rz(-2.4987614) q[0];
x q[1];
rz(0.68797942) q[2];
sx q[2];
rz(-1.844172) q[2];
sx q[2];
rz(2.2945987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0900967) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(-2.0983314) q[1];
rz(2.7111669) q[3];
sx q[3];
rz(-1.6820551) q[3];
sx q[3];
rz(2.984897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(2.1774489) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.490996) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
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