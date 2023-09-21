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
rz(1.3487863) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38219163) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(-0.6134183) q[0];
rz(-pi) q[1];
rz(-1.4853391) q[2];
sx q[2];
rz(-1.8329617) q[2];
sx q[2];
rz(2.0567577) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3218282) q[1];
sx q[1];
rz(-1.0294224) q[1];
sx q[1];
rz(-2.2775047) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9487713) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(-1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-0.34525004) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-0.96864831) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17125601) q[0];
sx q[0];
rz(-2.4617564) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3623747) q[0];
sx q[0];
rz(-1.172116) q[0];
sx q[0];
rz(-0.53074145) q[0];
rz(-pi) q[1];
rz(2.303896) q[2];
sx q[2];
rz(-2.838755) q[2];
sx q[2];
rz(-0.31485117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0888372) q[1];
sx q[1];
rz(-1.4375327) q[1];
sx q[1];
rz(-0.37957508) q[1];
rz(-pi) q[2];
rz(-2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(1.9259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9839086) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-0.18149158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719291) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(2.0757872) q[0];
rz(2.491465) q[2];
sx q[2];
rz(-2.1047154) q[2];
sx q[2];
rz(-0.18661737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1101802) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(-0.024755342) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39748945) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(-1.4504364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.7987569) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(-1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.36195155) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(-2.4848293) q[0];
rz(1.1991869) q[2];
sx q[2];
rz(-1.2174165) q[2];
sx q[2];
rz(1.6824739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(0.025440865) q[1];
x q[2];
rz(-0.36985107) q[3];
sx q[3];
rz(-2.5858013) q[3];
sx q[3];
rz(-2.7694447) q[3];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(0.68583268) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(-0.54840666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054664139) q[0];
sx q[0];
rz(-1.5318649) q[0];
sx q[0];
rz(-0.995308) q[0];
rz(-pi) q[1];
rz(-0.19210179) q[2];
sx q[2];
rz(-0.77741277) q[2];
sx q[2];
rz(-0.32595134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1141495) q[1];
sx q[1];
rz(-0.47361923) q[1];
sx q[1];
rz(-1.4373535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3324276) q[3];
sx q[3];
rz(-2.7276464) q[3];
sx q[3];
rz(0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(0.87695688) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(-0.99463314) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431865) q[0];
sx q[0];
rz(-0.61265677) q[0];
sx q[0];
rz(-1.1820656) q[0];
rz(0.0066130916) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(2.7343482) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.4410938) q[3];
sx q[3];
rz(-1.8860215) q[3];
sx q[3];
rz(-1.3198927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2991335) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-2.6966406) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(1.2114552) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.6192186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30949621) q[1];
sx q[1];
rz(-2.2239904) q[1];
sx q[1];
rz(-0.081835882) q[1];
rz(-pi) q[2];
rz(-0.92618561) q[3];
sx q[3];
rz(-1.1759967) q[3];
sx q[3];
rz(1.3046164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(0.04315367) q[3];
sx q[3];
rz(-1.4814601) q[3];
sx q[3];
rz(-2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1535004) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(-2.5832672) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(-2.4024898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16427134) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(3.0493899) q[0];
rz(-3.1355255) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(-2.9059682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6369789) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(-1.3330589) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3665175) q[3];
sx q[3];
rz(-0.26976997) q[3];
sx q[3];
rz(-1.1360053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(-0.39673355) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(-2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(-0.5272665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11589719) q[0];
sx q[0];
rz(-0.39931116) q[0];
sx q[0];
rz(-2.3621109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69006069) q[2];
sx q[2];
rz(-1.0933211) q[2];
sx q[2];
rz(0.32327393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3026915) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(2.332815) q[1];
rz(-pi) q[2];
rz(0.91068565) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(-1.1235352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-2.6623181) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(0.17818174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885926) q[0];
sx q[0];
rz(-1.0646001) q[0];
sx q[0];
rz(-0.83336713) q[0];
rz(-2.4536132) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(-2.2945987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61160117) q[1];
sx q[1];
rz(-1.0582663) q[1];
sx q[1];
rz(-0.26228735) q[1];
rz(1.6931157) q[3];
sx q[3];
rz(-1.1432075) q[3];
sx q[3];
rz(1.7784255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2330033) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(-1.2188777) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(-1.490996) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
rz(-2.8015295) q[3];
sx q[3];
rz(-0.9763413) q[3];
sx q[3];
rz(-1.7819596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];