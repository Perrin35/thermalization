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
rz(-2.1531818) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(1.3487863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006346) q[0];
sx q[0];
rz(-1.4057584) q[0];
sx q[0];
rz(1.4532695) q[0];
rz(-2.8785107) q[2];
sx q[2];
rz(-1.4882659) q[2];
sx q[2];
rz(-2.6334327) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(-2.3178046) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7705671) q[3];
sx q[3];
rz(-1.7598745) q[3];
sx q[3];
rz(0.13669361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(-0.96864831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43243877) q[0];
sx q[0];
rz(-2.0560987) q[0];
sx q[0];
rz(1.1164467) q[0];
rz(-pi) q[1];
rz(-2.935479) q[2];
sx q[2];
rz(-1.3473251) q[2];
sx q[2];
rz(-2.7001365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0527555) q[1];
sx q[1];
rz(-1.4375327) q[1];
sx q[1];
rz(2.7620176) q[1];
rz(-2.696051) q[3];
sx q[3];
rz(-1.1098776) q[3];
sx q[3];
rz(1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(2.0324198) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(-0.18149158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719291) q[0];
sx q[0];
rz(-2.0045067) q[0];
sx q[0];
rz(2.0757872) q[0];
rz(2.2096087) q[2];
sx q[2];
rz(-2.1188695) q[2];
sx q[2];
rz(1.0149479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.060170598) q[1];
sx q[1];
rz(-2.8675188) q[1];
sx q[1];
rz(1.4826135) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6604026) q[3];
sx q[3];
rz(-1.378873) q[3];
sx q[3];
rz(2.6709675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.7987569) q[2];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(-0.21903285) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(2.9825488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7796411) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(-0.65676332) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37695388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(-0.24573791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98475458) q[1];
sx q[1];
rz(-1.5459237) q[1];
sx q[1];
rz(-1.3590042) q[1];
rz(-pi) q[2];
rz(2.6166797) q[3];
sx q[3];
rz(-1.3789) q[3];
sx q[3];
rz(1.5167985) q[3];
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
rz(-0.019429026) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(-2.4470636) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6762125) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869285) q[0];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17715958) q[1];
sx q[1];
rz(-2.0398643) q[1];
sx q[1];
rz(0.068084929) q[1];
x q[2];
rz(1.1673101) q[3];
sx q[3];
rz(-1.6659123) q[3];
sx q[3];
rz(-1.1410037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8473062) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(2.2646358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(1.0992682) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.97904) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(-0.26035505) q[0];
rz(-pi) q[1];
rz(-3.1349796) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(2.7343482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1307581) q[1];
sx q[1];
rz(-2.8036183) q[1];
sx q[1];
rz(-1.5902751) q[1];
rz(-pi) q[2];
rz(-1.9739705) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(0.80727243) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(0.81594938) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-0.89609599) q[0];
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
rz(-2.1589926) q[0];
sx q[0];
rz(-1.9293961) q[0];
sx q[0];
rz(-0.067111777) q[0];
x q[1];
rz(2.5957017) q[2];
sx q[2];
rz(-2.2635169) q[2];
sx q[2];
rz(0.32165124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.311104) q[1];
sx q[1];
rz(-1.5058335) q[1];
sx q[1];
rz(-2.225609) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92618561) q[3];
sx q[3];
rz(-1.965596) q[3];
sx q[3];
rz(1.3046164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(0.55832541) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(-2.4024898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571592) q[0];
sx q[0];
rz(-1.5214349) q[0];
sx q[0];
rz(2.5779448) q[0];
x q[1];
rz(-3.1355255) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(0.23562442) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0308471) q[1];
sx q[1];
rz(-1.4111641) q[1];
sx q[1];
rz(2.2975132) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3063752) q[3];
sx q[3];
rz(-1.6248871) q[3];
sx q[3];
rz(-2.5097178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(2.0397662) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535764) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(-2.5226412) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(0.5272665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94811234) q[0];
sx q[0];
rz(-1.2939948) q[0];
sx q[0];
rz(-2.8500018) q[0];
rz(2.4591044) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(1.3860821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0638709) q[1];
sx q[1];
rz(-1.1472902) q[1];
sx q[1];
rz(2.0627756) q[1];
rz(-1.2654952) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(-0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
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
rz(-2.482174) q[1];
sx q[1];
rz(-0.17818174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3522211) q[0];
sx q[0];
rz(-2.2749315) q[0];
sx q[0];
rz(0.88130086) q[0];
rz(2.7257596) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(-2.7351565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0514959) q[1];
sx q[1];
rz(-1.7987383) q[1];
sx q[1];
rz(1.0432613) q[1];
rz(-pi) q[2];
rz(1.448477) q[3];
sx q[3];
rz(-1.9983851) q[3];
sx q[3];
rz(-1.3631671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(2.318312) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(1.9299097) q[2];
sx q[2];
rz(-1.5426987) q[2];
sx q[2];
rz(0.60088746) q[2];
rz(-0.34006313) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];