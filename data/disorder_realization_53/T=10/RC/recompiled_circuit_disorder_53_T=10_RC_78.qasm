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
rz(-0.53108162) q[1];
sx q[1];
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.792359) q[0];
sx q[0];
rz(-1.6867189) q[0];
sx q[0];
rz(-2.9754292) q[0];
rz(-pi) q[1];
rz(-0.30795745) q[2];
sx q[2];
rz(-2.8661558) q[2];
sx q[2];
rz(2.3759885) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(-0.82378806) q[1];
x q[2];
rz(1.3710255) q[3];
sx q[3];
rz(-1.3817182) q[3];
sx q[3];
rz(-3.004899) q[3];
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
rz(-2.7963426) q[2];
rz(2.9521862) q[3];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-0.41980699) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(2.125724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37594721) q[0];
sx q[0];
rz(-2.489466) q[0];
sx q[0];
rz(-2.447522) q[0];
rz(-pi) q[1];
rz(-1.7989356) q[2];
sx q[2];
rz(-1.7717138) q[2];
sx q[2];
rz(1.1756431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42900436) q[1];
sx q[1];
rz(-1.1947558) q[1];
sx q[1];
rz(1.4274548) q[1];
rz(-pi) q[2];
rz(-1.0677098) q[3];
sx q[3];
rz(-1.1745319) q[3];
sx q[3];
rz(0.5644507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(2.0324198) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.2976054) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(0.18149158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49849579) q[0];
sx q[0];
rz(-2.4884014) q[0];
sx q[0];
rz(-0.80723395) q[0];
x q[1];
rz(2.491465) q[2];
sx q[2];
rz(-2.1047154) q[2];
sx q[2];
rz(2.9549753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.060170598) q[1];
sx q[1];
rz(-2.8675188) q[1];
sx q[1];
rz(-1.6589792) q[1];
rz(-pi) q[2];
rz(0.39748945) q[3];
sx q[3];
rz(-0.51525138) q[3];
sx q[3];
rz(-1.6911563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(0.68853199) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(-2.9825488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36195155) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(2.4848293) q[0];
rz(-pi) q[1];
rz(0.7775457) q[2];
sx q[2];
rz(-2.6345207) q[2];
sx q[2];
rz(0.61447243) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58069431) q[1];
sx q[1];
rz(-1.3590707) q[1];
sx q[1];
rz(3.1161518) q[1];
x q[2];
rz(-0.52491297) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(1.6247941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(-3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46538019) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(2.45576) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6507051) q[0];
sx q[0];
rz(-2.1457932) q[0];
sx q[0];
rz(3.0951963) q[0];
rz(-2.9494909) q[2];
sx q[2];
rz(-2.3641799) q[2];
sx q[2];
rz(-0.32595134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4244528) q[1];
sx q[1];
rz(-1.6315178) q[1];
sx q[1];
rz(1.1007924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10336419) q[3];
sx q[3];
rz(-1.9723537) q[3];
sx q[3];
rz(-2.7523224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(-1.6748641) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(2.2646358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(0.99463314) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431865) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(1.1820656) q[0];
rz(0.50350952) q[2];
sx q[2];
rz(-1.5676055) q[2];
sx q[2];
rz(-1.9838331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1514046) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(3.1347472) q[1];
rz(-2.6732413) q[3];
sx q[3];
rz(-2.3845209) q[3];
sx q[3];
rz(-0.10145951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(2.6966406) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(2.2454967) q[0];
rz(2.232146) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(-3.0665841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98260005) q[0];
sx q[0];
rz(-1.2121965) q[0];
sx q[0];
rz(-3.0744809) q[0];
x q[1];
rz(0.80008614) q[2];
sx q[2];
rz(-1.1598088) q[2];
sx q[2];
rz(1.5223741) q[2];
rz(-pi) q[3];
x q[3];
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
x q[2];
rz(-2.215407) q[3];
sx q[3];
rz(-1.1759967) q[3];
sx q[3];
rz(1.8369762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-2.0090943) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(-2.7283227) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-0.84086001) q[1];
sx q[1];
rz(2.4024898) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16427134) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(-0.092202734) q[0];
rz(-pi) q[1];
rz(1.5584281) q[2];
sx q[2];
rz(-0.45607273) q[2];
sx q[2];
rz(0.2218483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11074556) q[1];
sx q[1];
rz(-1.4111641) q[1];
sx q[1];
rz(-0.84407945) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.056034485) q[3];
sx q[3];
rz(-1.8348215) q[3];
sx q[3];
rz(2.2173086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3070613) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(-2.0397662) q[2];
rz(-0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.5588201) q[1];
sx q[1];
rz(2.6143262) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11589719) q[0];
sx q[0];
rz(-2.7422815) q[0];
sx q[0];
rz(0.77948178) q[0];
x q[1];
rz(-2.4591044) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.3860821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4317324) q[1];
sx q[1];
rz(-2.0159971) q[1];
sx q[1];
rz(-0.47275895) q[1];
rz(-pi) q[2];
rz(2.230907) q[3];
sx q[3];
rz(-1.3798957) q[3];
sx q[3];
rz(1.1235352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28425372) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(-2.2100892) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
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
rz(0.60944027) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(2.6623181) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(0.17818174) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(0.88130086) q[0];
rz(-2.7257596) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(-0.40643613) q[2];
rz(-pi/2) q[3];
sx q[3];
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
x q[2];
rz(-1.6931157) q[3];
sx q[3];
rz(-1.1432075) q[3];
sx q[3];
rz(1.3631671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.9299097) q[2];
sx q[2];
rz(-1.598894) q[2];
sx q[2];
rz(-2.5407052) q[2];
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
