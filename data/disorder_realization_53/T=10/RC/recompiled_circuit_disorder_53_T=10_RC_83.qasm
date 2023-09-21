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
rz(-0.53108162) q[1];
sx q[1];
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.792359) q[0];
sx q[0];
rz(-1.6867189) q[0];
sx q[0];
rz(2.9754292) q[0];
rz(-pi) q[1];
rz(0.30795745) q[2];
sx q[2];
rz(-2.8661558) q[2];
sx q[2];
rz(0.76560417) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4781293) q[1];
sx q[1];
rz(-2.1609398) q[1];
sx q[1];
rz(-0.66901916) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19282135) q[3];
sx q[3];
rz(-1.76696) q[3];
sx q[3];
rz(1.7455268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-0.72221243) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(-2.125724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7091539) q[0];
sx q[0];
rz(-2.0560987) q[0];
sx q[0];
rz(-1.1164467) q[0];
rz(-pi) q[1];
rz(2.935479) q[2];
sx q[2];
rz(-1.7942675) q[2];
sx q[2];
rz(-2.7001365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42900436) q[1];
sx q[1];
rz(-1.9468369) q[1];
sx q[1];
rz(-1.4274548) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44554168) q[3];
sx q[3];
rz(-1.1098776) q[3];
sx q[3];
rz(1.9259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(0.88009673) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719291) q[0];
sx q[0];
rz(-2.0045067) q[0];
sx q[0];
rz(-1.0658054) q[0];
x q[1];
rz(-2.2096087) q[2];
sx q[2];
rz(-2.1188695) q[2];
sx q[2];
rz(-1.0149479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0814221) q[1];
sx q[1];
rz(-0.27407384) q[1];
sx q[1];
rz(-1.6589792) q[1];
rz(-pi) q[2];
rz(2.7441032) q[3];
sx q[3];
rz(-0.51525138) q[3];
sx q[3];
rz(1.6911563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.024753831) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.8163619) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.7929247) q[1];
sx q[1];
rz(0.15904388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36195155) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(0.65676332) q[0];
x q[1];
rz(-0.37695388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(0.24573791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70119584) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(1.6885944) q[1];
rz(-pi) q[2];
rz(-0.52491297) q[3];
sx q[3];
rz(-1.3789) q[3];
sx q[3];
rz(1.5167985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(-0.42567483) q[2];
rz(3.1221636) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(-2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-3.0131969) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(2.593186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5655366) q[0];
sx q[0];
rz(-2.5649374) q[0];
sx q[0];
rz(1.6422436) q[0];
rz(-pi) q[1];
rz(0.19210179) q[2];
sx q[2];
rz(-0.77741277) q[2];
sx q[2];
rz(0.32595134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9644331) q[1];
sx q[1];
rz(-1.1017283) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(-0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(2.1469595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69840616) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(-1.1820656) q[0];
rz(-3.1349796) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(2.7343482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1307581) q[1];
sx q[1];
rz(-2.8036183) q[1];
sx q[1];
rz(1.5513175) q[1];
x q[2];
rz(-1.9739705) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(-0.80727243) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(-0.81594938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(-2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(0.075008579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61178111) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(1.9301374) q[0];
rz(-pi) q[1];
rz(-0.80008614) q[2];
sx q[2];
rz(-1.1598088) q[2];
sx q[2];
rz(-1.5223741) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.311104) q[1];
sx q[1];
rz(-1.6357592) q[1];
sx q[1];
rz(-0.91598367) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48052629) q[3];
sx q[3];
rz(-0.98283813) q[3];
sx q[3];
rz(-0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16671495) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(1.2274851) q[2];
rz(0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(0.41326997) q[0];
rz(-2.5832672) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-0.73910284) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055187125) q[0];
sx q[0];
rz(-2.1336745) q[0];
sx q[0];
rz(1.5124209) q[0];
x q[1];
rz(-3.1355255) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(0.23562442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5046138) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(-1.3330589) q[1];
x q[2];
rz(1.8352175) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(2.5097178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(-2.0397662) q[2];
rz(0.39673355) q[3];
sx q[3];
rz(-1.4361897) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(-0.6189515) q[0];
rz(-2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(-2.6143262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11589719) q[0];
sx q[0];
rz(-2.7422815) q[0];
sx q[0];
rz(-2.3621109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68248827) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.7555106) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83890115) q[1];
sx q[1];
rz(-0.63759241) q[1];
sx q[1];
rz(-2.332815) q[1];
x q[2];
rz(1.2654952) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8573389) q[2];
sx q[2];
rz(-1.8062091) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(-2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(2.6623181) q[0];
rz(-0.88538623) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-0.17818174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885926) q[0];
sx q[0];
rz(-2.0769925) q[0];
sx q[0];
rz(0.83336713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(0.40643613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0308932) q[1];
sx q[1];
rz(-2.5712214) q[1];
sx q[1];
rz(2.0026026) q[1];
rz(-0.26161216) q[3];
sx q[3];
rz(-2.6978921) q[3];
sx q[3];
rz(1.651368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2330033) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(-1.490996) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
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
