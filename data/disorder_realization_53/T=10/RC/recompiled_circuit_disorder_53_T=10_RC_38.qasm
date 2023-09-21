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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.792359) q[0];
sx q[0];
rz(-1.4548737) q[0];
sx q[0];
rz(-0.16616343) q[0];
rz(2.8785107) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(-2.6334327) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(0.82378806) q[1];
x q[2];
rz(1.3710255) q[3];
sx q[3];
rz(-1.7598745) q[3];
sx q[3];
rz(3.004899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17125601) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(2.7217857) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(1.0158687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37594721) q[0];
sx q[0];
rz(-2.489466) q[0];
sx q[0];
rz(-0.69407065) q[0];
x q[1];
rz(-2.303896) q[2];
sx q[2];
rz(-0.3028377) q[2];
sx q[2];
rz(2.8267415) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3381172) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(0.34715279) q[1];
rz(-pi) q[2];
rz(-2.696051) q[3];
sx q[3];
rz(-2.0317151) q[3];
sx q[3];
rz(-1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-2.8217227) q[2];
sx q[2];
rz(2.0324198) q[2];
rz(2.7099113) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9839086) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(1.2616715) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7648375) q[0];
sx q[0];
rz(-1.1163045) q[0];
sx q[0];
rz(0.48671133) q[0];
x q[1];
rz(2.491465) q[2];
sx q[2];
rz(-2.1047154) q[2];
sx q[2];
rz(2.9549753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5955334) q[1];
sx q[1];
rz(-1.5946348) q[1];
sx q[1];
rz(1.2977352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4811901) q[3];
sx q[3];
rz(-1.378873) q[3];
sx q[3];
rz(-0.47062518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(1.7987569) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(-1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(2.9825488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7796411) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(-0.65676332) q[0];
x q[1];
rz(-0.37695388) q[2];
sx q[2];
rz(-1.2231584) q[2];
sx q[2];
rz(2.8958547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4403968) q[1];
sx q[1];
rz(-0.21322589) q[1];
sx q[1];
rz(1.6885944) q[1];
rz(-0.52491297) q[3];
sx q[3];
rz(-1.7626926) q[3];
sx q[3];
rz(1.6247941) q[3];
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
rz(-3.1221636) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6762125) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(-0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(2.593186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5655366) q[0];
sx q[0];
rz(-0.57665529) q[0];
sx q[0];
rz(1.6422436) q[0];
rz(-pi) q[1];
rz(1.3850645) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(0.059543691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9644331) q[1];
sx q[1];
rz(-2.0398643) q[1];
sx q[1];
rz(3.0735077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10336419) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(-2.7523224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-2.8319478) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10228957) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-2.4772494) q[1];
sx q[1];
rz(2.1469595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69840616) q[0];
sx q[0];
rz(-0.61265677) q[0];
sx q[0];
rz(1.959527) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50350952) q[2];
sx q[2];
rz(-1.5739872) q[2];
sx q[2];
rz(-1.1577595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0108346) q[1];
sx q[1];
rz(-0.33797435) q[1];
sx q[1];
rz(-1.5513175) q[1];
x q[2];
rz(-1.1676222) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(-2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(0.4449521) q[2];
rz(-0.80727243) q[3];
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
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31535661) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(-0.89609599) q[0];
rz(-2.232146) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(-3.0665841) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.6336332) q[0];
sx q[0];
rz(-1.2114552) q[0];
rz(-pi) q[1];
rz(2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.6192186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1753462) q[1];
sx q[1];
rz(-2.4840377) q[1];
sx q[1];
rz(-1.6772126) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(2.9748777) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.2274851) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-2.9149122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.84086001) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16427134) q[0];
sx q[0];
rz(-2.5760206) q[0];
sx q[0];
rz(-3.0493899) q[0];
x q[1];
rz(1.1147538) q[2];
sx q[2];
rz(-1.5653492) q[2];
sx q[2];
rz(-1.3378439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5046138) q[1];
sx q[1];
rz(-0.74090545) q[1];
sx q[1];
rz(-1.3330589) q[1];
rz(0.056034485) q[3];
sx q[3];
rz(-1.8348215) q[3];
sx q[3];
rz(0.92428401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(0.6189515) q[0];
rz(-1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(-2.6143262) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.4591044) q[2];
sx q[2];
rz(-2.3253369) q[2];
sx q[2];
rz(1.3860821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0638709) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(-1.0788171) q[1];
rz(-pi) q[2];
rz(1.2654952) q[3];
sx q[3];
rz(-2.4584208) q[3];
sx q[3];
rz(2.9340569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-2.5274726) q[3];
sx q[3];
rz(1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(0.17818174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(0.88130086) q[0];
rz(-pi) q[1];
rz(-0.68797942) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(-0.84699398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0514959) q[1];
sx q[1];
rz(-1.7987383) q[1];
sx q[1];
rz(2.0983314) q[1];
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
rz(-0.46897727) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(-3.0205884) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-3.1115816) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
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
