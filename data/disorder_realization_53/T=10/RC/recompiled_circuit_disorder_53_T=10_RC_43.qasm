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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38219163) q[0];
sx q[0];
rz(-0.2022976) q[0];
sx q[0];
rz(2.5281744) q[0];
rz(-1.4853391) q[2];
sx q[2];
rz(-1.308631) q[2];
sx q[2];
rz(-2.0567577) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4781293) q[1];
sx q[1];
rz(-0.98065286) q[1];
sx q[1];
rz(-0.66901916) q[1];
x q[2];
rz(-1.3710255) q[3];
sx q[3];
rz(-1.7598745) q[3];
sx q[3];
rz(-3.004899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7621883) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(-2.7963426) q[2];
rz(-0.18940645) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-2.4617564) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(2.7833815) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(-1.0158687) q[1];
rz(-pi) q[2];
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
rz(1.3426571) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(-1.1756431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7125883) q[1];
sx q[1];
rz(-1.9468369) q[1];
sx q[1];
rz(1.4274548) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2855371) q[3];
sx q[3];
rz(-0.62969172) q[3];
sx q[3];
rz(0.39470181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-2.0324198) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-0.3698529) q[3];
sx q[3];
rz(-3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576841) q[0];
sx q[0];
rz(-1.8439872) q[0];
sx q[0];
rz(-0.88009673) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42230168) q[0];
sx q[0];
rz(-1.137086) q[0];
sx q[0];
rz(2.0757872) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9319839) q[2];
sx q[2];
rz(-1.0227232) q[2];
sx q[2];
rz(1.0149479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.060170598) q[1];
sx q[1];
rz(-0.27407384) q[1];
sx q[1];
rz(1.6589792) q[1];
rz(-2.7441032) q[3];
sx q[3];
rz(-2.6263413) q[3];
sx q[3];
rz(1.6911563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1168388) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.8163619) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(-2.0480115) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4745859) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(-2.4530607) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(0.15904388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7796411) q[0];
sx q[0];
rz(-2.1719296) q[0];
sx q[0];
rz(2.4848293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1991869) q[2];
sx q[2];
rz(-1.2174165) q[2];
sx q[2];
rz(-1.6824739) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-0.025440865) q[1];
rz(-pi) q[2];
rz(-2.6166797) q[3];
sx q[3];
rz(-1.3789) q[3];
sx q[3];
rz(-1.5167985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8461385) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(2.7159178) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-0.21605505) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.54840666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6507051) q[0];
sx q[0];
rz(-2.1457932) q[0];
sx q[0];
rz(0.046396359) q[0];
x q[1];
rz(-0.76813278) q[2];
sx q[2];
rz(-1.4364725) q[2];
sx q[2];
rz(1.3825934) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4244528) q[1];
sx q[1];
rz(-1.5100749) q[1];
sx q[1];
rz(-2.0408003) q[1];
rz(-pi) q[2];
rz(-1.9742825) q[3];
sx q[3];
rz(-1.6659123) q[3];
sx q[3];
rz(2.000589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(0.99463314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1625527) q[0];
sx q[0];
rz(-1.0096706) q[0];
sx q[0];
rz(-2.8812376) q[0];
rz(-1.5744393) q[2];
sx q[2];
rz(-1.0672896) q[2];
sx q[2];
rz(2.7267981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1514046) q[1];
sx q[1];
rz(-1.908704) q[1];
sx q[1];
rz(-0.0068454725) q[1];
rz(-0.46835132) q[3];
sx q[3];
rz(-2.3845209) q[3];
sx q[3];
rz(-3.0401331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8424592) q[2];
sx q[2];
rz(-2.2392539) q[2];
sx q[2];
rz(2.6966406) q[2];
rz(-2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.8912631) q[1];
sx q[1];
rz(-3.0665841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1589926) q[0];
sx q[0];
rz(-1.9293961) q[0];
sx q[0];
rz(-0.067111777) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5957017) q[2];
sx q[2];
rz(-2.2635169) q[2];
sx q[2];
rz(2.8199414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1753462) q[1];
sx q[1];
rz(-0.65755492) q[1];
sx q[1];
rz(-1.46438) q[1];
x q[2];
rz(0.9644896) q[3];
sx q[3];
rz(-2.4006667) q[3];
sx q[3];
rz(0.20674202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3007326) q[1];
sx q[1];
rz(0.73910284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571592) q[0];
sx q[0];
rz(-1.6201577) q[0];
sx q[0];
rz(-2.5779448) q[0];
rz(-pi) q[1];
rz(-0.0060672005) q[2];
sx q[2];
rz(-2.0268315) q[2];
sx q[2];
rz(-2.9059682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5046138) q[1];
sx q[1];
rz(-2.4006872) q[1];
sx q[1];
rz(-1.3330589) q[1];
rz(-pi) q[2];
rz(1.8352175) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(2.5097178) q[3];
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
rz(-1.1018264) q[2];
rz(0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(0.81469369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(1.0194107) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(0.5272665) q[1];
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
x q[1];
rz(-2.4591044) q[2];
sx q[2];
rz(-0.81625578) q[2];
sx q[2];
rz(-1.3860821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83890115) q[1];
sx q[1];
rz(-2.5040002) q[1];
sx q[1];
rz(-2.332815) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2654952) q[3];
sx q[3];
rz(-2.4584208) q[3];
sx q[3];
rz(0.20753577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(0.81106097) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(-0.47927454) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(2.9634109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44706599) q[0];
sx q[0];
rz(-2.1994626) q[0];
sx q[0];
rz(0.64283128) q[0];
rz(-pi) q[1];
rz(0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(0.40643613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0900967) q[1];
sx q[1];
rz(-1.3428543) q[1];
sx q[1];
rz(-1.0432613) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7111669) q[3];
sx q[3];
rz(-1.6820551) q[3];
sx q[3];
rz(0.15669565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-2.4185833) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(3.0205884) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(0.96414375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.9299097) q[2];
sx q[2];
rz(-1.5426987) q[2];
sx q[2];
rz(0.60088746) q[2];
rz(-2.1929019) q[3];
sx q[3];
rz(-1.2908251) q[3];
sx q[3];
rz(2.7348107) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
