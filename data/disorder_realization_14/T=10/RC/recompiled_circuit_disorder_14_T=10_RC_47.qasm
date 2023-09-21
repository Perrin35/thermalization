OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(0.88357893) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1238026) q[0];
sx q[0];
rz(-1.666981) q[0];
sx q[0];
rz(-2.9642446) q[0];
rz(-1.892957) q[2];
sx q[2];
rz(-2.319921) q[2];
sx q[2];
rz(-1.8495714) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8625496) q[1];
sx q[1];
rz(-1.057813) q[1];
sx q[1];
rz(1.951745) q[1];
rz(-pi) q[2];
rz(-3.0953193) q[3];
sx q[3];
rz(-1.7414021) q[3];
sx q[3];
rz(1.5266649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1951695) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(2.9602489) q[2];
rz(0.26120734) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(0.38683495) q[0];
rz(-0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(-1.5418672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72774784) q[0];
sx q[0];
rz(-1.1607329) q[0];
sx q[0];
rz(-0.64322612) q[0];
rz(-pi) q[1];
rz(-2.9773211) q[2];
sx q[2];
rz(-2.0220244) q[2];
sx q[2];
rz(0.8115561) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35040417) q[1];
sx q[1];
rz(-1.5807307) q[1];
sx q[1];
rz(0.03573461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1003483) q[3];
sx q[3];
rz(-2.0083661) q[3];
sx q[3];
rz(-0.13089422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.4146457) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(1.8521076) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(-0.99197018) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6863166) q[0];
sx q[0];
rz(-1.0193829) q[0];
sx q[0];
rz(-0.21689143) q[0];
rz(-pi) q[1];
rz(-2.473258) q[2];
sx q[2];
rz(-1.4112345) q[2];
sx q[2];
rz(0.72232407) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2838973) q[1];
sx q[1];
rz(-1.1313492) q[1];
sx q[1];
rz(1.4717913) q[1];
rz(-pi) q[2];
rz(1.8157585) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(0.20416343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38301864) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(2.7475083) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2878993) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(-0.40507856) q[0];
rz(2.4515117) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(0.69782034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6676184) q[0];
sx q[0];
rz(-1.4809429) q[0];
sx q[0];
rz(3.1180179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5562416) q[2];
sx q[2];
rz(-2.5132837) q[2];
sx q[2];
rz(-1.6274239) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7792369) q[1];
sx q[1];
rz(-0.52473611) q[1];
sx q[1];
rz(1.1974105) q[1];
rz(-2.844595) q[3];
sx q[3];
rz(-1.9747435) q[3];
sx q[3];
rz(-1.7904074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71965376) q[2];
sx q[2];
rz(-1.7439338) q[2];
sx q[2];
rz(-1.5092124) q[2];
rz(2.737282) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25800911) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-1.0255381) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(-2.5122723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176661) q[0];
sx q[0];
rz(-1.9804945) q[0];
sx q[0];
rz(-1.2439338) q[0];
rz(-0.89485456) q[2];
sx q[2];
rz(-1.370315) q[2];
sx q[2];
rz(1.5948053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7484666) q[1];
sx q[1];
rz(-1.4561635) q[1];
sx q[1];
rz(2.6266891) q[1];
rz(-1.945433) q[3];
sx q[3];
rz(-2.5064838) q[3];
sx q[3];
rz(2.0600968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(2.9455744) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65258566) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(1.7350896) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(-1.3669744) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8946675) q[0];
sx q[0];
rz(-1.7917544) q[0];
sx q[0];
rz(0.59452438) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.235249) q[2];
sx q[2];
rz(-2.3830072) q[2];
sx q[2];
rz(2.2461265) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1176599) q[1];
sx q[1];
rz(-1.4791094) q[1];
sx q[1];
rz(-1.4183284) q[1];
x q[2];
rz(2.572445) q[3];
sx q[3];
rz(-1.2043118) q[3];
sx q[3];
rz(-1.6397427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(-2.2582167) q[2];
rz(-1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25154034) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(-1.2782156) q[0];
rz(-3.1037519) q[1];
sx q[1];
rz(-0.81532878) q[1];
sx q[1];
rz(1.3758804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2863231) q[0];
sx q[0];
rz(-2.3648242) q[0];
sx q[0];
rz(1.9957248) q[0];
rz(1.6872348) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(1.4637228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.059767698) q[1];
sx q[1];
rz(-1.4653112) q[1];
sx q[1];
rz(-3.0256773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4010299) q[3];
sx q[3];
rz(-0.3347291) q[3];
sx q[3];
rz(-1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4574796) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(2.4462162) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(2.4080283) q[0];
rz(0.60797524) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(2.9072993) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83345862) q[0];
sx q[0];
rz(-1.8292384) q[0];
sx q[0];
rz(-1.0779557) q[0];
x q[1];
rz(0.90227622) q[2];
sx q[2];
rz(-0.96792816) q[2];
sx q[2];
rz(0.55842802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1842321) q[1];
sx q[1];
rz(-0.82416526) q[1];
sx q[1];
rz(-2.8824473) q[1];
rz(0.33741823) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12604788) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.4253433) q[2];
rz(-1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(-1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54206806) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(0.9448005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244922) q[0];
sx q[0];
rz(-1.1860577) q[0];
sx q[0];
rz(2.4186633) q[0];
rz(-0.94061942) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(2.1911052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33465696) q[1];
sx q[1];
rz(-1.2922704) q[1];
sx q[1];
rz(-1.7486228) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9444669) q[3];
sx q[3];
rz(-1.2211495) q[3];
sx q[3];
rz(-2.8417348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(-2.8533868) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(1.5079927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11186803) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(-2.2286041) q[0];
rz(-0.37462014) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(0.8909117) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64542949) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(-1.8339001) q[0];
x q[1];
rz(1.1773445) q[2];
sx q[2];
rz(-0.65044728) q[2];
sx q[2];
rz(-2.8124867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76162321) q[1];
sx q[1];
rz(-2.7459811) q[1];
sx q[1];
rz(-0.53669866) q[1];
rz(0.18817801) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(2.1256413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9051819) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-0.51433993) q[2];
sx q[2];
rz(-0.30403501) q[2];
sx q[2];
rz(-2.4652849) q[2];
rz(2.1174259) q[3];
sx q[3];
rz(-1.6129941) q[3];
sx q[3];
rz(2.2254406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];