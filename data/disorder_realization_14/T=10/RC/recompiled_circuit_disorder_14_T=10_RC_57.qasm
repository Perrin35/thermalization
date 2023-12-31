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
rz(-2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(4.598773) q[1];
sx q[1];
rz(7.4814319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449013) q[0];
sx q[0];
rz(-2.9400819) q[0];
sx q[0];
rz(-2.6411396) q[0];
x q[1];
rz(1.892957) q[2];
sx q[2];
rz(-0.82167168) q[2];
sx q[2];
rz(1.2920213) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4858422) q[1];
sx q[1];
rz(-1.2409004) q[1];
sx q[1];
rz(0.54539036) q[1];
rz(-1.4000113) q[3];
sx q[3];
rz(-1.6163974) q[3];
sx q[3];
rz(0.036269773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1951695) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(0.18134376) q[2];
rz(-0.26120734) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-2.3852824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3023892) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4138448) q[0];
sx q[0];
rz(-1.1607329) q[0];
sx q[0];
rz(0.64322612) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2453169) q[2];
sx q[2];
rz(-2.6633334) q[2];
sx q[2];
rz(1.9667728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6501573) q[1];
sx q[1];
rz(-0.037089247) q[1];
sx q[1];
rz(-2.8703719) q[1];
x q[2];
rz(0.48331355) q[3];
sx q[3];
rz(-1.1477074) q[3];
sx q[3];
rz(1.9139293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60275045) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.4146457) q[2];
rz(-0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(-1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8787815) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(-1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(2.1496225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2881644) q[0];
sx q[0];
rz(-0.58840226) q[0];
sx q[0];
rz(-1.2342274) q[0];
x q[1];
rz(0.66833468) q[2];
sx q[2];
rz(-1.7303581) q[2];
sx q[2];
rz(2.4192686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5132644) q[1];
sx q[1];
rz(-2.6918415) q[1];
sx q[1];
rz(-0.20723923) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8157585) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(2.9374292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(0.26829159) q[2];
rz(-2.7475083) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85369337) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-0.40507856) q[0];
rz(2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(2.4437723) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0426548) q[0];
sx q[0];
rz(-1.5473167) q[0];
sx q[0];
rz(-1.6606746) q[0];
rz(-2.5970846) q[2];
sx q[2];
rz(-1.240057) q[2];
sx q[2];
rz(2.7059908) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2600425) q[1];
sx q[1];
rz(-1.7545732) q[1];
sx q[1];
rz(-2.0651414) q[1];
x q[2];
rz(0.29699765) q[3];
sx q[3];
rz(-1.1668491) q[3];
sx q[3];
rz(1.7904074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71965376) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.5092124) q[2];
rz(-2.737282) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(-1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(2.1160545) q[0];
rz(0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(0.62932032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176661) q[0];
sx q[0];
rz(-1.9804945) q[0];
sx q[0];
rz(-1.2439338) q[0];
rz(-pi) q[1];
rz(0.25482486) q[2];
sx q[2];
rz(-0.91081589) q[2];
sx q[2];
rz(-3.0072336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3931261) q[1];
sx q[1];
rz(-1.6854291) q[1];
sx q[1];
rz(2.6266891) q[1];
x q[2];
rz(-2.1719645) q[3];
sx q[3];
rz(-1.789635) q[3];
sx q[3];
rz(-2.9588483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59262529) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65258566) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(-0.18519369) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(-1.3669744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8946675) q[0];
sx q[0];
rz(-1.7917544) q[0];
sx q[0];
rz(-0.59452438) q[0];
x q[1];
rz(-2.8390859) q[2];
sx q[2];
rz(-0.86382804) q[2];
sx q[2];
rz(-0.44755852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1320912) q[1];
sx q[1];
rz(-2.9638634) q[1];
sx q[1];
rz(-2.1151665) q[1];
rz(-pi) q[2];
rz(-0.6188789) q[3];
sx q[3];
rz(-2.4757909) q[3];
sx q[3];
rz(-2.562079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25154034) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.7657123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0276887) q[0];
sx q[0];
rz(-1.8639599) q[0];
sx q[0];
rz(-2.3011074) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0443633) q[2];
sx q[2];
rz(-1.5175022) q[2];
sx q[2];
rz(3.1380944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3657276) q[1];
sx q[1];
rz(-0.15656808) q[1];
sx q[1];
rz(-0.74128976) q[1];
rz(-pi) q[2];
rz(3.0828956) q[3];
sx q[3];
rz(-1.241063) q[3];
sx q[3];
rz(1.2122455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4574796) q[2];
sx q[2];
rz(-2.423968) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(2.4080283) q[0];
rz(-2.5336174) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(2.9072993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5406571) q[0];
sx q[0];
rz(-2.0458851) q[0];
sx q[0];
rz(-0.29151543) q[0];
rz(2.2393164) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(-2.5831646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.329511) q[1];
sx q[1];
rz(-2.3595464) q[1];
sx q[1];
rz(-1.8409607) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0571363) q[3];
sx q[3];
rz(-0.64389766) q[3];
sx q[3];
rz(-0.38734303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.4253433) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.4956168) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(-1.9482127) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(-2.1967922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.955907) q[0];
sx q[0];
rz(-0.80230306) q[0];
sx q[0];
rz(-2.5923652) q[0];
x q[1];
rz(-2.8464727) q[2];
sx q[2];
rz(-1.1913538) q[2];
sx q[2];
rz(1.6377359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2855125) q[1];
sx q[1];
rz(-1.3998919) q[1];
sx q[1];
rz(0.28275615) q[1];
x q[2];
rz(0.19712574) q[3];
sx q[3];
rz(-1.2211495) q[3];
sx q[3];
rz(-0.2998578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21645674) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(-2.6618585) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11186803) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(-0.9129886) q[0];
rz(-0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(2.250681) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4961632) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(-1.8339001) q[0];
x q[1];
rz(1.1773445) q[2];
sx q[2];
rz(-2.4911454) q[2];
sx q[2];
rz(2.8124867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76162321) q[1];
sx q[1];
rz(-2.7459811) q[1];
sx q[1];
rz(2.604894) q[1];
rz(-pi) q[2];
rz(0.37836214) q[3];
sx q[3];
rz(-2.9394657) q[3];
sx q[3];
rz(-2.9581021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9051819) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(2.6399844) q[2];
rz(-1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.4176462) q[2];
sx q[2];
rz(-1.8344804) q[2];
sx q[2];
rz(-1.9305965) q[2];
rz(-0.049384762) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
