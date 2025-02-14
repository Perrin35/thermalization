OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29937509) q[0];
sx q[0];
rz(-2.8111281) q[0];
sx q[0];
rz(2.0781031) q[0];
rz(-0.039634135) q[1];
sx q[1];
rz(-0.57365817) q[1];
sx q[1];
rz(2.0613476) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40908694) q[0];
sx q[0];
rz(-1.5729781) q[0];
sx q[0];
rz(-1.1078784) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3990381) q[2];
sx q[2];
rz(-1.556201) q[2];
sx q[2];
rz(0.43806048) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0277522) q[1];
sx q[1];
rz(-0.58958399) q[1];
sx q[1];
rz(0.47926183) q[1];
x q[2];
rz(-2.4970666) q[3];
sx q[3];
rz(-1.9686507) q[3];
sx q[3];
rz(-0.81105876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0237191) q[2];
sx q[2];
rz(-1.7261852) q[2];
sx q[2];
rz(0.37303698) q[2];
rz(-2.5840664) q[3];
sx q[3];
rz(-1.6256465) q[3];
sx q[3];
rz(-0.47835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28906223) q[0];
sx q[0];
rz(-1.4339148) q[0];
sx q[0];
rz(1.0093932) q[0];
rz(-0.32762647) q[1];
sx q[1];
rz(-0.3265003) q[1];
sx q[1];
rz(2.0351298) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43531159) q[0];
sx q[0];
rz(-1.5967622) q[0];
sx q[0];
rz(-3.0089799) q[0];
rz(-pi) q[1];
rz(2.4566133) q[2];
sx q[2];
rz(-1.5826591) q[2];
sx q[2];
rz(1.4432743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4047186) q[1];
sx q[1];
rz(-1.1049184) q[1];
sx q[1];
rz(-0.89786713) q[1];
rz(-pi) q[2];
rz(-0.4119315) q[3];
sx q[3];
rz(-1.3917482) q[3];
sx q[3];
rz(1.826394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24836765) q[2];
sx q[2];
rz(-1.1838341) q[2];
sx q[2];
rz(-2.4719205) q[2];
rz(-1.9855965) q[3];
sx q[3];
rz(-2.7892734) q[3];
sx q[3];
rz(-0.34935752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30592331) q[0];
sx q[0];
rz(-2.7140129) q[0];
sx q[0];
rz(-1.8100354) q[0];
rz(1.2003027) q[1];
sx q[1];
rz(-0.2526865) q[1];
sx q[1];
rz(-2.4694064) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.975544) q[0];
sx q[0];
rz(-0.48915809) q[0];
sx q[0];
rz(0.25948317) q[0];
rz(-1.9750017) q[2];
sx q[2];
rz(-1.4160755) q[2];
sx q[2];
rz(0.2153309) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1609413) q[1];
sx q[1];
rz(-2.0066727) q[1];
sx q[1];
rz(1.7546045) q[1];
x q[2];
rz(1.6217885) q[3];
sx q[3];
rz(-1.4740406) q[3];
sx q[3];
rz(0.51550625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4868698) q[2];
sx q[2];
rz(-1.9630311) q[2];
sx q[2];
rz(-2.2910924) q[2];
rz(-2.1908098) q[3];
sx q[3];
rz(-0.50539223) q[3];
sx q[3];
rz(-2.2397485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.013407) q[0];
sx q[0];
rz(-0.81939092) q[0];
sx q[0];
rz(-2.4017781) q[0];
rz(0.20135227) q[1];
sx q[1];
rz(-1.9722152) q[1];
sx q[1];
rz(-2.9154725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660131) q[0];
sx q[0];
rz(-3.1034307) q[0];
sx q[0];
rz(1.1665795) q[0];
rz(-pi) q[1];
rz(-1.1925542) q[2];
sx q[2];
rz(-1.8762445) q[2];
sx q[2];
rz(-0.39337197) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12783957) q[1];
sx q[1];
rz(-2.3495449) q[1];
sx q[1];
rz(0.75551466) q[1];
x q[2];
rz(1.700541) q[3];
sx q[3];
rz(-1.1888388) q[3];
sx q[3];
rz(0.48473919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40253338) q[2];
sx q[2];
rz(-0.89503521) q[2];
sx q[2];
rz(-0.99739897) q[2];
rz(-0.41336695) q[3];
sx q[3];
rz(-1.2099384) q[3];
sx q[3];
rz(-1.0443784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1362374) q[0];
sx q[0];
rz(-2.4449466) q[0];
sx q[0];
rz(0.053939017) q[0];
rz(-2.9593762) q[1];
sx q[1];
rz(-0.91582623) q[1];
sx q[1];
rz(2.838476) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5526472) q[0];
sx q[0];
rz(-2.573138) q[0];
sx q[0];
rz(2.2442152) q[0];
rz(-pi) q[1];
rz(1.8218173) q[2];
sx q[2];
rz(-0.66749882) q[2];
sx q[2];
rz(-0.72161822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6408711) q[1];
sx q[1];
rz(-2.5815775) q[1];
sx q[1];
rz(1.9418094) q[1];
rz(-pi) q[2];
rz(-3.1040499) q[3];
sx q[3];
rz(-2.3830066) q[3];
sx q[3];
rz(-2.9428218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4431241) q[2];
sx q[2];
rz(-1.2681401) q[2];
sx q[2];
rz(-2.7847086) q[2];
rz(2.3667864) q[3];
sx q[3];
rz(-1.6322497) q[3];
sx q[3];
rz(2.7029412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88724941) q[0];
sx q[0];
rz(-1.9312504) q[0];
sx q[0];
rz(-2.3526225) q[0];
rz(-1.3428768) q[1];
sx q[1];
rz(-1.7306381) q[1];
sx q[1];
rz(2.3416187) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26316428) q[0];
sx q[0];
rz(-1.5845926) q[0];
sx q[0];
rz(0.8850125) q[0];
rz(-1.7943013) q[2];
sx q[2];
rz(-2.433521) q[2];
sx q[2];
rz(-2.023165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50417954) q[1];
sx q[1];
rz(-0.39200726) q[1];
sx q[1];
rz(2.8543695) q[1];
rz(-pi) q[2];
rz(0.81963934) q[3];
sx q[3];
rz(-2.5878518) q[3];
sx q[3];
rz(0.40938974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.051290123) q[2];
sx q[2];
rz(-2.3919969) q[2];
sx q[2];
rz(-1.7898111) q[2];
rz(2.6734062) q[3];
sx q[3];
rz(-2.0464094) q[3];
sx q[3];
rz(2.2075672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4139597) q[0];
sx q[0];
rz(-0.41828823) q[0];
sx q[0];
rz(-3.1391414) q[0];
rz(3.1353503) q[1];
sx q[1];
rz(-2.7793482) q[1];
sx q[1];
rz(-2.7856766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2438176) q[0];
sx q[0];
rz(-2.1718302) q[0];
sx q[0];
rz(2.7181781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32940253) q[2];
sx q[2];
rz(-0.74879941) q[2];
sx q[2];
rz(-1.8210653) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1568875) q[1];
sx q[1];
rz(-0.92159373) q[1];
sx q[1];
rz(0.32905103) q[1];
rz(2.9690817) q[3];
sx q[3];
rz(-0.54076946) q[3];
sx q[3];
rz(-1.4494165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1459085) q[2];
sx q[2];
rz(-1.9084787) q[2];
sx q[2];
rz(-2.0705059) q[2];
rz(2.5339825) q[3];
sx q[3];
rz(-1.1039609) q[3];
sx q[3];
rz(1.6149909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0869658) q[0];
sx q[0];
rz(-0.93769756) q[0];
sx q[0];
rz(-1.1338393) q[0];
rz(-1.0519823) q[1];
sx q[1];
rz(-1.6832422) q[1];
sx q[1];
rz(2.4005311) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9423964) q[0];
sx q[0];
rz(-0.62672296) q[0];
sx q[0];
rz(-2.6882307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95693077) q[2];
sx q[2];
rz(-1.9939878) q[2];
sx q[2];
rz(-2.500071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3723124) q[1];
sx q[1];
rz(-1.0456632) q[1];
sx q[1];
rz(1.5430081) q[1];
rz(-0.40048845) q[3];
sx q[3];
rz(-0.97972673) q[3];
sx q[3];
rz(1.0210291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97003254) q[2];
sx q[2];
rz(-0.67983183) q[2];
sx q[2];
rz(0.46530923) q[2];
rz(-2.4252637) q[3];
sx q[3];
rz(-1.497044) q[3];
sx q[3];
rz(1.959645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5301836) q[0];
sx q[0];
rz(-0.94785988) q[0];
sx q[0];
rz(1.9635669) q[0];
rz(2.1317962) q[1];
sx q[1];
rz(-1.806908) q[1];
sx q[1];
rz(-0.10733265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862122) q[0];
sx q[0];
rz(-0.79609603) q[0];
sx q[0];
rz(2.436321) q[0];
rz(-pi) q[1];
rz(-0.1520098) q[2];
sx q[2];
rz(-0.78897023) q[2];
sx q[2];
rz(-2.0308354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3299249) q[1];
sx q[1];
rz(-0.98333849) q[1];
sx q[1];
rz(-1.2932375) q[1];
x q[2];
rz(-0.97300164) q[3];
sx q[3];
rz(-1.3241555) q[3];
sx q[3];
rz(1.1688978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44059077) q[2];
sx q[2];
rz(-0.89459449) q[2];
sx q[2];
rz(-2.7015298) q[2];
rz(-2.9317686) q[3];
sx q[3];
rz(-1.194229) q[3];
sx q[3];
rz(2.8275209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1485727) q[0];
sx q[0];
rz(-0.58605376) q[0];
sx q[0];
rz(-1.3826189) q[0];
rz(0.67715174) q[1];
sx q[1];
rz(-1.4644198) q[1];
sx q[1];
rz(1.1322397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099927946) q[0];
sx q[0];
rz(-1.7986725) q[0];
sx q[0];
rz(-0.22875889) q[0];
rz(-0.79452823) q[2];
sx q[2];
rz(-1.3811888) q[2];
sx q[2];
rz(-0.42805782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1340227) q[1];
sx q[1];
rz(-1.5948199) q[1];
sx q[1];
rz(0.71716208) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1845334) q[3];
sx q[3];
rz(-0.99122916) q[3];
sx q[3];
rz(1.2141276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18572346) q[2];
sx q[2];
rz(-0.74803868) q[2];
sx q[2];
rz(1.424074) q[2];
rz(-0.87797034) q[3];
sx q[3];
rz(-0.22017559) q[3];
sx q[3];
rz(0.018208114) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672742) q[0];
sx q[0];
rz(-1.4602639) q[0];
sx q[0];
rz(-1.169745) q[0];
rz(0.36880233) q[1];
sx q[1];
rz(-1.7316876) q[1];
sx q[1];
rz(0.015451886) q[1];
rz(0.56024341) q[2];
sx q[2];
rz(-0.50749736) q[2];
sx q[2];
rz(0.56847405) q[2];
rz(2.5997509) q[3];
sx q[3];
rz(-0.65237696) q[3];
sx q[3];
rz(-0.44192016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
