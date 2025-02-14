OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.149481) q[0];
sx q[0];
rz(-0.7465201) q[0];
sx q[0];
rz(-0.56678766) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(2.221938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2956379) q[0];
sx q[0];
rz(-1.7311117) q[0];
sx q[0];
rz(2.3491067) q[0];
x q[1];
rz(-0.55962015) q[2];
sx q[2];
rz(-0.80794789) q[2];
sx q[2];
rz(-1.1839649) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7417042) q[1];
sx q[1];
rz(-2.2117858) q[1];
sx q[1];
rz(1.6287032) q[1];
rz(-pi) q[2];
rz(2.9018974) q[3];
sx q[3];
rz(-2.5928232) q[3];
sx q[3];
rz(-1.5769928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2016242) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(1.0533818) q[2];
rz(-0.71422226) q[3];
sx q[3];
rz(-2.768399) q[3];
sx q[3];
rz(1.2936973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193704) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(-2.569662) q[0];
rz(-1.2988623) q[1];
sx q[1];
rz(-2.3426901) q[1];
sx q[1];
rz(-2.3518708) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30366722) q[0];
sx q[0];
rz(-1.2695489) q[0];
sx q[0];
rz(3.0409854) q[0];
rz(-pi) q[1];
rz(-1.4105807) q[2];
sx q[2];
rz(-0.92276697) q[2];
sx q[2];
rz(-0.90479484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.870656) q[1];
sx q[1];
rz(-2.4693477) q[1];
sx q[1];
rz(0.45935528) q[1];
x q[2];
rz(2.0630453) q[3];
sx q[3];
rz(-2.1849602) q[3];
sx q[3];
rz(1.2224765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85084891) q[2];
sx q[2];
rz(-2.3727543) q[2];
sx q[2];
rz(1.7791465) q[2];
rz(0.29081523) q[3];
sx q[3];
rz(-1.3664061) q[3];
sx q[3];
rz(-0.10663685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0027851) q[0];
sx q[0];
rz(-2.0464351) q[0];
sx q[0];
rz(0.7005257) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-0.22451678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4493664) q[0];
sx q[0];
rz(-2.7590519) q[0];
sx q[0];
rz(-0.40479779) q[0];
rz(-0.098578886) q[2];
sx q[2];
rz(-2.1469627) q[2];
sx q[2];
rz(-2.2599932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91229328) q[1];
sx q[1];
rz(-2.5535431) q[1];
sx q[1];
rz(1.4347784) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6979917) q[3];
sx q[3];
rz(-2.4089097) q[3];
sx q[3];
rz(2.3310083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38480467) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(0.067042025) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(2.4750347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.1753569) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(-0.019158451) q[0];
rz(-1.3592023) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(-0.0044936831) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86165262) q[0];
sx q[0];
rz(-0.77095448) q[0];
sx q[0];
rz(-2.3907135) q[0];
rz(-pi) q[1];
rz(-1.121663) q[2];
sx q[2];
rz(-2.0259078) q[2];
sx q[2];
rz(-3.1246076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0862866) q[1];
sx q[1];
rz(-1.2125361) q[1];
sx q[1];
rz(-2.9614224) q[1];
rz(0.2035844) q[3];
sx q[3];
rz(-1.4798375) q[3];
sx q[3];
rz(-1.1435777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6898474) q[2];
sx q[2];
rz(-2.0514026) q[2];
sx q[2];
rz(1.4534265) q[2];
rz(-1.8870185) q[3];
sx q[3];
rz(-1.5213608) q[3];
sx q[3];
rz(0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314064) q[0];
sx q[0];
rz(-1.9047381) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(-2.7120178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21107656) q[0];
sx q[0];
rz(-1.9787702) q[0];
sx q[0];
rz(-3.1019449) q[0];
rz(2.863071) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(0.24519224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.030368806) q[1];
sx q[1];
rz(-0.75654034) q[1];
sx q[1];
rz(2.0732422) q[1];
x q[2];
rz(1.1384835) q[3];
sx q[3];
rz(-2.4412182) q[3];
sx q[3];
rz(0.56454235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9261711) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.8298836) q[2];
rz(2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(-3.0084012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8120414) q[0];
sx q[0];
rz(-0.1596182) q[0];
sx q[0];
rz(2.0352236) q[0];
rz(3.0270992) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(-3.0879424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0993529) q[0];
sx q[0];
rz(-1.5765549) q[0];
sx q[0];
rz(1.3151709) q[0];
rz(-pi) q[1];
rz(1.6752376) q[2];
sx q[2];
rz(-0.76600961) q[2];
sx q[2];
rz(-0.05482373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1785894) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(-1.4609171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63416449) q[3];
sx q[3];
rz(-2.0863219) q[3];
sx q[3];
rz(2.9662715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(0.50714058) q[2];
rz(1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(-1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(0.041286904) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(2.1521177) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7558003) q[0];
sx q[0];
rz(-0.72052252) q[0];
sx q[0];
rz(0.81554834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35571738) q[2];
sx q[2];
rz(-1.3928486) q[2];
sx q[2];
rz(-1.2426113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7478024) q[1];
sx q[1];
rz(-0.89266864) q[1];
sx q[1];
rz(-2.3403779) q[1];
rz(1.0999849) q[3];
sx q[3];
rz(-0.71094497) q[3];
sx q[3];
rz(-2.8266852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(-0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82666731) q[0];
sx q[0];
rz(-0.59190094) q[0];
sx q[0];
rz(2.9097606) q[0];
rz(-2.5595317) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(-1.2190762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1506596) q[0];
sx q[0];
rz(-2.5703589) q[0];
sx q[0];
rz(-2.9240963) q[0];
x q[1];
rz(-2.788576) q[2];
sx q[2];
rz(-1.6200049) q[2];
sx q[2];
rz(2.4491058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9207014) q[1];
sx q[1];
rz(-0.24020837) q[1];
sx q[1];
rz(2.2574053) q[1];
rz(1.0295415) q[3];
sx q[3];
rz(-0.89139056) q[3];
sx q[3];
rz(-2.816659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9570738) q[2];
sx q[2];
rz(-1.7902713) q[2];
sx q[2];
rz(-2.3547724) q[2];
rz(1.6400853) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(-1.7155581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96104923) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(-2.2013262) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-2.2856789) q[1];
sx q[1];
rz(2.99517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7112367) q[0];
sx q[0];
rz(-2.8161263) q[0];
sx q[0];
rz(0.40367608) q[0];
rz(-0.27250473) q[2];
sx q[2];
rz(-1.8159869) q[2];
sx q[2];
rz(-1.5281072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6633353) q[1];
sx q[1];
rz(-2.2406849) q[1];
sx q[1];
rz(-1.2557593) q[1];
x q[2];
rz(-1.6723456) q[3];
sx q[3];
rz(-2.3595105) q[3];
sx q[3];
rz(0.33354943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0545097) q[2];
sx q[2];
rz(-1.6528218) q[2];
sx q[2];
rz(-0.83756891) q[2];
rz(-2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.8785716) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.7735057) q[0];
rz(0.076016501) q[1];
sx q[1];
rz(-2.5612505) q[1];
sx q[1];
rz(2.0215633) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382384) q[0];
sx q[0];
rz(-1.3068259) q[0];
sx q[0];
rz(-1.2498943) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4330288) q[2];
sx q[2];
rz(-1.3363133) q[2];
sx q[2];
rz(-1.8013148) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6516797) q[1];
sx q[1];
rz(-1.8615684) q[1];
sx q[1];
rz(-0.099822961) q[1];
rz(-pi) q[2];
rz(1.4073652) q[3];
sx q[3];
rz(-1.4686958) q[3];
sx q[3];
rz(-0.40615505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(-0.24701992) q[2];
rz(-2.5714827) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(-1.1183687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433255) q[0];
sx q[0];
rz(-1.387431) q[0];
sx q[0];
rz(-0.36547216) q[0];
rz(-3.0430766) q[1];
sx q[1];
rz(-0.64058522) q[1];
sx q[1];
rz(-2.2596901) q[1];
rz(1.1578887) q[2];
sx q[2];
rz(-1.558424) q[2];
sx q[2];
rz(1.3267714) q[2];
rz(0.3803654) q[3];
sx q[3];
rz(-1.813253) q[3];
sx q[3];
rz(-1.062275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
