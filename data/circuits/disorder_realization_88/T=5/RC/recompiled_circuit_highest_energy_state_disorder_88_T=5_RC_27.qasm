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
rz(0.54230827) q[0];
sx q[0];
rz(-0.13442726) q[0];
sx q[0];
rz(-1.0472714) q[0];
rz(2.8174227) q[1];
sx q[1];
rz(-2.9919762) q[1];
sx q[1];
rz(-0.9447929) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36318521) q[0];
sx q[0];
rz(-1.249832) q[0];
sx q[0];
rz(-0.28807333) q[0];
rz(-0.23350291) q[2];
sx q[2];
rz(-2.2036607) q[2];
sx q[2];
rz(1.7720745) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8103247) q[1];
sx q[1];
rz(-0.70155084) q[1];
sx q[1];
rz(2.3980106) q[1];
x q[2];
rz(-1.7994491) q[3];
sx q[3];
rz(-2.6480669) q[3];
sx q[3];
rz(-2.760104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85113877) q[2];
sx q[2];
rz(-1.4522499) q[2];
sx q[2];
rz(-1.5770844) q[2];
rz(-2.5751298) q[3];
sx q[3];
rz(-0.62140673) q[3];
sx q[3];
rz(1.2982781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-0.7199921) q[0];
rz(-0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(-2.5879587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33687978) q[0];
sx q[0];
rz(-1.2186482) q[0];
sx q[0];
rz(-3.0180305) q[0];
rz(-2.4963412) q[2];
sx q[2];
rz(-0.49075365) q[2];
sx q[2];
rz(-3.0440999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1914042) q[1];
sx q[1];
rz(-2.2850415) q[1];
sx q[1];
rz(0.49226239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3838061) q[3];
sx q[3];
rz(-0.67000853) q[3];
sx q[3];
rz(-0.78950933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9686034) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(-2.7776862) q[2];
rz(0.11161741) q[3];
sx q[3];
rz(-2.2026187) q[3];
sx q[3];
rz(-0.76044559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97742057) q[0];
sx q[0];
rz(-0.82472473) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(2.3258356) q[1];
sx q[1];
rz(-2.0033671) q[1];
sx q[1];
rz(-2.7556748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49026981) q[0];
sx q[0];
rz(-1.2322591) q[0];
sx q[0];
rz(1.2744034) q[0];
rz(-pi) q[1];
rz(-2.1325941) q[2];
sx q[2];
rz(-1.2616065) q[2];
sx q[2];
rz(-2.5872292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2066127) q[1];
sx q[1];
rz(-2.7394208) q[1];
sx q[1];
rz(1.241317) q[1];
rz(-3.0073062) q[3];
sx q[3];
rz(-1.9818174) q[3];
sx q[3];
rz(0.81601303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58891121) q[2];
sx q[2];
rz(-2.3083355) q[2];
sx q[2];
rz(1.7595278) q[2];
rz(2.9803993) q[3];
sx q[3];
rz(-1.4313982) q[3];
sx q[3];
rz(2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2564119) q[0];
sx q[0];
rz(-1.8653402) q[0];
sx q[0];
rz(-1.8026344) q[0];
rz(1.8244052) q[1];
sx q[1];
rz(-0.97511292) q[1];
sx q[1];
rz(-2.8417974) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9343963) q[0];
sx q[0];
rz(-2.3863433) q[0];
sx q[0];
rz(2.62225) q[0];
x q[1];
rz(0.91641578) q[2];
sx q[2];
rz(-2.4407226) q[2];
sx q[2];
rz(0.13417164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7924031) q[1];
sx q[1];
rz(-1.9219256) q[1];
sx q[1];
rz(2.9955088) q[1];
rz(-pi) q[2];
rz(-0.13238975) q[3];
sx q[3];
rz(-0.93919884) q[3];
sx q[3];
rz(-1.6211444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51912159) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(-0.31912121) q[2];
rz(-0.090911344) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(1.3566141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(3.1392642) q[0];
rz(1.5706515) q[1];
sx q[1];
rz(-0.80051533) q[1];
sx q[1];
rz(1.194582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55694675) q[0];
sx q[0];
rz(-1.7255387) q[0];
sx q[0];
rz(-1.690288) q[0];
x q[1];
rz(0.04045893) q[2];
sx q[2];
rz(-0.85340696) q[2];
sx q[2];
rz(-0.77740146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0280605) q[1];
sx q[1];
rz(-2.781377) q[1];
sx q[1];
rz(-1.8647782) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8525328) q[3];
sx q[3];
rz(-2.1093371) q[3];
sx q[3];
rz(2.0592214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.883256) q[2];
sx q[2];
rz(-2.363435) q[2];
sx q[2];
rz(-0.077433132) q[2];
rz(-0.94610131) q[3];
sx q[3];
rz(-2.6587722) q[3];
sx q[3];
rz(-0.4536804) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9282064) q[0];
sx q[0];
rz(-0.086406924) q[0];
sx q[0];
rz(1.7267831) q[0];
rz(-0.66697031) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(0.40474969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20807438) q[0];
sx q[0];
rz(-2.7221823) q[0];
sx q[0];
rz(-2.389877) q[0];
rz(-pi) q[1];
rz(1.9046632) q[2];
sx q[2];
rz(-0.68511663) q[2];
sx q[2];
rz(1.8917055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25127801) q[1];
sx q[1];
rz(-2.0330486) q[1];
sx q[1];
rz(0.068271116) q[1];
rz(-pi) q[2];
rz(-2.7795198) q[3];
sx q[3];
rz(-1.1112461) q[3];
sx q[3];
rz(-1.5804231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82379597) q[2];
sx q[2];
rz(-0.78936374) q[2];
sx q[2];
rz(2.6981567) q[2];
rz(0.41380841) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(-0.88576353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7498748) q[0];
sx q[0];
rz(-1.7789919) q[0];
sx q[0];
rz(-0.66746563) q[0];
rz(-0.04774566) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(-1.1183636) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55505468) q[0];
sx q[0];
rz(-2.1004875) q[0];
sx q[0];
rz(1.0989957) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3862382) q[2];
sx q[2];
rz(-2.5552351) q[2];
sx q[2];
rz(-1.9826012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.075526) q[1];
sx q[1];
rz(-2.3956163) q[1];
sx q[1];
rz(2.4999208) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6947855) q[3];
sx q[3];
rz(-2.638804) q[3];
sx q[3];
rz(-0.10956746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.045791693) q[2];
sx q[2];
rz(-1.2597151) q[2];
sx q[2];
rz(-3.1328787) q[2];
rz(-1.2039315) q[3];
sx q[3];
rz(-2.7976076) q[3];
sx q[3];
rz(-3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0225723) q[0];
sx q[0];
rz(-0.38563269) q[0];
sx q[0];
rz(0.88985306) q[0];
rz(1.8564557) q[1];
sx q[1];
rz(-2.0882008) q[1];
sx q[1];
rz(-3.0056675) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724919) q[0];
sx q[0];
rz(-1.6181503) q[0];
sx q[0];
rz(-0.22914178) q[0];
rz(-pi) q[1];
rz(1.0338342) q[2];
sx q[2];
rz(-1.0091678) q[2];
sx q[2];
rz(-1.6800576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1358661) q[1];
sx q[1];
rz(-0.39038613) q[1];
sx q[1];
rz(0.54061215) q[1];
x q[2];
rz(1.9210757) q[3];
sx q[3];
rz(-1.4363343) q[3];
sx q[3];
rz(2.2148481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70302427) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(2.0795889) q[2];
rz(0.92987531) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(-2.3532531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36282614) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(2.7407001) q[0];
rz(-0.76155424) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(0.78899312) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61643314) q[0];
sx q[0];
rz(-0.1285006) q[0];
sx q[0];
rz(3.0165821) q[0];
x q[1];
rz(0.89575276) q[2];
sx q[2];
rz(-1.4491242) q[2];
sx q[2];
rz(-2.3020803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4915062) q[1];
sx q[1];
rz(-1.4515001) q[1];
sx q[1];
rz(-2.5507798) q[1];
rz(0.52481905) q[3];
sx q[3];
rz(-1.6136955) q[3];
sx q[3];
rz(0.08924291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5679428) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(0.72081494) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-0.78250116) q[3];
sx q[3];
rz(1.8318374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023130527) q[0];
sx q[0];
rz(-3.0266302) q[0];
sx q[0];
rz(-0.7777099) q[0];
rz(0.65981162) q[1];
sx q[1];
rz(-2.2751364) q[1];
sx q[1];
rz(0.39001098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3005053) q[0];
sx q[0];
rz(-1.705637) q[0];
sx q[0];
rz(-3.10411) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74441461) q[2];
sx q[2];
rz(-2.0995903) q[2];
sx q[2];
rz(0.86133682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9007787) q[1];
sx q[1];
rz(-1.1692746) q[1];
sx q[1];
rz(0.32380775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.635119) q[3];
sx q[3];
rz(-2.1586543) q[3];
sx q[3];
rz(-1.3660869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68221349) q[2];
sx q[2];
rz(-2.7028658) q[2];
sx q[2];
rz(0.48975804) q[2];
rz(-2.8343936) q[3];
sx q[3];
rz(-2.2178853) q[3];
sx q[3];
rz(-1.3908305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601892) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(1.0745984) q[1];
sx q[1];
rz(-2.2886724) q[1];
sx q[1];
rz(1.2235175) q[1];
rz(0.39911453) q[2];
sx q[2];
rz(-1.21798) q[2];
sx q[2];
rz(-1.1026364) q[2];
rz(0.60081595) q[3];
sx q[3];
rz(-1.0944585) q[3];
sx q[3];
rz(-0.16437358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
