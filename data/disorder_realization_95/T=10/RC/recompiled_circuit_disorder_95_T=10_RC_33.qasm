OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9165186) q[0];
sx q[0];
rz(-1.5009891) q[0];
sx q[0];
rz(-0.9401456) q[0];
rz(-2.2034982) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(0.11703581) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(0.27889241) q[1];
rz(-0.75818054) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71683305) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(0.99320937) q[0];
rz(-0.61848817) q[2];
sx q[2];
rz(-1.7881219) q[2];
sx q[2];
rz(2.6478812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2616927) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(-2.9628423) q[1];
rz(-pi) q[2];
rz(-1.3602123) q[3];
sx q[3];
rz(-0.31616351) q[3];
sx q[3];
rz(2.8663243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(2.0203967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7911581) q[0];
sx q[0];
rz(-2.3179623) q[0];
sx q[0];
rz(0.50823786) q[0];
rz(-pi) q[1];
rz(1.7240702) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(3.0561471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(2.4099039) q[1];
rz(-pi) q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(1.9672058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375835) q[0];
sx q[0];
rz(-1.596367) q[0];
sx q[0];
rz(0.63347647) q[0];
rz(-2.0747244) q[2];
sx q[2];
rz(-1.2328086) q[2];
sx q[2];
rz(-2.4525814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14028206) q[1];
sx q[1];
rz(-1.9754859) q[1];
sx q[1];
rz(0.54108041) q[1];
rz(1.6223525) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0557077) q[0];
sx q[0];
rz(-2.5143393) q[0];
sx q[0];
rz(1.6968615) q[0];
rz(-0.45102851) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(2.905068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83541218) q[1];
sx q[1];
rz(-2.7369943) q[1];
sx q[1];
rz(-2.6447891) q[1];
rz(-pi) q[2];
rz(-2.7025181) q[3];
sx q[3];
rz(-0.77507654) q[3];
sx q[3];
rz(2.4174158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(1.032069) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-0.2063624) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(-3.0138426) q[0];
rz(2.2588737) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(0.27872745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1604662) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(1.6756945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9456375) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241087) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(-0.80362513) q[0];
x q[1];
rz(1.4326101) q[2];
sx q[2];
rz(-1.0850564) q[2];
sx q[2];
rz(1.9463584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8163029) q[1];
sx q[1];
rz(-1.5029229) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(0.32254036) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.9205836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1655501) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(0.793215) q[0];
x q[1];
rz(0.17692716) q[2];
sx q[2];
rz(-1.4104341) q[2];
sx q[2];
rz(-0.71080506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.358566) q[1];
sx q[1];
rz(-1.4604124) q[1];
sx q[1];
rz(0.13659887) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2113167) q[3];
sx q[3];
rz(-2.3593138) q[3];
sx q[3];
rz(-0.6560916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(3.0158214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1544979) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(0.82018606) q[0];
x q[1];
rz(2.618082) q[2];
sx q[2];
rz(-2.7611809) q[2];
sx q[2];
rz(-1.3695804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8217433) q[1];
sx q[1];
rz(-2.2935767) q[1];
sx q[1];
rz(-2.7601932) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5500507) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(1.2951617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(-2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894293) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(-1.2055231) q[0];
x q[1];
rz(-2.8814949) q[2];
sx q[2];
rz(-2.0195228) q[2];
sx q[2];
rz(1.5425494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.6608095) q[1];
rz(-pi) q[2];
rz(-2.8285698) q[3];
sx q[3];
rz(-2.3232197) q[3];
sx q[3];
rz(-0.5648053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.3557419) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(1.981001) q[3];
sx q[3];
rz(-2.1174178) q[3];
sx q[3];
rz(-2.5934364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
