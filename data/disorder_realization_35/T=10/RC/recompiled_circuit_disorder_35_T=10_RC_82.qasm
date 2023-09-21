OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071872358) q[0];
sx q[0];
rz(-2.2315931) q[0];
sx q[0];
rz(1.1138492) q[0];
rz(0.65806234) q[2];
sx q[2];
rz(-2.1076638) q[2];
sx q[2];
rz(1.8368349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.642627) q[1];
sx q[1];
rz(-1.1945915) q[1];
sx q[1];
rz(-0.39019231) q[1];
x q[2];
rz(2.8291679) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.8189836) q[2];
rz(-2.8406075) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.422721) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(-1.3119112) q[0];
rz(0.53595397) q[2];
sx q[2];
rz(-2.0336656) q[2];
sx q[2];
rz(3.0218389) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71643752) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(1.1643216) q[1];
rz(-0.59229895) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(-0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66118801) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37704913) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(-0.069805108) q[0];
x q[1];
rz(-2.8163221) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(0.48649597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(-1.8797727) q[1];
rz(-2.0315941) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3804647) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(2.2855177) q[0];
rz(1.5393125) q[2];
sx q[2];
rz(-1.5891979) q[2];
sx q[2];
rz(-3.1061663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73211654) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(2.9684121) q[1];
rz(-pi) q[2];
rz(0.94621559) q[3];
sx q[3];
rz(-1.5120301) q[3];
sx q[3];
rz(2.0277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.4594706) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97809726) q[0];
sx q[0];
rz(-2.8601544) q[0];
sx q[0];
rz(-1.0174169) q[0];
rz(-pi) q[1];
rz(-0.8874138) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(-2.8487157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1369049) q[1];
sx q[1];
rz(-1.7675752) q[1];
sx q[1];
rz(0.23006567) q[1];
rz(-pi) q[2];
rz(-2.8206283) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(-2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(3.0029283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370699) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(1.5741332) q[0];
rz(1.4939098) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(-0.69918699) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87673346) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(2.369146) q[1];
rz(0.74216446) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(3.1165251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-0.62430635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(1.8486345) q[0];
rz(-pi) q[1];
rz(0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(2.6401273) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(-2.2570616) q[1];
rz(0.94654406) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(-0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(3.0338045) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-3.0126742) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3710204) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(1.3034986) q[0];
rz(2.244197) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(0.71133864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32313777) q[1];
sx q[1];
rz(-1.8433246) q[1];
sx q[1];
rz(0.44169359) q[1];
rz(-1.586732) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(-2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(1.0366584) q[0];
rz(0.57680723) q[2];
sx q[2];
rz(-1.2262605) q[2];
sx q[2];
rz(-2.92958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0235325) q[1];
sx q[1];
rz(-0.91731056) q[1];
sx q[1];
rz(2.1937624) q[1];
x q[2];
rz(-0.25165598) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-2.5949123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
rz(-0.47537739) q[2];
sx q[2];
rz(-1.4350841) q[2];
sx q[2];
rz(-3.140608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17391275) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(0.12853865) q[1];
rz(-pi) q[2];
rz(2.91594) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(-2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(2.6976363) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-1.2658723) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
rz(-2.0426345) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];