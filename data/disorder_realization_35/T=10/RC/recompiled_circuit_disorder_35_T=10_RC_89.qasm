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
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697203) q[0];
sx q[0];
rz(-2.2315931) q[0];
sx q[0];
rz(-1.1138492) q[0];
rz(-pi) q[1];
rz(-0.77180441) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(0.31847218) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4989657) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(0.39019231) q[1];
rz(2.911525) q[3];
sx q[3];
rz(-0.32031968) q[3];
sx q[3];
rz(2.7473118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26596507) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-2.1038726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99012016) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(0.56818509) q[0];
rz(-pi) q[1];
rz(-2.3677164) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-2.3351923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84320074) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(0.028224736) q[1];
rz(-2.188835) q[3];
sx q[3];
rz(-2.2817094) q[3];
sx q[3];
rz(-3.0666805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-0.57317615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37704913) q[0];
sx q[0];
rz(-1.1408313) q[0];
sx q[0];
rz(3.0717875) q[0];
rz(2.8163221) q[2];
sx q[2];
rz(-0.79954445) q[2];
sx q[2];
rz(-2.6550967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14904505) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(-0.64736127) q[1];
rz(-1.1099986) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.213221) q[0];
x q[1];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(-1.5359495) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73211654) q[1];
sx q[1];
rz(-0.64142694) q[1];
sx q[1];
rz(0.17318053) q[1];
rz(3.0691931) q[3];
sx q[3];
rz(-2.1941333) q[3];
sx q[3];
rz(0.49923957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5495816) q[0];
sx q[0];
rz(-1.3322543) q[0];
sx q[0];
rz(0.15079389) q[0];
x q[1];
rz(0.95894496) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-2.3252955) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(2.4232037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0110537) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(-1.2572446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(0.13866436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9753871) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(0.020394527) q[0];
rz(-pi) q[1];
rz(-2.9406592) q[2];
sx q[2];
rz(-1.6461419) q[2];
sx q[2];
rz(2.2546525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2648592) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(2.369146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1451374) q[3];
sx q[3];
rz(-0.87493757) q[3];
sx q[3];
rz(1.828572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1034575) q[0];
sx q[0];
rz(-0.75224829) q[0];
sx q[0];
rz(0.30970807) q[0];
x q[1];
rz(0.95958556) q[2];
sx q[2];
rz(-0.72003905) q[2];
sx q[2];
rz(0.50146539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99164167) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(2.409694) q[1];
x q[2];
rz(-0.94654406) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(3.0016622) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55420586) q[0];
sx q[0];
rz(-1.6763858) q[0];
sx q[0];
rz(-1.1734074) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89739563) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(-2.430254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7651556) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(-2.5625485) q[1];
rz(-pi) q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(-3.0144514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.3508266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5577561) q[0];
sx q[0];
rz(-1.6769874) q[0];
sx q[0];
rz(-1.7519959) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5647854) q[2];
sx q[2];
rz(-1.9153321) q[2];
sx q[2];
rz(-2.92958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0235325) q[1];
sx q[1];
rz(-0.91731056) q[1];
sx q[1];
rz(-2.1937624) q[1];
rz(1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(-2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(-1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4237758) q[0];
sx q[0];
rz(-2.2795838) q[0];
sx q[0];
rz(2.8481759) q[0];
rz(1.7231862) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(-1.5002804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17391275) q[1];
sx q[1];
rz(-2.850107) q[1];
sx q[1];
rz(-3.013054) q[1];
rz(-pi) q[2];
rz(-2.91594) q[3];
sx q[3];
rz(-1.435558) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-3.1267358) q[2];
sx q[2];
rz(-1.5236293) q[2];
sx q[2];
rz(-2.2894494) q[2];
rz(-1.0989582) q[3];
sx q[3];
rz(-1.8660587) q[3];
sx q[3];
rz(2.3793424) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
