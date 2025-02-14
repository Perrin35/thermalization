OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(-2.13621) q[0];
sx q[0];
rz(-0.96860743) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580509) q[0];
sx q[0];
rz(-0.59701118) q[0];
sx q[0];
rz(0.42177864) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81852976) q[2];
sx q[2];
rz(-1.6224183) q[2];
sx q[2];
rz(2.7843786) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1263784) q[1];
sx q[1];
rz(-0.83288899) q[1];
sx q[1];
rz(-0.53453858) q[1];
rz(-pi) q[2];
rz(1.6601059) q[3];
sx q[3];
rz(-1.4353283) q[3];
sx q[3];
rz(2.7394046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5030824) q[2];
sx q[2];
rz(-0.57380399) q[2];
sx q[2];
rz(-2.8199675) q[2];
rz(0.62421787) q[3];
sx q[3];
rz(-1.3892684) q[3];
sx q[3];
rz(1.6577087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44840789) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-0.93628991) q[0];
rz(-1.869465) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(-1.7535271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4988588) q[0];
sx q[0];
rz(-1.6279477) q[0];
sx q[0];
rz(-1.3512035) q[0];
rz(-pi) q[1];
rz(-0.77735591) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(-1.1071902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0213881) q[1];
sx q[1];
rz(-2.1600381) q[1];
sx q[1];
rz(0.91597775) q[1];
rz(-pi) q[2];
rz(-0.21556385) q[3];
sx q[3];
rz(-1.8290534) q[3];
sx q[3];
rz(2.7815931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0442514) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(-0.67330366) q[2];
rz(1.3168969) q[3];
sx q[3];
rz(-1.293106) q[3];
sx q[3];
rz(2.3015658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3299385) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(-0.04976186) q[0];
rz(-0.20637575) q[1];
sx q[1];
rz(-1.754909) q[1];
sx q[1];
rz(-1.6516364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7508036) q[0];
sx q[0];
rz(-0.62701462) q[0];
sx q[0];
rz(-0.20558447) q[0];
rz(-pi) q[1];
rz(2.4902053) q[2];
sx q[2];
rz(-0.66115253) q[2];
sx q[2];
rz(2.9686454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4088285) q[1];
sx q[1];
rz(-1.8547426) q[1];
sx q[1];
rz(2.8676621) q[1];
x q[2];
rz(-0.95155119) q[3];
sx q[3];
rz(-0.54446044) q[3];
sx q[3];
rz(-1.3109372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9516307) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(-3.087431) q[2];
rz(-1.9870029) q[3];
sx q[3];
rz(-1.2879939) q[3];
sx q[3];
rz(1.1156999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7390249) q[0];
sx q[0];
rz(-1.716528) q[0];
sx q[0];
rz(0.85790747) q[0];
rz(-2.4671381) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(-1.8106921) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.389849) q[0];
sx q[0];
rz(-1.3350335) q[0];
sx q[0];
rz(1.430473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0687549) q[2];
sx q[2];
rz(-0.34345657) q[2];
sx q[2];
rz(2.1844027) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34506306) q[1];
sx q[1];
rz(-1.3994675) q[1];
sx q[1];
rz(2.1354799) q[1];
rz(-pi) q[2];
rz(-1.5783806) q[3];
sx q[3];
rz(-0.41927734) q[3];
sx q[3];
rz(1.5602001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3746419) q[2];
sx q[2];
rz(-2.9810814) q[2];
sx q[2];
rz(3.0373419) q[2];
rz(2.5637964) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(0.92208636) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.47559) q[0];
sx q[0];
rz(0.9444899) q[0];
rz(-1.3635483) q[1];
sx q[1];
rz(-1.8133769) q[1];
sx q[1];
rz(-1.69467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3730225) q[0];
sx q[0];
rz(-1.3813586) q[0];
sx q[0];
rz(-0.85590881) q[0];
rz(-pi) q[1];
rz(3.1051272) q[2];
sx q[2];
rz(-1.3205832) q[2];
sx q[2];
rz(1.4831051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31476281) q[1];
sx q[1];
rz(-0.75195014) q[1];
sx q[1];
rz(-0.58677267) q[1];
x q[2];
rz(0.78131494) q[3];
sx q[3];
rz(-2.4014946) q[3];
sx q[3];
rz(0.61352713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.038534433) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(-0.72059694) q[2];
rz(0.12901148) q[3];
sx q[3];
rz(-1.3611662) q[3];
sx q[3];
rz(1.8425997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1593889) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(0.61656117) q[0];
rz(2.189134) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(1.8985101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12981249) q[0];
sx q[0];
rz(-1.8300548) q[0];
sx q[0];
rz(0.14136049) q[0];
rz(-2.7305349) q[2];
sx q[2];
rz(-0.90663821) q[2];
sx q[2];
rz(-2.1020087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4259498) q[1];
sx q[1];
rz(-0.851812) q[1];
sx q[1];
rz(-1.651161) q[1];
rz(-pi) q[2];
rz(-0.9493306) q[3];
sx q[3];
rz(-1.5104745) q[3];
sx q[3];
rz(-0.5746791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2359961) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(-0.074404152) q[2];
rz(2.060804) q[3];
sx q[3];
rz(-1.4932884) q[3];
sx q[3];
rz(1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3022795) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(3.0795414) q[0];
rz(-1.7265559) q[1];
sx q[1];
rz(-0.77295417) q[1];
sx q[1];
rz(-0.089990377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190529) q[0];
sx q[0];
rz(-0.74680674) q[0];
sx q[0];
rz(1.9487621) q[0];
rz(-pi) q[1];
rz(0.16570602) q[2];
sx q[2];
rz(-2.2519654) q[2];
sx q[2];
rz(2.990304) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5127232) q[1];
sx q[1];
rz(-0.97755331) q[1];
sx q[1];
rz(0.42494659) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59023436) q[3];
sx q[3];
rz(-2.1274421) q[3];
sx q[3];
rz(-2.7624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1354847) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(2.4583859) q[2];
rz(1.3225383) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(1.7542084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6739864) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(0.4678539) q[0];
rz(-2.4521008) q[1];
sx q[1];
rz(-0.71527022) q[1];
sx q[1];
rz(2.668344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9625585) q[0];
sx q[0];
rz(-2.9518572) q[0];
sx q[0];
rz(-1.1880072) q[0];
rz(-pi) q[1];
rz(2.6976349) q[2];
sx q[2];
rz(-1.7164162) q[2];
sx q[2];
rz(0.99010003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0508214) q[1];
sx q[1];
rz(-1.5552525) q[1];
sx q[1];
rz(1.5883716) q[1];
rz(2.8611818) q[3];
sx q[3];
rz(-1.8160928) q[3];
sx q[3];
rz(0.51285686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1271992) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(0.42352208) q[2];
rz(1.8262156) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(1.8486842) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92726436) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(1.8233874) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(3.1277025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0158169) q[0];
sx q[0];
rz(-0.36453521) q[0];
sx q[0];
rz(-1.4992072) q[0];
x q[1];
rz(-1.4133164) q[2];
sx q[2];
rz(-1.9183927) q[2];
sx q[2];
rz(1.9633788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.909643) q[1];
sx q[1];
rz(-1.0000537) q[1];
sx q[1];
rz(-0.29517031) q[1];
x q[2];
rz(0.84133103) q[3];
sx q[3];
rz(-1.6631261) q[3];
sx q[3];
rz(-1.1323557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58005971) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(2.0884936) q[2];
rz(-2.2140391) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(-0.48228669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574361) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(-0.60047737) q[0];
rz(-0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(1.762766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81299128) q[0];
sx q[0];
rz(-2.6909851) q[0];
sx q[0];
rz(-0.91186055) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37048933) q[2];
sx q[2];
rz(-1.8571808) q[2];
sx q[2];
rz(-3.1305673) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8104756) q[1];
sx q[1];
rz(-1.9974736) q[1];
sx q[1];
rz(1.3699378) q[1];
rz(0.61609488) q[3];
sx q[3];
rz(-1.2047676) q[3];
sx q[3];
rz(-0.6835365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.025734162) q[2];
sx q[2];
rz(-2.8074042) q[2];
sx q[2];
rz(2.1791229) q[2];
rz(-0.95364237) q[3];
sx q[3];
rz(-2.0680659) q[3];
sx q[3];
rz(1.3941221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3772603) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(1.7878905) q[1];
sx q[1];
rz(-1.1693015) q[1];
sx q[1];
rz(2.404626) q[1];
rz(-0.76112657) q[2];
sx q[2];
rz(-2.3115951) q[2];
sx q[2];
rz(-2.8106845) q[2];
rz(-2.8938724) q[3];
sx q[3];
rz(-2.6213624) q[3];
sx q[3];
rz(-2.6772671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
