OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89361184) q[0];
sx q[0];
rz(-2.3900744) q[0];
sx q[0];
rz(3.0889838) q[0];
x q[1];
rz(0.46484868) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(2.6236617) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6005046) q[1];
sx q[1];
rz(-1.0827912) q[1];
sx q[1];
rz(1.0469251) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(0.38744774) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(-3.1121029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55949975) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(-2.6452243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95784159) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(0.5387696) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4757231) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.385752) q[0];
sx q[0];
rz(-0.34163654) q[0];
sx q[0];
rz(1.2582448) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6443278) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(-0.49567859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8071825) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(0.73015726) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1242261) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-0.22180804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(0.25340432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2052106) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(3.0217231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5379982) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6654012) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(-2.979216) q[1];
rz(-pi) q[2];
rz(-2.4921791) q[3];
sx q[3];
rz(-1.489349) q[3];
sx q[3];
rz(-0.34929517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-0.45809349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375181) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(-3.1361561) q[0];
x q[1];
rz(0.96197084) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(-0.48987197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8372247) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.7423082) q[1];
rz(-pi) q[2];
rz(3.0466988) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(-2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-2.2221185) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4290659) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(-2.6119786) q[0];
rz(2.9847449) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(1.4175121) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3082723) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(1.0233378) q[1];
rz(0.24354981) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(1.5598284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(3.1183929) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(2.0297208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8060914) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(-2.0557687) q[0];
rz(-pi) q[1];
rz(0.53084897) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(-1.2209148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.062473) q[1];
sx q[1];
rz(-1.2812496) q[1];
sx q[1];
rz(-1.4909966) q[1];
x q[2];
rz(1.0681549) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(2.2964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-2.3866167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(-2.9478361) q[0];
rz(2.1515498) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(-2.5255447) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8440486) q[1];
sx q[1];
rz(-2.5522759) q[1];
sx q[1];
rz(0.20357666) q[1];
rz(-pi) q[2];
rz(-0.8927535) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(2.1219818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(3.1220904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7050539) q[0];
sx q[0];
rz(-0.22163135) q[0];
sx q[0];
rz(1.9428695) q[0];
x q[1];
rz(0.15140622) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(-0.99934794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86917415) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(1.8723349) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9762514) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(-0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0120221) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(2.9905031) q[0];
x q[1];
rz(1.7933153) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-0.36048181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9705829) q[1];
sx q[1];
rz(-0.90921558) q[1];
sx q[1];
rz(2.9243102) q[1];
rz(-pi) q[2];
rz(-1.6230574) q[3];
sx q[3];
rz(-1.2219547) q[3];
sx q[3];
rz(2.7319752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(3.1046042) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
