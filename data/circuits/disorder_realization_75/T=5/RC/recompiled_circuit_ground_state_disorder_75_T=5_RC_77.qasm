OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(0.2398332) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(-1.6286558) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1334907) q[0];
sx q[0];
rz(-1.9372131) q[0];
sx q[0];
rz(-0.13704637) q[0];
rz(-pi) q[1];
rz(2.0432908) q[2];
sx q[2];
rz(-1.8046981) q[2];
sx q[2];
rz(2.5462674) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4862772) q[1];
sx q[1];
rz(-1.4239171) q[1];
sx q[1];
rz(-2.6492061) q[1];
rz(-pi) q[2];
rz(-1.9001107) q[3];
sx q[3];
rz(-0.17767492) q[3];
sx q[3];
rz(0.8575646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98770398) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(-1.055701) q[2];
rz(-2.740247) q[3];
sx q[3];
rz(-1.515712) q[3];
sx q[3];
rz(1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6317247) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(-0.74613219) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(-0.064780898) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5279605) q[0];
sx q[0];
rz(-2.0191231) q[0];
sx q[0];
rz(-1.1290068) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6798052) q[2];
sx q[2];
rz(-1.7035988) q[2];
sx q[2];
rz(-2.6390569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8967817) q[1];
sx q[1];
rz(-2.9771388) q[1];
sx q[1];
rz(2.0816148) q[1];
x q[2];
rz(1.719789) q[3];
sx q[3];
rz(-2.9284366) q[3];
sx q[3];
rz(-0.085863559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-0.041291324) q[2];
rz(1.1193554) q[3];
sx q[3];
rz(-1.7740403) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(2.4422755) q[0];
rz(-0.62272561) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(2.8831388) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39614933) q[0];
sx q[0];
rz(-0.93702261) q[0];
sx q[0];
rz(-1.5556015) q[0];
x q[1];
rz(-2.223143) q[2];
sx q[2];
rz(-1.47965) q[2];
sx q[2];
rz(0.80207588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5194578) q[1];
sx q[1];
rz(-1.5506652) q[1];
sx q[1];
rz(-2.6926671) q[1];
rz(-pi) q[2];
rz(-1.650189) q[3];
sx q[3];
rz(-0.43118048) q[3];
sx q[3];
rz(3.0357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2059325) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(-0.70651954) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-0.8258515) q[3];
sx q[3];
rz(1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.127447) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(0.97417796) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(1.0008224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0865964) q[0];
sx q[0];
rz(-2.1657513) q[0];
sx q[0];
rz(2.9528098) q[0];
rz(-pi) q[1];
rz(-2.3712709) q[2];
sx q[2];
rz(-2.1674967) q[2];
sx q[2];
rz(-1.6092827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82881935) q[1];
sx q[1];
rz(-0.76550325) q[1];
sx q[1];
rz(2.440052) q[1];
rz(-1.521268) q[3];
sx q[3];
rz(-1.8911181) q[3];
sx q[3];
rz(1.76521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68690825) q[2];
sx q[2];
rz(-1.800622) q[2];
sx q[2];
rz(-0.81721133) q[2];
rz(-2.5456083) q[3];
sx q[3];
rz(-1.417336) q[3];
sx q[3];
rz(1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78438321) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(-2.0696409) q[0];
rz(1.1373854) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(0.17328182) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9715292) q[0];
sx q[0];
rz(-1.9636298) q[0];
sx q[0];
rz(2.6512572) q[0];
x q[1];
rz(-2.0614122) q[2];
sx q[2];
rz(-1.2254493) q[2];
sx q[2];
rz(2.5290348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90975033) q[1];
sx q[1];
rz(-0.76059231) q[1];
sx q[1];
rz(1.4666345) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9786879) q[3];
sx q[3];
rz(-1.4789733) q[3];
sx q[3];
rz(0.11549982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5002354) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(-0.87453169) q[2];
rz(-2.4747961) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(-2.3165406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85103971) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(-0.17459757) q[0];
rz(-0.54706508) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(-1.863716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9491315) q[0];
sx q[0];
rz(-1.8572079) q[0];
sx q[0];
rz(-3.0770296) q[0];
rz(-0.49680423) q[2];
sx q[2];
rz(-0.81648177) q[2];
sx q[2];
rz(-1.4085438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6656832) q[1];
sx q[1];
rz(-1.922125) q[1];
sx q[1];
rz(-2.372588) q[1];
rz(-pi) q[2];
rz(-0.90323351) q[3];
sx q[3];
rz(-2.0732911) q[3];
sx q[3];
rz(2.5046405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3829019) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(2.7721789) q[2];
rz(-2.3139125) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(-0.15414342) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5740042) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(-2.9045203) q[0];
rz(-1.7489307) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(1.0038092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8477551) q[0];
sx q[0];
rz(-1.7373573) q[0];
sx q[0];
rz(0.44035797) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4280274) q[2];
sx q[2];
rz(-1.2110365) q[2];
sx q[2];
rz(2.6204315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6663346) q[1];
sx q[1];
rz(-1.043519) q[1];
sx q[1];
rz(-1.9557245) q[1];
rz(-pi) q[2];
rz(-0.5587033) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(2.1897763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5281333) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(2.6252739) q[2];
rz(-2.3033219) q[3];
sx q[3];
rz(-1.4811938) q[3];
sx q[3];
rz(-2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-0.28165278) q[0];
sx q[0];
rz(-2.2273492) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.4062358) q[1];
sx q[1];
rz(2.2343238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4582918) q[0];
sx q[0];
rz(-2.2193546) q[0];
sx q[0];
rz(-3.1187727) q[0];
x q[1];
rz(-1.7444939) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(-0.71719682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4248391) q[1];
sx q[1];
rz(-2.3100634) q[1];
sx q[1];
rz(2.5118973) q[1];
rz(-pi) q[2];
rz(0.26772883) q[3];
sx q[3];
rz(-2.1526823) q[3];
sx q[3];
rz(-1.8393593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3153136) q[2];
sx q[2];
rz(-1.4297012) q[2];
sx q[2];
rz(2.2798174) q[2];
rz(-1.2365384) q[3];
sx q[3];
rz(-0.084241353) q[3];
sx q[3];
rz(-1.9305362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(2.1114517) q[0];
rz(-0.31632272) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(0.28824678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9194946) q[0];
sx q[0];
rz(-1.4256501) q[0];
sx q[0];
rz(-0.40444379) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4924303) q[2];
sx q[2];
rz(-1.9948975) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94142524) q[1];
sx q[1];
rz(-1.1516478) q[1];
sx q[1];
rz(2.2302385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28137286) q[3];
sx q[3];
rz(-2.4503539) q[3];
sx q[3];
rz(0.465525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2890702) q[2];
sx q[2];
rz(-1.728629) q[2];
sx q[2];
rz(2.9150325) q[2];
rz(1.6864927) q[3];
sx q[3];
rz(-2.2454567) q[3];
sx q[3];
rz(-0.69211012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914472) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(2.5164497) q[0];
rz(1.9236247) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(1.9295173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8344515) q[0];
sx q[0];
rz(-2.0491776) q[0];
sx q[0];
rz(-2.4605453) q[0];
x q[1];
rz(1.8292959) q[2];
sx q[2];
rz(-2.2257559) q[2];
sx q[2];
rz(1.0884681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.068262488) q[1];
sx q[1];
rz(-2.0113328) q[1];
sx q[1];
rz(0.063830062) q[1];
rz(0.30110601) q[3];
sx q[3];
rz(-1.8499377) q[3];
sx q[3];
rz(-0.10023403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43442279) q[2];
sx q[2];
rz(-1.7662798) q[2];
sx q[2];
rz(2.0909069) q[2];
rz(2.5896942) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(-1.2845854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(-2.3552409) q[1];
sx q[1];
rz(-2.4813589) q[1];
sx q[1];
rz(-2.9207041) q[1];
rz(-1.4715696) q[2];
sx q[2];
rz(-0.36505112) q[2];
sx q[2];
rz(0.65341841) q[2];
rz(-0.10765392) q[3];
sx q[3];
rz(-2.4526377) q[3];
sx q[3];
rz(-2.4124877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
