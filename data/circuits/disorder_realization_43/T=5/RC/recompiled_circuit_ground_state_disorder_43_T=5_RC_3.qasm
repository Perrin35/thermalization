OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8671829) q[0];
sx q[0];
rz(-2.6743968) q[0];
sx q[0];
rz(0.89611563) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28361904) q[0];
sx q[0];
rz(-1.9500226) q[0];
sx q[0];
rz(2.2573756) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3199086) q[2];
sx q[2];
rz(-2.0178621) q[2];
sx q[2];
rz(2.7607703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3442698) q[1];
sx q[1];
rz(-1.8472278) q[1];
sx q[1];
rz(1.7452507) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7475376) q[3];
sx q[3];
rz(-2.7480304) q[3];
sx q[3];
rz(1.8436028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2200615) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(-0.49911487) q[2];
rz(2.3257997) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(-0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4895184) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(-2.9079085) q[0];
rz(0.09985996) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(1.608009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312992) q[0];
sx q[0];
rz(-1.4742032) q[0];
sx q[0];
rz(-1.4486758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7031021) q[2];
sx q[2];
rz(-1.4027565) q[2];
sx q[2];
rz(0.7249671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3996904) q[1];
sx q[1];
rz(-0.61237017) q[1];
sx q[1];
rz(1.4111817) q[1];
rz(-2.1183999) q[3];
sx q[3];
rz(-0.55063081) q[3];
sx q[3];
rz(-0.44704244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0443772) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(0.054917939) q[2];
rz(-2.4954259) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(-3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9570479) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(2.3967337) q[0];
rz(-1.2767731) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(0.863711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.389457) q[0];
sx q[0];
rz(-2.6977718) q[0];
sx q[0];
rz(-3.0054166) q[0];
rz(-pi) q[1];
rz(0.50541877) q[2];
sx q[2];
rz(-0.8552455) q[2];
sx q[2];
rz(-2.9723013) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7924713) q[1];
sx q[1];
rz(-1.0949666) q[1];
sx q[1];
rz(0.70648593) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4808472) q[3];
sx q[3];
rz(-1.163394) q[3];
sx q[3];
rz(-1.416172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.471571) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(0.395533) q[2];
rz(1.0223355) q[3];
sx q[3];
rz(-0.46415713) q[3];
sx q[3];
rz(0.49183229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96524298) q[0];
sx q[0];
rz(-1.9597766) q[0];
sx q[0];
rz(0.055971948) q[0];
rz(-2.5296027) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(0.8459808) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2200945) q[0];
sx q[0];
rz(-1.8290231) q[0];
sx q[0];
rz(-2.1058215) q[0];
rz(2.0417287) q[2];
sx q[2];
rz(-2.467017) q[2];
sx q[2];
rz(1.2686307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.779102) q[1];
sx q[1];
rz(-1.7622708) q[1];
sx q[1];
rz(0.31656852) q[1];
rz(-pi) q[2];
rz(-0.9153644) q[3];
sx q[3];
rz(-0.5753606) q[3];
sx q[3];
rz(-1.8917454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0116288) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(2.4212627) q[2];
rz(-2.7885041) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(-0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406152) q[0];
sx q[0];
rz(-1.4530797) q[0];
sx q[0];
rz(2.154696) q[0];
rz(1.997021) q[1];
sx q[1];
rz(-2.3026376) q[1];
sx q[1];
rz(-2.1867337) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3026149) q[0];
sx q[0];
rz(-0.90388008) q[0];
sx q[0];
rz(0.033292183) q[0];
x q[1];
rz(-2.0667162) q[2];
sx q[2];
rz(-0.62037797) q[2];
sx q[2];
rz(1.7091144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8870626) q[1];
sx q[1];
rz(-1.945244) q[1];
sx q[1];
rz(-2.6545054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0480065) q[3];
sx q[3];
rz(-0.29055933) q[3];
sx q[3];
rz(2.235242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6285051) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(1.0465735) q[2];
rz(2.4322677) q[3];
sx q[3];
rz(-1.1449292) q[3];
sx q[3];
rz(2.5078702) q[3];
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
rz(2.3777622) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(2.8705226) q[0];
rz(0.079004869) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(-1.7105506) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56667152) q[0];
sx q[0];
rz(-1.3614166) q[0];
sx q[0];
rz(1.4230595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3777804) q[2];
sx q[2];
rz(-0.71893636) q[2];
sx q[2];
rz(2.1670408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8003098) q[1];
sx q[1];
rz(-0.73004111) q[1];
sx q[1];
rz(-0.37363877) q[1];
rz(-pi) q[2];
rz(-0.16057272) q[3];
sx q[3];
rz(-1.9844921) q[3];
sx q[3];
rz(0.90908066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0926823) q[2];
sx q[2];
rz(-1.1834669) q[2];
sx q[2];
rz(3.1177706) q[2];
rz(1.228099) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(2.9983799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1517076) q[0];
sx q[0];
rz(-2.8360974) q[0];
sx q[0];
rz(-3.049343) q[0];
rz(-0.30091885) q[1];
sx q[1];
rz(-1.2291127) q[1];
sx q[1];
rz(-2.0613964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75715059) q[0];
sx q[0];
rz(-0.35878962) q[0];
sx q[0];
rz(-1.8435988) q[0];
rz(-1.148573) q[2];
sx q[2];
rz(-1.2381089) q[2];
sx q[2];
rz(1.1115369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34970666) q[1];
sx q[1];
rz(-2.9936446) q[1];
sx q[1];
rz(0.011758864) q[1];
x q[2];
rz(0.37410847) q[3];
sx q[3];
rz(-2.4112933) q[3];
sx q[3];
rz(-0.41159831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0098308) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(0.62823137) q[2];
rz(0.96588165) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(2.8204744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1107165) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(1.4068756) q[0];
rz(0.12786099) q[1];
sx q[1];
rz(-1.4297994) q[1];
sx q[1];
rz(2.0157287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4425621) q[0];
sx q[0];
rz(-1.3872212) q[0];
sx q[0];
rz(0.72893127) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6091414) q[2];
sx q[2];
rz(-2.6895617) q[2];
sx q[2];
rz(-1.1029117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2919104) q[1];
sx q[1];
rz(-2.9314802) q[1];
sx q[1];
rz(-0.30847524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0142426) q[3];
sx q[3];
rz(-2.1428041) q[3];
sx q[3];
rz(-2.1406904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13059482) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(-3.0478743) q[2];
rz(-1.7081918) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(-1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.691064) q[0];
sx q[0];
rz(-1.0858303) q[0];
sx q[0];
rz(-1.0040671) q[0];
rz(1.1035236) q[1];
sx q[1];
rz(-0.48201489) q[1];
sx q[1];
rz(1.7335266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.133014) q[0];
sx q[0];
rz(-2.6264418) q[0];
sx q[0];
rz(3.1293082) q[0];
rz(1.4595217) q[2];
sx q[2];
rz(-1.2603052) q[2];
sx q[2];
rz(1.3729915) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1939331) q[1];
sx q[1];
rz(-1.1378189) q[1];
sx q[1];
rz(0.19755944) q[1];
rz(-pi) q[2];
rz(2.6226165) q[3];
sx q[3];
rz(-1.5277315) q[3];
sx q[3];
rz(1.2534634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39190009) q[2];
sx q[2];
rz(-1.6204648) q[2];
sx q[2];
rz(2.5938972) q[2];
rz(-2.8294166) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(-1.4148022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373229) q[0];
sx q[0];
rz(-1.4708568) q[0];
sx q[0];
rz(0.7789337) q[0];
rz(2.7670822) q[1];
sx q[1];
rz(-1.6284527) q[1];
sx q[1];
rz(-2.3096854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4171281) q[0];
sx q[0];
rz(-2.3848371) q[0];
sx q[0];
rz(-2.7555076) q[0];
rz(-pi) q[1];
rz(-2.358663) q[2];
sx q[2];
rz(-2.1179869) q[2];
sx q[2];
rz(-1.7683355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.191997) q[1];
sx q[1];
rz(-1.0284541) q[1];
sx q[1];
rz(-2.3344912) q[1];
x q[2];
rz(0.79672265) q[3];
sx q[3];
rz(-1.5925538) q[3];
sx q[3];
rz(0.12783229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-0.35934862) q[2];
sx q[2];
rz(-0.0084776004) q[2];
rz(-0.40555412) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(1.0943195) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(-3.0158184) q[2];
sx q[2];
rz(-1.3604506) q[2];
sx q[2];
rz(2.4764555) q[2];
rz(-0.052362818) q[3];
sx q[3];
rz(-3.047245) q[3];
sx q[3];
rz(0.17518763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
