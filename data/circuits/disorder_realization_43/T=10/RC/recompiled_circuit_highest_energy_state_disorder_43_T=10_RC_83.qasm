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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(-2.9414862) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057077335) q[0];
sx q[0];
rz(-0.50327089) q[0];
sx q[0];
rz(-1.0751192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0609786) q[2];
sx q[2];
rz(-2.1128383) q[2];
sx q[2];
rz(-1.8212183) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4551022) q[1];
sx q[1];
rz(-0.89340607) q[1];
sx q[1];
rz(-1.5090921) q[1];
rz(-pi) q[2];
x q[2];
rz(1.41134) q[3];
sx q[3];
rz(-0.48122367) q[3];
sx q[3];
rz(-0.10018292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.8270095) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44875328) q[0];
sx q[0];
rz(-1.3200545) q[0];
sx q[0];
rz(3.1186799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6131975) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(0.49371546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1545002) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(0.89548703) q[1];
x q[2];
rz(-1.3608906) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(-0.24649749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(-2.7454929) q[2];
rz(-1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-2.8025467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44797541) q[0];
sx q[0];
rz(-0.98784271) q[0];
sx q[0];
rz(-1.5459803) q[0];
rz(-1.6940368) q[2];
sx q[2];
rz(-1.7584821) q[2];
sx q[2];
rz(-1.0836243) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.590893) q[1];
sx q[1];
rz(-1.8807481) q[1];
sx q[1];
rz(-0.28642408) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1296785) q[3];
sx q[3];
rz(-0.77687009) q[3];
sx q[3];
rz(-2.7912031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622793) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(1.7904003) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(1.321235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6141339) q[0];
sx q[0];
rz(-1.2686994) q[0];
sx q[0];
rz(2.6149349) q[0];
rz(1.6826019) q[2];
sx q[2];
rz(-2.4833792) q[2];
sx q[2];
rz(1.3629701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3258354) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
rz(-1.9923237) q[3];
sx q[3];
rz(-1.1771923) q[3];
sx q[3];
rz(0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-0.72285405) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-0.29108873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81924651) q[0];
sx q[0];
rz(-0.3194321) q[0];
sx q[0];
rz(0.56630212) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5076324) q[2];
sx q[2];
rz(-0.77253714) q[2];
sx q[2];
rz(-0.26320266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-0.8833589) q[1];
sx q[1];
rz(2.8836125) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17125968) q[3];
sx q[3];
rz(-0.75920877) q[3];
sx q[3];
rz(1.1472792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(2.843294) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-2.1098302) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4435972) q[0];
sx q[0];
rz(-1.5866318) q[0];
sx q[0];
rz(-1.3965194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.736892) q[2];
sx q[2];
rz(-1.7621303) q[2];
sx q[2];
rz(1.9218685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25615869) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(-2.7536254) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86182819) q[3];
sx q[3];
rz(-2.3579881) q[3];
sx q[3];
rz(2.6605822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(1.3708699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2751728) q[0];
sx q[0];
rz(-0.93979561) q[0];
sx q[0];
rz(-1.5874528) q[0];
x q[1];
rz(-0.15372686) q[2];
sx q[2];
rz(-1.2187374) q[2];
sx q[2];
rz(-1.597126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.26105598) q[1];
sx q[1];
rz(-2.4789101) q[1];
sx q[1];
rz(-0.39775325) q[1];
x q[2];
rz(-2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(0.48279253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2844598) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-2.3585368) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.3979744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180041) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(0.82360928) q[0];
x q[1];
rz(-2.2538848) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(-0.53149022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74136342) q[1];
sx q[1];
rz(-1.5796164) q[1];
sx q[1];
rz(-2.2386902) q[1];
rz(-pi) q[2];
rz(-1.3574202) q[3];
sx q[3];
rz(-0.7979352) q[3];
sx q[3];
rz(-0.40003451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(1.4199055) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0592151) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(-1.5060471) q[0];
rz(-pi) q[1];
rz(-1.8516638) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(-0.40058595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7560278) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(1.0977618) q[1];
rz(-pi) q[2];
rz(-2.3929651) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(-0.77881294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-0.43711883) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(-0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8932874) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(0.020123578) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6548806) q[2];
sx q[2];
rz(-2.463474) q[2];
sx q[2];
rz(0.93914062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5367914) q[1];
sx q[1];
rz(-1.7030506) q[1];
sx q[1];
rz(2.3459646) q[1];
rz(-pi) q[2];
rz(2.203412) q[3];
sx q[3];
rz(-1.6611735) q[3];
sx q[3];
rz(-2.602102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71020469) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-2.6999264) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(-0.031671192) q[2];
sx q[2];
rz(-1.1902255) q[2];
sx q[2];
rz(-2.5524216) q[2];
rz(0.15968237) q[3];
sx q[3];
rz(-0.57542141) q[3];
sx q[3];
rz(0.66948359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
