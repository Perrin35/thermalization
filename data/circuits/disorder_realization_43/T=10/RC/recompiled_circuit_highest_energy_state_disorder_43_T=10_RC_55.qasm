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
rz(2.645675) q[0];
sx q[0];
rz(-1.1326651) q[0];
sx q[0];
rz(2.8854831) q[0];
rz(1.0806141) q[2];
sx q[2];
rz(-1.0287544) q[2];
sx q[2];
rz(1.8212183) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9871949) q[1];
sx q[1];
rz(-1.5227277) q[1];
sx q[1];
rz(-0.67832077) q[1];
x q[2];
rz(1.41134) q[3];
sx q[3];
rz(-0.48122367) q[3];
sx q[3];
rz(3.0414097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.8270095) q[2];
rz(-0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-2.2858009) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928394) q[0];
sx q[0];
rz(-1.3200545) q[0];
sx q[0];
rz(0.022912774) q[0];
x q[1];
rz(2.79736) q[2];
sx q[2];
rz(-1.5308799) q[2];
sx q[2];
rz(-1.0913864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1545002) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(2.2461056) q[1];
x q[2];
rz(0.26420748) q[3];
sx q[3];
rz(-2.4572861) q[3];
sx q[3];
rz(-0.58240376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(-0.68471471) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(2.8025467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1364839) q[0];
sx q[0];
rz(-1.5500796) q[0];
sx q[0];
rz(-2.5584975) q[0];
rz(2.9525063) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(0.4640641) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0319034) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(1.8930356) q[1];
rz(-0.77683461) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(-1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(-1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(1.3511924) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42962881) q[0];
sx q[0];
rz(-2.5416059) q[0];
sx q[0];
rz(0.55498755) q[0];
rz(-pi) q[1];
rz(3.0555326) q[2];
sx q[2];
rz(-0.91740184) q[2];
sx q[2];
rz(-1.5039521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8712052) q[1];
sx q[1];
rz(-0.14158881) q[1];
sx q[1];
rz(2.6006727) q[1];
rz(1.1492689) q[3];
sx q[3];
rz(-1.1771923) q[3];
sx q[3];
rz(0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(2.4187386) q[2];
rz(-0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(2.5816259) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-2.8505039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223461) q[0];
sx q[0];
rz(-2.8221606) q[0];
sx q[0];
rz(-0.56630212) q[0];
rz(-pi) q[1];
rz(0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(-1.8792626) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-2.2582338) q[1];
sx q[1];
rz(-2.8836125) q[1];
x q[2];
rz(-1.7311312) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(-1.7602518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(2.5391915) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(1.0317624) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2715797) q[0];
sx q[0];
rz(-1.3965415) q[0];
sx q[0];
rz(0.016079024) q[0];
rz(-pi) q[1];
rz(-0.40470064) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(1.9218685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.885434) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(-0.38796723) q[1];
rz(2.5661181) q[3];
sx q[3];
rz(-1.0053952) q[3];
sx q[3];
rz(-0.3994715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80456698) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(-0.8484146) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(-1.3708699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2751728) q[0];
sx q[0];
rz(-2.201797) q[0];
sx q[0];
rz(-1.5874528) q[0];
rz(-pi) q[1];
rz(1.2148803) q[2];
sx q[2];
rz(-1.7150262) q[2];
sx q[2];
rz(0.027050935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26105598) q[1];
sx q[1];
rz(-0.66268259) q[1];
sx q[1];
rz(-2.7438394) q[1];
rz(-pi) q[2];
rz(-2.4059345) q[3];
sx q[3];
rz(-0.94649411) q[3];
sx q[3];
rz(1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0520332) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(2.3585368) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0692978) q[0];
sx q[0];
rz(-1.1279449) q[0];
sx q[0];
rz(-0.45486562) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2538848) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(0.53149022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74136342) q[1];
sx q[1];
rz(-1.5796164) q[1];
sx q[1];
rz(0.90290248) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9277738) q[3];
sx q[3];
rz(-2.345746) q[3];
sx q[3];
rz(-2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.8208549) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(1.0844213) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0592151) q[0];
sx q[0];
rz(-0.70967662) q[0];
sx q[0];
rz(1.5060471) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25523369) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(-1.24025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.045001205) q[1];
sx q[1];
rz(-2.6634779) q[1];
sx q[1];
rz(-1.7284615) q[1];
rz(-0.19799216) q[3];
sx q[3];
rz(-0.75847236) q[3];
sx q[3];
rz(-0.93659479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(2.7044738) q[2];
rz(1.3433749) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.4105256) q[0];
rz(0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338884) q[0];
sx q[0];
rz(-1.5844463) q[0];
sx q[0];
rz(-0.74539124) q[0];
x q[1];
rz(2.5228517) q[2];
sx q[2];
rz(-1.2729984) q[2];
sx q[2];
rz(-2.1190475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9737967) q[1];
sx q[1];
rz(-2.3575511) q[1];
sx q[1];
rz(1.382949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7228863) q[3];
sx q[3];
rz(-0.63816164) q[3];
sx q[3];
rz(-1.1537976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
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
