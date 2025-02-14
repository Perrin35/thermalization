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
rz(0.2001065) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0845153) q[0];
sx q[0];
rz(-2.6383218) q[0];
sx q[0];
rz(2.0664735) q[0];
rz(2.4780689) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(0.51807846) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3568452) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(3.0650782) q[1];
x q[2];
rz(-1.0947919) q[3];
sx q[3];
rz(-1.6443569) q[3];
sx q[3];
rz(1.6122163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(1.3145831) q[2];
rz(-0.9032816) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(-0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(0.85579175) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44875328) q[0];
sx q[0];
rz(-1.8215382) q[0];
sx q[0];
rz(-0.022912774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11779479) q[2];
sx q[2];
rz(-0.34644768) q[2];
sx q[2];
rz(2.7730377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.012177906) q[1];
sx q[1];
rz(-1.2316634) q[1];
sx q[1];
rz(-1.1147092) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7807021) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(-2.8950952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(0.39609972) q[2];
rz(1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.7861231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(-0.060610108) q[0];
rz(-0.68471471) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(0.33904591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.402917) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(3.1039799) q[0];
rz(-pi) q[1];
rz(1.6940368) q[2];
sx q[2];
rz(-1.7584821) q[2];
sx q[2];
rz(-2.0579684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.590893) q[1];
sx q[1];
rz(-1.2608445) q[1];
sx q[1];
rz(-0.28642408) q[1];
x q[2];
rz(-1.5825082) q[3];
sx q[3];
rz(-0.79399601) q[3];
sx q[3];
rz(2.7744966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(-2.2792234) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(-1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.321235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(-0.55498755) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4589908) q[2];
sx q[2];
rz(-2.4833792) q[2];
sx q[2];
rz(-1.3629701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76391782) q[1];
sx q[1];
rz(-1.498068) q[1];
sx q[1];
rz(0.12158981) q[1];
x q[2];
rz(1.1492689) q[3];
sx q[3];
rz(-1.9644004) q[3];
sx q[3];
rz(-0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9185751) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-0.72285405) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-0.55996672) q[0];
rz(-0.2050744) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-0.29108873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2084853) q[0];
sx q[0];
rz(-1.7400842) q[0];
sx q[0];
rz(-2.8693958) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.061441378) q[2];
sx q[2];
rz(-0.80020088) q[2];
sx q[2];
rz(-0.35129181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-0.8833589) q[1];
sx q[1];
rz(2.8836125) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7311312) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(-1.7602518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2228955) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(-0.29829868) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(-1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(-0.6024012) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(2.1098302) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69799549) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(1.7450733) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7784836) q[2];
sx q[2];
rz(-1.9676932) q[2];
sx q[2];
rz(0.43235052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9578625) q[1];
sx q[1];
rz(-1.0938083) q[1];
sx q[1];
rz(1.3579923) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57547456) q[3];
sx q[3];
rz(-1.0053952) q[3];
sx q[3];
rz(0.3994715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80456698) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(-2.2931781) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(0.99789944) q[0];
rz(2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(-1.7707228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3034015) q[0];
sx q[0];
rz(-0.63119054) q[0];
sx q[0];
rz(3.1188008) q[0];
x q[1];
rz(2.9878658) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(-1.5444666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6295885) q[1];
sx q[1];
rz(-1.3301714) q[1];
sx q[1];
rz(-0.62368599) q[1];
rz(-pi) q[2];
rz(2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(2.7992547) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.3979744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358855) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(0.82360928) q[0];
x q[1];
rz(-0.2335642) q[2];
sx q[2];
rz(-1.2937045) q[2];
sx q[2];
rz(-1.2472926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82247558) q[1];
sx q[1];
rz(-2.2386595) q[1];
sx q[1];
rz(0.0112337) q[1];
rz(-1.7841724) q[3];
sx q[3];
rz(-2.3436574) q[3];
sx q[3];
rz(-0.40003451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.8208549) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19972292) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(1.0844213) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444839) q[0];
sx q[0];
rz(-0.86291828) q[0];
sx q[0];
rz(-3.0860712) q[0];
rz(-pi) q[1];
rz(-2.886359) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(-1.9013426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3855648) q[1];
sx q[1];
rz(-1.6431019) q[1];
sx q[1];
rz(-2.0438309) q[1];
rz(-pi) q[2];
rz(0.74862759) q[3];
sx q[3];
rz(-1.7065062) q[3];
sx q[3];
rz(0.77881294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(-1.7982177) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(-0.04537151) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338884) q[0];
sx q[0];
rz(-1.5844463) q[0];
sx q[0];
rz(-0.74539124) q[0];
rz(-1.2104697) q[2];
sx q[2];
rz(-0.98305741) q[2];
sx q[2];
rz(-0.34229842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1677959) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.7586437) q[1];
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
x q[1];
rz(2.431388) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(-0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(2.9001111) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63856335) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
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
rz(0.56959116) q[3];
sx q[3];
rz(-1.4841595) q[3];
sx q[3];
rz(-1.0356173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
