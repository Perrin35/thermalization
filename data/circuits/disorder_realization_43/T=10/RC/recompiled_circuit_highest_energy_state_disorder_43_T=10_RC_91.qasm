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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0845153) q[0];
sx q[0];
rz(-0.50327089) q[0];
sx q[0];
rz(1.0751192) q[0];
x q[1];
rz(1.0806141) q[2];
sx q[2];
rz(-2.1128383) q[2];
sx q[2];
rz(-1.8212183) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
x q[2];
rz(-1.0947919) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(1.5293763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(0.9032816) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(-2.7956853) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(-2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(0.78786293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928394) q[0];
sx q[0];
rz(-1.8215382) q[0];
sx q[0];
rz(-3.1186799) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0237979) q[2];
sx q[2];
rz(-2.795145) q[2];
sx q[2];
rz(0.36855498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98709244) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(-2.2461056) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7807021) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(-0.24649749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(-1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(1.7861231) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(0.060610108) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-2.8025467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7386757) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(-3.1039799) q[0];
rz(-0.18908638) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(-2.6775286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.590893) q[1];
sx q[1];
rz(-1.2608445) q[1];
sx q[1];
rz(-2.8551686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.364758) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(1.9296822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(0.40801868) q[0];
rz(1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.321235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21463284) q[0];
sx q[0];
rz(-1.0702735) q[0];
sx q[0];
rz(1.2248125) q[0];
rz(-pi) q[1];
rz(-1.4589908) q[2];
sx q[2];
rz(-0.65821346) q[2];
sx q[2];
rz(1.7786225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27038747) q[1];
sx q[1];
rz(-0.14158881) q[1];
sx q[1];
rz(-0.54091994) q[1];
rz(-2.3636171) q[3];
sx q[3];
rz(-0.56853349) q[3];
sx q[3];
rz(-1.83873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(0.72285405) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(2.8505039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223461) q[0];
sx q[0];
rz(-0.3194321) q[0];
sx q[0];
rz(0.56630212) q[0];
x q[1];
rz(0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(-1.8792626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8644789) q[1];
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
x q[1];
rz(-1.2228955) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(-1.1693303) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(1.5112618) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-2.1098302) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5866318) q[0];
sx q[0];
rz(1.7450733) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.736892) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(1.2197242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6556485) q[1];
sx q[1];
rz(-1.7595425) q[1];
sx q[1];
rz(2.6552378) q[1];
x q[2];
rz(-0.57547456) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(0.3994715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(-1.5186914) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(3.1019822) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(0.99789944) q[0];
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
rz(-1.8664198) q[0];
sx q[0];
rz(-0.93979561) q[0];
sx q[0];
rz(-1.5874528) q[0];
x q[1];
rz(1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(1.9667786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22860195) q[1];
sx q[1];
rz(-0.96768846) q[1];
sx q[1];
rz(1.8643537) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73565817) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(-1.566938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(0.029646309) q[2];
rz(-2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072294839) q[0];
sx q[0];
rz(-2.0136478) q[0];
sx q[0];
rz(-0.45486562) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88770788) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(-0.53149022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3009792) q[1];
sx q[1];
rz(-0.66794315) q[1];
sx q[1];
rz(1.5565558) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78432958) q[3];
sx q[3];
rz(-1.4186067) q[3];
sx q[3];
rz(-2.1209716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(-1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19972292) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-2.0571713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53756489) q[0];
sx q[0];
rz(-1.5286235) q[0];
sx q[0];
rz(2.2794363) q[0];
x q[1];
rz(-2.886359) q[2];
sx q[2];
rz(-1.8430147) q[2];
sx q[2];
rz(-1.24025) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.045001205) q[1];
sx q[1];
rz(-2.6634779) q[1];
sx q[1];
rz(1.4131312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74862759) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(0.77881294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(1.3433749) q[3];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042353543) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.7310671) q[0];
rz(-2.9208185) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(0.9332307) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30770424) q[0];
sx q[0];
rz(-1.5844463) q[0];
sx q[0];
rz(-0.74539124) q[0];
rz(-0.48671203) q[2];
sx q[2];
rz(-0.67811869) q[2];
sx q[2];
rz(-0.93914062) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9737967) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(-1.7586437) q[1];
rz(-pi) q[2];
rz(-3.0296928) q[3];
sx q[3];
rz(-2.2004232) q[3];
sx q[3];
rz(0.96523413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.491811) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
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
