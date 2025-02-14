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
rz(1.8097635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.645675) q[0];
sx q[0];
rz(-1.1326651) q[0];
sx q[0];
rz(2.8854831) q[0];
rz(-pi) q[1];
rz(-0.66352377) q[2];
sx q[2];
rz(-2.4276456) q[2];
sx q[2];
rz(2.6235142) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3568452) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(-0.076514449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0468007) q[3];
sx q[3];
rz(-1.6443569) q[3];
sx q[3];
rz(1.6122163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0693822) q[2];
sx q[2];
rz(-2.78077) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-2.2858009) q[0];
rz(-2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1163569) q[0];
sx q[0];
rz(-1.5486002) q[0];
sx q[0];
rz(1.3199914) q[0];
x q[1];
rz(-2.79736) q[2];
sx q[2];
rz(-1.5308799) q[2];
sx q[2];
rz(1.0913864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1545002) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(2.2461056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8773852) q[3];
sx q[3];
rz(-2.4572861) q[3];
sx q[3];
rz(0.58240376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(1.570545) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(1.7861231) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876672) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(0.33904591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0051088) q[0];
sx q[0];
rz(-1.5500796) q[0];
sx q[0];
rz(-2.5584975) q[0];
x q[1];
rz(-2.9525063) q[2];
sx q[2];
rz(-1.6918618) q[2];
sx q[2];
rz(-2.6775286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.590893) q[1];
sx q[1];
rz(-1.8807481) q[1];
sx q[1];
rz(-2.8551686) q[1];
x q[2];
rz(2.364758) q[3];
sx q[3];
rz(-1.5791487) q[3];
sx q[3];
rz(1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48245779) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(2.733574) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(-1.321235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21463284) q[0];
sx q[0];
rz(-1.0702735) q[0];
sx q[0];
rz(1.9167801) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6826019) q[2];
sx q[2];
rz(-0.65821346) q[2];
sx q[2];
rz(-1.7786225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8157573) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
rz(-pi) q[2];
rz(-2.7144954) q[3];
sx q[3];
rz(-1.9582886) q[3];
sx q[3];
rz(-2.1807293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(0.72285405) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-0.55996672) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-0.29108873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223461) q[0];
sx q[0];
rz(-2.8221606) q[0];
sx q[0];
rz(0.56630212) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6339602) q[2];
sx q[2];
rz(-0.77253714) q[2];
sx q[2];
rz(-0.26320266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4704125) q[1];
sx q[1];
rz(-0.7268097) q[1];
sx q[1];
rz(1.8720759) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7311312) q[3];
sx q[3];
rz(-0.82538) q[3];
sx q[3];
rz(1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(2.843294) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(-1.0317624) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5866318) q[0];
sx q[0];
rz(1.7450733) q[0];
rz(-pi) q[1];
rz(-0.40470064) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(-1.2197242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25615869) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(-2.7536254) q[1];
rz(-pi) q[2];
rz(2.2797645) q[3];
sx q[3];
rz(-2.3579881) q[3];
sx q[3];
rz(-0.48101048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80456698) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1161716) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(-1.7707228) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8557963) q[0];
sx q[0];
rz(-1.5573475) q[0];
sx q[0];
rz(0.6310668) q[0];
rz(-pi) q[1];
rz(-1.2148803) q[2];
sx q[2];
rz(-1.4265665) q[2];
sx q[2];
rz(-3.1145417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6295885) q[1];
sx q[1];
rz(-1.3301714) q[1];
sx q[1];
rz(-0.62368599) q[1];
rz(-pi) q[2];
rz(-0.79973914) q[3];
sx q[3];
rz(-0.99501401) q[3];
sx q[3];
rz(-2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.85713282) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(-2.7992547) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.7436183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358855) q[0];
sx q[0];
rz(-0.62380416) q[0];
sx q[0];
rz(0.82360928) q[0];
rz(-pi) q[1];
rz(-1.2863761) q[2];
sx q[2];
rz(-1.3462974) q[2];
sx q[2];
rz(0.38849354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74136342) q[1];
sx q[1];
rz(-1.5796164) q[1];
sx q[1];
rz(2.2386902) q[1];
rz(1.3574202) q[3];
sx q[3];
rz(-0.7979352) q[3];
sx q[3];
rz(-2.7415581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.8208549) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0592151) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(-1.6355455) q[0];
rz(-pi) q[1];
rz(0.83553548) q[2];
sx q[2];
rz(-0.37097574) q[2];
sx q[2];
rz(0.46984601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3855648) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(-2.0438309) q[1];
x q[2];
rz(1.3865269) q[3];
sx q[3];
rz(-0.83067719) q[3];
sx q[3];
rz(-2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(-3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2756903) q[0];
sx q[0];
rz(-0.82549107) q[0];
sx q[0];
rz(-1.5893713) q[0];
x q[1];
rz(-1.9311229) q[2];
sx q[2];
rz(-2.1585352) q[2];
sx q[2];
rz(2.7992942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9737967) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.7586437) q[1];
x q[2];
rz(-1.4187063) q[3];
sx q[3];
rz(-0.63816164) q[3];
sx q[3];
rz(-1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71020469) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030293) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(0.031671192) q[2];
sx q[2];
rz(-1.9513672) q[2];
sx q[2];
rz(0.58917108) q[2];
rz(1.6735703) q[3];
sx q[3];
rz(-2.1379875) q[3];
sx q[3];
rz(0.47982346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
