OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.31324759) q[0];
sx q[0];
rz(-2.3656443) q[0];
sx q[0];
rz(1.1021855) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(-1.1270539) q[1];
sx q[1];
rz(-2.9161646) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0505053) q[0];
sx q[0];
rz(-0.69982547) q[0];
sx q[0];
rz(-0.096629337) q[0];
rz(0.20689865) q[2];
sx q[2];
rz(-0.62600905) q[2];
sx q[2];
rz(-1.2956784) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.082583383) q[1];
sx q[1];
rz(-1.9349001) q[1];
sx q[1];
rz(2.8541982) q[1];
rz(-pi) q[2];
rz(2.647764) q[3];
sx q[3];
rz(-1.8803616) q[3];
sx q[3];
rz(2.1403596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6585091) q[2];
sx q[2];
rz(-0.67407125) q[2];
sx q[2];
rz(-3.1209514) q[2];
rz(2.9523197) q[3];
sx q[3];
rz(-0.34957379) q[3];
sx q[3];
rz(-1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610158) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(-2.1625157) q[0];
rz(-3.0304404) q[1];
sx q[1];
rz(-0.73768288) q[1];
sx q[1];
rz(2.0806064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2994591) q[0];
sx q[0];
rz(-2.4329022) q[0];
sx q[0];
rz(-0.75229074) q[0];
rz(-pi) q[1];
rz(1.5716288) q[2];
sx q[2];
rz(-1.7082518) q[2];
sx q[2];
rz(-0.47333131) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0506405) q[1];
sx q[1];
rz(-2.0555858) q[1];
sx q[1];
rz(1.3371972) q[1];
x q[2];
rz(-2.8021566) q[3];
sx q[3];
rz(-0.41012479) q[3];
sx q[3];
rz(-2.8462178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7079118) q[2];
sx q[2];
rz(-1.3206864) q[2];
sx q[2];
rz(0.39805463) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-2.8545696) q[3];
sx q[3];
rz(2.6557693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240391) q[0];
sx q[0];
rz(-2.6300639) q[0];
sx q[0];
rz(-1.0523354) q[0];
rz(-2.6984093) q[1];
sx q[1];
rz(-0.96142238) q[1];
sx q[1];
rz(2.4876432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0130413) q[0];
sx q[0];
rz(-1.4440026) q[0];
sx q[0];
rz(-2.3953505) q[0];
rz(-pi) q[1];
rz(-0.072418173) q[2];
sx q[2];
rz(-2.1927367) q[2];
sx q[2];
rz(2.3427013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7591651) q[1];
sx q[1];
rz(-0.8635206) q[1];
sx q[1];
rz(2.1334126) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4846333) q[3];
sx q[3];
rz(-2.8553319) q[3];
sx q[3];
rz(3.1050499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37086481) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(-0.23180836) q[2];
rz(-1.9306463) q[3];
sx q[3];
rz(-1.1823267) q[3];
sx q[3];
rz(-0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4732707) q[0];
sx q[0];
rz(-2.7372161) q[0];
sx q[0];
rz(-2.7647198) q[0];
rz(2.6301774) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(-2.2976141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91746861) q[0];
sx q[0];
rz(-0.11118764) q[0];
sx q[0];
rz(-1.412941) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9241946) q[2];
sx q[2];
rz(-1.0026284) q[2];
sx q[2];
rz(-1.8894486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5201903) q[1];
sx q[1];
rz(-2.108413) q[1];
sx q[1];
rz(-2.138184) q[1];
rz(-0.52167251) q[3];
sx q[3];
rz(-1.5893916) q[3];
sx q[3];
rz(-0.90584318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9722998) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(0.69804066) q[2];
rz(1.1997148) q[3];
sx q[3];
rz(-0.90568393) q[3];
sx q[3];
rz(2.5419295) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9508764) q[0];
sx q[0];
rz(-2.8671725) q[0];
sx q[0];
rz(3.0158667) q[0];
rz(-2.114864) q[1];
sx q[1];
rz(-1.7544361) q[1];
sx q[1];
rz(2.6153807) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3644117) q[0];
sx q[0];
rz(-0.48425774) q[0];
sx q[0];
rz(2.9817192) q[0];
rz(0.16321142) q[2];
sx q[2];
rz(-0.88055748) q[2];
sx q[2];
rz(2.9371967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3354644) q[1];
sx q[1];
rz(-1.6732235) q[1];
sx q[1];
rz(2.7598937) q[1];
rz(-pi) q[2];
rz(-2.4496042) q[3];
sx q[3];
rz(-0.92433911) q[3];
sx q[3];
rz(-0.54997589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0937664) q[2];
sx q[2];
rz(-2.3036849) q[2];
sx q[2];
rz(-0.31025904) q[2];
rz(-0.034916498) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0763615) q[0];
sx q[0];
rz(-1.4820453) q[0];
sx q[0];
rz(-0.37102997) q[0];
rz(-0.064844355) q[1];
sx q[1];
rz(-0.82227451) q[1];
sx q[1];
rz(-1.2670021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112087) q[0];
sx q[0];
rz(-0.74978854) q[0];
sx q[0];
rz(1.0553589) q[0];
rz(-2.7899988) q[2];
sx q[2];
rz(-2.1702655) q[2];
sx q[2];
rz(-0.64707398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5090086) q[1];
sx q[1];
rz(-1.1727837) q[1];
sx q[1];
rz(1.6081612) q[1];
rz(-pi) q[2];
rz(-1.2812595) q[3];
sx q[3];
rz(-2.0199593) q[3];
sx q[3];
rz(1.7790573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12396699) q[2];
sx q[2];
rz(-1.664868) q[2];
sx q[2];
rz(-1.5383447) q[2];
rz(0.41687837) q[3];
sx q[3];
rz(-0.50691253) q[3];
sx q[3];
rz(-2.9712334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040319547) q[0];
sx q[0];
rz(-2.9077001) q[0];
sx q[0];
rz(-2.7129569) q[0];
rz(-2.0173232) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(2.3720149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452827) q[0];
sx q[0];
rz(-0.0036024139) q[0];
sx q[0];
rz(1.035319) q[0];
rz(2.3891201) q[2];
sx q[2];
rz(-0.66900775) q[2];
sx q[2];
rz(-0.7758918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4691041) q[1];
sx q[1];
rz(-2.1616363) q[1];
sx q[1];
rz(1.2149151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0817716) q[3];
sx q[3];
rz(-1.8712412) q[3];
sx q[3];
rz(2.1238126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6577242) q[2];
sx q[2];
rz(-1.2843853) q[2];
sx q[2];
rz(-0.095495187) q[2];
rz(-0.5419845) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(-3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1456881) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(0.11181871) q[0];
rz(-1.3189141) q[1];
sx q[1];
rz(-1.2448064) q[1];
sx q[1];
rz(1.8359312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163286) q[0];
sx q[0];
rz(-0.56076196) q[0];
sx q[0];
rz(1.3365585) q[0];
rz(-pi) q[1];
rz(-1.7123004) q[2];
sx q[2];
rz(-1.8496528) q[2];
sx q[2];
rz(-0.57153801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9247465) q[1];
sx q[1];
rz(-1.1976744) q[1];
sx q[1];
rz(-1.2678964) q[1];
rz(-2.5091845) q[3];
sx q[3];
rz(-2.494209) q[3];
sx q[3];
rz(1.5349015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32885113) q[2];
sx q[2];
rz(-1.988827) q[2];
sx q[2];
rz(2.0218772) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(-1.9560811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.8609404) q[0];
sx q[0];
rz(-0.76186162) q[0];
sx q[0];
rz(2.5326488) q[0];
rz(2.2641585) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(0.61224365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5615879) q[0];
sx q[0];
rz(-1.9733493) q[0];
sx q[0];
rz(0.7640362) q[0];
x q[1];
rz(-0.52895542) q[2];
sx q[2];
rz(-0.63807708) q[2];
sx q[2];
rz(2.4412689) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5366) q[1];
sx q[1];
rz(-0.92291622) q[1];
sx q[1];
rz(-0.98439321) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5411392) q[3];
sx q[3];
rz(-1.5039467) q[3];
sx q[3];
rz(-0.56136405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.073254243) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(-0.18745984) q[2];
rz(-1.1120262) q[3];
sx q[3];
rz(-2.2234962) q[3];
sx q[3];
rz(-0.48346564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.40792313) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(-1.6581274) q[1];
sx q[1];
rz(-1.5189867) q[1];
sx q[1];
rz(-0.77267486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7917532) q[0];
sx q[0];
rz(-1.0966705) q[0];
sx q[0];
rz(-1.0081434) q[0];
x q[1];
rz(-1.4336606) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(1.1346045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0341943) q[1];
sx q[1];
rz(-1.6813283) q[1];
sx q[1];
rz(-2.8128036) q[1];
rz(-pi) q[2];
rz(2.3540007) q[3];
sx q[3];
rz(-0.95284546) q[3];
sx q[3];
rz(-1.7389115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3385758) q[2];
sx q[2];
rz(-2.1803941) q[2];
sx q[2];
rz(-0.72688603) q[2];
rz(-1.1072985) q[3];
sx q[3];
rz(-2.6328583) q[3];
sx q[3];
rz(-2.1342962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016639391) q[0];
sx q[0];
rz(-1.4062842) q[0];
sx q[0];
rz(1.9840065) q[0];
rz(-2.6955556) q[1];
sx q[1];
rz(-2.0560494) q[1];
sx q[1];
rz(2.5037419) q[1];
rz(0.28631306) q[2];
sx q[2];
rz(-0.41819093) q[2];
sx q[2];
rz(2.0750194) q[2];
rz(-0.42556006) q[3];
sx q[3];
rz(-1.4360518) q[3];
sx q[3];
rz(-0.19946972) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
