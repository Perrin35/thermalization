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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(1.1507432) q[0];
rz(-0.28252217) q[1];
sx q[1];
rz(-2.0158267) q[1];
sx q[1];
rz(1.1123302) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.756506) q[0];
sx q[0];
rz(-1.6876564) q[0];
sx q[0];
rz(2.8782842) q[0];
rz(-1.4006279) q[2];
sx q[2];
rz(-2.552816) q[2];
sx q[2];
rz(-2.6937304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9291275) q[1];
sx q[1];
rz(-1.0914789) q[1];
sx q[1];
rz(0.95770897) q[1];
x q[2];
rz(0.78948094) q[3];
sx q[3];
rz(-2.1095358) q[3];
sx q[3];
rz(-1.0305424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7925966) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(1.3953155) q[2];
rz(0.058517728) q[3];
sx q[3];
rz(-1.4165712) q[3];
sx q[3];
rz(1.3099366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7175452) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(2.6918217) q[0];
rz(-2.6768118) q[1];
sx q[1];
rz(-1.8315146) q[1];
sx q[1];
rz(-0.25258499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16134444) q[0];
sx q[0];
rz(-0.35156116) q[0];
sx q[0];
rz(1.5241966) q[0];
rz(-0.54754642) q[2];
sx q[2];
rz(-2.4035168) q[2];
sx q[2];
rz(-2.944327) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79073097) q[1];
sx q[1];
rz(-1.8523895) q[1];
sx q[1];
rz(-0.0029723309) q[1];
rz(-pi) q[2];
rz(-2.3303495) q[3];
sx q[3];
rz(-2.0245123) q[3];
sx q[3];
rz(-0.32920255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8978591) q[2];
sx q[2];
rz(-2.2458138) q[2];
sx q[2];
rz(3.1254357) q[2];
rz(2.2195418) q[3];
sx q[3];
rz(-0.99370304) q[3];
sx q[3];
rz(3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.09318) q[0];
sx q[0];
rz(-1.6572297) q[0];
sx q[0];
rz(-0.68159252) q[0];
rz(-0.59773481) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(-0.32424232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3675243) q[0];
sx q[0];
rz(-1.2307967) q[0];
sx q[0];
rz(-2.6794479) q[0];
rz(-pi) q[1];
rz(2.6347972) q[2];
sx q[2];
rz(-1.629771) q[2];
sx q[2];
rz(0.61117327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26962599) q[1];
sx q[1];
rz(-1.6128648) q[1];
sx q[1];
rz(-1.4021238) q[1];
rz(-2.1799654) q[3];
sx q[3];
rz(-1.7097843) q[3];
sx q[3];
rz(-2.9442996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2250259) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(1.0507091) q[2];
rz(0.74058908) q[3];
sx q[3];
rz(-1.2935484) q[3];
sx q[3];
rz(3.0802687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4968313) q[0];
sx q[0];
rz(-0.83503857) q[0];
sx q[0];
rz(-0.9541676) q[0];
rz(2.0046115) q[1];
sx q[1];
rz(-1.515772) q[1];
sx q[1];
rz(2.383393) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6369259) q[0];
sx q[0];
rz(-1.530679) q[0];
sx q[0];
rz(-1.2996496) q[0];
rz(-pi) q[1];
rz(-0.21110837) q[2];
sx q[2];
rz(-0.75497113) q[2];
sx q[2];
rz(-1.3364707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0727194) q[1];
sx q[1];
rz(-2.2617635) q[1];
sx q[1];
rz(3.0145538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.011298374) q[3];
sx q[3];
rz(-2.4027028) q[3];
sx q[3];
rz(2.9113462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0884462) q[2];
sx q[2];
rz(-0.55067486) q[2];
sx q[2];
rz(-1.8268364) q[2];
rz(-2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-2.4563167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984197) q[0];
sx q[0];
rz(-1.5564065) q[0];
sx q[0];
rz(-1.5013303) q[0];
rz(-0.38652626) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(1.5949257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3980557) q[0];
sx q[0];
rz(-2.4185582) q[0];
sx q[0];
rz(-2.8788752) q[0];
x q[1];
rz(0.96225454) q[2];
sx q[2];
rz(-1.8402765) q[2];
sx q[2];
rz(-2.7121921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23473492) q[1];
sx q[1];
rz(-1.4868007) q[1];
sx q[1];
rz(0.43387167) q[1];
rz(-pi) q[2];
rz(-2.5677979) q[3];
sx q[3];
rz(-1.6125896) q[3];
sx q[3];
rz(-0.52252636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74653643) q[2];
sx q[2];
rz(-2.655513) q[2];
sx q[2];
rz(1.0699832) q[2];
rz(-2.2561031) q[3];
sx q[3];
rz(-1.0642137) q[3];
sx q[3];
rz(-1.989367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2162061) q[0];
sx q[0];
rz(-2.8031271) q[0];
sx q[0];
rz(-2.0881407) q[0];
rz(-0.18928754) q[1];
sx q[1];
rz(-2.3129251) q[1];
sx q[1];
rz(-1.4922356) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5439592) q[0];
sx q[0];
rz(-1.9372371) q[0];
sx q[0];
rz(-2.9935915) q[0];
rz(0.21250658) q[2];
sx q[2];
rz(-1.8706609) q[2];
sx q[2];
rz(-0.99606711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63000667) q[1];
sx q[1];
rz(-1.6701856) q[1];
sx q[1];
rz(2.6376827) q[1];
x q[2];
rz(1.6118649) q[3];
sx q[3];
rz(-2.3784935) q[3];
sx q[3];
rz(-1.3535318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3377127) q[2];
sx q[2];
rz(-1.0151851) q[2];
sx q[2];
rz(-2.0886776) q[2];
rz(-3.0247011) q[3];
sx q[3];
rz(-1.0630853) q[3];
sx q[3];
rz(-1.6803928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18148024) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(-2.5970698) q[0];
rz(0.25360423) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(-0.25467083) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38886759) q[0];
sx q[0];
rz(-1.6427186) q[0];
sx q[0];
rz(-3.1353463) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7484863) q[2];
sx q[2];
rz(-1.4908893) q[2];
sx q[2];
rz(2.7167632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42159254) q[1];
sx q[1];
rz(-0.82029712) q[1];
sx q[1];
rz(3.0557219) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5768877) q[3];
sx q[3];
rz(-2.2484591) q[3];
sx q[3];
rz(1.1260179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2250681) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(-0.41898215) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.907828) q[3];
sx q[3];
rz(1.1121496) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58881775) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(-1.8137929) q[0];
rz(1.0325507) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(-2.5520777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96846554) q[0];
sx q[0];
rz(-1.0778946) q[0];
sx q[0];
rz(-2.1987134) q[0];
rz(-pi) q[1];
rz(-0.2915129) q[2];
sx q[2];
rz(-1.8290724) q[2];
sx q[2];
rz(1.8669389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8679598) q[1];
sx q[1];
rz(-0.74791779) q[1];
sx q[1];
rz(2.2272417) q[1];
rz(-2.2240713) q[3];
sx q[3];
rz(-0.60169501) q[3];
sx q[3];
rz(-1.2334005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0699658) q[2];
sx q[2];
rz(-0.67487851) q[2];
sx q[2];
rz(0.91342941) q[2];
rz(2.6895798) q[3];
sx q[3];
rz(-1.858523) q[3];
sx q[3];
rz(1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2695059) q[0];
sx q[0];
rz(-1.0259314) q[0];
sx q[0];
rz(1.6312697) q[0];
rz(1.8271081) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(-0.41183919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5740032) q[0];
sx q[0];
rz(-1.8471944) q[0];
sx q[0];
rz(2.2146691) q[0];
x q[1];
rz(-0.14526315) q[2];
sx q[2];
rz(-2.8194339) q[2];
sx q[2];
rz(-1.3346498) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6399314) q[1];
sx q[1];
rz(-1.2218214) q[1];
sx q[1];
rz(0.56569143) q[1];
x q[2];
rz(1.8985073) q[3];
sx q[3];
rz(-1.9372889) q[3];
sx q[3];
rz(-1.2372563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0044378) q[2];
sx q[2];
rz(-1.6496481) q[2];
sx q[2];
rz(2.906666) q[2];
rz(2.2233985) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(2.0535859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10201564) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(-1.8812195) q[0];
rz(-1.6873987) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(-1.4448602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65439502) q[0];
sx q[0];
rz(-2.5154468) q[0];
sx q[0];
rz(2.3930413) q[0];
rz(2.1865322) q[2];
sx q[2];
rz(-0.79074016) q[2];
sx q[2];
rz(3.0476211) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1989205) q[1];
sx q[1];
rz(-2.4148421) q[1];
sx q[1];
rz(-2.9533378) q[1];
rz(-pi) q[2];
x q[2];
rz(1.904018) q[3];
sx q[3];
rz(-1.929627) q[3];
sx q[3];
rz(-0.24397187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5789648) q[2];
sx q[2];
rz(-0.29890385) q[2];
sx q[2];
rz(-3.01801) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.8653423) q[3];
sx q[3];
rz(-2.6821274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4041483) q[0];
sx q[0];
rz(-1.3442208) q[0];
sx q[0];
rz(-1.5303045) q[0];
rz(-0.96724802) q[1];
sx q[1];
rz(-0.9728685) q[1];
sx q[1];
rz(-2.2546993) q[1];
rz(-1.6254454) q[2];
sx q[2];
rz(-2.1693632) q[2];
sx q[2];
rz(-2.7827435) q[2];
rz(-0.88913118) q[3];
sx q[3];
rz(-2.0778772) q[3];
sx q[3];
rz(-2.7419013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
