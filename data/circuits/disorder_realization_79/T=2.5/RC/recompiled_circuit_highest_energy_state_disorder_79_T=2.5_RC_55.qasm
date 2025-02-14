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
rz(0.6495629) q[0];
sx q[0];
rz(-0.49512884) q[0];
sx q[0];
rz(0.25547096) q[0];
rz(2.7893692) q[1];
sx q[1];
rz(-1.398634) q[1];
sx q[1];
rz(1.4056828) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79069239) q[0];
sx q[0];
rz(-3.1381099) q[0];
sx q[0];
rz(3.0538959) q[0];
rz(-pi) q[1];
rz(-0.1867577) q[2];
sx q[2];
rz(-2.6180912) q[2];
sx q[2];
rz(-1.40846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97326255) q[1];
sx q[1];
rz(-1.574206) q[1];
sx q[1];
rz(-1.9539321) q[1];
rz(0.051298012) q[3];
sx q[3];
rz(-0.56863943) q[3];
sx q[3];
rz(2.7647247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5070255) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(0.8325038) q[2];
rz(-2.3333874) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(2.8830849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7611564) q[0];
sx q[0];
rz(-0.34270898) q[0];
sx q[0];
rz(-2.9535182) q[0];
rz(3.0694118) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(1.6292876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3609027) q[0];
sx q[0];
rz(-1.6854428) q[0];
sx q[0];
rz(-2.8066382) q[0];
rz(-pi) q[1];
rz(-1.6167655) q[2];
sx q[2];
rz(-1.6013923) q[2];
sx q[2];
rz(-2.033288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20023274) q[1];
sx q[1];
rz(-1.6299452) q[1];
sx q[1];
rz(-0.099353921) q[1];
rz(-pi) q[2];
rz(0.54297437) q[3];
sx q[3];
rz(-0.88490769) q[3];
sx q[3];
rz(-2.7859774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(1.8485273) q[2];
rz(-0.34717789) q[3];
sx q[3];
rz(-2.8721589) q[3];
sx q[3];
rz(-0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8186571) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(-0.47054189) q[0];
rz(1.200354) q[1];
sx q[1];
rz(-0.73089868) q[1];
sx q[1];
rz(-1.8847195) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6171489) q[0];
sx q[0];
rz(-1.2192508) q[0];
sx q[0];
rz(1.0008414) q[0];
x q[1];
rz(1.435661) q[2];
sx q[2];
rz(-1.6622512) q[2];
sx q[2];
rz(0.29558638) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4670406) q[1];
sx q[1];
rz(-1.5853154) q[1];
sx q[1];
rz(-3.1402223) q[1];
rz(-pi) q[2];
rz(-1.9745578) q[3];
sx q[3];
rz(-0.65126538) q[3];
sx q[3];
rz(-2.2869956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5990126) q[2];
sx q[2];
rz(-0.97992367) q[2];
sx q[2];
rz(2.2194594) q[2];
rz(-1.7982091) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(-1.62742) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87329292) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(2.701395) q[0];
rz(-1.6104376) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(-0.24756113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0519686) q[0];
sx q[0];
rz(-1.0680833) q[0];
sx q[0];
rz(-2.6460656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6610959) q[2];
sx q[2];
rz(-1.4019499) q[2];
sx q[2];
rz(2.894998) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2046856) q[1];
sx q[1];
rz(-1.2720894) q[1];
sx q[1];
rz(-0.077110962) q[1];
rz(-pi) q[2];
rz(-0.64856802) q[3];
sx q[3];
rz(-1.5463136) q[3];
sx q[3];
rz(0.45059965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89746785) q[2];
sx q[2];
rz(-2.1110057) q[2];
sx q[2];
rz(-1.4424651) q[2];
rz(-2.9407732) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(-1.1403181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15631256) q[0];
sx q[0];
rz(-0.15287481) q[0];
sx q[0];
rz(0.54541624) q[0];
rz(-0.044005752) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(0.47637475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6355609) q[0];
sx q[0];
rz(-0.2579435) q[0];
sx q[0];
rz(2.8057465) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9721932) q[2];
sx q[2];
rz(-1.1938022) q[2];
sx q[2];
rz(2.4103386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4551983) q[1];
sx q[1];
rz(-2.2664811) q[1];
sx q[1];
rz(0.2376539) q[1];
x q[2];
rz(3.0073193) q[3];
sx q[3];
rz(-1.065514) q[3];
sx q[3];
rz(3.0485632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5215317) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(2.7812092) q[2];
rz(-0.22200577) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(-1.41159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5905269) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(-1.5565514) q[0];
rz(-3.0146154) q[1];
sx q[1];
rz(-1.8366837) q[1];
sx q[1];
rz(0.089687673) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68192775) q[0];
sx q[0];
rz(-1.0583504) q[0];
sx q[0];
rz(-2.1905023) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61059606) q[2];
sx q[2];
rz(-2.2572569) q[2];
sx q[2];
rz(-2.7629791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4597923) q[1];
sx q[1];
rz(-0.55306095) q[1];
sx q[1];
rz(-0.44981594) q[1];
rz(2.5406214) q[3];
sx q[3];
rz(-1.2699763) q[3];
sx q[3];
rz(-1.5301289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.108532) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(-0.25165558) q[2];
rz(-0.039994914) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(-2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43117487) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(2.726626) q[0];
rz(-1.6522495) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(0.32589486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.519442) q[0];
sx q[0];
rz(-1.9339193) q[0];
sx q[0];
rz(-2.0456929) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0517111) q[2];
sx q[2];
rz(-2.0873859) q[2];
sx q[2];
rz(-0.84026779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2532172) q[1];
sx q[1];
rz(-1.0723812) q[1];
sx q[1];
rz(1.2775757) q[1];
rz(-pi) q[2];
rz(-1.8741499) q[3];
sx q[3];
rz(-1.6432495) q[3];
sx q[3];
rz(-0.8006351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2275647) q[2];
sx q[2];
rz(-1.808337) q[2];
sx q[2];
rz(-2.0951159) q[2];
rz(2.922831) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(-1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662812) q[0];
sx q[0];
rz(-3.1240211) q[0];
sx q[0];
rz(-2.6783491) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(-1.499768) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69073286) q[0];
sx q[0];
rz(-1.5520649) q[0];
sx q[0];
rz(-2.0035726) q[0];
rz(-1.5157264) q[2];
sx q[2];
rz(-2.0158421) q[2];
sx q[2];
rz(2.9663756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25019885) q[1];
sx q[1];
rz(-1.3490744) q[1];
sx q[1];
rz(0.029832928) q[1];
rz(1.3070694) q[3];
sx q[3];
rz(-2.8805507) q[3];
sx q[3];
rz(1.3021951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90193343) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(1.8233914) q[2];
rz(-3.1384595) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(1.1656632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779293) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(2.928012) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.2108796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0830363) q[0];
sx q[0];
rz(-2.8175111) q[0];
sx q[0];
rz(-1.9084318) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69542747) q[2];
sx q[2];
rz(-1.8423242) q[2];
sx q[2];
rz(-2.5144387) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93608078) q[1];
sx q[1];
rz(-1.2911564) q[1];
sx q[1];
rz(2.1195009) q[1];
rz(-pi) q[2];
rz(-1.3812495) q[3];
sx q[3];
rz(-1.726103) q[3];
sx q[3];
rz(1.3819753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1750298) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(0.35061947) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(-1.143379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38458934) q[0];
sx q[0];
rz(-2.2191255) q[0];
sx q[0];
rz(-2.5272227) q[0];
rz(-0.059582926) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(-1.4245859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40648288) q[0];
sx q[0];
rz(-0.87672675) q[0];
sx q[0];
rz(0.31882341) q[0];
x q[1];
rz(0.13127998) q[2];
sx q[2];
rz(-1.5688217) q[2];
sx q[2];
rz(-1.9406089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9944133) q[1];
sx q[1];
rz(-2.96619) q[1];
sx q[1];
rz(-2.0270623) q[1];
rz(-pi) q[2];
rz(-2.0395117) q[3];
sx q[3];
rz(-2.6913683) q[3];
sx q[3];
rz(-2.1463822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0397348) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-0.52379918) q[2];
rz(2.0959181) q[3];
sx q[3];
rz(-1.2751251) q[3];
sx q[3];
rz(-0.4642134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47281784) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(-0.642943) q[2];
sx q[2];
rz(-1.3284773) q[2];
sx q[2];
rz(0.20198573) q[2];
rz(-0.37665357) q[3];
sx q[3];
rz(-0.49643825) q[3];
sx q[3];
rz(1.7986922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
