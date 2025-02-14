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
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(1.6311837) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(-0.55958334) q[1];
sx q[1];
rz(-0.72996563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529634) q[0];
sx q[0];
rz(-1.7173816) q[0];
sx q[0];
rz(-2.1781871) q[0];
rz(0.59487409) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(1.3485731) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3986721) q[1];
sx q[1];
rz(-1.5389331) q[1];
sx q[1];
rz(-0.58018654) q[1];
rz(-pi) q[2];
rz(-0.9040622) q[3];
sx q[3];
rz(-0.36012563) q[3];
sx q[3];
rz(-1.874543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(-0.81895858) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-0.5994125) q[0];
rz(-2.2741611) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(-2.3960466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2363385) q[0];
sx q[0];
rz(-2.8239692) q[0];
sx q[0];
rz(2.7228217) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33919427) q[2];
sx q[2];
rz(-2.3521947) q[2];
sx q[2];
rz(1.2123002) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88953253) q[1];
sx q[1];
rz(-0.65450689) q[1];
sx q[1];
rz(2.8783074) q[1];
rz(-pi) q[2];
rz(-0.87941951) q[3];
sx q[3];
rz(-1.169765) q[3];
sx q[3];
rz(-1.629231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(-1.7342742) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(-2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5582964) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(1.012828) q[0];
rz(-0.60802513) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(1.8720522) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758666) q[0];
sx q[0];
rz(-2.009543) q[0];
sx q[0];
rz(2.7164396) q[0];
x q[1];
rz(-0.89983268) q[2];
sx q[2];
rz(-1.6353893) q[2];
sx q[2];
rz(3.10499) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7439869) q[1];
sx q[1];
rz(-0.74938801) q[1];
sx q[1];
rz(-0.29008643) q[1];
x q[2];
rz(0.57165159) q[3];
sx q[3];
rz(-0.56662512) q[3];
sx q[3];
rz(2.1376615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(-0.21361175) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13852791) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(-1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-0.52070224) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10258612) q[0];
sx q[0];
rz(-1.163536) q[0];
sx q[0];
rz(-2.7558221) q[0];
rz(-pi) q[1];
rz(-2.7815656) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(-1.6619267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6286875) q[1];
sx q[1];
rz(-2.4081552) q[1];
sx q[1];
rz(-2.4705486) q[1];
x q[2];
rz(-0.52526955) q[3];
sx q[3];
rz(-0.97967463) q[3];
sx q[3];
rz(1.0159462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6167831) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(-1.7899803) q[2];
rz(-2.1221519) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(-1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(-0.098966448) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(2.6720572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5564559) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(1.742247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7623873) q[2];
sx q[2];
rz(-1.0673041) q[2];
sx q[2];
rz(-2.6852754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3217056) q[1];
sx q[1];
rz(-2.4155136) q[1];
sx q[1];
rz(-0.5989845) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23271493) q[3];
sx q[3];
rz(-0.56833) q[3];
sx q[3];
rz(-1.1956904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2609451) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(1.8355231) q[2];
rz(0.4246873) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(-0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.3223883) q[0];
rz(1.127683) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(-1.7599531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2974325) q[0];
sx q[0];
rz(-1.6571665) q[0];
sx q[0];
rz(2.3303836) q[0];
rz(-pi) q[1];
rz(0.83397978) q[2];
sx q[2];
rz(-2.4206941) q[2];
sx q[2];
rz(-2.3490459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4037212) q[1];
sx q[1];
rz(-0.779169) q[1];
sx q[1];
rz(0.26547758) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.840191) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(1.7168728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(-0.48119989) q[2];
rz(2.6324658) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(0.089381889) q[0];
rz(2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(0.99348974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29879323) q[0];
sx q[0];
rz(-1.6090819) q[0];
sx q[0];
rz(0.78193112) q[0];
rz(-0.87331302) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(0.26604929) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74718432) q[1];
sx q[1];
rz(-1.4239053) q[1];
sx q[1];
rz(-1.2042852) q[1];
x q[2];
rz(1.0519846) q[3];
sx q[3];
rz(-2.1362274) q[3];
sx q[3];
rz(-0.65549248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(-2.5636165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680229) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(1.3336257) q[0];
rz(-0.85583055) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(2.2241101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3429392) q[0];
sx q[0];
rz(-2.9173033) q[0];
sx q[0];
rz(2.0464315) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94709227) q[2];
sx q[2];
rz(-1.4741885) q[2];
sx q[2];
rz(1.7283224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6010103) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(2.90592) q[1];
x q[2];
rz(-0.64669426) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(-2.6961117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.8514006) q[2];
rz(0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(-1.2241036) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(-1.3714429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(-2.8369342) q[0];
rz(-pi) q[1];
rz(-0.65176378) q[2];
sx q[2];
rz(-0.65417505) q[2];
sx q[2];
rz(0.68133611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1094233) q[1];
sx q[1];
rz(-2.210683) q[1];
sx q[1];
rz(2.8830322) q[1];
rz(-pi) q[2];
rz(0.84884642) q[3];
sx q[3];
rz(-1.7023785) q[3];
sx q[3];
rz(-2.2717486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(-1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(0.57149291) q[0];
rz(0.68069619) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(-0.64819711) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.307319) q[0];
sx q[0];
rz(-1.5569485) q[0];
sx q[0];
rz(-1.5939006) q[0];
x q[1];
rz(-2.4785751) q[2];
sx q[2];
rz(-2.1923991) q[2];
sx q[2];
rz(2.4746975) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5773864) q[1];
sx q[1];
rz(-1.5485475) q[1];
sx q[1];
rz(-0.69907) q[1];
rz(-2.8043824) q[3];
sx q[3];
rz(-0.54223947) q[3];
sx q[3];
rz(-2.6553939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8615243) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(2.9649949) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8661154) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-1.2550884) q[2];
sx q[2];
rz(-1.9049302) q[2];
sx q[2];
rz(1.9096712) q[2];
rz(-1.0064784) q[3];
sx q[3];
rz(-1.1075533) q[3];
sx q[3];
rz(1.5816734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
