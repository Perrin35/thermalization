OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122983) q[0];
sx q[0];
rz(-2.0295143) q[0];
sx q[0];
rz(2.8192768) q[0];
x q[1];
rz(-1.5775561) q[2];
sx q[2];
rz(-1.5464371) q[2];
sx q[2];
rz(1.3855795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8569021) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(-2.2536709) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0715277) q[3];
sx q[3];
rz(-1.3032039) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(-2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(1.82812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.27539) q[0];
sx q[0];
rz(-1.7821527) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-pi) q[1];
rz(2.1140695) q[2];
sx q[2];
rz(-0.55210219) q[2];
sx q[2];
rz(3.0039624) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4310415) q[1];
sx q[1];
rz(-1.6574873) q[1];
sx q[1];
rz(-2.980742) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6766657) q[3];
sx q[3];
rz(-1.4956116) q[3];
sx q[3];
rz(-0.74010805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(2.9755039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27667339) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(2.9898781) q[0];
x q[1];
rz(1.0679929) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(-1.3082248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3661763) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(0.52904769) q[1];
rz(-pi) q[2];
rz(-2.283964) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(1.5426202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(2.3667483) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(-0.073547151) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5261126) q[2];
sx q[2];
rz(-1.9384906) q[2];
sx q[2];
rz(1.1553264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2312154) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(1.3740205) q[1];
rz(-pi) q[2];
rz(-2.9702441) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(-1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-3.0338874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340794) q[0];
sx q[0];
rz(-2.7764752) q[0];
sx q[0];
rz(-1.8318729) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92685076) q[2];
sx q[2];
rz(-2.0276208) q[2];
sx q[2];
rz(2.9608375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.50385034) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(-0.99650683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7961736) q[3];
sx q[3];
rz(-1.0343699) q[3];
sx q[3];
rz(2.9649343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(1.0169792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74422979) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(2.7683393) q[0];
rz(-1.5469656) q[2];
sx q[2];
rz(-1.7719291) q[2];
sx q[2];
rz(-0.80326524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6974666) q[1];
sx q[1];
rz(-2.0545469) q[1];
sx q[1];
rz(0.53375785) q[1];
x q[2];
rz(0.53652699) q[3];
sx q[3];
rz(-2.0201207) q[3];
sx q[3];
rz(0.41538737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-2.5069359) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780753) q[0];
sx q[0];
rz(-3.0419555) q[0];
sx q[0];
rz(2.7479322) q[0];
rz(-pi) q[1];
rz(-1.9856521) q[2];
sx q[2];
rz(-1.875669) q[2];
sx q[2];
rz(0.65326234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19983229) q[1];
sx q[1];
rz(-2.4344749) q[1];
sx q[1];
rz(-1.9553528) q[1];
rz(-pi) q[2];
rz(0.089902417) q[3];
sx q[3];
rz(-1.9780759) q[3];
sx q[3];
rz(-0.75059964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(-1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9650139) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(0.9126419) q[0];
rz(-pi) q[1];
rz(0.51034575) q[2];
sx q[2];
rz(-0.58594698) q[2];
sx q[2];
rz(2.3454391) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2891149) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(2.8647793) q[1];
rz(-2.7301634) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-0.065447741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2239969) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86635877) q[0];
sx q[0];
rz(-0.70952053) q[0];
sx q[0];
rz(1.3225681) q[0];
x q[1];
rz(-1.6998859) q[2];
sx q[2];
rz(-2.522905) q[2];
sx q[2];
rz(2.8512851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0867882) q[1];
sx q[1];
rz(-1.0541704) q[1];
sx q[1];
rz(-2.7727142) q[1];
rz(-pi) q[2];
rz(-0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(1.5411351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.3814829) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927267) q[0];
sx q[0];
rz(-2.6377569) q[0];
sx q[0];
rz(-3.0548274) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9276766) q[2];
sx q[2];
rz(-0.48366085) q[2];
sx q[2];
rz(-2.1733401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8379412) q[1];
sx q[1];
rz(-0.90180574) q[1];
sx q[1];
rz(2.6020357) q[1];
rz(-pi) q[2];
rz(-0.85716565) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(0.21881783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-1.3163153) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(-2.9968895) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
