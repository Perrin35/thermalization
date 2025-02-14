OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.705536) q[0];
sx q[0];
rz(-2.2826865) q[0];
sx q[0];
rz(-0.98281759) q[0];
rz(-1.9614027) q[1];
sx q[1];
rz(-0.72620121) q[1];
sx q[1];
rz(0.89816165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58098823) q[0];
sx q[0];
rz(-2.8074671) q[0];
sx q[0];
rz(0.32291205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9299081) q[2];
sx q[2];
rz(-1.058033) q[2];
sx q[2];
rz(0.36836926) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8232728) q[1];
sx q[1];
rz(-1.1538236) q[1];
sx q[1];
rz(-0.61064536) q[1];
rz(-0.14515669) q[3];
sx q[3];
rz(-1.6696128) q[3];
sx q[3];
rz(-2.1414808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0183705) q[2];
sx q[2];
rz(-1.3192588) q[2];
sx q[2];
rz(-0.74750626) q[2];
rz(-1.2627164) q[3];
sx q[3];
rz(-1.6044173) q[3];
sx q[3];
rz(-0.023716299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338407) q[0];
sx q[0];
rz(-3.096014) q[0];
sx q[0];
rz(1.0652834) q[0];
rz(-2.4899958) q[1];
sx q[1];
rz(-0.35531303) q[1];
sx q[1];
rz(1.6337055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9621443) q[0];
sx q[0];
rz(-2.1198258) q[0];
sx q[0];
rz(-2.7864576) q[0];
rz(-pi) q[1];
rz(-0.47453158) q[2];
sx q[2];
rz(-1.8900649) q[2];
sx q[2];
rz(-2.266707) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2245141) q[1];
sx q[1];
rz(-2.0339171) q[1];
sx q[1];
rz(0.45189894) q[1];
x q[2];
rz(-2.1820804) q[3];
sx q[3];
rz(-1.5942845) q[3];
sx q[3];
rz(1.5620205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6642586) q[2];
sx q[2];
rz(-1.6691672) q[2];
sx q[2];
rz(-0.092546917) q[2];
rz(2.7009098) q[3];
sx q[3];
rz(-0.40658545) q[3];
sx q[3];
rz(-1.996076) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6630845) q[0];
sx q[0];
rz(-2.9496851) q[0];
sx q[0];
rz(1.2364016) q[0];
rz(0.23794404) q[1];
sx q[1];
rz(-1.6704208) q[1];
sx q[1];
rz(-0.96756378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2264573) q[0];
sx q[0];
rz(-0.98513705) q[0];
sx q[0];
rz(-2.4047635) q[0];
rz(0.5401436) q[2];
sx q[2];
rz(-0.61068084) q[2];
sx q[2];
rz(1.7941098) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6131276) q[1];
sx q[1];
rz(-1.2674164) q[1];
sx q[1];
rz(0.19618285) q[1];
rz(-pi) q[2];
rz(1.7015412) q[3];
sx q[3];
rz(-2.7249911) q[3];
sx q[3];
rz(0.9804014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48270109) q[2];
sx q[2];
rz(-2.7693558) q[2];
sx q[2];
rz(-0.9185532) q[2];
rz(-2.3679768) q[3];
sx q[3];
rz(-0.88671237) q[3];
sx q[3];
rz(0.95019597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8755662) q[0];
sx q[0];
rz(-2.4022864) q[0];
sx q[0];
rz(0.29228041) q[0];
rz(0.87458163) q[1];
sx q[1];
rz(-2.5129109) q[1];
sx q[1];
rz(-0.94327092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2349512) q[0];
sx q[0];
rz(-2.4138193) q[0];
sx q[0];
rz(1.1207188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9246632) q[2];
sx q[2];
rz(-2.0299021) q[2];
sx q[2];
rz(0.87068671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0318289) q[1];
sx q[1];
rz(-2.2110143) q[1];
sx q[1];
rz(2.7494181) q[1];
rz(-pi) q[2];
rz(1.1499829) q[3];
sx q[3];
rz(-0.61987034) q[3];
sx q[3];
rz(3.0734143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6167407) q[2];
sx q[2];
rz(-2.0276232) q[2];
sx q[2];
rz(-0.88391602) q[2];
rz(2.0058477) q[3];
sx q[3];
rz(-2.2131049) q[3];
sx q[3];
rz(0.69042027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32490548) q[0];
sx q[0];
rz(-3.0563834) q[0];
sx q[0];
rz(-0.19164044) q[0];
rz(1.8457671) q[1];
sx q[1];
rz(-1.9486267) q[1];
sx q[1];
rz(0.4506909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2795581) q[0];
sx q[0];
rz(-2.4992895) q[0];
sx q[0];
rz(-1.4181644) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.500461) q[2];
sx q[2];
rz(-1.8425111) q[2];
sx q[2];
rz(2.4974785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.012659) q[1];
sx q[1];
rz(-1.6427354) q[1];
sx q[1];
rz(2.1484445) q[1];
rz(-pi) q[2];
rz(-0.29286082) q[3];
sx q[3];
rz(-1.5080875) q[3];
sx q[3];
rz(-0.69578275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6777163) q[2];
sx q[2];
rz(-1.8666942) q[2];
sx q[2];
rz(-0.95787588) q[2];
rz(-3.0887582) q[3];
sx q[3];
rz(-0.54532471) q[3];
sx q[3];
rz(-2.7132645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4868454) q[0];
sx q[0];
rz(-1.0423648) q[0];
sx q[0];
rz(-0.47252193) q[0];
rz(-2.3782702) q[1];
sx q[1];
rz(-0.7205874) q[1];
sx q[1];
rz(2.8514013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1520574) q[0];
sx q[0];
rz(-0.88071918) q[0];
sx q[0];
rz(1.2834653) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.824108) q[2];
sx q[2];
rz(-1.0429263) q[2];
sx q[2];
rz(1.9535397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4279033) q[1];
sx q[1];
rz(-0.91304243) q[1];
sx q[1];
rz(1.7005928) q[1];
rz(2.3004856) q[3];
sx q[3];
rz(-2.8770513) q[3];
sx q[3];
rz(0.11716784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7758238) q[2];
sx q[2];
rz(-0.86161986) q[2];
sx q[2];
rz(-2.5422868) q[2];
rz(-1.722909) q[3];
sx q[3];
rz(-1.320188) q[3];
sx q[3];
rz(-2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80239427) q[0];
sx q[0];
rz(-1.1856439) q[0];
sx q[0];
rz(-1.3108569) q[0];
rz(-1.5015548) q[1];
sx q[1];
rz(-0.40032598) q[1];
sx q[1];
rz(0.60360533) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5561605) q[0];
sx q[0];
rz(-1.7045665) q[0];
sx q[0];
rz(1.1473535) q[0];
x q[1];
rz(0.55159388) q[2];
sx q[2];
rz(-1.7579798) q[2];
sx q[2];
rz(2.1647705) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2767553) q[1];
sx q[1];
rz(-0.25028203) q[1];
sx q[1];
rz(1.8486057) q[1];
rz(1.3796174) q[3];
sx q[3];
rz(-0.83836765) q[3];
sx q[3];
rz(0.54938706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32460585) q[2];
sx q[2];
rz(-0.88674712) q[2];
sx q[2];
rz(0.19503197) q[2];
rz(0.65702355) q[3];
sx q[3];
rz(-1.7200836) q[3];
sx q[3];
rz(-3.0278964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8103771) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(3.069416) q[0];
rz(2.3585034) q[1];
sx q[1];
rz(-0.63242811) q[1];
sx q[1];
rz(2.2236688) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72361392) q[0];
sx q[0];
rz(-1.0999014) q[0];
sx q[0];
rz(1.8036475) q[0];
rz(-pi) q[1];
rz(-2.520944) q[2];
sx q[2];
rz(-1.7908515) q[2];
sx q[2];
rz(1.7372075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0705441) q[1];
sx q[1];
rz(-1.3335102) q[1];
sx q[1];
rz(-0.29233934) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9959171) q[3];
sx q[3];
rz(-1.8048058) q[3];
sx q[3];
rz(-2.0834415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5607295) q[2];
sx q[2];
rz(-1.9274638) q[2];
sx q[2];
rz(1.5011935) q[2];
rz(2.255127) q[3];
sx q[3];
rz(-1.3366046) q[3];
sx q[3];
rz(-3.0361573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5378872) q[0];
sx q[0];
rz(-2.7925346) q[0];
sx q[0];
rz(-0.50759298) q[0];
rz(-2.4137068) q[1];
sx q[1];
rz(-1.6517086) q[1];
sx q[1];
rz(-2.0354039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2685251) q[0];
sx q[0];
rz(-1.2374479) q[0];
sx q[0];
rz(-0.26645904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.699643) q[2];
sx q[2];
rz(-0.39296752) q[2];
sx q[2];
rz(0.134274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4708097) q[1];
sx q[1];
rz(-1.9341262) q[1];
sx q[1];
rz(0.35344191) q[1];
x q[2];
rz(-0.98420268) q[3];
sx q[3];
rz(-1.1569258) q[3];
sx q[3];
rz(-2.4936287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20405208) q[2];
sx q[2];
rz(-1.1037408) q[2];
sx q[2];
rz(-2.6943915) q[2];
rz(1.7043381) q[3];
sx q[3];
rz(-0.81366003) q[3];
sx q[3];
rz(-1.6566431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60780418) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(-0.5150038) q[0];
rz(-1.9683413) q[1];
sx q[1];
rz(-1.8517905) q[1];
sx q[1];
rz(0.44949964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4309301) q[0];
sx q[0];
rz(-1.8966513) q[0];
sx q[0];
rz(-0.023825721) q[0];
rz(-0.95566383) q[2];
sx q[2];
rz(-2.0086096) q[2];
sx q[2];
rz(-0.73305886) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3275309) q[1];
sx q[1];
rz(-1.1108634) q[1];
sx q[1];
rz(-0.67332585) q[1];
rz(-pi) q[2];
rz(-0.85136713) q[3];
sx q[3];
rz(-1.7568622) q[3];
sx q[3];
rz(-0.95209852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0939193) q[2];
sx q[2];
rz(-0.27457044) q[2];
sx q[2];
rz(-2.3982415) q[2];
rz(-0.36176935) q[3];
sx q[3];
rz(-1.0310562) q[3];
sx q[3];
rz(2.7638392) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7985228) q[0];
sx q[0];
rz(-1.9462076) q[0];
sx q[0];
rz(-2.5407347) q[0];
rz(2.9877904) q[1];
sx q[1];
rz(-1.3187131) q[1];
sx q[1];
rz(-2.9153894) q[1];
rz(1.2633688) q[2];
sx q[2];
rz(-1.6976327) q[2];
sx q[2];
rz(-0.007625811) q[2];
rz(-1.947851) q[3];
sx q[3];
rz(-0.98529639) q[3];
sx q[3];
rz(-2.1724971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
