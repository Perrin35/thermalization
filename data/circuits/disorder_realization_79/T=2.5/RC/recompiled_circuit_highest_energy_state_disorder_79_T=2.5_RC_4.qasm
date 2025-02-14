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
rz(-2.8861217) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(-1.7429587) q[1];
sx q[1];
rz(1.7359098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79069239) q[0];
sx q[0];
rz(-3.1381099) q[0];
sx q[0];
rz(-3.0538959) q[0];
rz(-pi) q[1];
rz(2.954835) q[2];
sx q[2];
rz(-0.52350145) q[2];
sx q[2];
rz(-1.7331327) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5454332) q[1];
sx q[1];
rz(-1.9539297) q[1];
sx q[1];
rz(0.0036762357) q[1];
rz(-0.56804232) q[3];
sx q[3];
rz(-1.5431817) q[3];
sx q[3];
rz(1.9909007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(0.8325038) q[2];
rz(2.3333874) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(0.25850779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.38043624) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(2.9535182) q[0];
rz(3.0694118) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(-1.5123051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.78069) q[0];
sx q[0];
rz(-1.6854428) q[0];
sx q[0];
rz(2.8066382) q[0];
x q[1];
rz(0.030628344) q[2];
sx q[2];
rz(-1.5248486) q[2];
sx q[2];
rz(0.46389889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9058414) q[1];
sx q[1];
rz(-3.0260147) q[1];
sx q[1];
rz(-0.53821941) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80792369) q[3];
sx q[3];
rz(-1.1594541) q[3];
sx q[3];
rz(-1.5613256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(1.2930653) q[2];
rz(0.34717789) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(-0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8186571) q[0];
sx q[0];
rz(-1.8523536) q[0];
sx q[0];
rz(-0.47054189) q[0];
rz(-1.200354) q[1];
sx q[1];
rz(-0.73089868) q[1];
sx q[1];
rz(1.8847195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82913933) q[0];
sx q[0];
rz(-1.0395673) q[0];
sx q[0];
rz(2.7307399) q[0];
x q[1];
rz(0.97312178) q[2];
sx q[2];
rz(-2.9785756) q[2];
sx q[2];
rz(2.4578641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4670406) q[1];
sx q[1];
rz(-1.5562773) q[1];
sx q[1];
rz(0.0013703636) q[1];
x q[2];
rz(-0.29095616) q[3];
sx q[3];
rz(-0.97948631) q[3];
sx q[3];
rz(-1.3475498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5990126) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-2.2194594) q[2];
rz(1.7982091) q[3];
sx q[3];
rz(-2.193439) q[3];
sx q[3];
rz(-1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87329292) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(-0.44019765) q[0];
rz(1.531155) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(-0.24756113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2092106) q[0];
sx q[0];
rz(-0.69057751) q[0];
sx q[0];
rz(0.85777046) q[0];
rz(-pi) q[1];
rz(-2.655022) q[2];
sx q[2];
rz(-2.9503194) q[2];
sx q[2];
rz(-0.74081206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1936613) q[1];
sx q[1];
rz(-0.30821092) q[1];
sx q[1];
rz(1.3256509) q[1];
rz(3.1010755) q[3];
sx q[3];
rz(-2.4926293) q[3];
sx q[3];
rz(-1.9891091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89746785) q[2];
sx q[2];
rz(-2.1110057) q[2];
sx q[2];
rz(-1.6991276) q[2];
rz(-0.20081946) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(-2.0012746) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9852801) q[0];
sx q[0];
rz(-0.15287481) q[0];
sx q[0];
rz(2.5961764) q[0];
rz(0.044005752) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(2.6652179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5060318) q[0];
sx q[0];
rz(-2.8836492) q[0];
sx q[0];
rz(-2.8057465) q[0];
x q[1];
rz(2.7353941) q[2];
sx q[2];
rz(-1.1990237) q[2];
sx q[2];
rz(0.68454725) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4551983) q[1];
sx q[1];
rz(-2.2664811) q[1];
sx q[1];
rz(-0.2376539) q[1];
rz(-pi) q[2];
rz(1.8082321) q[3];
sx q[3];
rz(-0.52132208) q[3];
sx q[3];
rz(-2.962474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62006092) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(-2.7812092) q[2];
rz(2.9195869) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5905269) q[0];
sx q[0];
rz(-1.5532302) q[0];
sx q[0];
rz(-1.5850413) q[0];
rz(-0.12697728) q[1];
sx q[1];
rz(-1.8366837) q[1];
sx q[1];
rz(3.051905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9162101) q[0];
sx q[0];
rz(-1.0400335) q[0];
sx q[0];
rz(0.60471524) q[0];
x q[1];
rz(2.181337) q[2];
sx q[2];
rz(-0.88443236) q[2];
sx q[2];
rz(-2.6853564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6818004) q[1];
sx q[1];
rz(-2.5885317) q[1];
sx q[1];
rz(-2.6917767) q[1];
rz(-pi) q[2];
rz(-1.9305658) q[3];
sx q[3];
rz(-2.1412951) q[3];
sx q[3];
rz(0.24107547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.108532) q[2];
sx q[2];
rz(-0.33691275) q[2];
sx q[2];
rz(-0.25165558) q[2];
rz(-0.039994914) q[3];
sx q[3];
rz(-2.7997041) q[3];
sx q[3];
rz(2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104178) q[0];
sx q[0];
rz(-3.0496821) q[0];
sx q[0];
rz(-0.41496667) q[0];
rz(1.4893432) q[1];
sx q[1];
rz(-0.057561189) q[1];
sx q[1];
rz(-0.32589486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6221507) q[0];
sx q[0];
rz(-1.9339193) q[0];
sx q[0];
rz(2.0456929) q[0];
x q[1];
rz(3.0517111) q[2];
sx q[2];
rz(-2.0873859) q[2];
sx q[2];
rz(-0.84026779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53908097) q[1];
sx q[1];
rz(-1.314114) q[1];
sx q[1];
rz(-0.51694444) q[1];
rz(-pi) q[2];
rz(0.075906673) q[3];
sx q[3];
rz(-1.2682639) q[3];
sx q[3];
rz(2.3940866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(-2.0951159) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(-1.3974765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662812) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(2.6783491) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-0.0015365096) q[1];
sx q[1];
rz(-1.6418246) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508598) q[0];
sx q[0];
rz(-1.5895278) q[0];
sx q[0];
rz(2.0035726) q[0];
x q[1];
rz(1.6258662) q[2];
sx q[2];
rz(-2.0158421) q[2];
sx q[2];
rz(-0.17521706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7565175) q[1];
sx q[1];
rz(-2.9179055) q[1];
sx q[1];
rz(-1.4392412) q[1];
x q[2];
rz(-1.3070694) q[3];
sx q[3];
rz(-0.26104195) q[3];
sx q[3];
rz(1.3021951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2396592) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(1.3182013) q[2];
rz(-3.1384595) q[3];
sx q[3];
rz(-1.9414709) q[3];
sx q[3];
rz(1.9759294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779293) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(0.21358061) q[1];
sx q[1];
rz(-0.029284632) q[1];
sx q[1];
rz(1.9307131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41325135) q[0];
sx q[0];
rz(-1.8759707) q[0];
sx q[0];
rz(0.11082173) q[0];
rz(-pi) q[1];
rz(0.40990746) q[2];
sx q[2];
rz(-0.73824595) q[2];
sx q[2];
rz(-1.8869836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3397202) q[1];
sx q[1];
rz(-2.0959217) q[1];
sx q[1];
rz(0.32466356) q[1];
x q[2];
rz(-1.7603432) q[3];
sx q[3];
rz(-1.726103) q[3];
sx q[3];
rz(-1.3819753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1750298) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(-0.19801298) q[2];
rz(2.7909732) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(1.143379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.38458934) q[0];
sx q[0];
rz(-2.2191255) q[0];
sx q[0];
rz(0.61436999) q[0];
rz(3.0820097) q[1];
sx q[1];
rz(-0.091484286) q[1];
sx q[1];
rz(-1.7170067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2587268) q[0];
sx q[0];
rz(-0.75260163) q[0];
sx q[0];
rz(-1.2105788) q[0];
rz(-pi) q[1];
rz(3.0103127) q[2];
sx q[2];
rz(-1.5688217) q[2];
sx q[2];
rz(-1.2009837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1471794) q[1];
sx q[1];
rz(-0.1754027) q[1];
sx q[1];
rz(1.1145304) q[1];
rz(-pi) q[2];
rz(1.9779102) q[3];
sx q[3];
rz(-1.3729248) q[3];
sx q[3];
rz(2.1382633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0397348) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(0.52379918) q[2];
rz(-1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(-2.6773793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47281784) q[0];
sx q[0];
rz(-1.6090809) q[0];
sx q[0];
rz(1.6627298) q[0];
rz(1.7300023) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(0.3910364) q[2];
sx q[2];
rz(-2.4606065) q[2];
sx q[2];
rz(1.4627964) q[2];
rz(-2.6749776) q[3];
sx q[3];
rz(-1.7468921) q[3];
sx q[3];
rz(0.56260059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
