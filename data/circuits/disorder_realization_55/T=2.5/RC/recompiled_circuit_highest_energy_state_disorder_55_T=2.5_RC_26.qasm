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
rz(-1.363938) q[0];
sx q[0];
rz(3.5600297) q[0];
sx q[0];
rz(10.726396) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(-1.6956704) q[1];
sx q[1];
rz(3.1248098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12110773) q[0];
sx q[0];
rz(-1.2099625) q[0];
sx q[0];
rz(-3.0617972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4799825) q[2];
sx q[2];
rz(-0.1802643) q[2];
sx q[2];
rz(1.6423051) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1396347) q[1];
sx q[1];
rz(-0.0070925709) q[1];
sx q[1];
rz(2.3435739) q[1];
x q[2];
rz(-2.7990325) q[3];
sx q[3];
rz(-0.16157074) q[3];
sx q[3];
rz(-1.2662966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(1.5793229) q[2];
rz(0.2114547) q[3];
sx q[3];
rz(-3.1410757) q[3];
sx q[3];
rz(-0.16099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(-1.3268693) q[0];
rz(-2.5621085) q[1];
sx q[1];
rz(-3.1377628) q[1];
sx q[1];
rz(-0.63900596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20637437) q[0];
sx q[0];
rz(-2.2274096) q[0];
sx q[0];
rz(0.73120631) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0206356) q[2];
sx q[2];
rz(-1.5874344) q[2];
sx q[2];
rz(1.5543092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9663869) q[1];
sx q[1];
rz(-1.5537973) q[1];
sx q[1];
rz(-1.5823889) q[1];
rz(0.08777471) q[3];
sx q[3];
rz(-2.3592161) q[3];
sx q[3];
rz(0.45004547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-0.13627626) q[2];
sx q[2];
rz(-1.603568) q[2];
rz(-1.5806574) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.2494217) q[0];
sx q[0];
rz(-2.6264661) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(-2.4341266) q[1];
sx q[1];
rz(-3.1222157) q[1];
sx q[1];
rz(-2.0170508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19605787) q[0];
sx q[0];
rz(-1.3302517) q[0];
sx q[0];
rz(1.3105767) q[0];
x q[1];
rz(-3.1156024) q[2];
sx q[2];
rz(-1.6857237) q[2];
sx q[2];
rz(2.9851802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4545817) q[1];
sx q[1];
rz(-3.0774766) q[1];
sx q[1];
rz(2.9518576) q[1];
rz(-pi) q[2];
rz(2.0865404) q[3];
sx q[3];
rz(-1.0760783) q[3];
sx q[3];
rz(2.5124541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4418929) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(-3.0921248) q[2];
rz(2.5345645) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98168755) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(3.1244151) q[0];
rz(-2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5944098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3348986) q[0];
sx q[0];
rz(-2.3133754) q[0];
sx q[0];
rz(-0.81636274) q[0];
rz(-pi) q[1];
rz(1.3457005) q[2];
sx q[2];
rz(-2.101311) q[2];
sx q[2];
rz(0.38166416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7878927) q[1];
sx q[1];
rz(-3.0038634) q[1];
sx q[1];
rz(2.6929123) q[1];
x q[2];
rz(0.64597102) q[3];
sx q[3];
rz(-1.5785909) q[3];
sx q[3];
rz(0.92168671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8881417) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(-2.4593501) q[2];
rz(-0.055179723) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(1.3050219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76160112) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(2.7165661) q[0];
rz(1.5399326) q[1];
sx q[1];
rz(-2.6585177) q[1];
sx q[1];
rz(-0.81533283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402712) q[0];
sx q[0];
rz(-3.0794332) q[0];
sx q[0];
rz(-2.9654853) q[0];
x q[1];
rz(-1.710395) q[2];
sx q[2];
rz(-0.023471467) q[2];
sx q[2];
rz(2.1144298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22742352) q[1];
sx q[1];
rz(-0.14761758) q[1];
sx q[1];
rz(-2.5116337) q[1];
rz(-2.8132503) q[3];
sx q[3];
rz(-2.7797065) q[3];
sx q[3];
rz(-2.7992424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1347947) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(-1.6689782) q[2];
rz(-2.2115479) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(-0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-0.023094026) q[0];
sx q[0];
rz(-1.4267138) q[0];
rz(-0.73000437) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(-1.1013365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520605) q[0];
sx q[0];
rz(-2.2666078) q[0];
sx q[0];
rz(0.3541644) q[0];
rz(-pi) q[1];
rz(1.8011212) q[2];
sx q[2];
rz(-1.4061951) q[2];
sx q[2];
rz(-0.71704656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8423398) q[1];
sx q[1];
rz(-1.4758759) q[1];
sx q[1];
rz(1.6948944) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6953129) q[3];
sx q[3];
rz(-1.5388371) q[3];
sx q[3];
rz(2.4886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71658984) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(-1.8343743) q[2];
rz(-2.763125) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3017479) q[0];
sx q[0];
rz(-1.8740338) q[0];
sx q[0];
rz(-2.2140455) q[0];
rz(1.357366) q[1];
sx q[1];
rz(-2.3097242) q[1];
sx q[1];
rz(1.5313139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5132234) q[0];
sx q[0];
rz(-1.5484865) q[0];
sx q[0];
rz(1.6091899) q[0];
rz(-pi) q[1];
rz(-0.53801914) q[2];
sx q[2];
rz(-2.6876861) q[2];
sx q[2];
rz(1.5858142) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.130333) q[1];
sx q[1];
rz(-1.5734871) q[1];
sx q[1];
rz(-1.6998768) q[1];
x q[2];
rz(2.0965936) q[3];
sx q[3];
rz(-0.57588314) q[3];
sx q[3];
rz(2.8623842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7808468) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(1.8878262) q[2];
rz(-2.4474261) q[3];
sx q[3];
rz(-2.4024051) q[3];
sx q[3];
rz(2.9085801) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6041782) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(-1.0373212) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(-1.6745837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5418422) q[0];
sx q[0];
rz(-1.5034202) q[0];
sx q[0];
rz(1.3753433) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0612512) q[2];
sx q[2];
rz(-0.27896491) q[2];
sx q[2];
rz(-2.7921048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8716177) q[1];
sx q[1];
rz(-1.570357) q[1];
sx q[1];
rz(-1.5719218) q[1];
x q[2];
rz(2.9503533) q[3];
sx q[3];
rz(-2.6789224) q[3];
sx q[3];
rz(-1.7709875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.011270114) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-2.6210426) q[3];
sx q[3];
rz(-3.136941) q[3];
sx q[3];
rz(-1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540102) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(1.822923) q[0];
rz(-1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(-1.5444548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8108687) q[0];
sx q[0];
rz(-1.4557342) q[0];
sx q[0];
rz(-0.48145357) q[0];
x q[1];
rz(2.804432) q[2];
sx q[2];
rz(-2.4558407) q[2];
sx q[2];
rz(-2.8711988) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.006598) q[1];
sx q[1];
rz(-1.8812351) q[1];
sx q[1];
rz(-2.9562692) q[1];
x q[2];
rz(2.5289152) q[3];
sx q[3];
rz(-0.78734382) q[3];
sx q[3];
rz(1.2932216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80457193) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(0.19680944) q[2];
rz(-1.192441) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(-0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(0.30109677) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(1.1916196) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(-1.5764538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3707664) q[0];
sx q[0];
rz(-1.8642555) q[0];
sx q[0];
rz(1.7459041) q[0];
rz(1.1249816) q[2];
sx q[2];
rz(-2.6458394) q[2];
sx q[2];
rz(-0.021136016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31101878) q[1];
sx q[1];
rz(-1.5717447) q[1];
sx q[1];
rz(-1.5702973) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2433167) q[3];
sx q[3];
rz(-0.18096033) q[3];
sx q[3];
rz(2.3955936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18371753) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(1.6831762) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-0.009549791) q[3];
sx q[3];
rz(-0.20235801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771582) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(1.5741813) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(1.4952343) q[2];
sx q[2];
rz(-1.5757558) q[2];
sx q[2];
rz(-1.4103665) q[2];
rz(-0.50441691) q[3];
sx q[3];
rz(-2.258959) q[3];
sx q[3];
rz(-2.7213982) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
