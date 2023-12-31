OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(2.4490693) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0723244) q[0];
sx q[0];
rz(-0.82075483) q[0];
sx q[0];
rz(-2.6988091) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5897883) q[2];
sx q[2];
rz(-1.3349512) q[2];
sx q[2];
rz(2.9896196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9618122) q[1];
sx q[1];
rz(-1.4084067) q[1];
sx q[1];
rz(-0.23602545) q[1];
x q[2];
rz(-2.5296506) q[3];
sx q[3];
rz(-2.3903333) q[3];
sx q[3];
rz(-2.2236168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(-1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40760621) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.227238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5204029) q[0];
sx q[0];
rz(-1.495201) q[0];
sx q[0];
rz(-1.441342) q[0];
rz(-pi) q[1];
rz(1.111755) q[2];
sx q[2];
rz(-2.0437721) q[2];
sx q[2];
rz(2.3943261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8933567) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-0.30171079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.482588) q[3];
sx q[3];
rz(-0.25203029) q[3];
sx q[3];
rz(-0.27740955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.041302117) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(0.89865249) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(0.333289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3086739) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(0.69899107) q[0];
rz(1.0684418) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(1.7114491) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8771613) q[1];
sx q[1];
rz(-0.96211551) q[1];
sx q[1];
rz(-0.18552893) q[1];
rz(-pi) q[2];
rz(-0.40368747) q[3];
sx q[3];
rz(-2.0848158) q[3];
sx q[3];
rz(-0.46883632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44160098) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.9365786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22574632) q[0];
sx q[0];
rz(-1.5374743) q[0];
sx q[0];
rz(-0.87945625) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92653805) q[2];
sx q[2];
rz(-1.7346003) q[2];
sx q[2];
rz(-1.0545727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9884831) q[1];
sx q[1];
rz(-1.8433237) q[1];
sx q[1];
rz(-2.2793798) q[1];
x q[2];
rz(-3.0699176) q[3];
sx q[3];
rz(-2.274401) q[3];
sx q[3];
rz(3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(2.7704346) q[2];
rz(1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(-2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054759653) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(2.5865119) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9585882) q[0];
sx q[0];
rz(-0.31561139) q[0];
sx q[0];
rz(1.6118227) q[0];
rz(-pi) q[1];
rz(1.1558756) q[2];
sx q[2];
rz(-1.6332111) q[2];
sx q[2];
rz(-1.8780564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2808387) q[1];
sx q[1];
rz(-0.44563952) q[1];
sx q[1];
rz(-1.526236) q[1];
x q[2];
rz(-3.0606907) q[3];
sx q[3];
rz(-1.5187851) q[3];
sx q[3];
rz(0.73976433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(-2.1991918) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-0.57975769) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.7738713) q[0];
sx q[0];
rz(-2.2905486) q[0];
rz(-2.8069205) q[2];
sx q[2];
rz(-0.13813189) q[2];
sx q[2];
rz(1.7931256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.786799) q[1];
sx q[1];
rz(-1.2878294) q[1];
sx q[1];
rz(2.016504) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88605373) q[3];
sx q[3];
rz(-1.7985117) q[3];
sx q[3];
rz(2.7453604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5876028) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(2.457298) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(2.6228242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2001901) q[0];
sx q[0];
rz(-1.7046283) q[0];
sx q[0];
rz(-0.83394136) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77063009) q[2];
sx q[2];
rz(-1.7297041) q[2];
sx q[2];
rz(-0.7030013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59473945) q[1];
sx q[1];
rz(-2.2097128) q[1];
sx q[1];
rz(0.98313318) q[1];
rz(-pi) q[2];
rz(-1.1701059) q[3];
sx q[3];
rz(-1.5466994) q[3];
sx q[3];
rz(1.3093684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1812487) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(1.7425591) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-0.74434892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8101013) q[0];
sx q[0];
rz(-1.6252675) q[0];
sx q[0];
rz(0.14793747) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0853945) q[2];
sx q[2];
rz(-1.7040164) q[2];
sx q[2];
rz(-1.0366057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5163755) q[1];
sx q[1];
rz(-0.55621925) q[1];
sx q[1];
rz(0.43362995) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5472502) q[3];
sx q[3];
rz(-2.0252725) q[3];
sx q[3];
rz(2.8794895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(2.5320833) q[2];
rz(2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(1.6375861) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(2.3666568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0477714) q[0];
sx q[0];
rz(-1.9592966) q[0];
sx q[0];
rz(2.2872778) q[0];
rz(-pi) q[1];
rz(-2.9853285) q[2];
sx q[2];
rz(-1.1322349) q[2];
sx q[2];
rz(-1.4807448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(-2.7369569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1433467) q[3];
sx q[3];
rz(-0.25902723) q[3];
sx q[3];
rz(1.6365901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(-2.9837218) q[2];
rz(-1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(2.0956031) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5770618) q[0];
sx q[0];
rz(-0.37527592) q[0];
sx q[0];
rz(-0.29348404) q[0];
x q[1];
rz(-1.5418391) q[2];
sx q[2];
rz(-2.1615897) q[2];
sx q[2];
rz(0.40320413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18098772) q[1];
sx q[1];
rz(-2.1121896) q[1];
sx q[1];
rz(-2.7156746) q[1];
rz(0.06185992) q[3];
sx q[3];
rz(-2.7354771) q[3];
sx q[3];
rz(1.5861685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(-2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8626704) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(0.40907787) q[2];
sx q[2];
rz(-1.2862051) q[2];
sx q[2];
rz(-0.4551879) q[2];
rz(-1.5076751) q[3];
sx q[3];
rz(-2.1199385) q[3];
sx q[3];
rz(-1.4169823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
