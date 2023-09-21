OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96903893) q[0];
sx q[0];
rz(-1.2533029) q[0];
sx q[0];
rz(-2.88455) q[0];
x q[1];
rz(-0.89262427) q[2];
sx q[2];
rz(-1.2783588) q[2];
sx q[2];
rz(-2.6543648) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.675128) q[1];
sx q[1];
rz(-0.29106859) q[1];
sx q[1];
rz(0.18578271) q[1];
x q[2];
rz(2.0831574) q[3];
sx q[3];
rz(-1.5844987) q[3];
sx q[3];
rz(2.0126359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.1738698) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(0.064963438) q[0];
rz(-0.57463542) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949462) q[0];
sx q[0];
rz(-1.6992237) q[0];
sx q[0];
rz(-0.24982474) q[0];
x q[1];
rz(-2.9827036) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(-1.5140669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17774432) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(2.6049155) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4317125) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(-2.3185454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(0.72845355) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-0.089086108) q[0];
sx q[0];
rz(0.41390093) q[0];
rz(0.74480199) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(-1.1643861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6779855) q[1];
sx q[1];
rz(-1.7079759) q[1];
sx q[1];
rz(2.3198818) q[1];
rz(1.4025027) q[3];
sx q[3];
rz(-0.28545359) q[3];
sx q[3];
rz(1.9608378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(-0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1082552) q[0];
sx q[0];
rz(-2.4318683) q[0];
sx q[0];
rz(1.8002585) q[0];
rz(-pi) q[1];
rz(-2.6631782) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(2.0699376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66407953) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(-1.627229) q[1];
rz(-2.187192) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(-2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9527234) q[0];
sx q[0];
rz(-1.6048604) q[0];
sx q[0];
rz(2.6012095) q[0];
rz(0.70154538) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(1.0813576) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69265428) q[1];
sx q[1];
rz(-0.9874953) q[1];
sx q[1];
rz(-1.7765691) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4928717) q[3];
sx q[3];
rz(-1.3789163) q[3];
sx q[3];
rz(-0.81024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(0.38875368) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42449441) q[0];
sx q[0];
rz(-1.6569123) q[0];
sx q[0];
rz(-0.17688246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7331946) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(-2.1085395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5387419) q[1];
sx q[1];
rz(-0.83487836) q[1];
sx q[1];
rz(-1.6470626) q[1];
rz(0.74495875) q[3];
sx q[3];
rz(-0.31674851) q[3];
sx q[3];
rz(1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(-2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(2.0986309) q[0];
rz(0.46328059) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8585513) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(1.5714684) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3473862) q[2];
sx q[2];
rz(-1.4457448) q[2];
sx q[2];
rz(1.320968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6003905) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(-2.3669835) q[1];
rz(-pi) q[2];
rz(-0.058810874) q[3];
sx q[3];
rz(-0.33274129) q[3];
sx q[3];
rz(-0.22418338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7730007) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(-2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8038717) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(-1.7659811) q[0];
x q[1];
rz(2.9647486) q[2];
sx q[2];
rz(-1.7272005) q[2];
sx q[2];
rz(-0.96207372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.265043) q[1];
sx q[1];
rz(-1.3339086) q[1];
sx q[1];
rz(-2.151728) q[1];
rz(-pi) q[2];
rz(2.3341228) q[3];
sx q[3];
rz(-1.8261357) q[3];
sx q[3];
rz(-2.1047999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8895421) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77712599) q[0];
sx q[0];
rz(-1.5317894) q[0];
sx q[0];
rz(-1.5968621) q[0];
rz(-pi) q[1];
rz(0.75195306) q[2];
sx q[2];
rz(-2.8402036) q[2];
sx q[2];
rz(1.4434659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5962332) q[1];
sx q[1];
rz(-1.4701478) q[1];
sx q[1];
rz(0.049493162) q[1];
rz(1.7248475) q[3];
sx q[3];
rz(-2.027958) q[3];
sx q[3];
rz(2.9726213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41708502) q[2];
sx q[2];
rz(-2.8736726) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(0.80741185) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080169454) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(0.28082401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.330864) q[0];
sx q[0];
rz(-0.64635902) q[0];
sx q[0];
rz(2.3848563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4091987) q[2];
sx q[2];
rz(-0.28389441) q[2];
sx q[2];
rz(0.64901272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7933958) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(-0.036199526) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9452768) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3165555) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(2.5785057) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6476718) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-2.3360844) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(0.35691805) q[3];
sx q[3];
rz(-1.9580943) q[3];
sx q[3];
rz(3.0860268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
