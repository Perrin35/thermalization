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
rz(-2.3075624) q[0];
sx q[0];
rz(-0.40606719) q[0];
sx q[0];
rz(1.1817482) q[0];
rz(-1.4015247) q[1];
sx q[1];
rz(-2.6670246) q[1];
sx q[1];
rz(-0.13303703) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849738) q[0];
sx q[0];
rz(-2.4388038) q[0];
sx q[0];
rz(1.982559) q[0];
rz(2.6288787) q[2];
sx q[2];
rz(-1.9641335) q[2];
sx q[2];
rz(2.8296997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.592748) q[1];
sx q[1];
rz(-1.8858375) q[1];
sx q[1];
rz(2.9058019) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99869149) q[3];
sx q[3];
rz(-1.8810728) q[3];
sx q[3];
rz(0.43768044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0066321) q[2];
sx q[2];
rz(-2.1217608) q[2];
sx q[2];
rz(0.68438619) q[2];
rz(2.0002174) q[3];
sx q[3];
rz(-0.26641521) q[3];
sx q[3];
rz(-0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.77966493) q[0];
sx q[0];
rz(-0.68988887) q[0];
sx q[0];
rz(0.46098125) q[0];
rz(-2.0224679) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(0.31203312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658995) q[0];
sx q[0];
rz(-1.5830432) q[0];
sx q[0];
rz(-0.027994305) q[0];
x q[1];
rz(3.0655601) q[2];
sx q[2];
rz(-2.3302148) q[2];
sx q[2];
rz(-1.8234314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4496135) q[1];
sx q[1];
rz(-2.9658554) q[1];
sx q[1];
rz(1.537561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9309405) q[3];
sx q[3];
rz(-0.68406287) q[3];
sx q[3];
rz(-1.6863914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1877039) q[2];
sx q[2];
rz(-1.1053332) q[2];
sx q[2];
rz(2.7375431) q[2];
rz(-2.7814501) q[3];
sx q[3];
rz(-1.4934243) q[3];
sx q[3];
rz(0.68296105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78524154) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(0.32401618) q[0];
rz(-2.0815381) q[1];
sx q[1];
rz(-0.29393229) q[1];
sx q[1];
rz(-0.70045984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16228904) q[0];
sx q[0];
rz(-1.7320181) q[0];
sx q[0];
rz(0.17031807) q[0];
rz(-pi) q[1];
rz(-1.9159928) q[2];
sx q[2];
rz(-1.7734062) q[2];
sx q[2];
rz(-2.4576296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4844651) q[1];
sx q[1];
rz(-1.8052237) q[1];
sx q[1];
rz(0.77204324) q[1];
rz(1.9869119) q[3];
sx q[3];
rz(-0.45588247) q[3];
sx q[3];
rz(-2.1385764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5141653) q[2];
sx q[2];
rz(-0.24719396) q[2];
sx q[2];
rz(-0.79666454) q[2];
rz(1.1967777) q[3];
sx q[3];
rz(-1.5408885) q[3];
sx q[3];
rz(2.5483907) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45348039) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(1.3230811) q[0];
rz(-2.6755013) q[1];
sx q[1];
rz(-0.90704647) q[1];
sx q[1];
rz(1.9922493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02657613) q[0];
sx q[0];
rz(-1.1372024) q[0];
sx q[0];
rz(-2.114891) q[0];
rz(-pi) q[1];
rz(1.0428814) q[2];
sx q[2];
rz(-2.4529874) q[2];
sx q[2];
rz(1.3618038) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6862029) q[1];
sx q[1];
rz(-1.7727994) q[1];
sx q[1];
rz(1.6167902) q[1];
x q[2];
rz(-2.9907865) q[3];
sx q[3];
rz(-1.9016087) q[3];
sx q[3];
rz(-0.88265935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.484802) q[2];
sx q[2];
rz(-2.4718086) q[2];
sx q[2];
rz(2.316324) q[2];
rz(1.8615865) q[3];
sx q[3];
rz(-1.620159) q[3];
sx q[3];
rz(2.4204204) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44329062) q[0];
sx q[0];
rz(-1.2606786) q[0];
sx q[0];
rz(-1.3662421) q[0];
rz(1.0900991) q[1];
sx q[1];
rz(-1.5049728) q[1];
sx q[1];
rz(1.8025788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3109717) q[0];
sx q[0];
rz(-1.3960103) q[0];
sx q[0];
rz(0.66077787) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62241222) q[2];
sx q[2];
rz(-0.96660766) q[2];
sx q[2];
rz(-2.6780918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73698649) q[1];
sx q[1];
rz(-1.6885719) q[1];
sx q[1];
rz(2.9280258) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2286148) q[3];
sx q[3];
rz(-1.2552731) q[3];
sx q[3];
rz(2.4339937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94685405) q[2];
sx q[2];
rz(-1.4118492) q[2];
sx q[2];
rz(0.33352387) q[2];
rz(-2.5528095) q[3];
sx q[3];
rz(-2.6684561) q[3];
sx q[3];
rz(0.80140448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3280585) q[0];
sx q[0];
rz(-1.7840339) q[0];
sx q[0];
rz(-1.300977) q[0];
rz(-1.0282372) q[1];
sx q[1];
rz(-0.49003777) q[1];
sx q[1];
rz(0.44233826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8632354) q[0];
sx q[0];
rz(-1.8148737) q[0];
sx q[0];
rz(-2.6574479) q[0];
rz(-pi) q[1];
rz(-2.6643941) q[2];
sx q[2];
rz(-1.0893359) q[2];
sx q[2];
rz(-1.19966) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79696733) q[1];
sx q[1];
rz(-1.4055924) q[1];
sx q[1];
rz(2.7112673) q[1];
rz(-0.30486561) q[3];
sx q[3];
rz(-0.33434048) q[3];
sx q[3];
rz(1.6337194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(2.6698574) q[2];
rz(-1.4419904) q[3];
sx q[3];
rz(-0.56157464) q[3];
sx q[3];
rz(-0.26427463) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1137375) q[0];
sx q[0];
rz(-1.6738482) q[0];
sx q[0];
rz(2.9748528) q[0];
rz(-0.79404229) q[1];
sx q[1];
rz(-1.9414976) q[1];
sx q[1];
rz(3.0418495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3152811) q[0];
sx q[0];
rz(-0.65412414) q[0];
sx q[0];
rz(-2.099206) q[0];
x q[1];
rz(0.95364596) q[2];
sx q[2];
rz(-1.8815478) q[2];
sx q[2];
rz(-2.1073532) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43370285) q[1];
sx q[1];
rz(-1.4447482) q[1];
sx q[1];
rz(-2.5688102) q[1];
rz(-pi) q[2];
rz(2.7230303) q[3];
sx q[3];
rz(-0.95223126) q[3];
sx q[3];
rz(1.5833861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2793067) q[2];
sx q[2];
rz(-1.0655094) q[2];
sx q[2];
rz(0.12969895) q[2];
rz(-1.8779523) q[3];
sx q[3];
rz(-1.2010776) q[3];
sx q[3];
rz(2.9960846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.170914) q[0];
sx q[0];
rz(-2.2190974) q[0];
sx q[0];
rz(0.42651919) q[0];
rz(-2.049394) q[1];
sx q[1];
rz(-1.1323606) q[1];
sx q[1];
rz(-1.0027286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.730091) q[0];
sx q[0];
rz(-2.1637332) q[0];
sx q[0];
rz(0.31296313) q[0];
rz(1.6636983) q[2];
sx q[2];
rz(-2.291541) q[2];
sx q[2];
rz(-2.0630974) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7510609) q[1];
sx q[1];
rz(-0.56274429) q[1];
sx q[1];
rz(-2.2319002) q[1];
rz(1.1137943) q[3];
sx q[3];
rz(-0.49906956) q[3];
sx q[3];
rz(2.6036865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87454522) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(-3.0181273) q[2];
rz(-2.5115783) q[3];
sx q[3];
rz(-1.8450626) q[3];
sx q[3];
rz(-2.4251078) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89710871) q[0];
sx q[0];
rz(-0.69863313) q[0];
sx q[0];
rz(-2.3022292) q[0];
rz(0.15946236) q[1];
sx q[1];
rz(-0.99226743) q[1];
sx q[1];
rz(-0.92963591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.000132) q[0];
sx q[0];
rz(-0.85171284) q[0];
sx q[0];
rz(1.9107463) q[0];
x q[1];
rz(0.18126296) q[2];
sx q[2];
rz(-2.5428704) q[2];
sx q[2];
rz(1.8229654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8690455) q[1];
sx q[1];
rz(-2.2490977) q[1];
sx q[1];
rz(-0.1609623) q[1];
x q[2];
rz(1.8760909) q[3];
sx q[3];
rz(-2.4809894) q[3];
sx q[3];
rz(1.5098287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9440426) q[2];
sx q[2];
rz(-1.1669179) q[2];
sx q[2];
rz(-2.3738532) q[2];
rz(0.36457101) q[3];
sx q[3];
rz(-1.3837827) q[3];
sx q[3];
rz(-3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.1345054) q[0];
sx q[0];
rz(-0.61408478) q[0];
sx q[0];
rz(2.2579204) q[0];
rz(-0.20835749) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(1.8066033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7430089) q[0];
sx q[0];
rz(-1.5129876) q[0];
sx q[0];
rz(-3.0763118) q[0];
x q[1];
rz(-1.7572549) q[2];
sx q[2];
rz(-2.0907195) q[2];
sx q[2];
rz(-2.7911719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22182759) q[1];
sx q[1];
rz(-1.034267) q[1];
sx q[1];
rz(2.6229658) q[1];
rz(-pi) q[2];
rz(2.4766165) q[3];
sx q[3];
rz(-2.1998365) q[3];
sx q[3];
rz(-1.4339989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9779382) q[2];
sx q[2];
rz(-1.6341219) q[2];
sx q[2];
rz(-3.0955691) q[2];
rz(-0.16511354) q[3];
sx q[3];
rz(-2.8486227) q[3];
sx q[3];
rz(1.9273812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56275) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(-2.9764755) q[1];
sx q[1];
rz(-1.6830291) q[1];
sx q[1];
rz(2.8370636) q[1];
rz(-1.8971741) q[2];
sx q[2];
rz(-1.1481297) q[2];
sx q[2];
rz(1.4884244) q[2];
rz(0.48043171) q[3];
sx q[3];
rz(-2.0365002) q[3];
sx q[3];
rz(-1.0040803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
