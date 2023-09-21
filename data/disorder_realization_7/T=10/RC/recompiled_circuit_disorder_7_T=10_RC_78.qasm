OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(-0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39591046) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(-1.0049055) q[0];
rz(-2.7896499) q[2];
sx q[2];
rz(-1.8004187) q[2];
sx q[2];
rz(-0.72598347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1498897) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(-1.0690772) q[1];
rz(2.9653373) q[3];
sx q[3];
rz(-1.8415383) q[3];
sx q[3];
rz(-2.3034629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-0.75479341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3283008) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(2.6585447) q[0];
x q[1];
rz(-2.7254843) q[2];
sx q[2];
rz(-2.5118718) q[2];
sx q[2];
rz(-2.1925418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.098298479) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(1.0272825) q[1];
rz(-pi) q[2];
rz(-1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(-2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(-0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.67702883) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-0.57998747) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45535116) q[0];
sx q[0];
rz(-0.65128122) q[0];
sx q[0];
rz(-0.63182802) q[0];
rz(-pi) q[1];
rz(1.9252831) q[2];
sx q[2];
rz(-1.4282994) q[2];
sx q[2];
rz(0.49152495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12990002) q[1];
sx q[1];
rz(-1.7923647) q[1];
sx q[1];
rz(1.8151912) q[1];
rz(-0.97614395) q[3];
sx q[3];
rz(-0.38446063) q[3];
sx q[3];
rz(-0.36851766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(0.48809537) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8433544) q[0];
sx q[0];
rz(-2.6697864) q[0];
sx q[0];
rz(-1.3993553) q[0];
x q[1];
rz(-0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(-0.70300245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6285703) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(3.0111074) q[1];
rz(-pi) q[2];
rz(2.5591764) q[3];
sx q[3];
rz(-1.9786668) q[3];
sx q[3];
rz(1.2553314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-0.50746894) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-2.7713293) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21393468) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(-1.8646122) q[0];
rz(0.032012149) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(-0.33547685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0177808) q[1];
sx q[1];
rz(-1.8291626) q[1];
sx q[1];
rz(-1.3892795) q[1];
rz(0.27627857) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.7477759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0575384) q[0];
sx q[0];
rz(-1.49465) q[0];
sx q[0];
rz(2.8094859) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1962842) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(1.9749157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5382232) q[1];
sx q[1];
rz(-1.9586316) q[1];
sx q[1];
rz(0.28660725) q[1];
rz(-pi) q[2];
rz(-2.0270258) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(1.2896468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-0.21025118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1649095) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(-3.0747736) q[0];
x q[1];
rz(-2.7619751) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0381895) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(1.7722539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5641714) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644972) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(-2.68481) q[0];
rz(-pi) q[1];
rz(1.1114242) q[2];
sx q[2];
rz(-2.9615235) q[2];
sx q[2];
rz(-1.1262696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8969438) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(1.1036554) q[1];
x q[2];
rz(-0.85985698) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(0.45339282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(-0.92064944) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(0.89909536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16222787) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(-2.2158932) q[0];
x q[1];
rz(2.8092438) q[2];
sx q[2];
rz(-1.0370478) q[2];
sx q[2];
rz(0.74408434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60914674) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(1.902641) q[1];
rz(-pi) q[2];
rz(-1.9409688) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(1.9716594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(2.7446279) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871646) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(0.88231477) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(0.099909401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46170235) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(0.093613503) q[1];
x q[2];
rz(-1.0730562) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(-2.4772714) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(1.9962911) q[3];
sx q[3];
rz(-2.8588061) q[3];
sx q[3];
rz(-2.2929946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];