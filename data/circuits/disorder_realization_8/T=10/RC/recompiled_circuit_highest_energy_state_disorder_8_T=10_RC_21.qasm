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
rz(2.3249792) q[0];
sx q[0];
rz(-2.6007574) q[0];
sx q[0];
rz(1.9833366) q[0];
rz(6.3732014) q[1];
sx q[1];
rz(6.7572588) q[1];
sx q[1];
rz(13.401539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0220563) q[0];
sx q[0];
rz(-1.2588663) q[0];
sx q[0];
rz(-2.3989912) q[0];
x q[1];
rz(-0.41423256) q[2];
sx q[2];
rz(-2.4576839) q[2];
sx q[2];
rz(-1.6773083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0660634) q[1];
sx q[1];
rz(-1.8154241) q[1];
sx q[1];
rz(-2.7615158) q[1];
x q[2];
rz(0.6583383) q[3];
sx q[3];
rz(-1.1000203) q[3];
sx q[3];
rz(-1.2623163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(-0.33207616) q[2];
rz(-0.24886985) q[3];
sx q[3];
rz(-1.1955465) q[3];
sx q[3];
rz(-0.88768774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87149549) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(-2.6089597) q[0];
rz(-0.10781413) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(2.8210988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1030514) q[0];
sx q[0];
rz(-1.5070931) q[0];
sx q[0];
rz(0.91821435) q[0];
x q[1];
rz(-1.8024496) q[2];
sx q[2];
rz(-1.7391053) q[2];
sx q[2];
rz(-1.2141808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1793666) q[1];
sx q[1];
rz(-1.7229862) q[1];
sx q[1];
rz(0.16140143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32933195) q[3];
sx q[3];
rz(-1.7747702) q[3];
sx q[3];
rz(-0.70131174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33114854) q[2];
sx q[2];
rz(-2.836477) q[2];
sx q[2];
rz(1.6395052) q[2];
rz(0.33809996) q[3];
sx q[3];
rz(-2.2564042) q[3];
sx q[3];
rz(0.52687183) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6239887) q[0];
sx q[0];
rz(-1.232134) q[0];
sx q[0];
rz(2.7841618) q[0];
rz(2.7492211) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(1.9270814) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8100963) q[0];
sx q[0];
rz(-1.4292681) q[0];
sx q[0];
rz(3.1241199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7262801) q[2];
sx q[2];
rz(-2.1124055) q[2];
sx q[2];
rz(-0.77237788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.743108) q[1];
sx q[1];
rz(-1.9581984) q[1];
sx q[1];
rz(0.0031664567) q[1];
rz(-pi) q[2];
rz(3.103996) q[3];
sx q[3];
rz(-2.1581833) q[3];
sx q[3];
rz(1.2560578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7963205) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(-2.7447682) q[2];
rz(1.2567629) q[3];
sx q[3];
rz(-1.4182914) q[3];
sx q[3];
rz(0.61029148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20763718) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(-1.9130094) q[0];
rz(1.928891) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(2.520715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1771496) q[0];
sx q[0];
rz(-1.4040134) q[0];
sx q[0];
rz(2.0509999) q[0];
x q[1];
rz(0.16571705) q[2];
sx q[2];
rz(-1.4172557) q[2];
sx q[2];
rz(2.4927947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8142227) q[1];
sx q[1];
rz(-1.3124221) q[1];
sx q[1];
rz(1.3034921) q[1];
x q[2];
rz(2.8534783) q[3];
sx q[3];
rz(-2.7441437) q[3];
sx q[3];
rz(1.1045052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.88317251) q[2];
sx q[2];
rz(-2.0032538) q[2];
sx q[2];
rz(0.15667285) q[2];
rz(0.91642085) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(0.55287439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857392) q[0];
sx q[0];
rz(-1.9802977) q[0];
sx q[0];
rz(-2.7440985) q[0];
rz(3.1187348) q[1];
sx q[1];
rz(-2.6409179) q[1];
sx q[1];
rz(2.4668677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0521654) q[0];
sx q[0];
rz(-2.8244655) q[0];
sx q[0];
rz(2.2990312) q[0];
x q[1];
rz(-0.65058913) q[2];
sx q[2];
rz(-2.9862767) q[2];
sx q[2];
rz(-2.3273205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6276363) q[1];
sx q[1];
rz(-1.4003203) q[1];
sx q[1];
rz(1.4161311) q[1];
x q[2];
rz(0.85002331) q[3];
sx q[3];
rz(-2.7394419) q[3];
sx q[3];
rz(-1.7274203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5270093) q[2];
sx q[2];
rz(-2.535227) q[2];
sx q[2];
rz(-2.442339) q[2];
rz(1.7953385) q[3];
sx q[3];
rz(-0.56763879) q[3];
sx q[3];
rz(2.4301372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4166477) q[0];
sx q[0];
rz(-0.41794932) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(-2.1669972) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(1.4385361) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33773003) q[0];
sx q[0];
rz(-0.32050214) q[0];
sx q[0];
rz(0.33945531) q[0];
x q[1];
rz(-0.010377093) q[2];
sx q[2];
rz(-1.8266018) q[2];
sx q[2];
rz(0.95358816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8827445) q[1];
sx q[1];
rz(-2.2387894) q[1];
sx q[1];
rz(-0.29752389) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4532279) q[3];
sx q[3];
rz(-2.6010644) q[3];
sx q[3];
rz(1.5557529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4633816) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(1.3055275) q[2];
rz(1.8259004) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(-2.2731884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0618133) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(1.4965936) q[0];
rz(-2.1553701) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(-0.49759069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4301028) q[0];
sx q[0];
rz(-0.27207366) q[0];
sx q[0];
rz(-1.2077232) q[0];
x q[1];
rz(2.1946218) q[2];
sx q[2];
rz(-2.0227891) q[2];
sx q[2];
rz(0.78540451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10863241) q[1];
sx q[1];
rz(-1.4147621) q[1];
sx q[1];
rz(1.969127) q[1];
x q[2];
rz(-0.31106205) q[3];
sx q[3];
rz(-2.8433244) q[3];
sx q[3];
rz(0.48613557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30279532) q[2];
sx q[2];
rz(-1.5747728) q[2];
sx q[2];
rz(-2.8509169) q[2];
rz(2.931328) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(2.1073585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57598376) q[0];
sx q[0];
rz(-1.1701595) q[0];
sx q[0];
rz(-1.1822816) q[0];
rz(0.15444175) q[1];
sx q[1];
rz(-1.7223822) q[1];
sx q[1];
rz(2.2363037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.824914) q[0];
sx q[0];
rz(-1.2980882) q[0];
sx q[0];
rz(3.1290359) q[0];
rz(-pi) q[1];
rz(-1.8087093) q[2];
sx q[2];
rz(-2.015997) q[2];
sx q[2];
rz(-0.92724909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0876079) q[1];
sx q[1];
rz(-1.7167257) q[1];
sx q[1];
rz(-0.75775679) q[1];
x q[2];
rz(-1.9117113) q[3];
sx q[3];
rz(-1.3671759) q[3];
sx q[3];
rz(-0.019817185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0901383) q[2];
sx q[2];
rz(-1.5807512) q[2];
sx q[2];
rz(-0.95019379) q[2];
rz(3.0692302) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(-2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224834) q[0];
sx q[0];
rz(-0.95471946) q[0];
sx q[0];
rz(-1.0954274) q[0];
rz(-1.0911881) q[1];
sx q[1];
rz(-2.2406816) q[1];
sx q[1];
rz(-0.99517623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54992095) q[0];
sx q[0];
rz(-0.66837817) q[0];
sx q[0];
rz(0.10381283) q[0];
x q[1];
rz(-1.9068933) q[2];
sx q[2];
rz(-1.8310412) q[2];
sx q[2];
rz(-1.1059831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1171451) q[1];
sx q[1];
rz(-0.096344171) q[1];
sx q[1];
rz(-0.49343719) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82251151) q[3];
sx q[3];
rz(-0.26897463) q[3];
sx q[3];
rz(1.176187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5114078) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(2.4616145) q[2];
rz(-1.6857111) q[3];
sx q[3];
rz(-2.2472491) q[3];
sx q[3];
rz(1.9623914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7162914) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(-2.5262078) q[0];
rz(-1.3353434) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(-1.7937484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772661) q[0];
sx q[0];
rz(-1.5249671) q[0];
sx q[0];
rz(0.35542506) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87839076) q[2];
sx q[2];
rz(-2.2579402) q[2];
sx q[2];
rz(2.2855482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9453958) q[1];
sx q[1];
rz(-2.7295503) q[1];
sx q[1];
rz(-2.7315188) q[1];
x q[2];
rz(0.12935454) q[3];
sx q[3];
rz(-0.760303) q[3];
sx q[3];
rz(3.1227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24796692) q[2];
sx q[2];
rz(-1.3004356) q[2];
sx q[2];
rz(0.51353961) q[2];
rz(-0.14704554) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(-1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87179398) q[0];
sx q[0];
rz(-0.88903058) q[0];
sx q[0];
rz(-0.24118184) q[0];
rz(-1.3006032) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(-1.2215963) q[2];
sx q[2];
rz(-0.39633718) q[2];
sx q[2];
rz(0.63971165) q[2];
rz(2.5592531) q[3];
sx q[3];
rz(-2.4262541) q[3];
sx q[3];
rz(-2.2179009) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
