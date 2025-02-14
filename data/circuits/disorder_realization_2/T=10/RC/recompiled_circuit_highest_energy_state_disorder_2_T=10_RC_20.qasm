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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(2.0896572) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5405226) q[0];
sx q[0];
rz(-0.94694505) q[0];
sx q[0];
rz(-0.076657587) q[0];
rz(-2.0016167) q[2];
sx q[2];
rz(-2.2099021) q[2];
sx q[2];
rz(-0.72102816) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6261648) q[1];
sx q[1];
rz(-2.0671131) q[1];
sx q[1];
rz(-1.4908916) q[1];
rz(-3.0188881) q[3];
sx q[3];
rz(-0.70022455) q[3];
sx q[3];
rz(-0.52470696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-2.0852883) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2317155) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(2.7114482) q[0];
rz(-0.12769708) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(1.4375623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970222) q[0];
sx q[0];
rz(-1.26957) q[0];
sx q[0];
rz(1.3006849) q[0];
rz(-pi) q[1];
rz(-3.0906648) q[2];
sx q[2];
rz(-1.5136592) q[2];
sx q[2];
rz(0.47540755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.592652) q[1];
sx q[1];
rz(-0.90528622) q[1];
sx q[1];
rz(-1.9926461) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7892465) q[3];
sx q[3];
rz(-0.56803136) q[3];
sx q[3];
rz(-1.8527135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(-0.44542584) q[2];
rz(0.40924254) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-0.87944952) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(1.0572761) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(-2.8143299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875558) q[0];
sx q[0];
rz(-1.7374514) q[0];
sx q[0];
rz(-1.0221857) q[0];
x q[1];
rz(-1.0017603) q[2];
sx q[2];
rz(-1.1998917) q[2];
sx q[2];
rz(-2.6628464) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79993805) q[1];
sx q[1];
rz(-0.1529049) q[1];
sx q[1];
rz(1.4351109) q[1];
rz(-0.25423519) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(-1.1439307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(-1.5039911) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(0.98201069) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686491) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.2329996) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208682) q[0];
sx q[0];
rz(-1.3107131) q[0];
sx q[0];
rz(-0.12518945) q[0];
rz(-pi) q[1];
rz(-2.2188051) q[2];
sx q[2];
rz(-1.4013259) q[2];
sx q[2];
rz(-0.10709912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8089227) q[1];
sx q[1];
rz(-0.70007174) q[1];
sx q[1];
rz(-1.1247404) q[1];
rz(-pi) q[2];
rz(1.5855012) q[3];
sx q[3];
rz(-1.3193519) q[3];
sx q[3];
rz(0.78212839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(1.5114463) q[2];
rz(1.1394507) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(-0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.9885709) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(0.26184729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42191045) q[0];
sx q[0];
rz(-1.6043264) q[0];
sx q[0];
rz(-1.6543341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0864429) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(-0.5400368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5763229) q[1];
sx q[1];
rz(-2.5082631) q[1];
sx q[1];
rz(0.38193462) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20528593) q[3];
sx q[3];
rz(-2.3488099) q[3];
sx q[3];
rz(2.1403811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5220945) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(1.7363133) q[2];
rz(0.74553472) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(0.94329992) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(-0.44031269) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(0.17131677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650314) q[0];
sx q[0];
rz(-2.2269571) q[0];
sx q[0];
rz(1.879346) q[0];
rz(-pi) q[1];
rz(-2.3883853) q[2];
sx q[2];
rz(-1.3454226) q[2];
sx q[2];
rz(-1.6373375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79119392) q[1];
sx q[1];
rz(-1.5499252) q[1];
sx q[1];
rz(-2.4430165) q[1];
x q[2];
rz(-2.4279057) q[3];
sx q[3];
rz(-2.1296576) q[3];
sx q[3];
rz(-1.9969459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(-1.9866379) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44523859) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(0.44662961) q[1];
sx q[1];
rz(-0.59372562) q[1];
sx q[1];
rz(-2.4536224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5410091) q[0];
sx q[0];
rz(-0.86555153) q[0];
sx q[0];
rz(1.3441104) q[0];
x q[1];
rz(-0.61882682) q[2];
sx q[2];
rz(-1.8558673) q[2];
sx q[2];
rz(3.0322591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7862305) q[1];
sx q[1];
rz(-1.3578102) q[1];
sx q[1];
rz(1.1287685) q[1];
rz(-pi) q[2];
rz(-2.4554208) q[3];
sx q[3];
rz(-0.99643512) q[3];
sx q[3];
rz(0.8818834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.530431) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(-2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.68194836) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(-1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(0.43509126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97532192) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(0.50384371) q[0];
rz(-2.9138759) q[2];
sx q[2];
rz(-0.68074742) q[2];
sx q[2];
rz(-1.5071703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7104946) q[1];
sx q[1];
rz(-0.63023797) q[1];
sx q[1];
rz(-1.789854) q[1];
x q[2];
rz(3.1396542) q[3];
sx q[3];
rz(-2.7172305) q[3];
sx q[3];
rz(-0.75943179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6626176) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-2.5153861) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58186746) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(-0.054280601) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(-1.8064226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30935449) q[0];
sx q[0];
rz(-2.1119499) q[0];
sx q[0];
rz(-0.65017976) q[0];
rz(-1.0112052) q[2];
sx q[2];
rz(-1.8768684) q[2];
sx q[2];
rz(-1.0572421) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74354913) q[1];
sx q[1];
rz(-0.7760007) q[1];
sx q[1];
rz(0.77568027) q[1];
rz(1.4328721) q[3];
sx q[3];
rz(-1.63878) q[3];
sx q[3];
rz(-1.2757511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6173031) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-2.7900901) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897301) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(-2.3824298) q[0];
rz(-2.0579386) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(-2.0738475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46684346) q[0];
sx q[0];
rz(-2.3567794) q[0];
sx q[0];
rz(-2.270257) q[0];
rz(-pi) q[1];
rz(-2.3876486) q[2];
sx q[2];
rz(-0.39421668) q[2];
sx q[2];
rz(-0.27335247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0800277) q[1];
sx q[1];
rz(-1.0288218) q[1];
sx q[1];
rz(-2.6656277) q[1];
rz(-0.87369793) q[3];
sx q[3];
rz(-1.3109968) q[3];
sx q[3];
rz(1.8099305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(2.9313226) q[2];
rz(0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-0.53019607) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982211) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.632985) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.2589113) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(2.0785594) q[3];
sx q[3];
rz(-2.8294143) q[3];
sx q[3];
rz(-2.1635273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
