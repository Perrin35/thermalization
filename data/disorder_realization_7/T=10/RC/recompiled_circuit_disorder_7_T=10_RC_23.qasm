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
rz(2.9404844) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62030828) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(-1.3814397) q[0];
x q[1];
rz(-2.7896499) q[2];
sx q[2];
rz(-1.3411739) q[2];
sx q[2];
rz(0.72598347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1605692) q[1];
sx q[1];
rz(-1.8796646) q[1];
sx q[1];
rz(-0.94998756) q[1];
rz(-0.17625531) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(2.3034629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24793808) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(2.1564116) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(-0.10736297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185703) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(2.1711711) q[0];
x q[1];
rz(0.58789247) q[2];
sx q[2];
rz(-1.3304454) q[2];
sx q[2];
rz(-2.8628778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.098298479) q[1];
sx q[1];
rz(-0.77273332) q[1];
sx q[1];
rz(-2.1143101) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7530305) q[3];
sx q[3];
rz(-1.4149168) q[3];
sx q[3];
rz(0.673783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(-2.4798933) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4989935) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(2.5901592) q[0];
rz(1.9252831) q[2];
sx q[2];
rz(-1.4282994) q[2];
sx q[2];
rz(0.49152495) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(-pi) q[2];
rz(-0.2228959) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(2.1427758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-0.55251399) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(-2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195388) q[0];
sx q[0];
rz(-1.4931803) q[0];
sx q[0];
rz(-1.104943) q[0];
rz(0.98056356) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(0.24888466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19255895) q[1];
sx q[1];
rz(-0.20201905) q[1];
sx q[1];
rz(-2.2662524) q[1];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-1.0415047) q[3];
sx q[3];
rz(-0.57100163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5421211) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(0.19651861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4465966) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(2.8269935) q[0];
x q[1];
rz(1.1097124) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(1.920514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6477485) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(2.8791048) q[1];
x q[2];
rz(-0.5444364) q[3];
sx q[3];
rz(-0.32005537) q[3];
sx q[3];
rz(0.65359945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(0.83459485) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4112339) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(-0.22986869) q[0];
x q[1];
rz(3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(0.44484777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-2.1761314) q[1];
rz(-2.4059541) q[3];
sx q[3];
rz(-1.2218352) q[3];
sx q[3];
rz(-2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94971925) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(-2.9313415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1649095) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(-0.066819013) q[0];
rz(-pi) q[1];
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
x q[0];
rz(1.5308295) q[1];
sx q[1];
rz(-1.7621343) q[1];
sx q[1];
rz(2.8192239) q[1];
rz(0.57742124) q[3];
sx q[3];
rz(-2.2038659) q[3];
sx q[3];
rz(0.046618332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92102345) q[0];
sx q[0];
rz(-1.1194293) q[0];
sx q[0];
rz(1.4054106) q[0];
rz(1.1114242) q[2];
sx q[2];
rz(-2.9615235) q[2];
sx q[2];
rz(-1.1262696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2446489) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(-2.0379373) q[1];
rz(2.5891853) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(0.92064944) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793648) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(-2.2158932) q[0];
x q[1];
rz(-2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(0.65288359) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2194849) q[1];
sx q[1];
rz(-1.2412294) q[1];
sx q[1];
rz(-0.12164128) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5713463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(-0.40035029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1966232) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(-2.1150132) q[0];
rz(-pi) q[1];
rz(2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(-0.099909401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6798903) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(3.0479792) q[1];
rz(-pi) q[2];
rz(-0.40542094) q[3];
sx q[3];
rz(-2.0339031) q[3];
sx q[3];
rz(-2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(3.0336822) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
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
