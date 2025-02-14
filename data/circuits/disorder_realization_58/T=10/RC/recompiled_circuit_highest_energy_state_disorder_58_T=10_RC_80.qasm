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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(8.6051056) q[1];
sx q[1];
rz(17.65045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1333251) q[0];
sx q[0];
rz(-3.016351) q[0];
sx q[0];
rz(2.6175346) q[0];
x q[1];
rz(-1.3875628) q[2];
sx q[2];
rz(-2.2361045) q[2];
sx q[2];
rz(2.5380425) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0173229) q[1];
sx q[1];
rz(-2.5376179) q[1];
sx q[1];
rz(1.9159957) q[1];
x q[2];
rz(-2.2584003) q[3];
sx q[3];
rz(-0.84005648) q[3];
sx q[3];
rz(-0.20945621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6370411) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(1.441628) q[2];
rz(1.0124892) q[3];
sx q[3];
rz(-1.9952521) q[3];
sx q[3];
rz(2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83990324) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(2.1121693) q[0];
rz(1.7816211) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(2.634868) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8520607) q[0];
sx q[0];
rz(-1.675005) q[0];
sx q[0];
rz(-3.122217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1400181) q[2];
sx q[2];
rz(-2.3990941) q[2];
sx q[2];
rz(3.1378821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3054333) q[1];
sx q[1];
rz(-0.99328128) q[1];
sx q[1];
rz(0.07825313) q[1];
x q[2];
rz(2.8414423) q[3];
sx q[3];
rz(-2.2734299) q[3];
sx q[3];
rz(-2.5182078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8257251) q[2];
sx q[2];
rz(-2.5434912) q[2];
sx q[2];
rz(-0.21710795) q[2];
rz(1.0202967) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(2.1609658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8751136) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-0.76817051) q[0];
rz(0.29113302) q[1];
sx q[1];
rz(-1.365064) q[1];
sx q[1];
rz(-1.1870144) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47531775) q[0];
sx q[0];
rz(-2.9143479) q[0];
sx q[0];
rz(-2.7623618) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0303346) q[2];
sx q[2];
rz(-3.0560962) q[2];
sx q[2];
rz(-2.4323815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95988454) q[1];
sx q[1];
rz(-1.2165637) q[1];
sx q[1];
rz(-1.541287) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1891425) q[3];
sx q[3];
rz(-0.72126167) q[3];
sx q[3];
rz(-0.14474584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27176303) q[2];
sx q[2];
rz(-2.2972079) q[2];
sx q[2];
rz(-0.96770206) q[2];
rz(0.23396954) q[3];
sx q[3];
rz(-2.440019) q[3];
sx q[3];
rz(2.7847737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135603) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(2.0618942) q[0];
rz(2.8375541) q[1];
sx q[1];
rz(-1.9875151) q[1];
sx q[1];
rz(1.8255723) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7181691) q[0];
sx q[0];
rz(-1.4470513) q[0];
sx q[0];
rz(2.607671) q[0];
x q[1];
rz(-1.448668) q[2];
sx q[2];
rz(-2.5425362) q[2];
sx q[2];
rz(1.0775104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5388064) q[1];
sx q[1];
rz(-0.99296184) q[1];
sx q[1];
rz(2.2972107) q[1];
rz(-1.997066) q[3];
sx q[3];
rz(-1.5139587) q[3];
sx q[3];
rz(2.6714966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5217846) q[2];
sx q[2];
rz(-1.9079756) q[2];
sx q[2];
rz(-3.0002777) q[2];
rz(2.5610793) q[3];
sx q[3];
rz(-1.4760735) q[3];
sx q[3];
rz(0.22835246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572131) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(1.6656531) q[0];
rz(-0.93302226) q[1];
sx q[1];
rz(-2.4304183) q[1];
sx q[1];
rz(1.3915541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45804241) q[0];
sx q[0];
rz(-1.6876966) q[0];
sx q[0];
rz(-2.3085672) q[0];
rz(-pi) q[1];
rz(-0.33220187) q[2];
sx q[2];
rz(-2.6966288) q[2];
sx q[2];
rz(1.6138168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36223026) q[1];
sx q[1];
rz(-1.0500776) q[1];
sx q[1];
rz(-3.1073276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0443893) q[3];
sx q[3];
rz(-2.5431051) q[3];
sx q[3];
rz(-2.583556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2453904) q[2];
sx q[2];
rz(-2.7707477) q[2];
sx q[2];
rz(0.24946269) q[2];
rz(1.4977411) q[3];
sx q[3];
rz(-1.4062873) q[3];
sx q[3];
rz(-2.7264285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6849218) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(-1.7556835) q[0];
rz(1.2617525) q[1];
sx q[1];
rz(-1.933814) q[1];
sx q[1];
rz(2.6127846) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3715293) q[0];
sx q[0];
rz(-1.125582) q[0];
sx q[0];
rz(-0.0078517954) q[0];
rz(-pi) q[1];
rz(0.98189378) q[2];
sx q[2];
rz(-1.0179276) q[2];
sx q[2];
rz(-0.58303631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13695194) q[1];
sx q[1];
rz(-1.6378504) q[1];
sx q[1];
rz(-0.93153684) q[1];
rz(2.4566023) q[3];
sx q[3];
rz(-1.6719975) q[3];
sx q[3];
rz(-1.2736959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1888107) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(1.3713651) q[2];
rz(0.2025226) q[3];
sx q[3];
rz(-1.4684418) q[3];
sx q[3];
rz(-1.9523581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093569) q[0];
sx q[0];
rz(-1.3487331) q[0];
sx q[0];
rz(-2.9877544) q[0];
rz(2.4065252) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(0.14911266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99054333) q[0];
sx q[0];
rz(-1.4874389) q[0];
sx q[0];
rz(-0.91611422) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17338578) q[2];
sx q[2];
rz(-1.9299704) q[2];
sx q[2];
rz(2.2196291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7478167) q[1];
sx q[1];
rz(-0.92499176) q[1];
sx q[1];
rz(-2.1834236) q[1];
rz(-pi) q[2];
rz(-1.0668236) q[3];
sx q[3];
rz(-1.6503346) q[3];
sx q[3];
rz(-0.37802896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4601124) q[2];
sx q[2];
rz(-1.6296547) q[2];
sx q[2];
rz(-0.45664772) q[2];
rz(-1.7234507) q[3];
sx q[3];
rz(-0.8816312) q[3];
sx q[3];
rz(-0.68896967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6316471) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(0.052635996) q[0];
rz(1.8317728) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(-1.2059258) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.238097) q[0];
sx q[0];
rz(-2.4795929) q[0];
sx q[0];
rz(2.2308626) q[0];
rz(-1.3383554) q[2];
sx q[2];
rz(-1.0224059) q[2];
sx q[2];
rz(-2.2218934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61319035) q[1];
sx q[1];
rz(-2.2397016) q[1];
sx q[1];
rz(1.7714798) q[1];
rz(-pi) q[2];
rz(1.2489592) q[3];
sx q[3];
rz(-2.494069) q[3];
sx q[3];
rz(-0.81784883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84411088) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(0.55997509) q[2];
rz(-1.0567793) q[3];
sx q[3];
rz(-2.4323075) q[3];
sx q[3];
rz(-0.81282508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63447222) q[0];
sx q[0];
rz(-1.7422603) q[0];
sx q[0];
rz(1.5647474) q[0];
rz(1.5693846) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-0.80894583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7943553) q[0];
sx q[0];
rz(-2.6449727) q[0];
sx q[0];
rz(-3.0625694) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39500631) q[2];
sx q[2];
rz(-0.85112903) q[2];
sx q[2];
rz(2.0499944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2613649) q[1];
sx q[1];
rz(-2.5305809) q[1];
sx q[1];
rz(-1.548012) q[1];
x q[2];
rz(-1.21502) q[3];
sx q[3];
rz(-2.9529533) q[3];
sx q[3];
rz(1.1805746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.481015) q[2];
sx q[2];
rz(-2.1889841) q[2];
sx q[2];
rz(0.048952254) q[2];
rz(-1.8598716) q[3];
sx q[3];
rz(-2.6661524) q[3];
sx q[3];
rz(1.6767282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.319741) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(-1.9191746) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(-0.45752057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8638141) q[0];
sx q[0];
rz(-0.84054986) q[0];
sx q[0];
rz(1.6919096) q[0];
x q[1];
rz(2.7399212) q[2];
sx q[2];
rz(-2.0435512) q[2];
sx q[2];
rz(0.061701802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9468675) q[1];
sx q[1];
rz(-1.3683142) q[1];
sx q[1];
rz(0.35670403) q[1];
rz(-2.6453552) q[3];
sx q[3];
rz(-0.75025815) q[3];
sx q[3];
rz(-1.6023146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7484625) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(1.040323) q[2];
rz(-2.6535502) q[3];
sx q[3];
rz(-1.4405684) q[3];
sx q[3];
rz(0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671065) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(-1.294301) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(-2.1637259) q[2];
sx q[2];
rz(-1.821283) q[2];
sx q[2];
rz(-1.4878426) q[2];
rz(-2.4518012) q[3];
sx q[3];
rz(-1.2648911) q[3];
sx q[3];
rz(-2.779724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
