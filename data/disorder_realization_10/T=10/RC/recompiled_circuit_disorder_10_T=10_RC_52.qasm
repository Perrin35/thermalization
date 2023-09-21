OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(-2.8432863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059271185) q[0];
sx q[0];
rz(-0.66723171) q[0];
sx q[0];
rz(3.0526403) q[0];
rz(-1.9036129) q[2];
sx q[2];
rz(-1.4822072) q[2];
sx q[2];
rz(-2.7804136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5391985) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(0.56166517) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.467642) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(-0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988929) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52662151) q[0];
sx q[0];
rz(-1.7978298) q[0];
sx q[0];
rz(-1.4814266) q[0];
rz(-pi) q[1];
rz(-0.50498982) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.1438952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3018803) q[1];
sx q[1];
rz(-1.0186968) q[1];
sx q[1];
rz(1.2549972) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15700335) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(-2.4412145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(-2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244814) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(-2.2729421) q[0];
rz(-pi) q[1];
rz(-1.9934898) q[2];
sx q[2];
rz(-0.87562497) q[2];
sx q[2];
rz(2.9425651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.599217) q[1];
sx q[1];
rz(-1.5812751) q[1];
sx q[1];
rz(3.0159365) q[1];
rz(-2.9946795) q[3];
sx q[3];
rz(-2.5111755) q[3];
sx q[3];
rz(2.279225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(-0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(-2.0137285) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(-1.9086054) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1283135) q[0];
sx q[0];
rz(-1.0767125) q[0];
sx q[0];
rz(-1.4923151) q[0];
rz(-pi) q[1];
rz(-1.0225251) q[2];
sx q[2];
rz(-2.0892482) q[2];
sx q[2];
rz(-3.0419635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2265046) q[1];
sx q[1];
rz(-0.18208948) q[1];
sx q[1];
rz(1.3703129) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22709417) q[3];
sx q[3];
rz(-2.2507651) q[3];
sx q[3];
rz(1.029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994613) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034886995) q[0];
sx q[0];
rz(-1.5505152) q[0];
sx q[0];
rz(-0.01625343) q[0];
x q[1];
rz(0.038482484) q[2];
sx q[2];
rz(-0.96652346) q[2];
sx q[2];
rz(-2.2272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37413874) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(2.2351082) q[1];
rz(2.6796474) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(2.7355821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(-1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-0.22512063) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(0.37757847) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7901944) q[0];
sx q[0];
rz(-1.5927918) q[0];
sx q[0];
rz(1.3047332) q[0];
rz(-2.8541366) q[2];
sx q[2];
rz(-1.524458) q[2];
sx q[2];
rz(-3.1099144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-2.6968323) q[1];
sx q[1];
rz(-1.8163535) q[1];
rz(-1.5402921) q[3];
sx q[3];
rz(-2.3654733) q[3];
sx q[3];
rz(3.0200849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(-2.8267982) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0525381) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(-0.25442466) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62691488) q[0];
sx q[0];
rz(-1.2509545) q[0];
sx q[0];
rz(-3.0266648) q[0];
rz(-pi) q[1];
rz(0.15433407) q[2];
sx q[2];
rz(-0.96417226) q[2];
sx q[2];
rz(0.91147214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88528819) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(-2.6726252) q[1];
rz(-0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(-1.3158201) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(-1.6865431) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5664068) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9109089) q[0];
sx q[0];
rz(-0.77373234) q[0];
sx q[0];
rz(-0.45549972) q[0];
x q[1];
rz(1.8640395) q[2];
sx q[2];
rz(-0.55985057) q[2];
sx q[2];
rz(-0.44550371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6187001) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(2.3431542) q[1];
rz(-2.4616562) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-0.27754647) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(2.5121571) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(2.004752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0658543) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(-2.0980741) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9929664) q[2];
sx q[2];
rz(-1.2264226) q[2];
sx q[2];
rz(0.070377199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5053619) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(1.2760217) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73590274) q[3];
sx q[3];
rz(-1.2888442) q[3];
sx q[3];
rz(1.6424996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-1.0409522) q[2];
rz(-0.0020290931) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(-2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.3025538) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(0.94605207) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-2.5295703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8004868) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(0.97408803) q[0];
rz(-0.34531784) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(-0.40030865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.260173) q[1];
sx q[1];
rz(-2.2248785) q[1];
sx q[1];
rz(-1.7263078) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1935812) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(-2.6138888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(-0.65336147) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(0.70070926) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(0.75469771) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(2.1379708) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(1.2470506) q[3];
sx q[3];
rz(-1.4907881) q[3];
sx q[3];
rz(-1.1992906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];