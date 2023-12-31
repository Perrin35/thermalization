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
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0878108) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(-1.5009297) q[0];
rz(-0.093703336) q[2];
sx q[2];
rz(-1.2393349) q[2];
sx q[2];
rz(-1.9625488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71516192) q[1];
sx q[1];
rz(-0.75725812) q[1];
sx q[1];
rz(-0.84233474) q[1];
x q[2];
rz(-0.88793036) q[3];
sx q[3];
rz(-3.008932) q[3];
sx q[3];
rz(0.56548972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-0.99386627) q[2];
rz(-0.99938756) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(-2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(0.47168628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0240078) q[0];
sx q[0];
rz(-1.6578668) q[0];
sx q[0];
rz(-2.9136806) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6623146) q[2];
sx q[2];
rz(-2.5844378) q[2];
sx q[2];
rz(-0.011475871) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43863152) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(0.57499927) q[1];
rz(-pi) q[2];
rz(-2.178533) q[3];
sx q[3];
rz(-1.7000323) q[3];
sx q[3];
rz(2.3605763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(-1.4556494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(2.2729421) q[0];
rz(-0.45708926) q[2];
sx q[2];
rz(-0.79490137) q[2];
sx q[2];
rz(2.7283816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1959343) q[1];
sx q[1];
rz(-3.0155026) q[1];
sx q[1];
rz(3.0581711) q[1];
x q[2];
rz(1.4643747) q[3];
sx q[3];
rz(-0.94821804) q[3];
sx q[3];
rz(0.6811844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-0.23770604) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9905332) q[0];
sx q[0];
rz(-0.49976832) q[0];
sx q[0];
rz(2.997056) q[0];
rz(-2.4013176) q[2];
sx q[2];
rz(-0.7358272) q[2];
sx q[2];
rz(0.98878132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9150881) q[1];
sx q[1];
rz(-0.18208948) q[1];
sx q[1];
rz(1.7712797) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2992371) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(1.3815051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(0.81400648) q[2];
rz(-1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(-1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-0.89170757) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(0.2125425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71056238) q[0];
sx q[0];
rz(-0.02598962) q[0];
sx q[0];
rz(0.89528577) q[0];
rz(3.1031102) q[2];
sx q[2];
rz(-0.96652346) q[2];
sx q[2];
rz(-0.91439263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7674539) q[1];
sx q[1];
rz(-2.1530188) q[1];
sx q[1];
rz(-0.9064845) q[1];
rz(-pi) q[2];
rz(-2.0810633) q[3];
sx q[3];
rz(-2.3465996) q[3];
sx q[3];
rz(1.0802964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-2.7640142) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7901944) q[0];
sx q[0];
rz(-1.5927918) q[0];
sx q[0];
rz(-1.3047332) q[0];
x q[1];
rz(-0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(-3.1099144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4695417) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(1.3252392) q[1];
x q[2];
rz(-1.5402921) q[3];
sx q[3];
rz(-2.3654733) q[3];
sx q[3];
rz(-0.12150773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1569415) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(2.8267982) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0525381) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(-2.9220707) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-0.25442466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62691488) q[0];
sx q[0];
rz(-1.8906381) q[0];
sx q[0];
rz(-0.11492782) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7888072) q[2];
sx q[2];
rz(-0.62354747) q[2];
sx q[2];
rz(-2.4965198) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9421778) q[1];
sx q[1];
rz(-2.605038) q[1];
sx q[1];
rz(0.55120991) q[1];
rz(-pi) q[2];
rz(2.0394578) q[3];
sx q[3];
rz(-1.161876) q[3];
sx q[3];
rz(-1.8917781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63021916) q[0];
sx q[0];
rz(-2.2492118) q[0];
sx q[0];
rz(1.1648965) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0303866) q[2];
sx q[2];
rz(-1.4166797) q[2];
sx q[2];
rz(2.266778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7745061) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(-2.7584502) q[1];
rz(-pi) q[2];
rz(-0.67993645) q[3];
sx q[3];
rz(-0.61938647) q[3];
sx q[3];
rz(0.030127545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-2.1267166) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(-2.004752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0658543) q[0];
sx q[0];
rz(-0.59959164) q[0];
sx q[0];
rz(1.0435186) q[0];
rz(-1.9929664) q[2];
sx q[2];
rz(-1.9151701) q[2];
sx q[2];
rz(-0.070377199) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6362308) q[1];
sx q[1];
rz(-2.8103235) q[1];
sx q[1];
rz(1.2760217) q[1];
rz(-2.4056899) q[3];
sx q[3];
rz(-1.8527485) q[3];
sx q[3];
rz(1.6424996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-2.1955406) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(-0.61202234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6291954) q[0];
sx q[0];
rz(-2.0266268) q[0];
sx q[0];
rz(2.3770611) q[0];
x q[1];
rz(0.34531784) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(0.40030865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.260173) q[1];
sx q[1];
rz(-2.2248785) q[1];
sx q[1];
rz(1.7263078) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58935921) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(-0.15976957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0845906) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-1.0036219) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.3238974) q[3];
sx q[3];
rz(-2.8084451) q[3];
sx q[3];
rz(0.13767903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
