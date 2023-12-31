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
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5814712) q[0];
sx q[0];
rz(-1.6257964) q[0];
sx q[0];
rz(2.4762857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3053427) q[2];
sx q[2];
rz(-0.34398088) q[2];
sx q[2];
rz(1.460182) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5391985) q[1];
sx q[1];
rz(-1.0326003) q[1];
sx q[1];
rz(0.56166517) q[1];
rz(1.467642) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(-2.8148357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(-0.47168628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90549201) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(2.7729176) q[0];
rz(1.291044) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(-0.56088698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3018803) q[1];
sx q[1];
rz(-2.1228959) q[1];
sx q[1];
rz(-1.2549972) q[1];
x q[2];
rz(-2.9845893) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(-0.70037819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.6859432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239379) q[0];
sx q[0];
rz(-1.0929937) q[0];
sx q[0];
rz(2.2295203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4007912) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(-1.4893116) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0297444) q[1];
sx q[1];
rz(-1.4451471) q[1];
sx q[1];
rz(1.5602342) q[1];
rz(-6/(13*pi)) q[3];
sx q[3];
rz(-2.5111755) q[3];
sx q[3];
rz(0.86236766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85686344) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-2.5033584) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9905332) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(-0.14453669) q[0];
rz(-2.5523283) q[2];
sx q[2];
rz(-1.1009842) q[2];
sx q[2];
rz(1.9643009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7113263) q[1];
sx q[1];
rz(-1.3923936) q[1];
sx q[1];
rz(-3.1049411) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22709417) q[3];
sx q[3];
rz(-0.89082754) q[3];
sx q[3];
rz(-2.1122776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-0.81400648) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-0.89170757) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4310303) q[0];
sx q[0];
rz(-0.02598962) q[0];
sx q[0];
rz(2.2463069) q[0];
x q[1];
rz(-1.5151305) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(2.2948613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58442851) q[1];
sx q[1];
rz(-0.85313988) q[1];
sx q[1];
rz(0.75259705) q[1];
rz(-0.84380031) q[3];
sx q[3];
rz(-1.926933) q[3];
sx q[3];
rz(0.86405495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-2.6203716) q[2];
rz(-1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.2683755) q[3];
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
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(-2.7640142) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9281884) q[0];
sx q[0];
rz(-1.3047991) q[0];
sx q[0];
rz(3.1187952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6191143) q[2];
sx q[2];
rz(-1.8579351) q[2];
sx q[2];
rz(1.5887807) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7403455) q[1];
sx q[1];
rz(-1.1402854) q[1];
sx q[1];
rz(0.11534782) q[1];
rz(-pi) q[2];
rz(-2.3466831) q[3];
sx q[3];
rz(-1.5494293) q[3];
sx q[3];
rz(-1.6705318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(2.8267982) q[3];
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
rz(-1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(-2.9220707) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(2.887168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97878153) q[0];
sx q[0];
rz(-2.8023976) q[0];
sx q[0];
rz(1.9041054) q[0];
x q[1];
rz(2.9872586) q[2];
sx q[2];
rz(-2.1774204) q[2];
sx q[2];
rz(-2.2301205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88528819) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(2.6726252) q[1];
x q[2];
rz(0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1640132) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(-2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106237) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5664068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63021916) q[0];
sx q[0];
rz(-0.89238088) q[0];
sx q[0];
rz(1.1648965) q[0];
x q[1];
rz(-0.17922108) q[2];
sx q[2];
rz(-1.0374829) q[2];
sx q[2];
rz(2.3537677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7944784) q[1];
sx q[1];
rz(-1.9273259) q[1];
sx q[1];
rz(-1.9630678) q[1];
rz(-pi) q[2];
rz(1.992222) q[3];
sx q[3];
rz(-2.0391658) q[3];
sx q[3];
rz(2.3295662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2660797) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(2.1267166) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(1.45654) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4515848) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(0.33126979) q[0];
x q[1];
rz(-2.7669737) q[2];
sx q[2];
rz(-1.9667452) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9469103) q[1];
sx q[1];
rz(-1.2543251) q[1];
sx q[1];
rz(0.099593347) q[1];
rz(-pi) q[2];
rz(-1.198248) q[3];
sx q[3];
rz(-0.87009831) q[3];
sx q[3];
rz(-2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(1.0409522) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(0.61202234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8004868) q[0];
sx q[0];
rz(-0.90011156) q[0];
sx q[0];
rz(0.97408803) q[0];
rz(-2.5435796) q[2];
sx q[2];
rz(-0.41053718) q[2];
sx q[2];
rz(1.7288127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(2.481639) q[1];
rz(-pi) q[2];
rz(2.6066783) q[3];
sx q[3];
rz(-1.8992918) q[3];
sx q[3];
rz(1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0570021) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4738612) q[2];
sx q[2];
rz(-0.79955352) q[2];
sx q[2];
rz(-0.59895589) q[2];
rz(0.084372088) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
