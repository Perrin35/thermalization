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
rz(-1.7582769) q[0];
sx q[0];
rz(4.5780616) q[0];
sx q[0];
rz(8.4599001) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(5.8536018) q[1];
sx q[1];
rz(11.516973) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222264) q[0];
sx q[0];
rz(-2.1458186) q[0];
sx q[0];
rz(0.40798835) q[0];
x q[1];
rz(-2.9994316) q[2];
sx q[2];
rz(-1.249525) q[2];
sx q[2];
rz(-2.4427593) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2856372) q[1];
sx q[1];
rz(-0.95889839) q[1];
sx q[1];
rz(2.6735071) q[1];
rz(-pi) q[2];
rz(2.562306) q[3];
sx q[3];
rz(-1.6484954) q[3];
sx q[3];
rz(0.9723817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45399484) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(3.122984) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740771) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(0.53502214) q[0];
rz(-0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.7417057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0947123) q[0];
sx q[0];
rz(-2.0167354) q[0];
sx q[0];
rz(-2.711854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56324701) q[2];
sx q[2];
rz(-1.0457102) q[2];
sx q[2];
rz(-1.6890749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8230086) q[1];
sx q[1];
rz(-0.46097791) q[1];
sx q[1];
rz(-2.5557842) q[1];
x q[2];
rz(-2.4055491) q[3];
sx q[3];
rz(-2.3793845) q[3];
sx q[3];
rz(-2.8413642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(-0.92607099) q[2];
rz(2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(-2.4721691) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7922908) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(1.5128304) q[0];
rz(-2.8577562) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(2.0294752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29942) q[0];
sx q[0];
rz(-1.5779691) q[0];
sx q[0];
rz(-0.0070863574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1826964) q[2];
sx q[2];
rz(-3.009353) q[2];
sx q[2];
rz(-0.43741495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97858799) q[1];
sx q[1];
rz(-1.7573865) q[1];
sx q[1];
rz(2.6546937) q[1];
rz(2.1792382) q[3];
sx q[3];
rz(-0.43399226) q[3];
sx q[3];
rz(2.2631462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6185559) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(-2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(-0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2570626) q[0];
sx q[0];
rz(-1.8842959) q[0];
sx q[0];
rz(-2.351981) q[0];
x q[1];
rz(2.7233549) q[2];
sx q[2];
rz(-0.23242885) q[2];
sx q[2];
rz(-1.4168036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9758983) q[1];
sx q[1];
rz(-1.1912575) q[1];
sx q[1];
rz(1.5089773) q[1];
rz(-2.9474218) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7249001) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(-2.5856384) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(-2.8506193) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81412643) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(-0.59447527) q[0];
rz(-0.56198436) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-2.1655703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.77895) q[0];
sx q[0];
rz(-2.5651826) q[0];
sx q[0];
rz(2.2767785) q[0];
rz(-0.75404928) q[2];
sx q[2];
rz(-1.591914) q[2];
sx q[2];
rz(-1.4625975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12078602) q[1];
sx q[1];
rz(-0.73985064) q[1];
sx q[1];
rz(-0.45005656) q[1];
rz(-pi) q[2];
rz(-1.7877949) q[3];
sx q[3];
rz(-1.2760538) q[3];
sx q[3];
rz(-2.6685696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(-2.5308934) q[2];
rz(2.8357909) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(-1.8122199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82764757) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(-0.77955359) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(0.73878845) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29863623) q[0];
sx q[0];
rz(-1.2901352) q[0];
sx q[0];
rz(-2.929351) q[0];
rz(-pi) q[1];
rz(2.8347978) q[2];
sx q[2];
rz(-1.580913) q[2];
sx q[2];
rz(-2.1515357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14911095) q[1];
sx q[1];
rz(-2.457239) q[1];
sx q[1];
rz(2.1696287) q[1];
x q[2];
rz(1.9442648) q[3];
sx q[3];
rz(-1.9091144) q[3];
sx q[3];
rz(-0.23972971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61591992) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(1.020485) q[0];
rz(0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(-2.3540156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0840698) q[0];
sx q[0];
rz(-2.2057187) q[0];
sx q[0];
rz(1.8674839) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8767005) q[2];
sx q[2];
rz(-0.26703003) q[2];
sx q[2];
rz(1.1970113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0246525) q[1];
sx q[1];
rz(-1.3924011) q[1];
sx q[1];
rz(-1.0712888) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.034463783) q[3];
sx q[3];
rz(-0.98317819) q[3];
sx q[3];
rz(1.8922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3202177) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(2.5569432) q[2];
rz(-0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(0.32824326) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(-1.6627056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2770689) q[0];
sx q[0];
rz(-1.3448937) q[0];
sx q[0];
rz(-0.27033131) q[0];
x q[1];
rz(2.4166862) q[2];
sx q[2];
rz(-0.6577684) q[2];
sx q[2];
rz(-2.7810514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11301576) q[1];
sx q[1];
rz(-1.8488819) q[1];
sx q[1];
rz(3.0427631) q[1];
rz(-pi) q[2];
rz(0.21027474) q[3];
sx q[3];
rz(-1.9326903) q[3];
sx q[3];
rz(-2.0023605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(2.1797144) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-2.403517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2797425) q[0];
sx q[0];
rz(-1.334001) q[0];
sx q[0];
rz(-2.7351232) q[0];
rz(-pi) q[1];
rz(2.5665087) q[2];
sx q[2];
rz(-0.7204537) q[2];
sx q[2];
rz(1.6937814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2710423) q[1];
sx q[1];
rz(-1.9890474) q[1];
sx q[1];
rz(-0.31914945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35154147) q[3];
sx q[3];
rz(-0.69720399) q[3];
sx q[3];
rz(0.54920025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(-2.2743478) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(2.0713461) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27425048) q[0];
sx q[0];
rz(-0.19495067) q[0];
sx q[0];
rz(2.250953) q[0];
rz(-pi) q[1];
rz(-0.81999166) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(-1.9209678) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.874843) q[1];
sx q[1];
rz(-2.0834181) q[1];
sx q[1];
rz(-1.1260518) q[1];
x q[2];
rz(1.460288) q[3];
sx q[3];
rz(-1.432087) q[3];
sx q[3];
rz(-1.8294301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(-1.2316068) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3043542) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-0.41481836) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(-1.6258705) q[2];
sx q[2];
rz(-2.016042) q[2];
sx q[2];
rz(2.945937) q[2];
rz(2.5255193) q[3];
sx q[3];
rz(-2.2884634) q[3];
sx q[3];
rz(-0.4175755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
