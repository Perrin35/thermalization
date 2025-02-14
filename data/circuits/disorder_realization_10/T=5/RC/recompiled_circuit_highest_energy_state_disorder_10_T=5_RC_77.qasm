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
rz(-1.7614814) q[0];
sx q[0];
rz(-1.5505646) q[0];
sx q[0];
rz(-2.759759) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782717) q[0];
sx q[0];
rz(-0.98211432) q[0];
sx q[0];
rz(-0.6529385) q[0];
rz(-pi) q[1];
rz(2.7713791) q[2];
sx q[2];
rz(-1.5523124) q[2];
sx q[2];
rz(-3.0188675) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1273444) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(-2.8365652) q[1];
rz(0.96278874) q[3];
sx q[3];
rz(-2.5217524) q[3];
sx q[3];
rz(-2.7328797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(-1.5748242) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.1760733) q[0];
rz(-0.067497079) q[1];
sx q[1];
rz(-0.57756966) q[1];
sx q[1];
rz(2.8228021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8628937) q[0];
sx q[0];
rz(-2.3434533) q[0];
sx q[0];
rz(2.7006335) q[0];
rz(-pi) q[1];
rz(-1.8068878) q[2];
sx q[2];
rz(-1.1523048) q[2];
sx q[2];
rz(2.921517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6401854) q[1];
sx q[1];
rz(-2.1245271) q[1];
sx q[1];
rz(-2.6928765) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0158668) q[3];
sx q[3];
rz(-1.2866402) q[3];
sx q[3];
rz(-1.3181669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0789644) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(0.47743615) q[2];
rz(-1.3580953) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(0.0099946578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.39182144) q[0];
sx q[0];
rz(-1.8056159) q[0];
sx q[0];
rz(-1.1392449) q[0];
rz(-2.7353752) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11230532) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(-1.4112006) q[0];
rz(-2.1966372) q[2];
sx q[2];
rz(-1.2765108) q[2];
sx q[2];
rz(-1.4465894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.660608) q[1];
sx q[1];
rz(-1.0502195) q[1];
sx q[1];
rz(1.03517) q[1];
rz(-pi) q[2];
rz(1.0615248) q[3];
sx q[3];
rz(-1.8647927) q[3];
sx q[3];
rz(-1.7727838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-0.25585678) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(2.381109) q[0];
rz(2.7898232) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(-0.15028353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4486769) q[0];
sx q[0];
rz(-1.5741072) q[0];
sx q[0];
rz(0.0019372367) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1885095) q[2];
sx q[2];
rz(-1.7312382) q[2];
sx q[2];
rz(0.40942243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1898415) q[1];
sx q[1];
rz(-0.69165666) q[1];
sx q[1];
rz(1.0361978) q[1];
rz(-pi) q[2];
rz(1.3129381) q[3];
sx q[3];
rz(-2.0153196) q[3];
sx q[3];
rz(-0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1111697) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(0.58491582) q[2];
rz(0.065718204) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(-1.0544624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718662) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(0.72750339) q[0];
rz(-1.8456521) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(2.4268699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172404) q[0];
sx q[0];
rz(-2.9091798) q[0];
sx q[0];
rz(2.8281141) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6863912) q[2];
sx q[2];
rz(-1.8697479) q[2];
sx q[2];
rz(1.7345127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2421869) q[1];
sx q[1];
rz(-1.9493628) q[1];
sx q[1];
rz(2.942571) q[1];
rz(-pi) q[2];
rz(-2.1470137) q[3];
sx q[3];
rz(-1.2354697) q[3];
sx q[3];
rz(1.9063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99594816) q[2];
sx q[2];
rz(-1.5278634) q[2];
sx q[2];
rz(-1.145251) q[2];
rz(-2.9389985) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(-2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2957434) q[0];
sx q[0];
rz(-2.8491617) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(-1.3289183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5707626) q[0];
sx q[0];
rz(-1.098363) q[0];
sx q[0];
rz(-2.5357312) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0014056) q[2];
sx q[2];
rz(-1.5221688) q[2];
sx q[2];
rz(1.0076154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3420978) q[1];
sx q[1];
rz(-1.8256542) q[1];
sx q[1];
rz(-1.8138769) q[1];
rz(-pi) q[2];
rz(0.82002016) q[3];
sx q[3];
rz(-2.8201172) q[3];
sx q[3];
rz(-3.0177651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38587511) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(-2.7937549) q[2];
rz(-0.31473413) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1494074) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(2.8490207) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(2.1975885) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2601449) q[0];
sx q[0];
rz(-1.8788188) q[0];
sx q[0];
rz(-0.016168895) q[0];
x q[1];
rz(0.47796504) q[2];
sx q[2];
rz(-2.0066119) q[2];
sx q[2];
rz(-0.27406853) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71783644) q[1];
sx q[1];
rz(-1.5482117) q[1];
sx q[1];
rz(-1.3976239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0349501) q[3];
sx q[3];
rz(-0.45848819) q[3];
sx q[3];
rz(1.5036086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(1.7000807) q[2];
rz(3.1169543) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(-0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330419) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(0.1380052) q[0];
rz(1.0778614) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-2.2053351) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753006) q[0];
sx q[0];
rz(-1.202824) q[0];
sx q[0];
rz(1.8938246) q[0];
x q[1];
rz(-2.594844) q[2];
sx q[2];
rz(-2.2348352) q[2];
sx q[2];
rz(-0.67686096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.028658) q[1];
sx q[1];
rz(-1.7535661) q[1];
sx q[1];
rz(0.014928047) q[1];
rz(0.68558399) q[3];
sx q[3];
rz(-1.3635677) q[3];
sx q[3];
rz(-2.384546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.998698) q[2];
sx q[2];
rz(-1.2724718) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-0.43012777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97813022) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(0.36460707) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(2.2273831) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5152138) q[0];
sx q[0];
rz(-0.45013407) q[0];
sx q[0];
rz(-1.2594957) q[0];
x q[1];
rz(-3.0905452) q[2];
sx q[2];
rz(-2.406359) q[2];
sx q[2];
rz(0.027054199) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57371512) q[1];
sx q[1];
rz(-2.2063631) q[1];
sx q[1];
rz(-1.741841) q[1];
rz(2.6007341) q[3];
sx q[3];
rz(-1.4424994) q[3];
sx q[3];
rz(3.1311664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8387575) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(1.5128822) q[2];
rz(-0.58879876) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673679) q[0];
sx q[0];
rz(-2.5741757) q[0];
sx q[0];
rz(-0.2844511) q[0];
rz(2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(-0.42025748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3271844) q[0];
sx q[0];
rz(-1.1564009) q[0];
sx q[0];
rz(0.23420686) q[0];
rz(-pi) q[1];
rz(-0.45966123) q[2];
sx q[2];
rz(-0.84520413) q[2];
sx q[2];
rz(2.9265253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6984562) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(3.0522947) q[1];
rz(2.612667) q[3];
sx q[3];
rz(-2.5074337) q[3];
sx q[3];
rz(-0.24833873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9922716) q[2];
sx q[2];
rz(-0.48506609) q[2];
sx q[2];
rz(2.1743656) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(3.0178487) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(0.32044784) q[2];
sx q[2];
rz(-1.4437699) q[2];
sx q[2];
rz(-0.52158471) q[2];
rz(-3.0322778) q[3];
sx q[3];
rz(-2.8277665) q[3];
sx q[3];
rz(1.2841429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
