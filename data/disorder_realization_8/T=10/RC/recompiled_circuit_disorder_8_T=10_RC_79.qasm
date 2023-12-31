OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(-0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68569505) q[0];
sx q[0];
rz(-1.0766317) q[0];
sx q[0];
rz(-2.6463638) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4279847) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(1.7507391) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20476725) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(0.25524615) q[1];
rz(-pi) q[2];
rz(-2.9526887) q[3];
sx q[3];
rz(-1.2652664) q[3];
sx q[3];
rz(2.4019935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(2.105666) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(2.3764215) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457263) q[0];
sx q[0];
rz(-1.6329137) q[0];
sx q[0];
rz(1.6371884) q[0];
rz(-pi) q[1];
rz(-2.2234369) q[2];
sx q[2];
rz(-1.2456206) q[2];
sx q[2];
rz(0.52415372) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6389097) q[1];
sx q[1];
rz(-1.1306292) q[1];
sx q[1];
rz(1.4795951) q[1];
x q[2];
rz(1.4088267) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0654304) q[0];
sx q[0];
rz(-0.87834529) q[0];
sx q[0];
rz(2.6718219) q[0];
x q[1];
rz(-0.33883314) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(-1.1178521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0573404) q[1];
sx q[1];
rz(-0.80576128) q[1];
sx q[1];
rz(-1.8112103) q[1];
x q[2];
rz(1.5699925) q[3];
sx q[3];
rz(-1.6254079) q[3];
sx q[3];
rz(0.87517232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(0.51309103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95272321) q[0];
sx q[0];
rz(-0.9523305) q[0];
sx q[0];
rz(-1.6678715) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80715837) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(1.6023139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0780371) q[1];
sx q[1];
rz(-0.47680285) q[1];
sx q[1];
rz(2.8366691) q[1];
rz(-pi) q[2];
rz(1.7170834) q[3];
sx q[3];
rz(-2.9330367) q[3];
sx q[3];
rz(-1.8640765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-2.9117057) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992633) q[0];
sx q[0];
rz(-3.0780601) q[0];
sx q[0];
rz(-0.97139831) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4511756) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(2.5809443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7580326) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(-3.1314965) q[1];
rz(-1.945799) q[3];
sx q[3];
rz(-0.71205322) q[3];
sx q[3];
rz(-0.28111162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-2.8862254) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(-0.88476673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038267604) q[0];
sx q[0];
rz(-1.9511576) q[0];
sx q[0];
rz(0.079770712) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0753291) q[2];
sx q[2];
rz(-0.21050669) q[2];
sx q[2];
rz(1.9312242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(-2.8273696) q[1];
x q[2];
rz(2.1045477) q[3];
sx q[3];
rz(-1.3579206) q[3];
sx q[3];
rz(2.6900929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(-1.3686251) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(3.080522) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-0.73289245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31156763) q[0];
sx q[0];
rz(-0.68651474) q[0];
sx q[0];
rz(-2.4870883) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1204719) q[2];
sx q[2];
rz(-2.3602544) q[2];
sx q[2];
rz(2.6598425) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.43436189) q[1];
sx q[1];
rz(-1.5642628) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6781428) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(-2.6220208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-0.27871305) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3364209) q[0];
sx q[0];
rz(-1.3609017) q[0];
sx q[0];
rz(1.7091442) q[0];
rz(-pi) q[1];
rz(2.7630373) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(2.2373667) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(-1.8885683) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7940984) q[3];
sx q[3];
rz(-1.6468862) q[3];
sx q[3];
rz(3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.999058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0363077) q[0];
sx q[0];
rz(-1.6616016) q[0];
sx q[0];
rz(1.8663976) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70334401) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(-1.4700996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.090230378) q[1];
sx q[1];
rz(-1.4453332) q[1];
sx q[1];
rz(-2.5999703) q[1];
rz(-pi) q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(-2.4329894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(-2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.2385626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0093875) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(1.042004) q[0];
rz(3.086834) q[2];
sx q[2];
rz(-0.75373947) q[2];
sx q[2];
rz(-2.4278305) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.726767) q[1];
sx q[1];
rz(-1.4531724) q[1];
sx q[1];
rz(-2.2284501) q[1];
rz(-pi) q[2];
x q[2];
rz(0.978312) q[3];
sx q[3];
rz(-1.4761792) q[3];
sx q[3];
rz(0.30944165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(3.0977541) q[2];
rz(-1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713521) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(0.025370601) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.2123443) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
rz(-0.17037114) q[3];
sx q[3];
rz(-0.80633612) q[3];
sx q[3];
rz(2.7663305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
