OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3026128) q[0];
sx q[0];
rz(4.8470654) q[0];
sx q[0];
rz(9.6416311) q[0];
rz(2.6537553) q[1];
sx q[1];
rz(-2.2626329) q[1];
sx q[1];
rz(0.15329696) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76591208) q[0];
sx q[0];
rz(-1.3049676) q[0];
sx q[0];
rz(-1.8701118) q[0];
rz(-2.3011123) q[2];
sx q[2];
rz(-2.4908713) q[2];
sx q[2];
rz(0.67091984) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5342321) q[1];
sx q[1];
rz(-2.5776754) q[1];
sx q[1];
rz(2.6644071) q[1];
rz(2.9847095) q[3];
sx q[3];
rz(-0.81251729) q[3];
sx q[3];
rz(-2.7650327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6903901) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(0.85980493) q[2];
rz(-2.9016923) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(2.8587604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4848223) q[0];
sx q[0];
rz(-1.7253933) q[0];
sx q[0];
rz(2.2221478) q[0];
rz(2.2942309) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(1.1712317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396512) q[0];
sx q[0];
rz(-1.8293132) q[0];
sx q[0];
rz(-0.052019428) q[0];
rz(-0.33185766) q[2];
sx q[2];
rz(-0.88311354) q[2];
sx q[2];
rz(-0.0092384641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50368631) q[1];
sx q[1];
rz(-1.3136547) q[1];
sx q[1];
rz(-1.9288344) q[1];
x q[2];
rz(1.592335) q[3];
sx q[3];
rz(-1.0327368) q[3];
sx q[3];
rz(-3.1137385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0688811) q[2];
sx q[2];
rz(-2.2525807) q[2];
sx q[2];
rz(1.9286801) q[2];
rz(2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(-2.4685278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251959) q[0];
sx q[0];
rz(-1.238751) q[0];
sx q[0];
rz(0.58309251) q[0];
rz(-1.8816226) q[1];
sx q[1];
rz(-0.59695736) q[1];
sx q[1];
rz(-1.3139542) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2221916) q[0];
sx q[0];
rz(-1.254651) q[0];
sx q[0];
rz(-1.1911946) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5384244) q[2];
sx q[2];
rz(-0.83725032) q[2];
sx q[2];
rz(-0.59226474) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3594234) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(-1.729639) q[1];
rz(-pi) q[2];
rz(-1.2168808) q[3];
sx q[3];
rz(-2.6887213) q[3];
sx q[3];
rz(-2.4726506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6446357) q[2];
sx q[2];
rz(-2.3778215) q[2];
sx q[2];
rz(2.606707) q[2];
rz(2.6584117) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(3.0218637) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508761) q[0];
sx q[0];
rz(-3.0829939) q[0];
sx q[0];
rz(1.4709877) q[0];
rz(-1.2141256) q[1];
sx q[1];
rz(-1.43707) q[1];
sx q[1];
rz(-0.4462744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2504172) q[0];
sx q[0];
rz(-1.5022523) q[0];
sx q[0];
rz(1.7467996) q[0];
rz(0.52881119) q[2];
sx q[2];
rz(-1.5870278) q[2];
sx q[2];
rz(1.123614) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41149263) q[1];
sx q[1];
rz(-2.2124083) q[1];
sx q[1];
rz(1.8909251) q[1];
x q[2];
rz(-2/(11*pi)) q[3];
sx q[3];
rz(-2.8984275) q[3];
sx q[3];
rz(-1.1806219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2176167) q[2];
sx q[2];
rz(-0.46102229) q[2];
sx q[2];
rz(0.32611845) q[2];
rz(2.7255132) q[3];
sx q[3];
rz(-1.6385498) q[3];
sx q[3];
rz(2.3638341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.07051) q[0];
sx q[0];
rz(-1.5776881) q[0];
sx q[0];
rz(-0.34410205) q[0];
rz(-0.35585078) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(2.8867302) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6159812) q[0];
sx q[0];
rz(-0.92824844) q[0];
sx q[0];
rz(-1.8985974) q[0];
rz(-1.7698257) q[2];
sx q[2];
rz(-0.47037273) q[2];
sx q[2];
rz(-0.34188893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50824158) q[1];
sx q[1];
rz(-2.3741747) q[1];
sx q[1];
rz(-0.91752802) q[1];
x q[2];
rz(0.42307573) q[3];
sx q[3];
rz(-2.3300749) q[3];
sx q[3];
rz(-2.047416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43763375) q[2];
sx q[2];
rz(-0.86486977) q[2];
sx q[2];
rz(-0.09058365) q[2];
rz(-0.80638742) q[3];
sx q[3];
rz(-1.839829) q[3];
sx q[3];
rz(-1.8346132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059747132) q[0];
sx q[0];
rz(-1.9192001) q[0];
sx q[0];
rz(0.61862373) q[0];
rz(0.6126569) q[1];
sx q[1];
rz(-0.63059348) q[1];
sx q[1];
rz(0.15288606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02462834) q[0];
sx q[0];
rz(-1.3354339) q[0];
sx q[0];
rz(-0.23925608) q[0];
x q[1];
rz(-2.9866508) q[2];
sx q[2];
rz(-2.5563588) q[2];
sx q[2];
rz(0.80807782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1807993) q[1];
sx q[1];
rz(-1.5174148) q[1];
sx q[1];
rz(-1.175715) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.116947) q[3];
sx q[3];
rz(-0.99089115) q[3];
sx q[3];
rz(-0.49296311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20554081) q[2];
sx q[2];
rz(-2.3369393) q[2];
sx q[2];
rz(-1.12977) q[2];
rz(2.1352844) q[3];
sx q[3];
rz(-1.8215424) q[3];
sx q[3];
rz(1.4973076) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0928918) q[0];
sx q[0];
rz(-2.0755656) q[0];
sx q[0];
rz(2.997828) q[0];
rz(0.20214209) q[1];
sx q[1];
rz(-2.6136687) q[1];
sx q[1];
rz(1.082083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7162299) q[0];
sx q[0];
rz(-0.83304616) q[0];
sx q[0];
rz(2.0701755) q[0];
rz(-pi) q[1];
rz(-2.4074102) q[2];
sx q[2];
rz(-2.1819644) q[2];
sx q[2];
rz(-2.2847459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94596568) q[1];
sx q[1];
rz(-2.2711427) q[1];
sx q[1];
rz(-0.27865748) q[1];
rz(-0.96636678) q[3];
sx q[3];
rz(-2.2141738) q[3];
sx q[3];
rz(-0.21277393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14273345) q[2];
sx q[2];
rz(-1.961668) q[2];
sx q[2];
rz(0.72706968) q[2];
rz(-1.8266228) q[3];
sx q[3];
rz(-1.8948137) q[3];
sx q[3];
rz(1.6453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975248) q[0];
sx q[0];
rz(-2.0638564) q[0];
sx q[0];
rz(1.9986073) q[0];
rz(0.22363981) q[1];
sx q[1];
rz(-0.63839212) q[1];
sx q[1];
rz(-1.9167831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1849104) q[0];
sx q[0];
rz(-1.3089797) q[0];
sx q[0];
rz(1.8641406) q[0];
x q[1];
rz(-1.9838422) q[2];
sx q[2];
rz(-2.2487679) q[2];
sx q[2];
rz(1.3785481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30149192) q[1];
sx q[1];
rz(-0.21243851) q[1];
sx q[1];
rz(1.849624) q[1];
rz(-1.4644044) q[3];
sx q[3];
rz(-1.0118183) q[3];
sx q[3];
rz(-0.91181741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58330047) q[2];
sx q[2];
rz(-1.4536828) q[2];
sx q[2];
rz(2.3395786) q[2];
rz(-1.924104) q[3];
sx q[3];
rz(-0.13974443) q[3];
sx q[3];
rz(1.8961934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8121346) q[0];
sx q[0];
rz(-1.7725002) q[0];
sx q[0];
rz(1.5611956) q[0];
rz(0.87248692) q[1];
sx q[1];
rz(-2.8552738) q[1];
sx q[1];
rz(0.40503851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0319388) q[0];
sx q[0];
rz(-0.73470107) q[0];
sx q[0];
rz(0.74560179) q[0];
rz(-2.9895859) q[2];
sx q[2];
rz(-2.5436253) q[2];
sx q[2];
rz(-2.1961308) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26012684) q[1];
sx q[1];
rz(-2.2592789) q[1];
sx q[1];
rz(0.75551178) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14696481) q[3];
sx q[3];
rz(-1.5515447) q[3];
sx q[3];
rz(-1.4052773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27434719) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(-0.61152968) q[2];
rz(-1.3036171) q[3];
sx q[3];
rz(-2.7795064) q[3];
sx q[3];
rz(2.005579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63854727) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(-3.0882623) q[0];
rz(2.5697925) q[1];
sx q[1];
rz(-1.5170393) q[1];
sx q[1];
rz(-1.2791876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17634468) q[0];
sx q[0];
rz(-1.2851479) q[0];
sx q[0];
rz(1.833451) q[0];
rz(-pi) q[1];
rz(2.9933571) q[2];
sx q[2];
rz(-1.094813) q[2];
sx q[2];
rz(-2.2479284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51579976) q[1];
sx q[1];
rz(-0.98550057) q[1];
sx q[1];
rz(-0.093823508) q[1];
rz(-pi) q[2];
rz(-1.199715) q[3];
sx q[3];
rz(-1.367192) q[3];
sx q[3];
rz(2.9534087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6751487) q[2];
sx q[2];
rz(-1.6089336) q[2];
sx q[2];
rz(-2.4565728) q[2];
rz(-2.3576665) q[3];
sx q[3];
rz(-2.7214366) q[3];
sx q[3];
rz(2.4043731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6017799) q[0];
sx q[0];
rz(-1.2417326) q[0];
sx q[0];
rz(-1.7057521) q[0];
rz(-2.336179) q[1];
sx q[1];
rz(-2.6782811) q[1];
sx q[1];
rz(-2.4334999) q[1];
rz(-0.28975363) q[2];
sx q[2];
rz(-1.3633681) q[2];
sx q[2];
rz(1.0492556) q[2];
rz(-0.41498923) q[3];
sx q[3];
rz(-1.6798974) q[3];
sx q[3];
rz(0.85174042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
