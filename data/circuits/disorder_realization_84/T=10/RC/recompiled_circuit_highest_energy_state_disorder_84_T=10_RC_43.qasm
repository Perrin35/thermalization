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
rz(0.52047211) q[0];
sx q[0];
rz(-2.1055129) q[0];
sx q[0];
rz(2.5007024) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(-3.0566888) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215243) q[0];
sx q[0];
rz(-0.86948538) q[0];
sx q[0];
rz(-3.0140956) q[0];
rz(-1.5395959) q[2];
sx q[2];
rz(-1.5123431) q[2];
sx q[2];
rz(-1.8606869) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7218923) q[1];
sx q[1];
rz(-1.2806935) q[1];
sx q[1];
rz(1.4516524) q[1];
rz(-2.1617684) q[3];
sx q[3];
rz(-1.238046) q[3];
sx q[3];
rz(2.7036288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.889633) q[2];
sx q[2];
rz(-1.2026938) q[2];
sx q[2];
rz(-1.6331875) q[2];
rz(-0.83186045) q[3];
sx q[3];
rz(-2.8155477) q[3];
sx q[3];
rz(1.7970201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93083301) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(-0.26500901) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(-1.4211242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9321971) q[0];
sx q[0];
rz(-1.8350112) q[0];
sx q[0];
rz(2.2351859) q[0];
rz(-pi) q[1];
rz(-2.6735953) q[2];
sx q[2];
rz(-2.1350522) q[2];
sx q[2];
rz(-2.6411438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4231734) q[1];
sx q[1];
rz(-1.8224026) q[1];
sx q[1];
rz(-3.0965641) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1116074) q[3];
sx q[3];
rz(-1.4032149) q[3];
sx q[3];
rz(0.94887892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77218324) q[2];
sx q[2];
rz(-0.99516827) q[2];
sx q[2];
rz(-0.65917242) q[2];
rz(2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(1.5825745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.42940656) q[0];
sx q[0];
rz(-0.61539188) q[0];
sx q[0];
rz(-1.8581101) q[0];
rz(-0.42690024) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.9740419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1438863) q[0];
sx q[0];
rz(-2.6397815) q[0];
sx q[0];
rz(0.83117475) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0658355) q[2];
sx q[2];
rz(-1.0965818) q[2];
sx q[2];
rz(1.1467288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0223654) q[1];
sx q[1];
rz(-1.4767547) q[1];
sx q[1];
rz(-0.29081197) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5646788) q[3];
sx q[3];
rz(-0.31539279) q[3];
sx q[3];
rz(-1.2030935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9350962) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(0.11779724) q[2];
rz(1.3055034) q[3];
sx q[3];
rz(-2.2820303) q[3];
sx q[3];
rz(2.1507202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369668) q[0];
sx q[0];
rz(-1.0517629) q[0];
sx q[0];
rz(-0.032531746) q[0];
rz(0.98371983) q[1];
sx q[1];
rz(-2.3301221) q[1];
sx q[1];
rz(-2.5987015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8642917) q[0];
sx q[0];
rz(-1.8435209) q[0];
sx q[0];
rz(2.4433873) q[0];
rz(-3.0383598) q[2];
sx q[2];
rz(-1.2599753) q[2];
sx q[2];
rz(-0.23676591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9223843) q[1];
sx q[1];
rz(-0.83679779) q[1];
sx q[1];
rz(0.34511415) q[1];
rz(-pi) q[2];
rz(0.832295) q[3];
sx q[3];
rz(-1.6820388) q[3];
sx q[3];
rz(2.5865366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5618318) q[2];
sx q[2];
rz(-2.0115972) q[2];
sx q[2];
rz(-2.180991) q[2];
rz(-2.1660755) q[3];
sx q[3];
rz(-2.2701008) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3112711) q[0];
sx q[0];
rz(-0.68205849) q[0];
sx q[0];
rz(-0.3325381) q[0];
rz(-0.63811103) q[1];
sx q[1];
rz(-2.5159914) q[1];
sx q[1];
rz(-0.45613751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35121954) q[0];
sx q[0];
rz(-2.995092) q[0];
sx q[0];
rz(1.1266493) q[0];
rz(0.47686968) q[2];
sx q[2];
rz(-0.84706351) q[2];
sx q[2];
rz(0.70666955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8673381) q[1];
sx q[1];
rz(-0.98877866) q[1];
sx q[1];
rz(3.0349005) q[1];
rz(-1.1061323) q[3];
sx q[3];
rz(-2.439508) q[3];
sx q[3];
rz(-0.15248577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0629603) q[2];
sx q[2];
rz(-0.76138622) q[2];
sx q[2];
rz(0.45219335) q[2];
rz(2.8771628) q[3];
sx q[3];
rz(-0.14614883) q[3];
sx q[3];
rz(1.9047846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0662769) q[0];
sx q[0];
rz(-2.3105268) q[0];
sx q[0];
rz(2.4237295) q[0];
rz(1.5526937) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(0.20725651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2882746) q[0];
sx q[0];
rz(-2.69876) q[0];
sx q[0];
rz(-0.063837039) q[0];
rz(1.7278094) q[2];
sx q[2];
rz(-0.42745379) q[2];
sx q[2];
rz(-2.4326716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0644857) q[1];
sx q[1];
rz(-2.047309) q[1];
sx q[1];
rz(-0.48397343) q[1];
x q[2];
rz(-1.1795711) q[3];
sx q[3];
rz(-0.28124725) q[3];
sx q[3];
rz(0.96042552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3805716) q[2];
sx q[2];
rz(-1.0552768) q[2];
sx q[2];
rz(2.5804248) q[2];
rz(-1.0593972) q[3];
sx q[3];
rz(-2.0249517) q[3];
sx q[3];
rz(1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(-0.84681502) q[0];
rz(-2.1532586) q[1];
sx q[1];
rz(-1.3761995) q[1];
sx q[1];
rz(2.3086937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4109548) q[0];
sx q[0];
rz(-2.9175076) q[0];
sx q[0];
rz(0.48416324) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6360699) q[2];
sx q[2];
rz(-2.9794783) q[2];
sx q[2];
rz(-2.7270728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0918232) q[1];
sx q[1];
rz(-2.1668129) q[1];
sx q[1];
rz(0.45097643) q[1];
x q[2];
rz(-3.0924721) q[3];
sx q[3];
rz(-1.4572772) q[3];
sx q[3];
rz(0.037179557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46574584) q[2];
sx q[2];
rz(-2.4387359) q[2];
sx q[2];
rz(0.77009002) q[2];
rz(-0.13104023) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(1.2469163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.34762621) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(2.6499709) q[0];
rz(-2.8288016) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(-0.082399592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51615825) q[0];
sx q[0];
rz(-1.5580728) q[0];
sx q[0];
rz(1.660343) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5748207) q[2];
sx q[2];
rz(-0.9563891) q[2];
sx q[2];
rz(-2.4943309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56697733) q[1];
sx q[1];
rz(-1.4094556) q[1];
sx q[1];
rz(-1.7477504) q[1];
x q[2];
rz(1.8776603) q[3];
sx q[3];
rz(-0.32684775) q[3];
sx q[3];
rz(-0.89034789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4440492) q[2];
sx q[2];
rz(-1.2080344) q[2];
sx q[2];
rz(3.1235798) q[2];
rz(-0.87120122) q[3];
sx q[3];
rz(-0.33659354) q[3];
sx q[3];
rz(-0.66756162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910601) q[0];
sx q[0];
rz(-2.58707) q[0];
sx q[0];
rz(0.46448034) q[0];
rz(-0.49098268) q[1];
sx q[1];
rz(-1.6698488) q[1];
sx q[1];
rz(1.6741265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6086162) q[0];
sx q[0];
rz(-2.5166844) q[0];
sx q[0];
rz(-3.141538) q[0];
x q[1];
rz(2.6355714) q[2];
sx q[2];
rz(-1.8401166) q[2];
sx q[2];
rz(1.276818) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0027568) q[1];
sx q[1];
rz(-1.3228387) q[1];
sx q[1];
rz(0.062422189) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4201066) q[3];
sx q[3];
rz(-2.2389447) q[3];
sx q[3];
rz(0.20340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.014650194) q[2];
sx q[2];
rz(-1.6501959) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(-0.93975449) q[3];
sx q[3];
rz(-2.0846114) q[3];
sx q[3];
rz(-1.1872928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9068271) q[0];
sx q[0];
rz(-2.6563788) q[0];
sx q[0];
rz(2.1038798) q[0];
rz(-2.0872133) q[1];
sx q[1];
rz(-1.3155921) q[1];
sx q[1];
rz(1.0494999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18436156) q[0];
sx q[0];
rz(-1.6205233) q[0];
sx q[0];
rz(-3.141409) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56042258) q[2];
sx q[2];
rz(-1.9935521) q[2];
sx q[2];
rz(2.1566856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8423374) q[1];
sx q[1];
rz(-1.5399057) q[1];
sx q[1];
rz(2.1717291) q[1];
x q[2];
rz(-0.97941937) q[3];
sx q[3];
rz(-2.0775177) q[3];
sx q[3];
rz(-3.0924071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99027571) q[2];
sx q[2];
rz(-1.3070725) q[2];
sx q[2];
rz(-1.446373) q[2];
rz(-2.9826048) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(-0.75001636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0010407) q[0];
sx q[0];
rz(-2.1646071) q[0];
sx q[0];
rz(0.80143308) q[0];
rz(1.6544381) q[1];
sx q[1];
rz(-0.14688891) q[1];
sx q[1];
rz(2.6085703) q[1];
rz(0.85114807) q[2];
sx q[2];
rz(-1.0590886) q[2];
sx q[2];
rz(1.5589489) q[2];
rz(0.98900908) q[3];
sx q[3];
rz(-1.0186356) q[3];
sx q[3];
rz(2.9531425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
