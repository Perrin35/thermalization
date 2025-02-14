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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(-2.5007024) q[0];
rz(2.9042397) q[1];
sx q[1];
rz(-0.64164716) q[1];
sx q[1];
rz(-3.0566888) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083402) q[0];
sx q[0];
rz(-1.6680934) q[0];
sx q[0];
rz(0.86546524) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48975421) q[2];
sx q[2];
rz(-0.066250525) q[2];
sx q[2];
rz(0.79023933) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0247268) q[1];
sx q[1];
rz(-1.6849396) q[1];
sx q[1];
rz(0.29205871) q[1];
x q[2];
rz(-2.1259948) q[3];
sx q[3];
rz(-2.4732504) q[3];
sx q[3];
rz(-0.67977842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.889633) q[2];
sx q[2];
rz(-1.9388988) q[2];
sx q[2];
rz(1.6331875) q[2];
rz(0.83186045) q[3];
sx q[3];
rz(-2.8155477) q[3];
sx q[3];
rz(-1.7970201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93083301) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(0.26500901) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(1.7204684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20939556) q[0];
sx q[0];
rz(-1.3065814) q[0];
sx q[0];
rz(-0.90640675) q[0];
rz(-pi) q[1];
rz(2.1900329) q[2];
sx q[2];
rz(-2.4251221) q[2];
sx q[2];
rz(-1.2576511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4231734) q[1];
sx q[1];
rz(-1.3191901) q[1];
sx q[1];
rz(-3.0965641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3953929) q[3];
sx q[3];
rz(-2.9713745) q[3];
sx q[3];
rz(-1.1268009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77218324) q[2];
sx q[2];
rz(-2.1464244) q[2];
sx q[2];
rz(0.65917242) q[2];
rz(0.92165056) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(-1.5825745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42940656) q[0];
sx q[0];
rz(-0.61539188) q[0];
sx q[0];
rz(1.8581101) q[0];
rz(-0.42690024) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.9740419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1438863) q[0];
sx q[0];
rz(-2.6397815) q[0];
sx q[0];
rz(-0.83117475) q[0];
rz(-pi) q[1];
rz(-2.6111772) q[2];
sx q[2];
rz(-1.1258719) q[2];
sx q[2];
rz(-2.4702768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0223654) q[1];
sx q[1];
rz(-1.664838) q[1];
sx q[1];
rz(-2.8507807) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5646788) q[3];
sx q[3];
rz(-2.8261999) q[3];
sx q[3];
rz(-1.9384991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9350962) q[2];
sx q[2];
rz(-1.6556211) q[2];
sx q[2];
rz(3.0237954) q[2];
rz(-1.3055034) q[3];
sx q[3];
rz(-0.85956231) q[3];
sx q[3];
rz(-0.99087244) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30462581) q[0];
sx q[0];
rz(-2.0898297) q[0];
sx q[0];
rz(0.032531746) q[0];
rz(-2.1578728) q[1];
sx q[1];
rz(-2.3301221) q[1];
sx q[1];
rz(0.54289114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8642917) q[0];
sx q[0];
rz(-1.2980718) q[0];
sx q[0];
rz(2.4433873) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2603733) q[2];
sx q[2];
rz(-0.32698787) q[2];
sx q[2];
rz(0.56337683) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.429814) q[1];
sx q[1];
rz(-2.3444054) q[1];
sx q[1];
rz(-1.2120257) q[1];
rz(1.4063605) q[3];
sx q[3];
rz(-0.74526605) q[3];
sx q[3];
rz(2.0045053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5618318) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(-0.96060166) q[2];
rz(-0.97551712) q[3];
sx q[3];
rz(-0.87149182) q[3];
sx q[3];
rz(0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8303216) q[0];
sx q[0];
rz(-2.4595342) q[0];
sx q[0];
rz(0.3325381) q[0];
rz(-0.63811103) q[1];
sx q[1];
rz(-0.6256012) q[1];
sx q[1];
rz(0.45613751) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77959427) q[0];
sx q[0];
rz(-1.5080305) q[0];
sx q[0];
rz(-1.4383352) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0917408) q[2];
sx q[2];
rz(-2.2992813) q[2];
sx q[2];
rz(3.0974744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0597555) q[1];
sx q[1];
rz(-2.55099) q[1];
sx q[1];
rz(-1.4103622) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1061323) q[3];
sx q[3];
rz(-2.439508) q[3];
sx q[3];
rz(2.9891069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0786324) q[2];
sx q[2];
rz(-2.3802064) q[2];
sx q[2];
rz(2.6893993) q[2];
rz(-0.26442987) q[3];
sx q[3];
rz(-0.14614883) q[3];
sx q[3];
rz(-1.2368081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0753157) q[0];
sx q[0];
rz(-2.3105268) q[0];
sx q[0];
rz(0.7178632) q[0];
rz(1.588899) q[1];
sx q[1];
rz(-1.6009067) q[1];
sx q[1];
rz(-0.20725651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22482797) q[0];
sx q[0];
rz(-1.5434573) q[0];
sx q[0];
rz(0.44204373) q[0];
rz(-1.147993) q[2];
sx q[2];
rz(-1.5059274) q[2];
sx q[2];
rz(2.1366304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9178324) q[1];
sx q[1];
rz(-0.66559911) q[1];
sx q[1];
rz(-0.83719801) q[1];
rz(-pi) q[2];
rz(1.8317811) q[3];
sx q[3];
rz(-1.4647604) q[3];
sx q[3];
rz(2.1539254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76102105) q[2];
sx q[2];
rz(-1.0552768) q[2];
sx q[2];
rz(-2.5804248) q[2];
rz(-1.0593972) q[3];
sx q[3];
rz(-1.1166409) q[3];
sx q[3];
rz(-1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.2137432) q[0];
sx q[0];
rz(-2.1708467) q[0];
sx q[0];
rz(-2.2947776) q[0];
rz(2.1532586) q[1];
sx q[1];
rz(-1.3761995) q[1];
sx q[1];
rz(0.83289897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5079437) q[0];
sx q[0];
rz(-1.4671774) q[0];
sx q[0];
rz(-0.1990464) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9994644) q[2];
sx q[2];
rz(-1.649039) q[2];
sx q[2];
rz(2.4852666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78645241) q[1];
sx q[1];
rz(-1.9398098) q[1];
sx q[1];
rz(-0.92496101) q[1];
x q[2];
rz(3.0924721) q[3];
sx q[3];
rz(-1.4572772) q[3];
sx q[3];
rz(3.1044131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6758468) q[2];
sx q[2];
rz(-0.70285672) q[2];
sx q[2];
rz(-2.3715026) q[2];
rz(0.13104023) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(1.8946764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34762621) q[0];
sx q[0];
rz(-0.93526953) q[0];
sx q[0];
rz(0.4916218) q[0];
rz(-2.8288016) q[1];
sx q[1];
rz(-0.28952315) q[1];
sx q[1];
rz(3.0591931) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9461878) q[0];
sx q[0];
rz(-0.090443693) q[0];
sx q[0];
rz(-1.7121332) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5748207) q[2];
sx q[2];
rz(-0.9563891) q[2];
sx q[2];
rz(2.4943309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1090549) q[1];
sx q[1];
rz(-1.7454285) q[1];
sx q[1];
rz(-2.9777378) q[1];
x q[2];
rz(1.2639323) q[3];
sx q[3];
rz(-0.32684775) q[3];
sx q[3];
rz(0.89034789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6975434) q[2];
sx q[2];
rz(-1.9335582) q[2];
sx q[2];
rz(0.018012878) q[2];
rz(-2.2703914) q[3];
sx q[3];
rz(-2.8049991) q[3];
sx q[3];
rz(-0.66756162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1910601) q[0];
sx q[0];
rz(-2.58707) q[0];
sx q[0];
rz(-0.46448034) q[0];
rz(-0.49098268) q[1];
sx q[1];
rz(-1.6698488) q[1];
sx q[1];
rz(1.6741265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6086836) q[0];
sx q[0];
rz(-0.94588806) q[0];
sx q[0];
rz(1.5707569) q[0];
x q[1];
rz(1.2651132) q[2];
sx q[2];
rz(-1.0846429) q[2];
sx q[2];
rz(-0.4403688) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75338335) q[1];
sx q[1];
rz(-0.25553726) q[1];
sx q[1];
rz(1.3292043) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76045658) q[3];
sx q[3];
rz(-2.1158614) q[3];
sx q[3];
rz(2.2732002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1269425) q[2];
sx q[2];
rz(-1.4913968) q[2];
sx q[2];
rz(-1.8102185) q[2];
rz(2.2018382) q[3];
sx q[3];
rz(-1.0569812) q[3];
sx q[3];
rz(1.1872928) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2347655) q[0];
sx q[0];
rz(-2.6563788) q[0];
sx q[0];
rz(2.1038798) q[0];
rz(1.0543793) q[1];
sx q[1];
rz(-1.3155921) q[1];
sx q[1];
rz(1.0494999) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18436156) q[0];
sx q[0];
rz(-1.6205233) q[0];
sx q[0];
rz(0.00018361365) q[0];
rz(2.5811701) q[2];
sx q[2];
rz(-1.9935521) q[2];
sx q[2];
rz(2.1566856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8912208) q[1];
sx q[1];
rz(-2.1714021) q[1];
sx q[1];
rz(-0.037446219) q[1];
rz(-0.78759362) q[3];
sx q[3];
rz(-2.3830722) q[3];
sx q[3];
rz(-2.1473928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1513169) q[2];
sx q[2];
rz(-1.8345202) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1405519) q[0];
sx q[0];
rz(-0.97698553) q[0];
sx q[0];
rz(-2.3401596) q[0];
rz(-1.6544381) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(-0.64143388) q[2];
sx q[2];
rz(-2.1830255) q[2];
sx q[2];
rz(0.39354702) q[2];
rz(0.63538649) q[3];
sx q[3];
rz(-1.0839331) q[3];
sx q[3];
rz(-1.4270368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
