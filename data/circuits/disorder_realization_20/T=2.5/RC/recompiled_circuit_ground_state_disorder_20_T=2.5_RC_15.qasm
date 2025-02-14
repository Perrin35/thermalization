OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29937509) q[0];
sx q[0];
rz(-2.8111281) q[0];
sx q[0];
rz(2.0781031) q[0];
rz(-0.039634135) q[1];
sx q[1];
rz(-0.57365817) q[1];
sx q[1];
rz(-1.080245) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40908694) q[0];
sx q[0];
rz(-1.5729781) q[0];
sx q[0];
rz(2.0337142) q[0];
rz(-pi) q[1];
rz(-1.3990381) q[2];
sx q[2];
rz(-1.556201) q[2];
sx q[2];
rz(0.43806048) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4689881) q[1];
sx q[1];
rz(-2.0867587) q[1];
sx q[1];
rz(-1.8700061) q[1];
rz(-0.61038252) q[3];
sx q[3];
rz(-0.74225194) q[3];
sx q[3];
rz(0.28377747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1178736) q[2];
sx q[2];
rz(-1.4154075) q[2];
sx q[2];
rz(0.37303698) q[2];
rz(-2.5840664) q[3];
sx q[3];
rz(-1.5159461) q[3];
sx q[3];
rz(-2.663234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8525304) q[0];
sx q[0];
rz(-1.7076778) q[0];
sx q[0];
rz(1.0093932) q[0];
rz(-0.32762647) q[1];
sx q[1];
rz(-2.8150924) q[1];
sx q[1];
rz(1.1064628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43531159) q[0];
sx q[0];
rz(-1.5967622) q[0];
sx q[0];
rz(3.0089799) q[0];
x q[1];
rz(-2.4566133) q[2];
sx q[2];
rz(-1.5589336) q[2];
sx q[2];
rz(1.4432743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49010897) q[1];
sx q[1];
rz(-2.1613908) q[1];
sx q[1];
rz(-2.5701853) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42454607) q[3];
sx q[3];
rz(-2.694482) q[3];
sx q[3];
rz(-3.0101484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24836765) q[2];
sx q[2];
rz(-1.1838341) q[2];
sx q[2];
rz(0.66967213) q[2];
rz(-1.1559961) q[3];
sx q[3];
rz(-0.35231927) q[3];
sx q[3];
rz(2.7922351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30592331) q[0];
sx q[0];
rz(-0.42757973) q[0];
sx q[0];
rz(-1.3315573) q[0];
rz(1.2003027) q[1];
sx q[1];
rz(-0.2526865) q[1];
sx q[1];
rz(0.67218626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5066619) q[0];
sx q[0];
rz(-1.4499393) q[0];
sx q[0];
rz(-0.47516993) q[0];
rz(-pi) q[1];
rz(-1.1665909) q[2];
sx q[2];
rz(-1.4160755) q[2];
sx q[2];
rz(-0.2153309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3954288) q[1];
sx q[1];
rz(-2.6708467) q[1];
sx q[1];
rz(-2.7676299) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5198042) q[3];
sx q[3];
rz(-1.4740406) q[3];
sx q[3];
rz(-0.51550625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4868698) q[2];
sx q[2];
rz(-1.1785616) q[2];
sx q[2];
rz(-0.85050026) q[2];
rz(2.1908098) q[3];
sx q[3];
rz(-2.6362004) q[3];
sx q[3];
rz(0.90184414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.013407) q[0];
sx q[0];
rz(-2.3222017) q[0];
sx q[0];
rz(0.73981458) q[0];
rz(0.20135227) q[1];
sx q[1];
rz(-1.9722152) q[1];
sx q[1];
rz(0.22612017) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126919) q[0];
sx q[0];
rz(-1.5357112) q[0];
sx q[0];
rz(3.1265774) q[0];
rz(0.32711012) q[2];
sx q[2];
rz(-1.2108821) q[2];
sx q[2];
rz(1.2963595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0579083) q[1];
sx q[1];
rz(-1.0261361) q[1];
sx q[1];
rz(-0.96353957) q[1];
x q[2];
rz(1.4410517) q[3];
sx q[3];
rz(-1.1888388) q[3];
sx q[3];
rz(2.6568535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40253338) q[2];
sx q[2];
rz(-0.89503521) q[2];
sx q[2];
rz(-2.1441937) q[2];
rz(-2.7282257) q[3];
sx q[3];
rz(-1.9316542) q[3];
sx q[3];
rz(2.0972142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1362374) q[0];
sx q[0];
rz(-2.4449466) q[0];
sx q[0];
rz(-0.053939017) q[0];
rz(0.18221642) q[1];
sx q[1];
rz(-2.2257664) q[1];
sx q[1];
rz(0.30311662) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5678068) q[0];
sx q[0];
rz(-1.2284096) q[0];
sx q[0];
rz(-2.0339147) q[0];
rz(-pi) q[1];
rz(1.8218173) q[2];
sx q[2];
rz(-0.66749882) q[2];
sx q[2];
rz(2.4199744) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.388465) q[1];
sx q[1];
rz(-1.7645986) q[1];
sx q[1];
rz(1.0419921) q[1];
x q[2];
rz(-2.3833586) q[3];
sx q[3];
rz(-1.5966187) q[3];
sx q[3];
rz(1.7968221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4431241) q[2];
sx q[2];
rz(-1.8734525) q[2];
sx q[2];
rz(0.35688409) q[2];
rz(2.3667864) q[3];
sx q[3];
rz(-1.509343) q[3];
sx q[3];
rz(0.43865144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.2543432) q[0];
sx q[0];
rz(-1.9312504) q[0];
sx q[0];
rz(-0.78897011) q[0];
rz(1.7987159) q[1];
sx q[1];
rz(-1.7306381) q[1];
sx q[1];
rz(-0.79997396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8784284) q[0];
sx q[0];
rz(-1.5845926) q[0];
sx q[0];
rz(-0.8850125) q[0];
rz(-pi) q[1];
rz(0.87514295) q[2];
sx q[2];
rz(-1.4261386) q[2];
sx q[2];
rz(0.62334594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9467873) q[1];
sx q[1];
rz(-1.9459263) q[1];
sx q[1];
rz(-1.4542143) q[1];
rz(-0.39928945) q[3];
sx q[3];
rz(-1.1762816) q[3];
sx q[3];
rz(2.7189915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0903025) q[2];
sx q[2];
rz(-2.3919969) q[2];
sx q[2];
rz(1.7898111) q[2];
rz(-0.46818647) q[3];
sx q[3];
rz(-2.0464094) q[3];
sx q[3];
rz(-0.93402544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139597) q[0];
sx q[0];
rz(-0.41828823) q[0];
sx q[0];
rz(3.1391414) q[0];
rz(3.1353503) q[1];
sx q[1];
rz(-0.36224449) q[1];
sx q[1];
rz(-0.35591602) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2438176) q[0];
sx q[0];
rz(-0.96976244) q[0];
sx q[0];
rz(0.42341451) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4202815) q[2];
sx q[2];
rz(-1.7928267) q[2];
sx q[2];
rz(-0.49563615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1568875) q[1];
sx q[1];
rz(-2.2199989) q[1];
sx q[1];
rz(-2.8125416) q[1];
rz(-pi) q[2];
rz(-1.6735092) q[3];
sx q[3];
rz(-1.0389162) q[3];
sx q[3];
rz(-1.892688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9956841) q[2];
sx q[2];
rz(-1.2331139) q[2];
sx q[2];
rz(1.0710867) q[2];
rz(0.60761014) q[3];
sx q[3];
rz(-1.1039609) q[3];
sx q[3];
rz(-1.6149909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0869658) q[0];
sx q[0];
rz(-2.2038951) q[0];
sx q[0];
rz(2.0077534) q[0];
rz(1.0519823) q[1];
sx q[1];
rz(-1.6832422) q[1];
sx q[1];
rz(0.7410616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0042564226) q[0];
sx q[0];
rz(-1.8305873) q[0];
sx q[0];
rz(0.57704837) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95693077) q[2];
sx q[2];
rz(-1.1476048) q[2];
sx q[2];
rz(0.6415216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7692803) q[1];
sx q[1];
rz(-1.0456632) q[1];
sx q[1];
rz(1.5985846) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2005733) q[3];
sx q[3];
rz(-1.9004603) q[3];
sx q[3];
rz(-0.31808149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97003254) q[2];
sx q[2];
rz(-0.67983183) q[2];
sx q[2];
rz(2.6762834) q[2];
rz(0.71632898) q[3];
sx q[3];
rz(-1.497044) q[3];
sx q[3];
rz(-1.1819476) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61140907) q[0];
sx q[0];
rz(-2.1937328) q[0];
sx q[0];
rz(1.9635669) q[0];
rz(-2.1317962) q[1];
sx q[1];
rz(-1.3346846) q[1];
sx q[1];
rz(3.03426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5553805) q[0];
sx q[0];
rz(-0.79609603) q[0];
sx q[0];
rz(-0.70527161) q[0];
x q[1];
rz(-0.1520098) q[2];
sx q[2];
rz(-0.78897023) q[2];
sx q[2];
rz(1.1107572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91574749) q[1];
sx q[1];
rz(-1.3406995) q[1];
sx q[1];
rz(0.6055931) q[1];
rz(-pi) q[2];
rz(2.168591) q[3];
sx q[3];
rz(-1.8174371) q[3];
sx q[3];
rz(-1.1688978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44059077) q[2];
sx q[2];
rz(-0.89459449) q[2];
sx q[2];
rz(0.44006285) q[2];
rz(2.9317686) q[3];
sx q[3];
rz(-1.9473636) q[3];
sx q[3];
rz(-0.31407174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99301991) q[0];
sx q[0];
rz(-2.5555389) q[0];
sx q[0];
rz(-1.7589737) q[0];
rz(2.4644409) q[1];
sx q[1];
rz(-1.6771728) q[1];
sx q[1];
rz(1.1322397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234194) q[0];
sx q[0];
rz(-1.7935402) q[0];
sx q[0];
rz(1.3370418) q[0];
x q[1];
rz(-1.8381565) q[2];
sx q[2];
rz(-2.3472381) q[2];
sx q[2];
rz(-1.8092138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59076485) q[1];
sx q[1];
rz(-2.4240997) q[1];
sx q[1];
rz(3.1050502) q[1];
rz(-pi) q[2];
rz(-2.4664425) q[3];
sx q[3];
rz(-1.0680001) q[3];
sx q[3];
rz(-2.4166963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18572346) q[2];
sx q[2];
rz(-0.74803868) q[2];
sx q[2];
rz(-1.7175187) q[2];
rz(0.87797034) q[3];
sx q[3];
rz(-2.9214171) q[3];
sx q[3];
rz(-3.1233845) q[3];
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
rz(-1.8743185) q[0];
sx q[0];
rz(-1.6813288) q[0];
sx q[0];
rz(1.9718476) q[0];
rz(0.36880233) q[1];
sx q[1];
rz(-1.7316876) q[1];
sx q[1];
rz(0.015451886) q[1];
rz(-0.56024341) q[2];
sx q[2];
rz(-2.6340953) q[2];
sx q[2];
rz(-2.5731186) q[2];
rz(0.54184171) q[3];
sx q[3];
rz(-2.4892157) q[3];
sx q[3];
rz(2.6996725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
