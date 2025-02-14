OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3684664) q[0];
sx q[0];
rz(-2.3946895) q[0];
sx q[0];
rz(-0.84063831) q[0];
rz(-3.0199938) q[1];
sx q[1];
rz(-1.8687948) q[1];
sx q[1];
rz(2.8425541) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6071607) q[0];
sx q[0];
rz(-0.91437712) q[0];
sx q[0];
rz(-0.54404152) q[0];
rz(-2.0571124) q[2];
sx q[2];
rz(-1.63967) q[2];
sx q[2];
rz(1.8191847) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4579276) q[1];
sx q[1];
rz(-1.8103765) q[1];
sx q[1];
rz(-1.7182106) q[1];
rz(2.0255346) q[3];
sx q[3];
rz(-1.699243) q[3];
sx q[3];
rz(-1.0806395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76089871) q[2];
sx q[2];
rz(-1.0502522) q[2];
sx q[2];
rz(-2.8139581) q[2];
rz(1.3753752) q[3];
sx q[3];
rz(-1.7211434) q[3];
sx q[3];
rz(0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37680092) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(-2.2862527) q[0];
rz(-0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(-0.45113742) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83734918) q[0];
sx q[0];
rz(-0.74241246) q[0];
sx q[0];
rz(-1.8038595) q[0];
rz(-1.3786267) q[2];
sx q[2];
rz(-0.8079257) q[2];
sx q[2];
rz(-2.5469123) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0939744) q[1];
sx q[1];
rz(-0.4741569) q[1];
sx q[1];
rz(2.7535901) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1633918) q[3];
sx q[3];
rz(-1.1278858) q[3];
sx q[3];
rz(2.9748084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1796639) q[2];
sx q[2];
rz(-1.1397866) q[2];
sx q[2];
rz(-0.12953225) q[2];
rz(-0.1772964) q[3];
sx q[3];
rz(-2.5640021) q[3];
sx q[3];
rz(0.052848335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.391908) q[0];
sx q[0];
rz(-2.7295697) q[0];
sx q[0];
rz(0.55927292) q[0];
rz(-3.0468805) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(-0.47725484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1806948) q[0];
sx q[0];
rz(-1.8261693) q[0];
sx q[0];
rz(1.9810505) q[0];
x q[1];
rz(2.7626286) q[2];
sx q[2];
rz(-0.38658374) q[2];
sx q[2];
rz(-1.5075589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86384799) q[1];
sx q[1];
rz(-1.7949634) q[1];
sx q[1];
rz(0.78296354) q[1];
rz(-2.6179254) q[3];
sx q[3];
rz(-0.91991495) q[3];
sx q[3];
rz(-1.4613348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(1.8837455) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(-3.0959082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78087085) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(-2.1029396) q[0];
rz(0.61141283) q[1];
sx q[1];
rz(-0.88714209) q[1];
sx q[1];
rz(-0.92322737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34133615) q[0];
sx q[0];
rz(-1.0025327) q[0];
sx q[0];
rz(0.41109127) q[0];
rz(0.011767894) q[2];
sx q[2];
rz(-1.9885049) q[2];
sx q[2];
rz(-2.2209446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23455305) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(-1.6978463) q[1];
x q[2];
rz(1.3096894) q[3];
sx q[3];
rz(-0.93702836) q[3];
sx q[3];
rz(-1.5719617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-0.46773657) q[2];
sx q[2];
rz(2.5941217) q[2];
rz(0.064420961) q[3];
sx q[3];
rz(-1.3037325) q[3];
sx q[3];
rz(2.2678383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15544686) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(1.4601532) q[0];
rz(-2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(-0.11988457) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.368705) q[0];
sx q[0];
rz(-1.3471706) q[0];
sx q[0];
rz(-2.2740433) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5566543) q[2];
sx q[2];
rz(-2.1368933) q[2];
sx q[2];
rz(-1.0370129) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2561349) q[1];
sx q[1];
rz(-1.3559582) q[1];
sx q[1];
rz(-1.4900521) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0644887) q[3];
sx q[3];
rz(-2.1112006) q[3];
sx q[3];
rz(-0.04549724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.039006058) q[2];
sx q[2];
rz(-2.234499) q[2];
sx q[2];
rz(-0.26068035) q[2];
rz(2.2634704) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(-0.78316435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53073019) q[0];
sx q[0];
rz(-3.0513638) q[0];
sx q[0];
rz(-1.0020142) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(0.9446876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1169918) q[0];
sx q[0];
rz(-0.62378609) q[0];
sx q[0];
rz(-0.37895112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5414921) q[2];
sx q[2];
rz(-0.99008152) q[2];
sx q[2];
rz(2.3490259) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73026555) q[1];
sx q[1];
rz(-1.7217727) q[1];
sx q[1];
rz(0.32457268) q[1];
x q[2];
rz(-1.764749) q[3];
sx q[3];
rz(-2.5911281) q[3];
sx q[3];
rz(2.6725519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79823309) q[2];
sx q[2];
rz(-1.5340021) q[2];
sx q[2];
rz(0.043924335) q[2];
rz(1.515306) q[3];
sx q[3];
rz(-1.1224727) q[3];
sx q[3];
rz(1.4656434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2783022) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(-2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-2.6712766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3669339) q[0];
sx q[0];
rz(-2.1485188) q[0];
sx q[0];
rz(-1.7455186) q[0];
rz(-0.83018556) q[2];
sx q[2];
rz(-0.52415327) q[2];
sx q[2];
rz(-0.94607991) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1086169) q[1];
sx q[1];
rz(-1.2101296) q[1];
sx q[1];
rz(0.74957871) q[1];
rz(-0.11300762) q[3];
sx q[3];
rz(-1.5586434) q[3];
sx q[3];
rz(-0.069610217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(-0.31141591) q[2];
rz(-2.9733859) q[3];
sx q[3];
rz(-1.5752537) q[3];
sx q[3];
rz(2.8581207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300221) q[0];
sx q[0];
rz(-0.011307414) q[0];
sx q[0];
rz(2.2054963) q[0];
rz(-2.6452737) q[1];
sx q[1];
rz(-2.4744108) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054809991) q[0];
sx q[0];
rz(-2.5544689) q[0];
sx q[0];
rz(2.5820288) q[0];
x q[1];
rz(-1.5580721) q[2];
sx q[2];
rz(-0.1977405) q[2];
sx q[2];
rz(1.5418574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5230323) q[1];
sx q[1];
rz(-1.9123239) q[1];
sx q[1];
rz(-0.073445436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27630287) q[3];
sx q[3];
rz(-1.7401764) q[3];
sx q[3];
rz(-0.15514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5101667) q[2];
sx q[2];
rz(-1.7743856) q[2];
sx q[2];
rz(1.4754971) q[2];
rz(0.79536074) q[3];
sx q[3];
rz(-2.9795591) q[3];
sx q[3];
rz(-1.8719261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0610166) q[0];
sx q[0];
rz(-0.59122714) q[0];
sx q[0];
rz(-3.0812145) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(-0.9789595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87099052) q[0];
sx q[0];
rz(-0.63828642) q[0];
sx q[0];
rz(2.9459475) q[0];
rz(-pi) q[1];
rz(-0.27411119) q[2];
sx q[2];
rz(-2.2345671) q[2];
sx q[2];
rz(-1.2520777) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1711848) q[1];
sx q[1];
rz(-1.1029585) q[1];
sx q[1];
rz(-2.4572608) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99154226) q[3];
sx q[3];
rz(-1.5994306) q[3];
sx q[3];
rz(0.34285173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0673361) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(-0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(0.0028751956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6947967) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(0.52892518) q[0];
rz(-0.18813285) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(0.58427748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463294) q[0];
sx q[0];
rz(-1.8285311) q[0];
sx q[0];
rz(1.3168174) q[0];
rz(-pi) q[1];
rz(-2.7033349) q[2];
sx q[2];
rz(-1.9263679) q[2];
sx q[2];
rz(2.8388765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1288695) q[1];
sx q[1];
rz(-1.1616316) q[1];
sx q[1];
rz(1.6175458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1891865) q[3];
sx q[3];
rz(-1.3460396) q[3];
sx q[3];
rz(0.99638961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3820485) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(0.071852597) q[2];
rz(2.1182649) q[3];
sx q[3];
rz(-1.5691248) q[3];
sx q[3];
rz(0.48880997) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857984) q[0];
sx q[0];
rz(-0.6577984) q[0];
sx q[0];
rz(-0.26554769) q[0];
rz(2.4664948) q[1];
sx q[1];
rz(-1.5470807) q[1];
sx q[1];
rz(3.0463228) q[1];
rz(3.0037389) q[2];
sx q[2];
rz(-2.6225435) q[2];
sx q[2];
rz(-0.23487716) q[2];
rz(1.4102546) q[3];
sx q[3];
rz(-1.7358801) q[3];
sx q[3];
rz(0.87270234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
