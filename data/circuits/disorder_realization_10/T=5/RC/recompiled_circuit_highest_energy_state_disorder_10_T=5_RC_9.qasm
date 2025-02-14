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
rz(7.8742134) q[0];
sx q[0];
rz(5.9013517) q[0];
rz(-0.39185697) q[1];
sx q[1];
rz(4.691603) q[1];
sx q[1];
rz(7.2075972) q[1];
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
rz(2.4886542) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37021356) q[2];
sx q[2];
rz(-1.5892803) q[2];
sx q[2];
rz(-0.12272515) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0142483) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(-0.30502747) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0409058) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.6774981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(1.5667685) q[2];
rz(1.9751679) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(-1.9655193) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(2.8228021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27869895) q[0];
sx q[0];
rz(-2.3434533) q[0];
sx q[0];
rz(-0.44095914) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(2.3815439) q[1];
sx q[1];
rz(-2.4439619) q[1];
sx q[1];
rz(-2.182644) q[1];
rz(-pi) q[2];
rz(3.0158668) q[3];
sx q[3];
rz(-1.2866402) q[3];
sx q[3];
rz(1.3181669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0789644) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(2.6641565) q[2];
rz(-1.7834974) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7497712) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(0.40621743) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11230532) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(1.7303921) q[0];
rz(-pi) q[1];
rz(0.35786104) q[2];
sx q[2];
rz(-2.1659019) q[2];
sx q[2];
rz(2.8107363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4809847) q[1];
sx q[1];
rz(-2.0913731) q[1];
sx q[1];
rz(-2.1064227) q[1];
rz(-pi) q[2];
rz(2.0800679) q[3];
sx q[3];
rz(-1.8647927) q[3];
sx q[3];
rz(1.7727838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.5072301) q[2];
sx q[2];
rz(-2.0646084) q[2];
rz(-1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-2.8857359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(2.381109) q[0];
rz(-0.35176945) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(-0.15028353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0194797) q[0];
sx q[0];
rz(-1.5727336) q[0];
sx q[0];
rz(1.5674855) q[0];
rz(-0.17268659) q[2];
sx q[2];
rz(-1.9479247) q[2];
sx q[2];
rz(-1.0972301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80881883) q[1];
sx q[1];
rz(-1.2398232) q[1];
sx q[1];
rz(2.1899151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6500449) q[3];
sx q[3];
rz(-2.632049) q[3];
sx q[3];
rz(2.5423637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1111697) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-2.5566768) q[2];
rz(0.065718204) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.86972648) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(0.72750339) q[0];
rz(1.8456521) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(0.71472275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097209771) q[0];
sx q[0];
rz(-1.3499027) q[0];
sx q[0];
rz(-1.4979375) q[0];
rz(-0.3008414) q[2];
sx q[2];
rz(-1.46035) q[2];
sx q[2];
rz(0.12953239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89940577) q[1];
sx q[1];
rz(-1.9493628) q[1];
sx q[1];
rz(2.942571) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1470137) q[3];
sx q[3];
rz(-1.9061229) q[3];
sx q[3];
rz(-1.2352342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99594816) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(1.9963416) q[2];
rz(-2.9389985) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(-2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-1.8458493) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(-0.16786815) q[0];
rz(-0.56495086) q[1];
sx q[1];
rz(-2.0940557) q[1];
sx q[1];
rz(-1.3289183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58121366) q[0];
sx q[0];
rz(-0.74958505) q[0];
sx q[0];
rz(-2.4102274) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8064427) q[2];
sx q[2];
rz(-0.14832917) q[2];
sx q[2];
rz(2.2466765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.708864) q[1];
sx q[1];
rz(-1.335718) q[1];
sx q[1];
rz(-0.26223305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8096446) q[3];
sx q[3];
rz(-1.3535366) q[3];
sx q[3];
rz(-2.1715733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7557175) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(-0.34783777) q[2];
rz(2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(-1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.6353761) q[0];
sx q[0];
rz(0.95112479) q[0];
rz(-0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(2.1975885) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4571465) q[0];
sx q[0];
rz(-1.5553885) q[0];
sx q[0];
rz(-1.2627361) q[0];
rz(-pi) q[1];
rz(1.0877785) q[2];
sx q[2];
rz(-1.1406787) q[2];
sx q[2];
rz(-1.5120235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4170226) q[1];
sx q[1];
rz(-0.17462433) q[1];
sx q[1];
rz(-1.440446) q[1];
x q[2];
rz(1.0349501) q[3];
sx q[3];
rz(-0.45848819) q[3];
sx q[3];
rz(-1.5036086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6501179) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(-3.1169543) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(-0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80855075) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(3.0035875) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.4638008) q[1];
sx q[1];
rz(-0.93625751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27565378) q[0];
sx q[0];
rz(-1.8714974) q[0];
sx q[0];
rz(-2.7554437) q[0];
x q[1];
rz(-0.9844043) q[2];
sx q[2];
rz(-2.3086562) q[2];
sx q[2];
rz(3.0393103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4605751) q[1];
sx q[1];
rz(-1.5854757) q[1];
sx q[1];
rz(-1.753586) q[1];
rz(-0.3206203) q[3];
sx q[3];
rz(-0.71133876) q[3];
sx q[3];
rz(-2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.998698) q[2];
sx q[2];
rz(-1.2724718) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(3.0190492) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97813022) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(2.2273831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5152138) q[0];
sx q[0];
rz(-0.45013407) q[0];
sx q[0];
rz(-1.882097) q[0];
rz(1.5246824) q[2];
sx q[2];
rz(-0.83674016) q[2];
sx q[2];
rz(-0.095832713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57371512) q[1];
sx q[1];
rz(-0.93522954) q[1];
sx q[1];
rz(-1.3997517) q[1];
x q[2];
rz(1.4214324) q[3];
sx q[3];
rz(-1.0348667) q[3];
sx q[3];
rz(1.5045297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(-1.5128822) q[2];
rz(0.58879876) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(-0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.47422472) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(2.4868763) q[1];
sx q[1];
rz(-1.3755211) q[1];
sx q[1];
rz(-2.7213352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3271844) q[0];
sx q[0];
rz(-1.9851917) q[0];
sx q[0];
rz(2.9073858) q[0];
rz(-pi) q[1];
rz(2.6819314) q[2];
sx q[2];
rz(-0.84520413) q[2];
sx q[2];
rz(-0.21506735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4431364) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(0.08929792) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2154142) q[3];
sx q[3];
rz(-1.0338262) q[3];
sx q[3];
rz(-0.87600183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(-0.9672271) q[2];
rz(-2.1238756) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(-2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9781072) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(-3.0178487) q[1];
sx q[1];
rz(-2.1950304) q[1];
sx q[1];
rz(-1.310941) q[1];
rz(0.38519771) q[2];
sx q[2];
rz(-0.34389797) q[2];
sx q[2];
rz(1.4138538) q[2];
rz(-1.5354034) q[3];
sx q[3];
rz(-1.2589068) q[3];
sx q[3];
rz(1.3990228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
