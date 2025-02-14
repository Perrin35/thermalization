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
rz(2.2171807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.763321) q[0];
sx q[0];
rz(-2.1594783) q[0];
sx q[0];
rz(0.6529385) q[0];
rz(0.37021356) q[2];
sx q[2];
rz(-1.5523124) q[2];
sx q[2];
rz(-0.12272515) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3081835) q[1];
sx q[1];
rz(-1.4410958) q[1];
sx q[1];
rz(0.42713487) q[1];
x q[2];
rz(0.96278874) q[3];
sx q[3];
rz(-0.6198403) q[3];
sx q[3];
rz(2.7328797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9367289) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(1.5748242) q[2];
rz(-1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(-0.69275698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(-0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(0.31879058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87321157) q[0];
sx q[0];
rz(-0.86641524) q[0];
sx q[0];
rz(1.1581139) q[0];
x q[1];
rz(-2.6574357) q[2];
sx q[2];
rz(-0.47704298) q[2];
sx q[2];
rz(-2.8271528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8242256) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(-0.96940689) q[1];
x q[2];
rz(-1.2845006) q[3];
sx q[3];
rz(-1.4501374) q[3];
sx q[3];
rz(2.9243822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(0.47743615) q[2];
rz(1.3580953) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(0.0099946578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7497712) q[0];
sx q[0];
rz(-1.8056159) q[0];
sx q[0];
rz(-1.1392449) q[0];
rz(-2.7353752) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(-2.4635945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0292873) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(1.4112006) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7837316) q[2];
sx q[2];
rz(-0.97569078) q[2];
sx q[2];
rz(0.33085631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78732508) q[1];
sx q[1];
rz(-0.72871043) q[1];
sx q[1];
rz(-2.4142152) q[1];
rz(-pi) q[2];
rz(-0.33379995) q[3];
sx q[3];
rz(-1.085328) q[3];
sx q[3];
rz(2.7791821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(-1.8814603) q[3];
sx q[3];
rz(-2.1144805) q[3];
sx q[3];
rz(0.25585678) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28427163) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(-2.381109) q[0];
rz(-0.35176945) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(0.15028353) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6929157) q[0];
sx q[0];
rz(-1.5674855) q[0];
sx q[0];
rz(3.1396554) q[0];
x q[1];
rz(1.9801122) q[2];
sx q[2];
rz(-0.41305734) q[2];
sx q[2];
rz(1.5395791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3327738) q[1];
sx q[1];
rz(-1.2398232) q[1];
sx q[1];
rz(0.95167758) q[1];
rz(1.8286546) q[3];
sx q[3];
rz(-1.1262731) q[3];
sx q[3];
rz(-0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1111697) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-2.5566768) q[2];
rz(-3.0758744) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2718662) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(-2.4140893) q[0];
rz(-1.8456521) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(2.4268699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172404) q[0];
sx q[0];
rz(-2.9091798) q[0];
sx q[0];
rz(-0.31347855) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7834846) q[2];
sx q[2];
rz(-0.3198959) q[2];
sx q[2];
rz(1.7826155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59698518) q[1];
sx q[1];
rz(-1.3860354) q[1];
sx q[1];
rz(-1.1853335) q[1];
x q[2];
rz(2.1470137) q[3];
sx q[3];
rz(-1.2354697) q[3];
sx q[3];
rz(1.2352342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(-1.145251) q[2];
rz(2.9389985) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(2.8128459) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458493) q[0];
sx q[0];
rz(-2.8491617) q[0];
sx q[0];
rz(-0.16786815) q[0];
rz(2.5766418) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(1.3289183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.57083) q[0];
sx q[0];
rz(-2.0432297) q[0];
sx q[0];
rz(-2.5357312) q[0];
x q[1];
rz(-1.5216878) q[2];
sx q[2];
rz(-1.7108166) q[2];
sx q[2];
rz(-2.5852709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3420978) q[1];
sx q[1];
rz(-1.3159385) q[1];
sx q[1];
rz(1.8138769) q[1];
x q[2];
rz(-2.9181913) q[3];
sx q[3];
rz(-1.8039244) q[3];
sx q[3];
rz(-0.65321556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38587511) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(-2.7937549) q[2];
rz(2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(-1.1381963) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(-0.95112479) q[0];
rz(0.29257193) q[1];
sx q[1];
rz(-0.89076275) q[1];
sx q[1];
rz(2.1975885) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88144775) q[0];
sx q[0];
rz(-1.8788188) q[0];
sx q[0];
rz(-0.016168895) q[0];
rz(1.0877785) q[2];
sx q[2];
rz(-1.1406787) q[2];
sx q[2];
rz(1.6295691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2846825) q[1];
sx q[1];
rz(-1.7439242) q[1];
sx q[1];
rz(0.022927479) q[1];
rz(-pi) q[2];
rz(-1.1694447) q[3];
sx q[3];
rz(-1.7987393) q[3];
sx q[3];
rz(-2.5851188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(0.024638351) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(0.98954454) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330419) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(3.0035875) q[0];
rz(2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(2.2053351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8659389) q[0];
sx q[0];
rz(-1.8714974) q[0];
sx q[0];
rz(0.38614892) q[0];
rz(-pi) q[1];
rz(-0.9844043) q[2];
sx q[2];
rz(-2.3086562) q[2];
sx q[2];
rz(-0.10228233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11293462) q[1];
sx q[1];
rz(-1.3880265) q[1];
sx q[1];
rz(-3.1266646) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8209723) q[3];
sx q[3];
rz(-2.4302539) q[3];
sx q[3];
rz(2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.3188837) q[2];
rz(0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97813022) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(-0.36460707) q[0];
rz(-1.3847903) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(2.2273831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62637882) q[0];
sx q[0];
rz(-2.6914586) q[0];
sx q[0];
rz(1.882097) q[0];
x q[1];
rz(1.5246824) q[2];
sx q[2];
rz(-0.83674016) q[2];
sx q[2];
rz(-0.095832713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89489854) q[1];
sx q[1];
rz(-1.7082038) q[1];
sx q[1];
rz(-0.64260428) q[1];
x q[2];
rz(-0.24550415) q[3];
sx q[3];
rz(-0.55439083) q[3];
sx q[3];
rz(1.3504775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(1.5128822) q[2];
rz(2.5527939) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(-2.485386) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
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
rz(pi/2) q[2];
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
rz(-2.6819314) q[2];
sx q[2];
rz(-0.84520413) q[2];
sx q[2];
rz(-2.9265253) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0033231) q[1];
sx q[1];
rz(-1.6594634) q[1];
sx q[1];
rz(1.4515463) q[1];
x q[2];
rz(-2.612667) q[3];
sx q[3];
rz(-0.63415895) q[3];
sx q[3];
rz(-0.24833873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(-2.1743656) q[2];
rz(-2.1238756) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9781072) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(3.0178487) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(0.38519771) q[2];
sx q[2];
rz(-0.34389797) q[2];
sx q[2];
rz(1.4138538) q[2];
rz(0.10931482) q[3];
sx q[3];
rz(-2.8277665) q[3];
sx q[3];
rz(1.2841429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
