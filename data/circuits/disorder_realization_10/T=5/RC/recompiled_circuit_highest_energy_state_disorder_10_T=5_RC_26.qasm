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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782717) q[0];
sx q[0];
rz(-0.98211432) q[0];
sx q[0];
rz(2.4886542) q[0];
rz(-pi) q[1];
rz(-3.0905444) q[2];
sx q[2];
rz(-2.7709392) q[2];
sx q[2];
rz(1.6459344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1273444) q[1];
sx q[1];
rz(-0.44522983) q[1];
sx q[1];
rz(2.8365652) q[1];
x q[2];
rz(-2.754491) q[3];
sx q[3];
rz(-2.0678068) q[3];
sx q[3];
rz(-0.29868515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(-1.5748242) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(-0.69275698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.9655193) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-0.57756966) q[1];
sx q[1];
rz(0.31879058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8628937) q[0];
sx q[0];
rz(-2.3434533) q[0];
sx q[0];
rz(-0.44095914) q[0];
rz(1.8068878) q[2];
sx q[2];
rz(-1.9892879) q[2];
sx q[2];
rz(-0.22007569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31736703) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(0.96940689) q[1];
x q[2];
rz(1.2845006) q[3];
sx q[3];
rz(-1.4501374) q[3];
sx q[3];
rz(0.21721043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(-0.47743615) q[2];
rz(1.7834974) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(-3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39182144) q[0];
sx q[0];
rz(-1.8056159) q[0];
sx q[0];
rz(-2.0023477) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(-2.4635945) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8157522) q[0];
sx q[0];
rz(-1.4817955) q[0];
sx q[0];
rz(-0.98317105) q[0];
rz(-pi) q[1];
rz(-0.35786104) q[2];
sx q[2];
rz(-2.1659019) q[2];
sx q[2];
rz(-2.8107363) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4809847) q[1];
sx q[1];
rz(-1.0502195) q[1];
sx q[1];
rz(-1.03517) q[1];
x q[2];
rz(-0.33379995) q[3];
sx q[3];
rz(-2.0562647) q[3];
sx q[3];
rz(0.36241058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(1.8814603) q[3];
sx q[3];
rz(-2.1144805) q[3];
sx q[3];
rz(-0.25585678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.857321) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(0.76048365) q[0];
rz(2.7898232) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(0.15028353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9780804) q[0];
sx q[0];
rz(-0.0038359782) q[0];
sx q[0];
rz(-1.0413961) q[0];
x q[1];
rz(-1.1885095) q[2];
sx q[2];
rz(-1.4103544) q[2];
sx q[2];
rz(2.7321702) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6071668) q[1];
sx q[1];
rz(-2.1517506) q[1];
sx q[1];
rz(-0.39926417) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.683879) q[3];
sx q[3];
rz(-1.8031075) q[3];
sx q[3];
rz(1.4087698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.030423) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-0.58491582) q[2];
rz(-3.0758744) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.86972648) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(0.72750339) q[0];
rz(1.2959405) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(2.4268699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.489577) q[0];
sx q[0];
rz(-1.4997109) q[0];
sx q[0];
rz(-2.9201304) q[0];
rz(2.8407513) q[2];
sx q[2];
rz(-1.6812426) q[2];
sx q[2];
rz(3.0120603) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59698518) q[1];
sx q[1];
rz(-1.7555573) q[1];
sx q[1];
rz(-1.1853335) q[1];
rz(-pi) q[2];
rz(-2.7477164) q[3];
sx q[3];
rz(-2.1112313) q[3];
sx q[3];
rz(-0.5462164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-0.069772094) q[3];
sx q[3];
rz(-2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.2957434) q[0];
sx q[0];
rz(-2.8491617) q[0];
sx q[0];
rz(0.16786815) q[0];
rz(2.5766418) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(1.3289183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.57083) q[0];
sx q[0];
rz(-2.0432297) q[0];
sx q[0];
rz(-0.60586141) q[0];
rz(-pi) q[1];
rz(1.5216878) q[2];
sx q[2];
rz(-1.430776) q[2];
sx q[2];
rz(-2.5852709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5766474) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(2.3957361) q[1];
rz(-pi) q[2];
rz(1.331948) q[3];
sx q[3];
rz(-1.788056) q[3];
sx q[3];
rz(0.9700194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7557175) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(-0.34783777) q[2];
rz(2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(-0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(2.1975885) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2601449) q[0];
sx q[0];
rz(-1.8788188) q[0];
sx q[0];
rz(-0.016168895) q[0];
rz(-pi) q[1];
rz(2.6636276) q[2];
sx q[2];
rz(-1.1349808) q[2];
sx q[2];
rz(-0.27406853) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72457006) q[1];
sx q[1];
rz(-2.9669683) q[1];
sx q[1];
rz(-1.440446) q[1];
x q[2];
rz(2.8947325) q[3];
sx q[3];
rz(-1.9611957) q[3];
sx q[3];
rz(-0.91871432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(-3.1169543) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(2.1520481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330419) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(3.0035875) q[0];
rz(1.0778614) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-2.2053351) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753006) q[0];
sx q[0];
rz(-1.202824) q[0];
sx q[0];
rz(1.247768) q[0];
rz(-pi) q[1];
rz(0.54674863) q[2];
sx q[2];
rz(-2.2348352) q[2];
sx q[2];
rz(-0.67686096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.028658) q[1];
sx q[1];
rz(-1.3880265) q[1];
sx q[1];
rz(-3.1266646) q[1];
x q[2];
rz(-2.4560087) q[3];
sx q[3];
rz(-1.778025) q[3];
sx q[3];
rz(2.384546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.3188837) q[2];
rz(3.0190492) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97813022) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(-0.36460707) q[0];
rz(-1.7568024) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(-0.91420954) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4791661) q[0];
sx q[0];
rz(-1.7044596) q[0];
sx q[0];
rz(-1.1396618) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6169102) q[2];
sx q[2];
rz(-2.3048525) q[2];
sx q[2];
rz(-0.095832713) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5678775) q[1];
sx q[1];
rz(-0.93522954) q[1];
sx q[1];
rz(1.3997517) q[1];
rz(-pi) q[2];
rz(-1.7201603) q[3];
sx q[3];
rz(-1.0348667) q[3];
sx q[3];
rz(1.5045297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30283516) q[2];
sx q[2];
rz(-2.5573533) q[2];
sx q[2];
rz(1.6287104) q[2];
rz(-2.5527939) q[3];
sx q[3];
rz(-1.9353119) q[3];
sx q[3];
rz(-2.485386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6673679) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(2.8571416) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.3755211) q[1];
sx q[1];
rz(2.7213352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3493372) q[0];
sx q[0];
rz(-0.47264034) q[0];
sx q[0];
rz(-2.0562754) q[0];
rz(-pi) q[1];
rz(-2.6819314) q[2];
sx q[2];
rz(-0.84520413) q[2];
sx q[2];
rz(-2.9265253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4431364) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(0.08929792) q[1];
x q[2];
rz(-2.612667) q[3];
sx q[3];
rz(-0.63415895) q[3];
sx q[3];
rz(2.8932539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(-2.1238756) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(-0.12374395) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(1.704557) q[2];
sx q[2];
rz(-1.8885713) q[2];
sx q[2];
rz(1.0071913) q[2];
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
