OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(0.41391882) q[0];
rz(3.4105372) q[1];
sx q[1];
rz(6.6528448) q[1];
sx q[1];
rz(11.289968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9577173) q[0];
sx q[0];
rz(-1.4801637) q[0];
sx q[0];
rz(-0.091770016) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90246713) q[2];
sx q[2];
rz(-1.5215708) q[2];
sx q[2];
rz(2.7572981) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40117404) q[1];
sx q[1];
rz(-0.83073069) q[1];
sx q[1];
rz(0.38895815) q[1];
x q[2];
rz(1.7539361) q[3];
sx q[3];
rz(-0.98455849) q[3];
sx q[3];
rz(2.8318015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74497491) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(1.940894) q[2];
rz(-2.4602304) q[3];
sx q[3];
rz(-1.6187637) q[3];
sx q[3];
rz(2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6354527) q[0];
sx q[0];
rz(-0.30158392) q[0];
sx q[0];
rz(0.26148456) q[0];
rz(-1.6188949) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(-1.2759298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44156333) q[0];
sx q[0];
rz(-1.6671379) q[0];
sx q[0];
rz(-0.55192134) q[0];
x q[1];
rz(1.0307113) q[2];
sx q[2];
rz(-0.99100494) q[2];
sx q[2];
rz(-2.3183398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1497685) q[1];
sx q[1];
rz(-0.64095488) q[1];
sx q[1];
rz(-0.561552) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7959782) q[3];
sx q[3];
rz(-2.9888267) q[3];
sx q[3];
rz(-1.6880715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9246989) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(0.72803289) q[2];
rz(1.3701471) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(3.106015) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1074693) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(0.6849826) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(-1.0122976) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489126) q[0];
sx q[0];
rz(-1.9143701) q[0];
sx q[0];
rz(1.9846041) q[0];
rz(-pi) q[1];
rz(1.3755765) q[2];
sx q[2];
rz(-1.3998271) q[2];
sx q[2];
rz(2.7107216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60437459) q[1];
sx q[1];
rz(-0.24362016) q[1];
sx q[1];
rz(-2.591406) q[1];
rz(-pi) q[2];
rz(-0.52148444) q[3];
sx q[3];
rz(-2.3064724) q[3];
sx q[3];
rz(-2.2948752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.181902) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(2.336179) q[2];
rz(-0.81625932) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(-0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212946) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(-0.62776172) q[0];
rz(-0.10920564) q[1];
sx q[1];
rz(-2.2716378) q[1];
sx q[1];
rz(2.3039718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7895486) q[0];
sx q[0];
rz(-1.0800011) q[0];
sx q[0];
rz(-1.646203) q[0];
rz(2.1659746) q[2];
sx q[2];
rz(-1.383923) q[2];
sx q[2];
rz(2.486791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1091812) q[1];
sx q[1];
rz(-2.4285168) q[1];
sx q[1];
rz(-2.9362657) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37552278) q[3];
sx q[3];
rz(-1.2831472) q[3];
sx q[3];
rz(0.2660397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3179021) q[2];
sx q[2];
rz(-0.48419848) q[2];
sx q[2];
rz(-0.34453264) q[2];
rz(-1.4065929) q[3];
sx q[3];
rz(-0.35511261) q[3];
sx q[3];
rz(-3.0936892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59871167) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(0.5058381) q[0];
rz(-2.089031) q[1];
sx q[1];
rz(-0.86741766) q[1];
sx q[1];
rz(0.94598407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3179683) q[0];
sx q[0];
rz(-1.8232947) q[0];
sx q[0];
rz(-0.28244762) q[0];
x q[1];
rz(0.37028124) q[2];
sx q[2];
rz(-0.38114377) q[2];
sx q[2];
rz(-0.76424341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9717488) q[1];
sx q[1];
rz(-1.9439812) q[1];
sx q[1];
rz(0.62974522) q[1];
rz(-pi) q[2];
rz(-1.4220174) q[3];
sx q[3];
rz(-1.6384203) q[3];
sx q[3];
rz(0.54883445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0537009) q[2];
sx q[2];
rz(-2.044007) q[2];
sx q[2];
rz(-1.5778479) q[2];
rz(1.0846694) q[3];
sx q[3];
rz(-0.8756777) q[3];
sx q[3];
rz(2.5904371) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797193) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(-0.24000034) q[0];
rz(-2.2198246) q[1];
sx q[1];
rz(-2.0926937) q[1];
sx q[1];
rz(-1.9711432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76239785) q[0];
sx q[0];
rz(-2.1403011) q[0];
sx q[0];
rz(-1.5194917) q[0];
rz(-pi) q[1];
rz(0.35640772) q[2];
sx q[2];
rz(-1.4517759) q[2];
sx q[2];
rz(1.0072034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41783782) q[1];
sx q[1];
rz(-2.0643294) q[1];
sx q[1];
rz(0.51512169) q[1];
rz(-pi) q[2];
rz(1.7975054) q[3];
sx q[3];
rz(-0.41019687) q[3];
sx q[3];
rz(2.6997379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3971098) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(0.58103713) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-2.5751028) q[3];
sx q[3];
rz(0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.13643232) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(-0.50049472) q[0];
rz(-2.9035134) q[1];
sx q[1];
rz(-1.4512738) q[1];
sx q[1];
rz(-0.64661017) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85486551) q[0];
sx q[0];
rz(-1.6056132) q[0];
sx q[0];
rz(-1.5725732) q[0];
x q[1];
rz(1.6853203) q[2];
sx q[2];
rz(-1.537064) q[2];
sx q[2];
rz(-2.9303868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.096488086) q[1];
sx q[1];
rz(-2.641171) q[1];
sx q[1];
rz(-1.3184271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.929333) q[3];
sx q[3];
rz(-1.6286116) q[3];
sx q[3];
rz(-1.1392347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87767345) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(0.70362299) q[2];
rz(1.2879114) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(0.50584403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25381655) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(1.8164841) q[0];
rz(2.3690986) q[1];
sx q[1];
rz(-2.019181) q[1];
sx q[1];
rz(-2.966029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85925519) q[0];
sx q[0];
rz(-1.757375) q[0];
sx q[0];
rz(-0.76754359) q[0];
rz(-pi) q[1];
rz(2.2261691) q[2];
sx q[2];
rz(-2.2594514) q[2];
sx q[2];
rz(-1.6885533) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9055134) q[1];
sx q[1];
rz(-1.4314974) q[1];
sx q[1];
rz(2.0361316) q[1];
rz(2.4252093) q[3];
sx q[3];
rz(-1.4759403) q[3];
sx q[3];
rz(1.5945292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0385711) q[2];
sx q[2];
rz(-2.0671637) q[2];
sx q[2];
rz(1.0938905) q[2];
rz(-1.1755747) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(-0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232101) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(1.3257931) q[0];
rz(-2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(-2.3525499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51368749) q[0];
sx q[0];
rz(-0.18842489) q[0];
sx q[0];
rz(-1.0303251) q[0];
rz(-2.5297935) q[2];
sx q[2];
rz(-0.42610301) q[2];
sx q[2];
rz(-0.95624051) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4261158) q[1];
sx q[1];
rz(-1.9144692) q[1];
sx q[1];
rz(-2.2261376) q[1];
x q[2];
rz(-0.77603839) q[3];
sx q[3];
rz(-1.0641085) q[3];
sx q[3];
rz(-2.6069178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21251707) q[2];
sx q[2];
rz(-2.2079461) q[2];
sx q[2];
rz(2.2361501) q[2];
rz(0.91600156) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(-1.0441095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020029) q[0];
sx q[0];
rz(-0.8332533) q[0];
sx q[0];
rz(0.92779094) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-2.5826726) q[1];
sx q[1];
rz(-0.13339001) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095208406) q[0];
sx q[0];
rz(-0.36050561) q[0];
sx q[0];
rz(0.62490873) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8850408) q[2];
sx q[2];
rz(-0.5945328) q[2];
sx q[2];
rz(2.9192215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9815753) q[1];
sx q[1];
rz(-2.3381553) q[1];
sx q[1];
rz(-2.395242) q[1];
rz(-2.9319256) q[3];
sx q[3];
rz(-1.121796) q[3];
sx q[3];
rz(-1.0202788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(-1.8782328) q[2];
rz(-1.0697621) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(0.076816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99061154) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(1.6928584) q[1];
sx q[1];
rz(-0.81304638) q[1];
sx q[1];
rz(-0.1437694) q[1];
rz(-1.0461367) q[2];
sx q[2];
rz(-1.8423648) q[2];
sx q[2];
rz(0.8348196) q[2];
rz(-1.9142022) q[3];
sx q[3];
rz(-1.8029984) q[3];
sx q[3];
rz(2.3912994) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
