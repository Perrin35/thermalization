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
rz(2.8528557) q[0];
sx q[0];
rz(5.5878162) q[0];
sx q[0];
rz(6.5509808) q[0];
rz(-2.7195622) q[1];
sx q[1];
rz(-0.92075092) q[1];
sx q[1];
rz(1.86778) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86912495) q[0];
sx q[0];
rz(-2.8960138) q[0];
sx q[0];
rz(-1.7564943) q[0];
x q[1];
rz(2.2681178) q[2];
sx q[2];
rz(-1.0823853) q[2];
sx q[2];
rz(-1.4174051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85292683) q[1];
sx q[1];
rz(-1.5024606) q[1];
sx q[1];
rz(2.4803376) q[1];
rz(-pi) q[2];
rz(-3.1056728) q[3];
sx q[3];
rz(-1.714141) q[3];
sx q[3];
rz(-2.8570321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6196809) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(3.0989975) q[2];
rz(0.28111449) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(-2.5458096) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5113145) q[0];
sx q[0];
rz(-1.1742641) q[0];
sx q[0];
rz(-2.0654772) q[0];
rz(-1.0307182) q[1];
sx q[1];
rz(-2.0298256) q[1];
sx q[1];
rz(-1.3105185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56633184) q[0];
sx q[0];
rz(-1.5916414) q[0];
sx q[0];
rz(0.91253176) q[0];
rz(2.7380472) q[2];
sx q[2];
rz(-2.2650044) q[2];
sx q[2];
rz(-0.36118868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84634631) q[1];
sx q[1];
rz(-2.935886) q[1];
sx q[1];
rz(-0.61726112) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53383975) q[3];
sx q[3];
rz(-1.2152142) q[3];
sx q[3];
rz(1.3188286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92404667) q[2];
sx q[2];
rz(-0.7520389) q[2];
sx q[2];
rz(-2.1853866) q[2];
rz(1.8337967) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(3.0202878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8420551) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(-0.036238413) q[0];
rz(0.80129519) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5755441) q[0];
sx q[0];
rz(-1.4641718) q[0];
sx q[0];
rz(1.7315277) q[0];
rz(-pi) q[1];
rz(3.0984584) q[2];
sx q[2];
rz(-0.97402527) q[2];
sx q[2];
rz(1.5297001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8138995) q[1];
sx q[1];
rz(-2.3997953) q[1];
sx q[1];
rz(2.2375536) q[1];
x q[2];
rz(-0.43576305) q[3];
sx q[3];
rz(-0.73832694) q[3];
sx q[3];
rz(-2.3070525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7654045) q[2];
sx q[2];
rz(-1.7690965) q[2];
sx q[2];
rz(-0.21793951) q[2];
rz(2.229522) q[3];
sx q[3];
rz(-1.6488766) q[3];
sx q[3];
rz(-2.2899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(2.3260314) q[0];
rz(1.0527481) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(1.7900593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4690875) q[0];
sx q[0];
rz(-2.6658305) q[0];
sx q[0];
rz(-0.16938727) q[0];
rz(-pi) q[1];
rz(2.4250373) q[2];
sx q[2];
rz(-1.637023) q[2];
sx q[2];
rz(2.9481681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38511577) q[1];
sx q[1];
rz(-2.6653892) q[1];
sx q[1];
rz(-1.6333196) q[1];
x q[2];
rz(-2.4249486) q[3];
sx q[3];
rz(-1.1782681) q[3];
sx q[3];
rz(1.5115304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1286596) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(-2.2991119) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(-2.1459818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.9812444) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(0.55111849) q[0];
rz(2.5235858) q[1];
sx q[1];
rz(-1.2851241) q[1];
sx q[1];
rz(-1.0689703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7525714) q[0];
sx q[0];
rz(-2.4271936) q[0];
sx q[0];
rz(2.8117489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0238012) q[2];
sx q[2];
rz(-2.0108622) q[2];
sx q[2];
rz(-2.5534782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.080090962) q[1];
sx q[1];
rz(-1.4451761) q[1];
sx q[1];
rz(2.3088423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8884747) q[3];
sx q[3];
rz(-2.6147644) q[3];
sx q[3];
rz(1.4344858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7908287) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(-2.186415) q[2];
rz(0.59349924) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.3957297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7668358) q[0];
sx q[0];
rz(-0.042348472) q[0];
sx q[0];
rz(1.7497077) q[0];
rz(-1.998924) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(1.3124189) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5757933) q[0];
sx q[0];
rz(-2.2941219) q[0];
sx q[0];
rz(0.91200836) q[0];
x q[1];
rz(2.4644971) q[2];
sx q[2];
rz(-2.1506718) q[2];
sx q[2];
rz(3.0486272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10656938) q[1];
sx q[1];
rz(-0.5173389) q[1];
sx q[1];
rz(-0.10915233) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3950609) q[3];
sx q[3];
rz(-0.32009691) q[3];
sx q[3];
rz(1.8863581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7890847) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(-1.4136723) q[2];
rz(0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(2.5889034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308052) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(-0.68156534) q[0];
rz(-1.956578) q[1];
sx q[1];
rz(-2.3763035) q[1];
sx q[1];
rz(-2.3550745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241572) q[0];
sx q[0];
rz(-0.52612129) q[0];
sx q[0];
rz(-2.4301162) q[0];
x q[1];
rz(0.925073) q[2];
sx q[2];
rz(-2.4074915) q[2];
sx q[2];
rz(-0.9303329) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29132641) q[1];
sx q[1];
rz(-0.89610142) q[1];
sx q[1];
rz(3.0043601) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1004518) q[3];
sx q[3];
rz(-2.3015071) q[3];
sx q[3];
rz(-1.6139489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6098392) q[2];
sx q[2];
rz(-0.24007758) q[2];
sx q[2];
rz(1.5717724) q[2];
rz(1.2872559) q[3];
sx q[3];
rz(-1.5232892) q[3];
sx q[3];
rz(-0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59931961) q[0];
sx q[0];
rz(-0.71745187) q[0];
sx q[0];
rz(1.188311) q[0];
rz(-3.1207454) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(-2.5426224) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99020236) q[0];
sx q[0];
rz(-1.6693475) q[0];
sx q[0];
rz(2.0771785) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39283237) q[2];
sx q[2];
rz(-2.246703) q[2];
sx q[2];
rz(2.3623938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.54810337) q[1];
sx q[1];
rz(-2.2240891) q[1];
sx q[1];
rz(-0.95667019) q[1];
rz(0.4407626) q[3];
sx q[3];
rz(-2.2513933) q[3];
sx q[3];
rz(-1.1632533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.282436) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(-2.6089) q[3];
sx q[3];
rz(-2.0566514) q[3];
sx q[3];
rz(2.2044619) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027997967) q[0];
sx q[0];
rz(-0.6830712) q[0];
sx q[0];
rz(-0.37044507) q[0];
rz(2.0619552) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(0.058578514) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97998842) q[0];
sx q[0];
rz(-1.7287041) q[0];
sx q[0];
rz(0.61451332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3434161) q[2];
sx q[2];
rz(-1.5274962) q[2];
sx q[2];
rz(-2.3169278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1540268) q[1];
sx q[1];
rz(-0.79001629) q[1];
sx q[1];
rz(-1.7055737) q[1];
rz(-pi) q[2];
rz(-2.305916) q[3];
sx q[3];
rz(-1.1181076) q[3];
sx q[3];
rz(1.2293466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2728682) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(0.33099428) q[2];
rz(-3.1281779) q[3];
sx q[3];
rz(-2.1881073) q[3];
sx q[3];
rz(1.0564055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5633504) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(2.6934534) q[0];
rz(-0.093712417) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.9261446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5543723) q[0];
sx q[0];
rz(-1.1379114) q[0];
sx q[0];
rz(-0.032399633) q[0];
x q[1];
rz(-2.7461461) q[2];
sx q[2];
rz(-0.85265358) q[2];
sx q[2];
rz(1.8888231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7772953) q[1];
sx q[1];
rz(-0.60579311) q[1];
sx q[1];
rz(-1.08849) q[1];
x q[2];
rz(-0.13588174) q[3];
sx q[3];
rz(-1.1757506) q[3];
sx q[3];
rz(-0.53077936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.965968) q[2];
sx q[2];
rz(-1.5886687) q[2];
sx q[2];
rz(0.69941163) q[2];
rz(-1.2753963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(-1.2709966) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324343) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(-1.746183) q[1];
sx q[1];
rz(-1.9231053) q[1];
sx q[1];
rz(-1.3833192) q[1];
rz(0.66580843) q[2];
sx q[2];
rz(-1.6338909) q[2];
sx q[2];
rz(0.8798107) q[2];
rz(0.46635177) q[3];
sx q[3];
rz(-1.4644571) q[3];
sx q[3];
rz(-1.1639948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
