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
rz(0.056869153) q[0];
sx q[0];
rz(2.9480204) q[0];
sx q[0];
rz(9.8326346) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9303794) q[0];
sx q[0];
rz(-1.6193074) q[0];
sx q[0];
rz(-1.1332953) q[0];
rz(-pi) q[1];
rz(0.27215927) q[2];
sx q[2];
rz(-2.6725884) q[2];
sx q[2];
rz(-0.16527612) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.08886628) q[1];
sx q[1];
rz(-1.7824934) q[1];
sx q[1];
rz(-0.39525169) q[1];
x q[2];
rz(-2.5276189) q[3];
sx q[3];
rz(-2.2212265) q[3];
sx q[3];
rz(-0.87343317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5539598) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(2.9239192) q[2];
rz(-0.31146464) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(-1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4199453) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(0.69197792) q[0];
rz(-0.75633374) q[1];
sx q[1];
rz(-0.5223918) q[1];
sx q[1];
rz(1.5914894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1082008) q[0];
sx q[0];
rz(-2.4188359) q[0];
sx q[0];
rz(-0.18527822) q[0];
x q[1];
rz(-1.3266356) q[2];
sx q[2];
rz(-2.7327644) q[2];
sx q[2];
rz(0.93016184) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8988425) q[1];
sx q[1];
rz(-1.4489343) q[1];
sx q[1];
rz(-0.16074462) q[1];
rz(-2.2859073) q[3];
sx q[3];
rz(-2.4484903) q[3];
sx q[3];
rz(2.2732796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456566) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(-1.650882) q[3];
sx q[3];
rz(-1.8473293) q[3];
sx q[3];
rz(1.2283121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.8656798) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(3.1269585) q[0];
rz(-0.89302653) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-2.9258974) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1954262) q[0];
sx q[0];
rz(-3.0901244) q[0];
sx q[0];
rz(1.6610751) q[0];
x q[1];
rz(1.2273618) q[2];
sx q[2];
rz(-2.2015155) q[2];
sx q[2];
rz(1.7523927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5657422) q[1];
sx q[1];
rz(-1.8025959) q[1];
sx q[1];
rz(1.7017468) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87743296) q[3];
sx q[3];
rz(-2.139353) q[3];
sx q[3];
rz(-1.3662149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0594242) q[2];
sx q[2];
rz(-1.9614204) q[2];
sx q[2];
rz(2.1480985) q[2];
rz(-2.4332186) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(3.0181001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74743903) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(-2.4476449) q[0];
rz(2.8548062) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(2.0538816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74831731) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(1.0626777) q[0];
rz(-0.55107815) q[2];
sx q[2];
rz(-2.0281396) q[2];
sx q[2];
rz(1.5440905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3301311) q[1];
sx q[1];
rz(-0.68492678) q[1];
sx q[1];
rz(-1.0792988) q[1];
x q[2];
rz(-0.052549683) q[3];
sx q[3];
rz(-0.71501117) q[3];
sx q[3];
rz(-2.8831567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22471681) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(-1.3678331) q[2];
rz(1.2380838) q[3];
sx q[3];
rz(-1.5725458) q[3];
sx q[3];
rz(0.21627538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5991768) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(-0.61781484) q[0];
rz(2.0333596) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(-2.2584426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9882433) q[0];
sx q[0];
rz(-2.3267965) q[0];
sx q[0];
rz(1.5809466) q[0];
rz(-1.4083715) q[2];
sx q[2];
rz(-1.139763) q[2];
sx q[2];
rz(1.169432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1146436) q[1];
sx q[1];
rz(-2.0609239) q[1];
sx q[1];
rz(-3.0672706) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1843076) q[3];
sx q[3];
rz(-2.4924879) q[3];
sx q[3];
rz(1.7136991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21657476) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(-1.7803171) q[2];
rz(-1.2733634) q[3];
sx q[3];
rz(-1.4342156) q[3];
sx q[3];
rz(2.1586965) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044416044) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(0.62028766) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(1.1995859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.967554) q[0];
sx q[0];
rz(-3.0332374) q[0];
sx q[0];
rz(1.3757785) q[0];
x q[1];
rz(2.994147) q[2];
sx q[2];
rz(-1.3612399) q[2];
sx q[2];
rz(3.0990019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6288249) q[1];
sx q[1];
rz(-0.56124748) q[1];
sx q[1];
rz(-2.7732549) q[1];
x q[2];
rz(2.7971917) q[3];
sx q[3];
rz(-1.4592429) q[3];
sx q[3];
rz(-1.2408181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(-0.29117584) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(0.62180716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156242) q[0];
sx q[0];
rz(-1.7142897) q[0];
sx q[0];
rz(0.17844644) q[0];
rz(-2.2715691) q[1];
sx q[1];
rz(-2.2585637) q[1];
sx q[1];
rz(-0.80284405) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3960796) q[0];
sx q[0];
rz(-0.95142309) q[0];
sx q[0];
rz(-1.6016763) q[0];
x q[1];
rz(-0.58391352) q[2];
sx q[2];
rz(-1.3656934) q[2];
sx q[2];
rz(2.5896304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1969779) q[1];
sx q[1];
rz(-0.77003819) q[1];
sx q[1];
rz(0.47702392) q[1];
x q[2];
rz(-0.25814806) q[3];
sx q[3];
rz(-0.95883639) q[3];
sx q[3];
rz(1.1560837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37658438) q[2];
sx q[2];
rz(-1.8439801) q[2];
sx q[2];
rz(-0.89361781) q[2];
rz(1.5417967) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(-0.60683513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1641418) q[0];
sx q[0];
rz(-0.41005382) q[0];
sx q[0];
rz(2.9264911) q[0];
rz(-2.2970301) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(-2.3473158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85958033) q[0];
sx q[0];
rz(-1.5589801) q[0];
sx q[0];
rz(-3.1150453) q[0];
rz(0.5735917) q[2];
sx q[2];
rz(-2.6124138) q[2];
sx q[2];
rz(1.0768138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8848443) q[1];
sx q[1];
rz(-1.0326385) q[1];
sx q[1];
rz(0.10579188) q[1];
rz(-pi) q[2];
rz(0.78171697) q[3];
sx q[3];
rz(-1.8060883) q[3];
sx q[3];
rz(2.5552251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8636785) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(1.7318783) q[2];
rz(0.75285161) q[3];
sx q[3];
rz(-1.1466305) q[3];
sx q[3];
rz(-2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3442605) q[0];
sx q[0];
rz(-2.7334038) q[0];
sx q[0];
rz(0.97287384) q[0];
rz(1.6196039) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(2.5111228) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3469543) q[0];
sx q[0];
rz(-1.7525867) q[0];
sx q[0];
rz(1.7111045) q[0];
rz(1.5233598) q[2];
sx q[2];
rz(-0.15711297) q[2];
sx q[2];
rz(-2.4542798) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18955293) q[1];
sx q[1];
rz(-1.5209043) q[1];
sx q[1];
rz(2.3734809) q[1];
rz(-pi) q[2];
rz(-2.7961066) q[3];
sx q[3];
rz(-2.957323) q[3];
sx q[3];
rz(-2.1994182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.015532739) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(1.3014334) q[2];
rz(-1.5850916) q[3];
sx q[3];
rz(-0.82587487) q[3];
sx q[3];
rz(2.7263156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62366098) q[0];
sx q[0];
rz(-0.15075891) q[0];
sx q[0];
rz(-0.49322042) q[0];
rz(1.8863691) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(0.36466041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0756149) q[0];
sx q[0];
rz(-1.9034874) q[0];
sx q[0];
rz(-0.14603024) q[0];
rz(-0.45367806) q[2];
sx q[2];
rz(-0.96900207) q[2];
sx q[2];
rz(1.2860677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11110567) q[1];
sx q[1];
rz(-2.3685622) q[1];
sx q[1];
rz(1.6513325) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8871898) q[3];
sx q[3];
rz(-1.5324161) q[3];
sx q[3];
rz(-0.61172133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6898592) q[2];
sx q[2];
rz(-1.006459) q[2];
sx q[2];
rz(0.9922007) q[2];
rz(0.36710468) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(-2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51919666) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(0.92338152) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(-1.7708764) q[2];
sx q[2];
rz(-1.4309819) q[2];
sx q[2];
rz(-0.33071721) q[2];
rz(0.70960511) q[3];
sx q[3];
rz(-2.3313739) q[3];
sx q[3];
rz(-0.45562638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
