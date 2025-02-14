OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(2.9958041) q[0];
sx q[0];
rz(11.267405) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(4.4823449) q[1];
sx q[1];
rz(9.525099) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31831384) q[0];
sx q[0];
rz(-0.66368503) q[0];
sx q[0];
rz(0.084245988) q[0];
rz(-pi) q[1];
rz(0.60007976) q[2];
sx q[2];
rz(-2.588495) q[2];
sx q[2];
rz(1.2002522) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2130174) q[1];
sx q[1];
rz(-1.8452106) q[1];
sx q[1];
rz(-1.3540512) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8842949) q[3];
sx q[3];
rz(-2.1380205) q[3];
sx q[3];
rz(0.49916609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7413062) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(2.3347704) q[2];
rz(-2.412879) q[3];
sx q[3];
rz(-0.7586793) q[3];
sx q[3];
rz(1.9369283) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5173986) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(2.8345795) q[0];
rz(-1.864805) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(1.101864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80230882) q[0];
sx q[0];
rz(-1.3972939) q[0];
sx q[0];
rz(1.0109148) q[0];
rz(1.5069783) q[2];
sx q[2];
rz(-1.7498657) q[2];
sx q[2];
rz(3.0916391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15436048) q[1];
sx q[1];
rz(-1.7061966) q[1];
sx q[1];
rz(-1.7495278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3526221) q[3];
sx q[3];
rz(-1.0010825) q[3];
sx q[3];
rz(2.8017442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1133984) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(-0.096435189) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(0.69795394) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(-1.6290889) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(-0.2581183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982957) q[0];
sx q[0];
rz(-2.8683897) q[0];
sx q[0];
rz(1.0323204) q[0];
x q[1];
rz(-2.4674008) q[2];
sx q[2];
rz(-3.0140244) q[2];
sx q[2];
rz(-1.3889165) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8363179) q[1];
sx q[1];
rz(-1.5658251) q[1];
sx q[1];
rz(-3.1286376) q[1];
rz(-1.3985004) q[3];
sx q[3];
rz(-1.3329643) q[3];
sx q[3];
rz(-2.3034277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(-3.1094587) q[2];
rz(-0.26767996) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5699919) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(3.1056504) q[1];
sx q[1];
rz(-0.033999559) q[1];
sx q[1];
rz(-0.34119225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0787028) q[0];
sx q[0];
rz(-2.4475532) q[0];
sx q[0];
rz(2.7258123) q[0];
x q[1];
rz(-1.8282537) q[2];
sx q[2];
rz(-2.1419868) q[2];
sx q[2];
rz(1.5866304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0387151) q[1];
sx q[1];
rz(-2.2897118) q[1];
sx q[1];
rz(0.9714367) q[1];
rz(1.1197508) q[3];
sx q[3];
rz(-1.0195135) q[3];
sx q[3];
rz(2.0752843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(1.535447) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75761211) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(-1.0773995) q[0];
rz(-2.7491838) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(-2.1108625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01973132) q[0];
sx q[0];
rz(-1.1543589) q[0];
sx q[0];
rz(2.5520578) q[0];
x q[1];
rz(1.0723128) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(-2.0778401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2490152) q[1];
sx q[1];
rz(-1.4059781) q[1];
sx q[1];
rz(1.7958162) q[1];
rz(2.7983973) q[3];
sx q[3];
rz(-2.0569909) q[3];
sx q[3];
rz(0.23517683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82421676) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(2.3023494) q[2];
rz(2.2281036) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.4588673) q[0];
sx q[0];
rz(-2.904628) q[0];
sx q[0];
rz(1.4759901) q[0];
rz(-0.39235517) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(-2.5501693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05706035) q[0];
sx q[0];
rz(-2.1520237) q[0];
sx q[0];
rz(0.44598647) q[0];
rz(-1.7514958) q[2];
sx q[2];
rz(-1.2998253) q[2];
sx q[2];
rz(2.9621552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.082669584) q[1];
sx q[1];
rz(-1.6367404) q[1];
sx q[1];
rz(1.0643427) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2675681) q[3];
sx q[3];
rz(-1.2748147) q[3];
sx q[3];
rz(-0.10914739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8412987) q[2];
sx q[2];
rz(-0.66606194) q[2];
sx q[2];
rz(2.5692614) q[2];
rz(0.22219292) q[3];
sx q[3];
rz(-0.43155813) q[3];
sx q[3];
rz(-0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38729024) q[0];
sx q[0];
rz(-0.13986762) q[0];
sx q[0];
rz(2.733316) q[0];
rz(-0.73221842) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(0.29762038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8430804) q[0];
sx q[0];
rz(-1.1491927) q[0];
sx q[0];
rz(-1.1888191) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6082798) q[2];
sx q[2];
rz(-2.6465694) q[2];
sx q[2];
rz(1.9062454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.630328) q[1];
sx q[1];
rz(-1.9115834) q[1];
sx q[1];
rz(-2.9326669) q[1];
rz(-0.17642658) q[3];
sx q[3];
rz(-2.4301612) q[3];
sx q[3];
rz(2.3030858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1709661) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(2.4453898) q[2];
rz(-2.2512186) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(-1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-2.9288375) q[0];
rz(-0.46956024) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(0.75417095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096587688) q[0];
sx q[0];
rz(-1.8593328) q[0];
sx q[0];
rz(2.9456784) q[0];
rz(0.32525678) q[2];
sx q[2];
rz(-1.7647247) q[2];
sx q[2];
rz(3.0162899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2060125) q[1];
sx q[1];
rz(-1.4569062) q[1];
sx q[1];
rz(-2.4920032) q[1];
x q[2];
rz(0.76393083) q[3];
sx q[3];
rz(-2.3124933) q[3];
sx q[3];
rz(1.0111077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(0.78224409) q[2];
rz(-1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049147216) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(2.1771722) q[0];
rz(1.8655221) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(1.6395578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3827189) q[0];
sx q[0];
rz(-1.3387464) q[0];
sx q[0];
rz(2.9358572) q[0];
rz(-pi) q[1];
rz(2.4110498) q[2];
sx q[2];
rz(-1.2337832) q[2];
sx q[2];
rz(2.4180129) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.137845) q[1];
sx q[1];
rz(-1.5937576) q[1];
sx q[1];
rz(-2.0883043) q[1];
rz(0.067599452) q[3];
sx q[3];
rz(-1.7177594) q[3];
sx q[3];
rz(0.33892469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8883349) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(-0.9453195) q[2];
rz(2.3156598) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(-0.53502423) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076040529) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(2.3384576) q[0];
rz(-1.5638634) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(2.8520083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591551) q[0];
sx q[0];
rz(-1.4652325) q[0];
sx q[0];
rz(-3.1292874) q[0];
rz(-pi) q[1];
rz(-2.824823) q[2];
sx q[2];
rz(-1.8837591) q[2];
sx q[2];
rz(0.53887689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9601396) q[1];
sx q[1];
rz(-1.795869) q[1];
sx q[1];
rz(-1.9378661) q[1];
rz(-1.2373459) q[3];
sx q[3];
rz(-1.8419918) q[3];
sx q[3];
rz(1.1734133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-2.181459) q[2];
rz(0.58297408) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.690602) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(3.1008537) q[1];
sx q[1];
rz(-2.4650885) q[1];
sx q[1];
rz(-3.0104641) q[1];
rz(-2.0410983) q[2];
sx q[2];
rz(-1.294877) q[2];
sx q[2];
rz(-0.69653947) q[2];
rz(-1.6419646) q[3];
sx q[3];
rz(-1.6215848) q[3];
sx q[3];
rz(3.0231089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
