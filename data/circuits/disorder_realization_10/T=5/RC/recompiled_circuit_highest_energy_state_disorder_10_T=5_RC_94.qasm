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
rz(0.43542433) q[0];
sx q[0];
rz(-0.84917899) q[0];
sx q[0];
rz(0.8325141) q[0];
rz(-0.37021356) q[2];
sx q[2];
rz(-1.5892803) q[2];
sx q[2];
rz(-0.12272515) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0142483) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(-2.8365652) q[1];
rz(-pi) q[2];
rz(2.1006868) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.4640946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(-1.5667685) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(2.4488357) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-0.80802149) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(-0.31879058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2683811) q[0];
sx q[0];
rz(-2.2751774) q[0];
sx q[0];
rz(1.1581139) q[0];
rz(-pi) q[1];
rz(0.484157) q[2];
sx q[2];
rz(-2.6645497) q[2];
sx q[2];
rz(2.8271528) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8242256) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(-2.1721858) q[1];
rz(1.2845006) q[3];
sx q[3];
rz(-1.6914552) q[3];
sx q[3];
rz(2.9243822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0626283) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(2.6641565) q[2];
rz(-1.7834974) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(-0.0099946578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7497712) q[0];
sx q[0];
rz(-1.8056159) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8157522) q[0];
sx q[0];
rz(-1.6597972) q[0];
sx q[0];
rz(0.98317105) q[0];
x q[1];
rz(-1.0933206) q[2];
sx q[2];
rz(-0.68308631) q[2];
sx q[2];
rz(0.25743279) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3542676) q[1];
sx q[1];
rz(-0.72871043) q[1];
sx q[1];
rz(0.72737741) q[1];
rz(0.33379995) q[3];
sx q[3];
rz(-1.085328) q[3];
sx q[3];
rz(-2.7791821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.5072301) q[2];
sx q[2];
rz(-2.0646084) q[2];
rz(-1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(0.25585678) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28427163) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(2.381109) q[0];
rz(-0.35176945) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(0.15028353) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4486769) q[0];
sx q[0];
rz(-1.5741072) q[0];
sx q[0];
rz(-0.0019372367) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9689061) q[2];
sx q[2];
rz(-1.9479247) q[2];
sx q[2];
rz(-2.0443626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6071668) q[1];
sx q[1];
rz(-2.1517506) q[1];
sx q[1];
rz(-0.39926417) q[1];
x q[2];
rz(1.3129381) q[3];
sx q[3];
rz(-2.0153196) q[3];
sx q[3];
rz(-0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.030423) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(0.58491582) q[2];
rz(-3.0758744) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(-1.0544624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2718662) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(0.72750339) q[0];
rz(1.8456521) q[1];
sx q[1];
rz(-1.0546874) q[1];
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
rz(0.22146225) q[0];
rz(-pi) q[1];
rz(-0.35810808) q[2];
sx q[2];
rz(-0.3198959) q[2];
sx q[2];
rz(1.3589771) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3989563) q[1];
sx q[1];
rz(-0.42544895) q[1];
sx q[1];
rz(2.0320973) q[1];
rz(-pi) q[2];
rz(-2.7477164) q[3];
sx q[3];
rz(-1.0303613) q[3];
sx q[3];
rz(0.5462164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(-1.9963416) q[2];
rz(0.20259419) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(-2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458493) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(-1.3289183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8362373) q[0];
sx q[0];
rz(-2.1026045) q[0];
sx q[0];
rz(-1.0145857) q[0];
rz(-2.8064427) q[2];
sx q[2];
rz(-2.9932635) q[2];
sx q[2];
rz(2.2466765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5649453) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(-0.74585657) q[1];
rz(0.22340138) q[3];
sx q[3];
rz(-1.3376682) q[3];
sx q[3];
rz(-2.4883771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7557175) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(-2.7937549) q[2];
rz(-2.8268585) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(0.95112479) q[0];
rz(-0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(-0.94400418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68444618) q[0];
sx q[0];
rz(-1.5553885) q[0];
sx q[0];
rz(-1.2627361) q[0];
x q[1];
rz(-2.6636276) q[2];
sx q[2];
rz(-1.1349808) q[2];
sx q[2];
rz(-2.8675241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72457006) q[1];
sx q[1];
rz(-2.9669683) q[1];
sx q[1];
rz(-1.440446) q[1];
rz(-2.1066426) q[3];
sx q[3];
rz(-0.45848819) q[3];
sx q[3];
rz(1.637984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(-1.7000807) q[2];
rz(0.024638351) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(-0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330419) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(0.1380052) q[0];
rz(-1.0778614) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(2.2053351) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9662921) q[0];
sx q[0];
rz(-1.9387687) q[0];
sx q[0];
rz(-1.247768) q[0];
rz(-pi) q[1];
rz(2.594844) q[2];
sx q[2];
rz(-0.90675747) q[2];
sx q[2];
rz(-0.67686096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1106134) q[1];
sx q[1];
rz(-2.958221) q[1];
sx q[1];
rz(-1.6513837) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3055757) q[3];
sx q[3];
rz(-0.90258963) q[3];
sx q[3];
rz(0.98047719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.998698) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.3188837) q[2];
rz(-0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634624) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(-2.7769856) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(2.2273831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171998) q[0];
sx q[0];
rz(-1.1437609) q[0];
sx q[0];
rz(2.9946505) q[0];
rz(3.0905452) q[2];
sx q[2];
rz(-2.406359) q[2];
sx q[2];
rz(-0.027054199) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57371512) q[1];
sx q[1];
rz(-0.93522954) q[1];
sx q[1];
rz(-1.3997517) q[1];
rz(-pi) q[2];
rz(2.6007341) q[3];
sx q[3];
rz(-1.6990933) q[3];
sx q[3];
rz(0.010426253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8387575) q[2];
sx q[2];
rz(-2.5573533) q[2];
sx q[2];
rz(1.6287104) q[2];
rz(-2.5527939) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(-0.65620667) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673679) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(2.8571416) q[0];
rz(2.4868763) q[1];
sx q[1];
rz(-1.3755211) q[1];
sx q[1];
rz(0.42025748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66061879) q[0];
sx q[0];
rz(-1.7848564) q[0];
sx q[0];
rz(-1.995489) q[0];
rz(-1.1070232) q[2];
sx q[2];
rz(-2.3056185) q[2];
sx q[2];
rz(-0.42586621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3455167) q[1];
sx q[1];
rz(-0.14847595) q[1];
sx q[1];
rz(-0.92904894) q[1];
rz(2.612667) q[3];
sx q[3];
rz(-2.5074337) q[3];
sx q[3];
rz(2.8932539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9922716) q[2];
sx q[2];
rz(-0.48506609) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(-0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(3.0178487) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(-1.704557) q[2];
sx q[2];
rz(-1.2530213) q[2];
sx q[2];
rz(-2.1344013) q[2];
rz(3.0322778) q[3];
sx q[3];
rz(-0.31382618) q[3];
sx q[3];
rz(-1.8574497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
