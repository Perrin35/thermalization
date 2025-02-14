OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.33136) q[0];
sx q[0];
rz(-1.2547837) q[0];
sx q[0];
rz(0.64594185) q[0];
rz(-1.4589925) q[1];
sx q[1];
rz(-0.12812935) q[1];
sx q[1];
rz(2.473414) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0872535) q[0];
sx q[0];
rz(-0.26832661) q[0];
sx q[0];
rz(0.4383723) q[0];
rz(-2.8142848) q[2];
sx q[2];
rz(-2.1919554) q[2];
sx q[2];
rz(-2.6747764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1937374) q[1];
sx q[1];
rz(-1.9661709) q[1];
sx q[1];
rz(2.1389524) q[1];
rz(2.1124437) q[3];
sx q[3];
rz(-2.3437269) q[3];
sx q[3];
rz(-0.11395479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6231923) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(-2.0584959) q[2];
rz(-2.8862503) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(2.0032739) q[1];
sx q[1];
rz(-2.1470224) q[1];
sx q[1];
rz(2.9046362) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4073674) q[0];
sx q[0];
rz(-1.7091284) q[0];
sx q[0];
rz(0.23081918) q[0];
rz(-2.9541624) q[2];
sx q[2];
rz(-2.0399562) q[2];
sx q[2];
rz(0.7721061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8023876) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(-1.3783546) q[1];
rz(-0.93795012) q[3];
sx q[3];
rz(-0.68162912) q[3];
sx q[3];
rz(1.5555895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(3.0212413) q[2];
rz(0.58593166) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2633857) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(-2.1534488) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-0.71690503) q[1];
sx q[1];
rz(-1.5536701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9016759) q[0];
sx q[0];
rz(-1.5696948) q[0];
sx q[0];
rz(0.015908054) q[0];
rz(-1.5018671) q[2];
sx q[2];
rz(-1.7048827) q[2];
sx q[2];
rz(-0.20251911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9389438) q[1];
sx q[1];
rz(-1.5420785) q[1];
sx q[1];
rz(1.8320626) q[1];
x q[2];
rz(-2.5293143) q[3];
sx q[3];
rz(-2.4796072) q[3];
sx q[3];
rz(0.16890165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3042018) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(-0.070579441) q[2];
rz(-2.6523759) q[3];
sx q[3];
rz(-0.990812) q[3];
sx q[3];
rz(-0.14744559) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(1.7406933) q[0];
rz(-0.1700302) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.9084575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2182541) q[0];
sx q[0];
rz(-1.4460576) q[0];
sx q[0];
rz(-1.1472923) q[0];
x q[1];
rz(-0.64575671) q[2];
sx q[2];
rz(-0.82953549) q[2];
sx q[2];
rz(-0.84596201) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.199536) q[1];
sx q[1];
rz(-0.28008533) q[1];
sx q[1];
rz(-1.5166609) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6771116) q[3];
sx q[3];
rz(-1.5815524) q[3];
sx q[3];
rz(-1.2594204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3418545) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(0.62961659) q[2];
rz(-3.0047505) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(-1.4153882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(-3.0392905) q[0];
rz(-1.4981859) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(2.7844875) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461991) q[0];
sx q[0];
rz(-1.7850189) q[0];
sx q[0];
rz(-0.22837436) q[0];
rz(-pi) q[1];
rz(1.5979366) q[2];
sx q[2];
rz(-1.6939591) q[2];
sx q[2];
rz(-0.20448286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50430027) q[1];
sx q[1];
rz(-1.2724814) q[1];
sx q[1];
rz(-2.6105085) q[1];
rz(0.83729845) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(1.069998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47306481) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(1.4332917) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(0.1639666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54356164) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(-1.438197) q[0];
rz(-1.9012798) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(-1.7887438) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082093246) q[0];
sx q[0];
rz(-1.664986) q[0];
sx q[0];
rz(1.9777498) q[0];
rz(-pi) q[1];
rz(-2.6154269) q[2];
sx q[2];
rz(-1.7169098) q[2];
sx q[2];
rz(2.7046159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8554012) q[1];
sx q[1];
rz(-1.0637861) q[1];
sx q[1];
rz(0.42393522) q[1];
rz(1.4574962) q[3];
sx q[3];
rz(-1.1778514) q[3];
sx q[3];
rz(1.9696152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.214434) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(-1.8292142) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44535962) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(-1.4129289) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.5556346) q[1];
sx q[1];
rz(0.63953343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7752658) q[0];
sx q[0];
rz(-1.1970988) q[0];
sx q[0];
rz(-0.92832066) q[0];
x q[1];
rz(2.9908871) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(3.0091803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6653135) q[1];
sx q[1];
rz(-1.8355882) q[1];
sx q[1];
rz(0.65515838) q[1];
x q[2];
rz(-1.9652548) q[3];
sx q[3];
rz(-1.6812075) q[3];
sx q[3];
rz(-0.73964707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21706906) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-0.45219839) q[2];
rz(-2.9186115) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.2740998) q[1];
sx q[1];
rz(2.9071992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680555) q[0];
sx q[0];
rz(-2.148365) q[0];
sx q[0];
rz(-0.67649354) q[0];
rz(0.35104819) q[2];
sx q[2];
rz(-0.75134885) q[2];
sx q[2];
rz(-0.66354942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.80658) q[1];
sx q[1];
rz(-1.3836821) q[1];
sx q[1];
rz(0.92430618) q[1];
rz(-pi) q[2];
rz(1.5041483) q[3];
sx q[3];
rz(-1.8469583) q[3];
sx q[3];
rz(-0.23279382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648935) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(2.0546761) q[0];
rz(1.2326321) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.8392275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8933181) q[0];
sx q[0];
rz(-0.23856197) q[0];
sx q[0];
rz(2.6435411) q[0];
x q[1];
rz(-0.97424284) q[2];
sx q[2];
rz(-1.8678209) q[2];
sx q[2];
rz(0.048962083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.088868695) q[1];
sx q[1];
rz(-2.0985328) q[1];
sx q[1];
rz(2.0305968) q[1];
rz(-0.62303657) q[3];
sx q[3];
rz(-1.5861694) q[3];
sx q[3];
rz(-2.8492209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46633259) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(0.41435286) q[2];
rz(2.3710592) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-2.2079302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82938021) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(2.6729551) q[0];
rz(1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-0.31148568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5632919) q[0];
sx q[0];
rz(-1.6285163) q[0];
sx q[0];
rz(-0.85203926) q[0];
x q[1];
rz(0.1172322) q[2];
sx q[2];
rz(-1.6571317) q[2];
sx q[2];
rz(0.59938474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7817318) q[1];
sx q[1];
rz(-1.6174498) q[1];
sx q[1];
rz(-1.0830888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0823147) q[3];
sx q[3];
rz(-1.0955819) q[3];
sx q[3];
rz(0.39840523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0353388) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(-1.7485471) q[2];
rz(-2.2140908) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.8504234) q[0];
sx q[0];
rz(-0.32857729) q[0];
sx q[0];
rz(1.487442) q[0];
rz(0.18118478) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(-1.3109315) q[2];
sx q[2];
rz(-2.5928072) q[2];
sx q[2];
rz(-0.44347632) q[2];
rz(2.2929706) q[3];
sx q[3];
rz(-0.80905882) q[3];
sx q[3];
rz(-0.73425135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
