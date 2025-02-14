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
rz(1.6826001) q[1];
sx q[1];
rz(3.269722) q[1];
sx q[1];
rz(10.092957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39820078) q[0];
sx q[0];
rz(-1.8132134) q[0];
sx q[0];
rz(1.4546118) q[0];
rz(-pi) q[1];
rz(-2.2179888) q[2];
sx q[2];
rz(-1.8353113) q[2];
sx q[2];
rz(2.2326927) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9757802) q[1];
sx q[1];
rz(-2.4620612) q[1];
sx q[1];
rz(-2.2305016) q[1];
rz(-2.6553538) q[3];
sx q[3];
rz(-2.231153) q[3];
sx q[3];
rz(-0.59729353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51840034) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(2.0584959) q[2];
rz(-0.25534233) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(1.9820836) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(0.033578385) q[0];
rz(-2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-0.23695645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9457696) q[0];
sx q[0];
rz(-1.3422215) q[0];
sx q[0];
rz(-1.7128471) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.218538) q[2];
sx q[2];
rz(-2.6389803) q[2];
sx q[2];
rz(-1.1692804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0396871) q[1];
sx q[1];
rz(-2.9475005) q[1];
sx q[1];
rz(-1.7029352) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1501447) q[3];
sx q[3];
rz(-1.1889403) q[3];
sx q[3];
rz(-2.6086311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.529155) q[2];
sx q[2];
rz(-1.0470904) q[2];
sx q[2];
rz(-0.1203514) q[2];
rz(-0.58593166) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878207) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(2.1534488) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(-1.5879226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2399167) q[0];
sx q[0];
rz(-1.5696948) q[0];
sx q[0];
rz(-3.1256846) q[0];
rz(-pi) q[1];
rz(-3.0071909) q[2];
sx q[2];
rz(-1.6391058) q[2];
sx q[2];
rz(-1.3775064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37582477) q[1];
sx q[1];
rz(-1.3096403) q[1];
sx q[1];
rz(3.1118666) q[1];
rz(-pi) q[2];
rz(-2.5738945) q[3];
sx q[3];
rz(-1.9318707) q[3];
sx q[3];
rz(-0.89601024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83739088) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(-1.7406933) q[0];
rz(0.1700302) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(1.9084575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62183731) q[0];
sx q[0];
rz(-0.4404236) q[0];
sx q[0];
rz(-1.8669403) q[0];
rz(-pi) q[1];
rz(2.4242086) q[2];
sx q[2];
rz(-2.0307433) q[2];
sx q[2];
rz(-1.1954952) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94205663) q[1];
sx q[1];
rz(-0.28008533) q[1];
sx q[1];
rz(1.6249318) q[1];
rz(2.6771116) q[3];
sx q[3];
rz(-1.5600403) q[3];
sx q[3];
rz(1.8821723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79973811) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(-0.62961659) q[2];
rz(0.13684212) q[3];
sx q[3];
rz(-2.5119731) q[3];
sx q[3];
rz(-1.7262044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28302309) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(3.0392905) q[0];
rz(-1.4981859) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(-0.35710517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84182318) q[0];
sx q[0];
rz(-1.3477316) q[0];
sx q[0];
rz(-1.3510432) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12320766) q[2];
sx q[2];
rz(-1.5977309) q[2];
sx q[2];
rz(-1.3696485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50430027) q[1];
sx q[1];
rz(-1.8691113) q[1];
sx q[1];
rz(2.6105085) q[1];
x q[2];
rz(2.3042942) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(2.0715947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6685278) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(-1.3225887) q[2];
rz(-1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(0.1639666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598031) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(-1.7033956) q[0];
rz(-1.2403129) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(1.3528489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082093246) q[0];
sx q[0];
rz(-1.664986) q[0];
sx q[0];
rz(-1.9777498) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4022302) q[2];
sx q[2];
rz(-1.050808) q[2];
sx q[2];
rz(-1.2181768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53735087) q[1];
sx q[1];
rz(-2.4927995) q[1];
sx q[1];
rz(2.2082445) q[1];
rz(-pi) q[2];
rz(2.7463673) q[3];
sx q[3];
rz(-1.6754284) q[3];
sx q[3];
rz(0.44236174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9271586) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44535962) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(1.4129289) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-0.63953343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93776751) q[0];
sx q[0];
rz(-0.97903189) q[0];
sx q[0];
rz(2.6861565) q[0];
x q[1];
rz(-1.6463747) q[2];
sx q[2];
rz(-1.7210782) q[2];
sx q[2];
rz(1.7145451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4762791) q[1];
sx q[1];
rz(-1.3060044) q[1];
sx q[1];
rz(-0.65515838) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1763379) q[3];
sx q[3];
rz(-1.4603851) q[3];
sx q[3];
rz(2.4019456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21706906) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(0.2229812) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(2.9862459) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0105791) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(2.9796694) q[0];
rz(2.7936392) q[1];
sx q[1];
rz(-1.2740998) q[1];
sx q[1];
rz(0.23439342) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18412515) q[0];
sx q[0];
rz(-2.1229366) q[0];
sx q[0];
rz(2.2669621) q[0];
x q[1];
rz(0.35104819) q[2];
sx q[2];
rz(-0.75134885) q[2];
sx q[2];
rz(2.4780432) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.80658) q[1];
sx q[1];
rz(-1.3836821) q[1];
sx q[1];
rz(0.92430618) q[1];
x q[2];
rz(-1.5041483) q[3];
sx q[3];
rz(-1.8469583) q[3];
sx q[3];
rz(0.23279382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(0.74328077) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(-2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648935) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(-2.0546761) q[0];
rz(-1.9089606) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(1.8392275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24827458) q[0];
sx q[0];
rz(-2.9030307) q[0];
sx q[0];
rz(-2.6435411) q[0];
x q[1];
rz(-0.35436578) q[2];
sx q[2];
rz(-2.1379037) q[2];
sx q[2];
rz(1.815955) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2757146) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(-0.65071836) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5185561) q[3];
sx q[3];
rz(-1.5554232) q[3];
sx q[3];
rz(-2.8492209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(2.7272398) q[2];
rz(2.3710592) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-2.2079302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82938021) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(2.6729551) q[0];
rz(-1.2736443) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(-2.830107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5783008) q[0];
sx q[0];
rz(-1.5130763) q[0];
sx q[0];
rz(-0.85203926) q[0];
rz(-pi) q[1];
rz(-2.5045372) q[2];
sx q[2];
rz(-0.14547507) q[2];
sx q[2];
rz(-1.5381952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2356687) q[1];
sx q[1];
rz(-1.0836658) q[1];
sx q[1];
rz(3.0887927) q[1];
rz(-1.0823147) q[3];
sx q[3];
rz(-2.0460108) q[3];
sx q[3];
rz(-2.7431874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(-1.7485471) q[2];
rz(-2.2140908) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(0.86021304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(2.9604079) q[1];
sx q[1];
rz(-2.1962427) q[1];
sx q[1];
rz(-0.93999351) q[1];
rz(-1.8306611) q[2];
sx q[2];
rz(-0.54878546) q[2];
sx q[2];
rz(2.6981163) q[2];
rz(-2.2374033) q[3];
sx q[3];
rz(-2.0695569) q[3];
sx q[3];
rz(0.29026779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
