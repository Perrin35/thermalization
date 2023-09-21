OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8960436) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(-1.8401237) q[0];
rz(-pi) q[1];
rz(0.98923367) q[2];
sx q[2];
rz(-2.6510112) q[2];
sx q[2];
rz(-1.4197592) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45935985) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-3.1205936) q[1];
rz(-pi) q[2];
rz(-2.5979795) q[3];
sx q[3];
rz(-1.3136778) q[3];
sx q[3];
rz(1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(-2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(1.1308934) q[0];
rz(-pi) q[1];
rz(1.2752091) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(0.57397599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-0.54547711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3470207) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(-1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910599) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(0.022457794) q[0];
rz(2.7408319) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(2.1378627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61422435) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(0.36348344) q[1];
rz(-pi) q[2];
rz(0.73266352) q[3];
sx q[3];
rz(-1.817173) q[3];
sx q[3];
rz(-2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(0.8262659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6986615) q[0];
sx q[0];
rz(-0.94416617) q[0];
sx q[0];
rz(0.71039623) q[0];
x q[1];
rz(-0.063886558) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(2.137616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5570045) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(-1.2332548) q[1];
rz(2.6310001) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(2.5523394) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782398) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(-0.99476238) q[0];
rz(-pi) q[1];
rz(2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(-0.15304676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6286271) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(0.97475027) q[1];
rz(-pi) q[2];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(-2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20443944) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(0.8888437) q[0];
x q[1];
rz(-0.18896582) q[2];
sx q[2];
rz(-0.97550387) q[2];
sx q[2];
rz(-1.8408066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99299586) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(2.6246043) q[1];
rz(0.018499231) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(0.040239008) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38567625) q[0];
sx q[0];
rz(-1.809281) q[0];
sx q[0];
rz(-3.0576865) q[0];
rz(0.046182403) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(-1.501776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.28133) q[1];
sx q[1];
rz(-1.4308235) q[1];
sx q[1];
rz(2.1920188) q[1];
rz(-1.2651029) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7751986) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.9154857) q[0];
rz(-0.6263528) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(0.49382892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-3.0909096) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8137915) q[3];
sx q[3];
rz(-2.1468688) q[3];
sx q[3];
rz(-2.2097261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9547573) q[0];
sx q[0];
rz(-1.9592013) q[0];
sx q[0];
rz(0.79083058) q[0];
rz(-pi) q[1];
rz(0.64847704) q[2];
sx q[2];
rz(-2.4816374) q[2];
sx q[2];
rz(2.6476268) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7281108) q[1];
sx q[1];
rz(-0.15004798) q[1];
sx q[1];
rz(-2.8996182) q[1];
rz(-pi) q[2];
rz(-0.42616578) q[3];
sx q[3];
rz(-2.6821972) q[3];
sx q[3];
rz(3.1118432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.4642749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2357764) q[0];
sx q[0];
rz(-2.896744) q[0];
sx q[0];
rz(-0.4723627) q[0];
x q[1];
rz(-1.2400869) q[2];
sx q[2];
rz(-0.61243528) q[2];
sx q[2];
rz(1.9233179) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97800868) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-0.34983695) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3947992) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(-1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-2.0586661) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(-1.9477378) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
