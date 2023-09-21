OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8960436) q[0];
sx q[0];
rz(-0.25078068) q[0];
sx q[0];
rz(1.8401237) q[0];
x q[1];
rz(2.8561864) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(-2.0602496) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5073587) q[1];
sx q[1];
rz(-2.6886352) q[1];
sx q[1];
rz(-1.5276315) q[1];
rz(-0.54361312) q[3];
sx q[3];
rz(-1.3136778) q[3];
sx q[3];
rz(-1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(0.84709644) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0711489) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8424938) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-2.9849844) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-0.54547711) q[1];
x q[2];
rz(1.133425) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(-0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(-0.5273231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.880078) q[0];
sx q[0];
rz(-1.563619) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2450561) q[2];
sx q[2];
rz(-0.92391652) q[2];
sx q[2];
rz(1.6270454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1151162) q[1];
sx q[1];
rz(-2.3079702) q[1];
sx q[1];
rz(-1.2181746) q[1];
rz(-pi) q[2];
rz(0.35964386) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(-1.5776724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-0.8262659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59506455) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(2.3331649) q[0];
rz(2.4414908) q[2];
sx q[2];
rz(-1.5296017) q[2];
sx q[2];
rz(0.61566478) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5570045) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(1.9083379) q[1];
x q[2];
rz(2.6310001) q[3];
sx q[3];
rz(-0.69805745) q[3];
sx q[3];
rz(-1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78695801) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.5251072) q[0];
rz(-2.4992906) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(1.9833267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51296556) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(-2.1668424) q[1];
rz(-pi) q[2];
rz(1.5836309) q[3];
sx q[3];
rz(-1.7463297) q[3];
sx q[3];
rz(2.3042391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(2.9615013) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20443944) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(-2.252749) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9526268) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(-1.8408066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.084701531) q[1];
sx q[1];
rz(-1.4086205) q[1];
sx q[1];
rz(2.8497036) q[1];
x q[2];
rz(-1.5363541) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(0.29355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95773762) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(-3.1013536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043645488) q[0];
sx q[0];
rz(-2.8890434) q[0];
sx q[0];
rz(-1.2388242) q[0];
x q[1];
rz(-1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(1.6183311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.096958728) q[1];
sx q[1];
rz(-2.5068336) q[1];
sx q[1];
rz(-1.8083014) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(-0.039507341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818427) q[0];
sx q[0];
rz(-1.2568226) q[0];
sx q[0];
rz(-2.7008675) q[0];
x q[1];
rz(2.2340441) q[2];
sx q[2];
rz(-0.86576701) q[2];
sx q[2];
rz(0.384534) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2124894) q[1];
sx q[1];
rz(-2.266532) q[1];
sx q[1];
rz(1.5284258) q[1];
rz(1.1106311) q[3];
sx q[3];
rz(-2.4880829) q[3];
sx q[3];
rz(-1.6516407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.8766778) q[1];
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
x q[1];
rz(2.5876797) q[2];
sx q[2];
rz(-1.950112) q[2];
sx q[2];
rz(-1.6162789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1688655) q[1];
sx q[1];
rz(-1.4251514) q[1];
sx q[1];
rz(-1.5345854) q[1];
rz(-pi) q[2];
rz(1.7725138) q[3];
sx q[3];
rz(-1.1551876) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.5309155) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9368162) q[0];
sx q[0];
rz(-1.6813155) q[0];
sx q[0];
rz(0.21893455) q[0];
rz(-2.1572838) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(-0.078660065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5058532) q[1];
sx q[1];
rz(-2.3931599) q[1];
sx q[1];
rz(1.9745419) q[1];
rz(-1.3947992) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(2.396092) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-2.5901569) q[2];
sx q[2];
rz(-2.3497992) q[2];
sx q[2];
rz(-2.1247911) q[2];
rz(0.54995723) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];