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
rz(2.0776318) q[0];
sx q[0];
rz(-1.504717) q[0];
sx q[0];
rz(1.328701) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.152359) q[2];
sx q[2];
rz(-0.49058149) q[2];
sx q[2];
rz(1.4197592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.634234) q[1];
sx q[1];
rz(-2.6886352) q[1];
sx q[1];
rz(1.6139612) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47031109) q[3];
sx q[3];
rz(-2.5458126) q[3];
sx q[3];
rz(-2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(2.2944962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9053099) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(-0.73487868) q[0];
rz(-pi) q[1];
rz(-1.8663835) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(0.57397599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-0.54547711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0081677) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(-0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.2615146) q[0];
sx q[0];
rz(-1.5779737) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7408319) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(-2.1378627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(0.36348344) q[1];
x q[2];
rz(1.2445883) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(1.0969485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(-0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-0.8262659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47047939) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(2.3040422) q[0];
rz(0.70010186) q[2];
sx q[2];
rz(-1.5296017) q[2];
sx q[2];
rz(-0.61566478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0127718) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(-3.0613042) q[1];
x q[2];
rz(-2.6310001) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(-1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3591101) q[0];
sx q[0];
rz(-1.5251274) q[0];
sx q[0];
rz(-0.029687667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(2.2861779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2684106) q[1];
sx q[1];
rz(-1.0007443) q[1];
sx q[1];
rz(-2.8084055) q[1];
x q[2];
rz(-2.9660452) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(3.0517975) q[0];
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
rz(-0.20443944) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(-0.8888437) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96713709) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(2.9784163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5345726) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(-1.7400017) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0635707) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95773762) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.4181597) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(-0.040239008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7559164) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(0.083906108) q[0];
rz(-pi) q[1];
rz(0.046182403) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(-1.501776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0446339) q[1];
sx q[1];
rz(-2.5068336) q[1];
sx q[1];
rz(-1.3332913) q[1];
rz(0.31733613) q[3];
sx q[3];
rz(-1.279497) q[3];
sx q[3];
rz(-2.3288162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0597499) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(-0.44072515) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6263528) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(-0.49382892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-3.0909096) q[1];
rz(-pi) q[2];
rz(-1.1106311) q[3];
sx q[3];
rz(-0.65350973) q[3];
sx q[3];
rz(1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(-1.8766778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919968) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(-1.0438265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5876797) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(1.5253138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39667323) q[1];
sx q[1];
rz(-1.6066237) q[1];
sx q[1];
rz(2.9958535) q[1];
x q[2];
rz(0.42616578) q[3];
sx q[3];
rz(-2.6821972) q[3];
sx q[3];
rz(-3.1118432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(0.35486737) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(-0.24542228) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2357764) q[0];
sx q[0];
rz(-2.896744) q[0];
sx q[0];
rz(-2.66923) q[0];
x q[1];
rz(-1.2400869) q[2];
sx q[2];
rz(-2.5291574) q[2];
sx q[2];
rz(-1.9233179) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36832419) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(2.2776105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7467935) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(0.71184288) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(0.57337702) q[3];
sx q[3];
rz(-1.8918512) q[3];
sx q[3];
rz(0.20791114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
