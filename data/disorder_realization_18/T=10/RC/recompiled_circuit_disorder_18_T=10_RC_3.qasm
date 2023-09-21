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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8960436) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.301469) q[0];
rz(-0.28540622) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(-1.0813431) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5073587) q[1];
sx q[1];
rz(-2.6886352) q[1];
sx q[1];
rz(1.5276315) q[1];
x q[2];
rz(-2.6712816) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(-2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(-2.3117075) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9053099) q[0];
sx q[0];
rz(-2.5295527) q[0];
sx q[0];
rz(0.73487868) q[0];
x q[1];
rz(2.3624624) q[2];
sx q[2];
rz(-1.3574294) q[2];
sx q[2];
rz(-1.2029635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8250679) q[1];
sx q[1];
rz(-2.1148588) q[1];
sx q[1];
rz(1.6506509) q[1];
rz(2.2291525) q[3];
sx q[3];
rz(-0.53386253) q[3];
sx q[3];
rz(2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8427211) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(-0.30928916) q[0];
rz(-pi) q[1];
rz(1.8965365) q[2];
sx q[2];
rz(-0.92391652) q[2];
sx q[2];
rz(1.5145472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0264764) q[1];
sx q[1];
rz(-2.3079702) q[1];
sx q[1];
rz(-1.923418) q[1];
rz(-pi) q[2];
rz(0.73266352) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(-0.25783595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-2.3153268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59506455) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(-2.3331649) q[0];
rz(1.6246396) q[2];
sx q[2];
rz(-0.87140897) q[2];
sx q[2];
rz(2.2211423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.367825) q[1];
sx q[1];
rz(-0.34562472) q[1];
sx q[1];
rz(-1.7961545) q[1];
rz(-0.51059254) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(0.61757913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782398) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(2.1468303) q[0];
rz(-pi) q[1];
rz(-2.4992906) q[2];
sx q[2];
rz(-1.9335434) q[2];
sx q[2];
rz(-1.9833267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(2.0425668) q[1];
rz(-pi) q[2];
rz(1.5579617) q[3];
sx q[3];
rz(-1.395263) q[3];
sx q[3];
rz(-0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-2.5036507) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20443944) q[0];
sx q[0];
rz(-2.3817872) q[0];
sx q[0];
rz(-2.252749) q[0];
rz(-1.300235) q[2];
sx q[2];
rz(-2.5205043) q[2];
sx q[2];
rz(-1.6294711) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5345726) q[1];
sx q[1];
rz(-1.2828476) q[1];
sx q[1];
rz(-1.401591) q[1];
rz(1.0780219) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.4181597) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(0.040239008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(-1.8100912) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0954103) q[2];
sx q[2];
rz(-1.9765515) q[2];
sx q[2];
rz(-1.501776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.096958728) q[1];
sx q[1];
rz(-2.5068336) q[1];
sx q[1];
rz(-1.3332913) q[1];
rz(-pi) q[2];
rz(-1.8764898) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.5054024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(1.9154857) q[0];
x q[1];
rz(-0.90754857) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(-0.384534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-0.69617747) q[1];
rz(-pi) q[2];
rz(2.8137915) q[3];
sx q[3];
rz(-2.1468688) q[3];
sx q[3];
rz(2.2097261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1160693) q[0];
sx q[0];
rz(-0.86206305) q[0];
sx q[0];
rz(-0.52225964) q[0];
rz(-pi) q[1];
rz(-0.64847704) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(2.6476268) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1688655) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(-1.5345854) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7725138) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(-0.49858529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(-0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9058162) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(-2.66923) q[0];
rz(-0.98430888) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(0.078660065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5058532) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(1.1670508) q[1];
x q[2];
rz(2.0189832) q[3];
sx q[3];
rz(-2.9467694) q[3];
sx q[3];
rz(-0.71988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(-2.0132813) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.55143572) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(-2.5916354) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];