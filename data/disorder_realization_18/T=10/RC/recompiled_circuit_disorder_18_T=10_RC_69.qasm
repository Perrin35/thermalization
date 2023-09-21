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
rz(1.7083038) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5231397) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(-0.068058204) q[0];
rz(2.8561864) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(-2.0602496) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1022542) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(1.1182055) q[1];
rz(-pi) q[2];
rz(0.47031109) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(-2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.9887259) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23628274) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(-2.406714) q[0];
rz(2.8424938) q[2];
sx q[2];
rz(-0.80183376) q[2];
sx q[2];
rz(-2.9849844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9287195) q[1];
sx q[1];
rz(-1.5024912) q[1];
sx q[1];
rz(-0.54547711) q[1];
x q[2];
rz(-2.794572) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.2615146) q[0];
sx q[0];
rz(-1.5779737) q[0];
rz(-2.7408319) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(-1.00373) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4434102) q[1];
sx q[1];
rz(-1.3122307) q[1];
sx q[1];
rz(0.76885389) q[1];
rz(-1.8970044) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(1.0969485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5465281) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(2.3331649) q[0];
rz(0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(-2.5259279) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77376765) q[1];
sx q[1];
rz(-0.34562472) q[1];
sx q[1];
rz(-1.3454382) q[1];
x q[2];
rz(-2.6310001) q[3];
sx q[3];
rz(-0.69805745) q[3];
sx q[3];
rz(-1.5654246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546346) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.5251072) q[0];
rz(-pi) q[1];
rz(2.5768443) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(-2.2861779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6981814) q[1];
sx q[1];
rz(-2.4907506) q[1];
sx q[1];
rz(-1.0990259) q[1];
x q[2];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-2.9615013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0463379) q[0];
sx q[0];
rz(-1.006608) q[0];
sx q[0];
rz(0.5395704) q[0];
rz(-2.1744556) q[2];
sx q[2];
rz(-1.7269616) q[2];
sx q[2];
rz(2.9784163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(0.51698835) q[1];
rz(-pi) q[2];
x q[2];
rz(0.018499231) q[3];
sx q[3];
rz(-2.0634994) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(3.1013536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9366074) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.3315014) q[0];
rz(-pi) q[1];
x q[1];
rz(0.046182403) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(-1.501776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(-2.9700301) q[1];
rz(-pi) q[2];
rz(1.8764898) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0030901) q[2];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.5054024) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36639402) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(-1.226107) q[0];
x q[1];
rz(0.90754857) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(0.384534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2124894) q[1];
sx q[1];
rz(-0.87506064) q[1];
sx q[1];
rz(1.5284258) q[1];
x q[2];
rz(0.96945854) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9547573) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(-2.3507621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5876797) q[2];
sx q[2];
rz(-1.950112) q[2];
sx q[2];
rz(-1.6162789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4134819) q[1];
sx q[1];
rz(-0.15004798) q[1];
sx q[1];
rz(0.24197443) q[1];
rz(-1.3690788) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(-0.49858529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
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
rz(2.2357764) q[0];
sx q[0];
rz(-2.896744) q[0];
sx q[0];
rz(2.66923) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9015058) q[2];
sx q[2];
rz(-0.61243528) q[2];
sx q[2];
rz(1.9233179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97800868) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(-0.34983695) q[1];
x q[2];
rz(-2.0189832) q[3];
sx q[3];
rz(-0.19482329) q[3];
sx q[3];
rz(-0.71988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-1.0829265) q[2];
sx q[2];
rz(-0.91960533) q[2];
sx q[2];
rz(-1.4056924) q[2];
rz(1.1938548) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];