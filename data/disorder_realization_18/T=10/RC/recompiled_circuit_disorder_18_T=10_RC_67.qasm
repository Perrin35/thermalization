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
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24554907) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.301469) q[0];
rz(-pi) q[1];
rz(-2.8561864) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(-2.0602496) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45935985) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(3.1205936) q[1];
rz(-2.5979795) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(1.4527814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.2343181) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0704437) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(2.0106993) q[0];
rz(-pi) q[1];
rz(2.8424938) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-0.15660827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(2.5961155) q[1];
rz(2.2291525) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(0.57953366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505328) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(-3.1191349) q[0];
rz(-pi) q[1];
rz(-2.7408319) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(-2.1378627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61422435) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(0.36348344) q[1];
x q[2];
rz(1.8970044) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(-1.5184901) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-0.8262659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6711133) q[0];
sx q[0];
rz(-2.2320036) q[0];
sx q[0];
rz(2.3040422) q[0];
rz(-pi) q[1];
x q[1];
rz(0.063886558) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(1.0039767) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0127718) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(0.080288447) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5097333) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5740009) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(-0.15304676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6981814) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(2.0425668) q[1];
rz(-pi) q[2];
rz(2.9660452) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(-2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-2.9615013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9371532) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(-0.8888437) q[0];
rz(-pi) q[1];
rz(0.18896582) q[2];
sx q[2];
rz(-0.97550387) q[2];
sx q[2];
rz(-1.300786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99299586) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-0.51698835) q[1];
rz(-1.6052386) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(0.040239008) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7559164) q[0];
sx q[0];
rz(-1.809281) q[0];
sx q[0];
rz(-0.083906108) q[0];
rz(-pi) q[1];
x q[1];
rz(0.046182403) q[2];
sx q[2];
rz(-1.9765515) q[2];
sx q[2];
rz(-1.6398167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8602627) q[1];
sx q[1];
rz(-1.7107692) q[1];
sx q[1];
rz(0.94957385) q[1];
rz(-pi) q[2];
rz(-0.31733613) q[3];
sx q[3];
rz(-1.279497) q[3];
sx q[3];
rz(-0.81277646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.5054024) q[1];
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
rz(-0.90754857) q[2];
sx q[2];
rz(-0.86576701) q[2];
sx q[2];
rz(-2.7570587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-0.69617747) q[1];
rz(-1.1106311) q[3];
sx q[3];
rz(-0.65350973) q[3];
sx q[3];
rz(1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
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
rz(-1.7130518) q[1];
sx q[1];
rz(1.8766778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919968) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(-2.0977661) q[0];
x q[1];
rz(-0.55391295) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(-1.5253138) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39667323) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(0.14573914) q[1];
x q[2];
rz(-2.7183652) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(1.9870027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-0.35486737) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(2.5852809) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.3532191) q[0];
sx q[0];
rz(1.4575973) q[0];
rz(-2.9173031) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(-1.6155417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36832419) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(2.2776105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3947992) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(0.71184288) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(-2.5682156) q[3];
sx q[3];
rz(-1.8918512) q[3];
sx q[3];
rz(0.20791114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
