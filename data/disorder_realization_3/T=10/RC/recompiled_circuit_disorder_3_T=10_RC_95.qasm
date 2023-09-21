OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(-2.3303633) q[0];
x q[1];
rz(2.1404999) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(-2.9542838) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58798446) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(1.9967805) q[1];
x q[2];
rz(1.6730509) q[3];
sx q[3];
rz(-1.3073321) q[3];
sx q[3];
rz(-1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8115494) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-0.051785843) q[0];
rz(-0.45692921) q[2];
sx q[2];
rz(-2.3655342) q[2];
sx q[2];
rz(-1.4571783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0454228) q[1];
sx q[1];
rz(-1.3816557) q[1];
sx q[1];
rz(-1.063709) q[1];
rz(2.5856421) q[3];
sx q[3];
rz(-2.2009938) q[3];
sx q[3];
rz(0.48660183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(-2.8853436) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8921593) q[0];
sx q[0];
rz(-1.306635) q[0];
sx q[0];
rz(2.507472) q[0];
rz(-0.44552866) q[2];
sx q[2];
rz(-2.0085213) q[2];
sx q[2];
rz(-2.90403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5445404) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(-0.87045963) q[1];
rz(-pi) q[2];
rz(2.9647397) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(0.73192275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703092) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(-1.0259823) q[0];
x q[1];
rz(-0.133693) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(2.0267817) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54938984) q[1];
sx q[1];
rz(-1.7648186) q[1];
sx q[1];
rz(-0.26280304) q[1];
rz(-1.6963523) q[3];
sx q[3];
rz(-1.351965) q[3];
sx q[3];
rz(-2.9343176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7933465) q[0];
sx q[0];
rz(-3.0214546) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-2.5227929) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(0.4180846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.080575374) q[1];
sx q[1];
rz(-1.6643545) q[1];
sx q[1];
rz(-2.1993125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6990715) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(2.4181441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(3.0715122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228186) q[0];
sx q[0];
rz(-2.9111324) q[0];
sx q[0];
rz(2.3019058) q[0];
rz(-pi) q[1];
rz(-2.1224535) q[2];
sx q[2];
rz(-0.86691228) q[2];
sx q[2];
rz(0.64955074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3190223) q[1];
sx q[1];
rz(-2.550107) q[1];
sx q[1];
rz(-2.849008) q[1];
x q[2];
rz(0.035404215) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(-0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(2.5202259) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(-3.133657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82848362) q[0];
sx q[0];
rz(-0.61921739) q[0];
sx q[0];
rz(2.6909268) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0646348) q[2];
sx q[2];
rz(-0.75876615) q[2];
sx q[2];
rz(-2.2366692) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4850033) q[1];
sx q[1];
rz(-0.90066972) q[1];
sx q[1];
rz(-0.90799241) q[1];
rz(-pi) q[2];
rz(0.31422024) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(-0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.5023124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973328) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(-pi) q[1];
rz(2.084923) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(-2.4664997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9099721) q[1];
sx q[1];
rz(-2.2044704) q[1];
sx q[1];
rz(-2.2909067) q[1];
rz(-pi) q[2];
rz(0.092076093) q[3];
sx q[3];
rz(-1.7756988) q[3];
sx q[3];
rz(-0.74842194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1736261) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(3.0941512) q[0];
rz(-pi) q[1];
x q[1];
rz(2.142799) q[2];
sx q[2];
rz(-2.1527094) q[2];
sx q[2];
rz(-1.7172161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7855362) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(-0.15798012) q[1];
x q[2];
rz(1.7008408) q[3];
sx q[3];
rz(-1.7107043) q[3];
sx q[3];
rz(3.1062982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(-0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7693217) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-3.0110714) q[0];
rz(-pi) q[1];
rz(2.0595423) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-2.3479455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10510124) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(0.75863104) q[1];
rz(0.15354746) q[3];
sx q[3];
rz(-0.921075) q[3];
sx q[3];
rz(-2.6447907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-0.50921847) q[2];
sx q[2];
rz(-1.5621395) q[2];
sx q[2];
rz(-0.12315673) q[2];
rz(-1.8018467) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
