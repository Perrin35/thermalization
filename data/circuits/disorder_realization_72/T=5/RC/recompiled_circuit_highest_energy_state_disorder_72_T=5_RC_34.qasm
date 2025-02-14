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
rz(0.25207818) q[0];
sx q[0];
rz(3.7171465) q[0];
sx q[0];
rz(9.3024749) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487028) q[0];
sx q[0];
rz(-2.1445591) q[0];
sx q[0];
rz(2.696373) q[0];
x q[1];
rz(-0.78661175) q[2];
sx q[2];
rz(-2.9742138) q[2];
sx q[2];
rz(-2.6712863) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9895175) q[1];
sx q[1];
rz(-0.65836009) q[1];
sx q[1];
rz(-2.1170298) q[1];
rz(0.24738048) q[3];
sx q[3];
rz(-1.5906963) q[3];
sx q[3];
rz(-0.54075356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4965839) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(1.4487779) q[2];
rz(-0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(2.9922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377624) q[0];
sx q[0];
rz(-1.2685403) q[0];
sx q[0];
rz(-0.25800905) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(2.9749427) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7014871) q[0];
sx q[0];
rz(-2.1512554) q[0];
sx q[0];
rz(-1.6850182) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7562253) q[2];
sx q[2];
rz(-1.8358942) q[2];
sx q[2];
rz(2.6814658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38063654) q[1];
sx q[1];
rz(-1.6957307) q[1];
sx q[1];
rz(-0.71604587) q[1];
rz(2.1364501) q[3];
sx q[3];
rz(-2.0453302) q[3];
sx q[3];
rz(0.33250721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0717281) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(-1.2321164) q[2];
rz(-2.2169436) q[3];
sx q[3];
rz(-0.93622127) q[3];
sx q[3];
rz(-2.0360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6919747) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(-3.1396507) q[0];
rz(0.034189668) q[1];
sx q[1];
rz(-1.2119774) q[1];
sx q[1];
rz(-1.5984104) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5078093) q[0];
sx q[0];
rz(-2.6144321) q[0];
sx q[0];
rz(0.5306923) q[0];
rz(1.9667186) q[2];
sx q[2];
rz(-2.2668419) q[2];
sx q[2];
rz(3.0484704) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.10101) q[1];
sx q[1];
rz(-1.2173182) q[1];
sx q[1];
rz(-1.5291052) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.218762) q[3];
sx q[3];
rz(-1.7668955) q[3];
sx q[3];
rz(-0.01250532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2446642) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-1.0373235) q[2];
rz(-1.1050998) q[3];
sx q[3];
rz(-1.6048071) q[3];
sx q[3];
rz(2.0028116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26745519) q[0];
sx q[0];
rz(-2.656811) q[0];
sx q[0];
rz(1.7304035) q[0];
rz(0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(-3.0944518) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706076) q[0];
sx q[0];
rz(-1.8520675) q[0];
sx q[0];
rz(-0.60312834) q[0];
rz(-pi) q[1];
rz(1.3762952) q[2];
sx q[2];
rz(-1.10656) q[2];
sx q[2];
rz(0.08438202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.667283) q[1];
sx q[1];
rz(-1.2172926) q[1];
sx q[1];
rz(-2.1992285) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8465038) q[3];
sx q[3];
rz(-2.2722244) q[3];
sx q[3];
rz(-1.6158582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.822927) q[2];
sx q[2];
rz(-1.9959799) q[2];
sx q[2];
rz(-1.218943) q[2];
rz(0.31560358) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(-2.992673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3217992) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(0.34859443) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(-1.3016275) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4322223) q[0];
sx q[0];
rz(-0.8101058) q[0];
sx q[0];
rz(0.10929426) q[0];
x q[1];
rz(1.9867861) q[2];
sx q[2];
rz(-2.092053) q[2];
sx q[2];
rz(0.74300569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5036936) q[1];
sx q[1];
rz(-2.854373) q[1];
sx q[1];
rz(2.2566608) q[1];
x q[2];
rz(-2.1215277) q[3];
sx q[3];
rz(-1.8740843) q[3];
sx q[3];
rz(1.0098977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.18315135) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(-1.6661673) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(-1.6492122) q[0];
rz(-2.1624508) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(-1.917256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164152) q[0];
sx q[0];
rz(-1.5402216) q[0];
sx q[0];
rz(2.9435817) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3516828) q[2];
sx q[2];
rz(-2.1093528) q[2];
sx q[2];
rz(1.6679273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8414014) q[1];
sx q[1];
rz(-0.61392861) q[1];
sx q[1];
rz(-1.5924686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1093418) q[3];
sx q[3];
rz(-0.68285817) q[3];
sx q[3];
rz(-0.40753713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(0.9160308) q[2];
rz(1.3845059) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0093805669) q[0];
sx q[0];
rz(-0.051932422) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-2.7047899) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(-0.11016914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5075275) q[0];
sx q[0];
rz(-1.8095836) q[0];
sx q[0];
rz(0.33383835) q[0];
rz(-pi) q[1];
rz(-2.4442441) q[2];
sx q[2];
rz(-1.9983872) q[2];
sx q[2];
rz(0.22874895) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4136849) q[1];
sx q[1];
rz(-1.2251108) q[1];
sx q[1];
rz(0.0045846049) q[1];
x q[2];
rz(-0.25400587) q[3];
sx q[3];
rz(-2.0585103) q[3];
sx q[3];
rz(-2.8539477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7097912) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(0.94376454) q[2];
rz(0.56378311) q[3];
sx q[3];
rz(-2.0457334) q[3];
sx q[3];
rz(-1.110466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540045) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(2.7960844) q[0];
rz(-2.0461138) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(-1.5798205) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876707) q[0];
sx q[0];
rz(-2.3527288) q[0];
sx q[0];
rz(-0.044290941) q[0];
x q[1];
rz(0.061959668) q[2];
sx q[2];
rz(-0.91455787) q[2];
sx q[2];
rz(-1.3861173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4862389) q[1];
sx q[1];
rz(-1.9353239) q[1];
sx q[1];
rz(-2.2097387) q[1];
x q[2];
rz(-2.6437279) q[3];
sx q[3];
rz(-1.5988873) q[3];
sx q[3];
rz(1.3974578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(-0.53449574) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6119824) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.3105422) q[0];
rz(0.88690859) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(-1.8720253) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55063215) q[0];
sx q[0];
rz(-2.7661588) q[0];
sx q[0];
rz(1.2143137) q[0];
rz(1.1940895) q[2];
sx q[2];
rz(-2.7393638) q[2];
sx q[2];
rz(0.62709432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15358217) q[1];
sx q[1];
rz(-1.8021823) q[1];
sx q[1];
rz(1.2910976) q[1];
rz(-pi) q[2];
rz(-0.2036152) q[3];
sx q[3];
rz(-1.609966) q[3];
sx q[3];
rz(-0.43653442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.612192) q[2];
sx q[2];
rz(-1.3357013) q[2];
sx q[2];
rz(2.9086435) q[2];
rz(-1.8565146) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-0.3860093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(0.097271517) q[0];
rz(2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(-2.9248617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6410206) q[0];
sx q[0];
rz(-1.7413229) q[0];
sx q[0];
rz(0.056030355) q[0];
rz(-pi) q[1];
rz(0.31483908) q[2];
sx q[2];
rz(-0.14893571) q[2];
sx q[2];
rz(0.51477942) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8179124) q[1];
sx q[1];
rz(-0.92880946) q[1];
sx q[1];
rz(-2.2873408) q[1];
rz(-1.2249464) q[3];
sx q[3];
rz(-0.93805635) q[3];
sx q[3];
rz(0.94055292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29297605) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5826591) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(0.16286477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.0804629) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(-2.0551266) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(-1.3475781) q[2];
sx q[2];
rz(-0.032608727) q[2];
sx q[2];
rz(-0.41667117) q[2];
rz(-1.6752406) q[3];
sx q[3];
rz(-2.4267195) q[3];
sx q[3];
rz(0.60703312) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
