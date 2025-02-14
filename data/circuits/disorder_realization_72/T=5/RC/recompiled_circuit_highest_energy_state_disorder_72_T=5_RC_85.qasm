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
rz(-2.8895145) q[0];
sx q[0];
rz(-0.57555389) q[0];
sx q[0];
rz(0.12230305) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3275546) q[0];
sx q[0];
rz(-0.71056847) q[0];
sx q[0];
rz(-0.98301218) q[0];
rz(-pi) q[1];
rz(0.11876583) q[2];
sx q[2];
rz(-1.4525754) q[2];
sx q[2];
rz(1.2615276) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.029509228) q[1];
sx q[1];
rz(-1.8942299) q[1];
sx q[1];
rz(0.98677633) q[1];
x q[2];
rz(-2.8942122) q[3];
sx q[3];
rz(-1.5508964) q[3];
sx q[3];
rz(-2.6008391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64500874) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(-1.6928147) q[2];
rz(0.19041666) q[3];
sx q[3];
rz(-2.0864291) q[3];
sx q[3];
rz(2.9922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038302) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(2.8835836) q[0];
rz(1.5500655) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(2.9749427) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23388966) q[0];
sx q[0];
rz(-2.5512716) q[0];
sx q[0];
rz(0.17206828) q[0];
x q[1];
rz(-2.5160997) q[2];
sx q[2];
rz(-0.46395597) q[2];
sx q[2];
rz(2.6044012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7609561) q[1];
sx q[1];
rz(-1.6957307) q[1];
sx q[1];
rz(0.71604587) q[1];
rz(-pi) q[2];
rz(-2.1364501) q[3];
sx q[3];
rz(-1.0962624) q[3];
sx q[3];
rz(-2.8090854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0698645) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(-1.9094763) q[2];
rz(0.92464906) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(2.0360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.449618) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(3.1396507) q[0];
rz(-3.107403) q[1];
sx q[1];
rz(-1.2119774) q[1];
sx q[1];
rz(1.5431822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5078093) q[0];
sx q[0];
rz(-2.6144321) q[0];
sx q[0];
rz(2.6109004) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7091647) q[2];
sx q[2];
rz(-2.3574867) q[2];
sx q[2];
rz(-2.6570005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0405827) q[1];
sx q[1];
rz(-1.9242745) q[1];
sx q[1];
rz(-1.5291052) q[1];
rz(-pi) q[2];
rz(0.9228306) q[3];
sx q[3];
rz(-1.7668955) q[3];
sx q[3];
rz(-0.01250532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2446642) q[2];
sx q[2];
rz(-0.43411532) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741375) q[0];
sx q[0];
rz(-2.656811) q[0];
sx q[0];
rz(-1.7304035) q[0];
rz(2.5022068) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(0.04714084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9586724) q[0];
sx q[0];
rz(-2.483568) q[0];
sx q[0];
rz(0.47112314) q[0];
rz(-pi) q[1];
rz(2.669692) q[2];
sx q[2];
rz(-1.7444897) q[2];
sx q[2];
rz(-1.398441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47430965) q[1];
sx q[1];
rz(-1.2172926) q[1];
sx q[1];
rz(0.94236417) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2940801) q[3];
sx q[3];
rz(-1.7948331) q[3];
sx q[3];
rz(0.23875313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.822927) q[2];
sx q[2];
rz(-1.9959799) q[2];
sx q[2];
rz(-1.9226496) q[2];
rz(2.8259891) q[3];
sx q[3];
rz(-0.054840755) q[3];
sx q[3];
rz(-2.992673) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8197935) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(-2.7929982) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.8399651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5900629) q[0];
sx q[0];
rz(-0.7669391) q[0];
sx q[0];
rz(1.6849031) q[0];
rz(-pi) q[1];
rz(-0.56058863) q[2];
sx q[2];
rz(-1.9288256) q[2];
sx q[2];
rz(2.5303417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5983513) q[1];
sx q[1];
rz(-1.390402) q[1];
sx q[1];
rz(1.3460657) q[1];
rz(-pi) q[2];
rz(2.7896406) q[3];
sx q[3];
rz(-1.0478596) q[3];
sx q[3];
rz(2.3992996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(-1.4754254) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(-0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(-1.4923805) q[0];
rz(0.97914186) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(-1.2243366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92517744) q[0];
sx q[0];
rz(-1.5402216) q[0];
sx q[0];
rz(2.9435817) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7899098) q[2];
sx q[2];
rz(-1.0322399) q[2];
sx q[2];
rz(-1.4736654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30019125) q[1];
sx q[1];
rz(-0.61392861) q[1];
sx q[1];
rz(1.5924686) q[1];
x q[2];
rz(-2.200279) q[3];
sx q[3];
rz(-1.8555897) q[3];
sx q[3];
rz(-2.3465057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41902038) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(-2.2255619) q[2];
rz(-1.7570868) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(-0.70703435) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1322121) q[0];
sx q[0];
rz(-0.051932422) q[0];
sx q[0];
rz(-0.95426553) q[0];
rz(0.43680278) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(3.0314235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5075275) q[0];
sx q[0];
rz(-1.8095836) q[0];
sx q[0];
rz(2.8077543) q[0];
rz(-pi) q[1];
rz(2.5244401) q[2];
sx q[2];
rz(-0.79884702) q[2];
sx q[2];
rz(1.3399194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71437826) q[1];
sx q[1];
rz(-0.34571474) q[1];
sx q[1];
rz(-1.5835254) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0131936) q[3];
sx q[3];
rz(-2.5964649) q[3];
sx q[3];
rz(0.21827182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7097912) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(-0.94376454) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(-0.34550825) q[0];
rz(2.0461138) q[1];
sx q[1];
rz(-1.6637207) q[1];
sx q[1];
rz(-1.5798205) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0559383) q[0];
sx q[0];
rz(-1.5393747) q[0];
sx q[0];
rz(0.78837331) q[0];
x q[1];
rz(-1.4905632) q[2];
sx q[2];
rz(-2.4828665) q[2];
sx q[2];
rz(1.6541437) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53198481) q[1];
sx q[1];
rz(-0.72276211) q[1];
sx q[1];
rz(-1.0015798) q[1];
rz(-pi) q[2];
rz(2.6437279) q[3];
sx q[3];
rz(-1.5988873) q[3];
sx q[3];
rz(-1.3974578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-0.26143917) q[2];
sx q[2];
rz(-1.4307107) q[2];
rz(-2.6070969) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(-0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5296103) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(1.8310504) q[0];
rz(-2.2546841) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(1.2695674) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9717012) q[0];
sx q[0];
rz(-1.2200238) q[0];
sx q[0];
rz(3.0049075) q[0];
rz(-pi) q[1];
rz(1.9475031) q[2];
sx q[2];
rz(-0.40222886) q[2];
sx q[2];
rz(-2.5144983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9880105) q[1];
sx q[1];
rz(-1.8021823) q[1];
sx q[1];
rz(1.2910976) q[1];
x q[2];
rz(-1.6107913) q[3];
sx q[3];
rz(-1.3673395) q[3];
sx q[3];
rz(-2.0154161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5294007) q[2];
sx q[2];
rz(-1.3357013) q[2];
sx q[2];
rz(-2.9086435) q[2];
rz(1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-0.3860093) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(2.9248617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9602338) q[0];
sx q[0];
rz(-2.9621819) q[0];
sx q[0];
rz(1.8852194) q[0];
x q[1];
rz(2.8267536) q[2];
sx q[2];
rz(-2.9926569) q[2];
sx q[2];
rz(-2.6268132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4137553) q[1];
sx q[1];
rz(-1.0168795) q[1];
sx q[1];
rz(-2.3604849) q[1];
x q[2];
rz(1.2249464) q[3];
sx q[3];
rz(-0.93805635) q[3];
sx q[3];
rz(2.2010397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29297605) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5589335) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804629) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(-1.0864661) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(1.6025966) q[2];
sx q[2];
rz(-1.563579) q[2];
sx q[2];
rz(-2.210571) q[2];
rz(-2.2829655) q[3];
sx q[3];
rz(-1.6391907) q[3];
sx q[3];
rz(-1.0427604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
