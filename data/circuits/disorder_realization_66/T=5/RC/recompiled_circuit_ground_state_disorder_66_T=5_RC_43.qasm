OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(-2.4165805) q[0];
sx q[0];
rz(-0.83591914) q[0];
rz(-2.3140276) q[1];
sx q[1];
rz(-0.27048549) q[1];
sx q[1];
rz(-2.7343813) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044687957) q[0];
sx q[0];
rz(-2.7735307) q[0];
sx q[0];
rz(-0.9442455) q[0];
rz(-0.99518203) q[2];
sx q[2];
rz(-0.42127171) q[2];
sx q[2];
rz(-2.5025299) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29556698) q[1];
sx q[1];
rz(-0.4163308) q[1];
sx q[1];
rz(-2.8041072) q[1];
rz(-2.148197) q[3];
sx q[3];
rz(-0.90809408) q[3];
sx q[3];
rz(-2.4277594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.779458) q[2];
sx q[2];
rz(-0.3388277) q[2];
sx q[2];
rz(-2.1753878) q[2];
rz(-2.2400014) q[3];
sx q[3];
rz(-2.2880771) q[3];
sx q[3];
rz(1.0454267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7411984) q[0];
sx q[0];
rz(-0.61120954) q[0];
sx q[0];
rz(-3.0389767) q[0];
rz(-1.1688894) q[1];
sx q[1];
rz(-0.97016197) q[1];
sx q[1];
rz(0.51365596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1875604) q[0];
sx q[0];
rz(-0.49278773) q[0];
sx q[0];
rz(1.7634747) q[0];
rz(-pi) q[1];
rz(-2.6046126) q[2];
sx q[2];
rz(-1.6936079) q[2];
sx q[2];
rz(-1.7857743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2741435) q[1];
sx q[1];
rz(-1.871487) q[1];
sx q[1];
rz(0.97598981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97764303) q[3];
sx q[3];
rz(-0.50948516) q[3];
sx q[3];
rz(2.2397985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1932842) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(2.9330758) q[2];
rz(-0.68566132) q[3];
sx q[3];
rz(-2.0833569) q[3];
sx q[3];
rz(0.16511495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0971646) q[0];
sx q[0];
rz(-1.8655638) q[0];
sx q[0];
rz(1.2318508) q[0];
rz(0.34230289) q[1];
sx q[1];
rz(-2.0399703) q[1];
sx q[1];
rz(-1.5493578) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7881675) q[0];
sx q[0];
rz(-0.36398712) q[0];
sx q[0];
rz(-0.53382477) q[0];
rz(-2.9804055) q[2];
sx q[2];
rz(-1.4518132) q[2];
sx q[2];
rz(-2.7143402) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0682837) q[1];
sx q[1];
rz(-1.4335634) q[1];
sx q[1];
rz(0.84899606) q[1];
rz(-pi) q[2];
rz(-0.54247673) q[3];
sx q[3];
rz(-2.0869792) q[3];
sx q[3];
rz(-2.8621648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0937097) q[2];
sx q[2];
rz(-1.1992531) q[2];
sx q[2];
rz(2.1236146) q[2];
rz(-1.4114981) q[3];
sx q[3];
rz(-2.394815) q[3];
sx q[3];
rz(1.5311034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133076) q[0];
sx q[0];
rz(-1.7787373) q[0];
sx q[0];
rz(0.12411975) q[0];
rz(0.92503754) q[1];
sx q[1];
rz(-1.8412453) q[1];
sx q[1];
rz(-2.4031924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644965) q[0];
sx q[0];
rz(-2.3214046) q[0];
sx q[0];
rz(-3.0565681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52780789) q[2];
sx q[2];
rz(-2.4619812) q[2];
sx q[2];
rz(-2.2371593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4324146) q[1];
sx q[1];
rz(-1.3644618) q[1];
sx q[1];
rz(-0.56992759) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4703512) q[3];
sx q[3];
rz(-0.99712205) q[3];
sx q[3];
rz(-0.84000444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063252123) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(-1.1032907) q[2];
rz(0.25857806) q[3];
sx q[3];
rz(-1.0753814) q[3];
sx q[3];
rz(-2.8598089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2461808) q[0];
sx q[0];
rz(-2.2639332) q[0];
sx q[0];
rz(-1.8038764) q[0];
rz(2.5216263) q[1];
sx q[1];
rz(-1.4620616) q[1];
sx q[1];
rz(1.5043129) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02649826) q[0];
sx q[0];
rz(-2.8015602) q[0];
sx q[0];
rz(-0.97397255) q[0];
rz(2.4392534) q[2];
sx q[2];
rz(-1.6026671) q[2];
sx q[2];
rz(1.7447217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6474939) q[1];
sx q[1];
rz(-1.1859736) q[1];
sx q[1];
rz(-0.61487756) q[1];
rz(-3.0359836) q[3];
sx q[3];
rz(-2.3007959) q[3];
sx q[3];
rz(0.13326204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52516809) q[2];
sx q[2];
rz(-1.1526356) q[2];
sx q[2];
rz(1.1241414) q[2];
rz(-0.21640402) q[3];
sx q[3];
rz(-0.83306044) q[3];
sx q[3];
rz(1.0696577) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6935317) q[0];
sx q[0];
rz(-2.3275284) q[0];
sx q[0];
rz(0.46367177) q[0];
rz(0.13557735) q[1];
sx q[1];
rz(-0.99628535) q[1];
sx q[1];
rz(-0.3840951) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6679055) q[0];
sx q[0];
rz(-1.1039005) q[0];
sx q[0];
rz(1.4422732) q[0];
rz(-pi) q[1];
rz(-0.76445396) q[2];
sx q[2];
rz(-1.7181509) q[2];
sx q[2];
rz(0.58576983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0277703) q[1];
sx q[1];
rz(-2.7384187) q[1];
sx q[1];
rz(1.1330963) q[1];
rz(1.4469524) q[3];
sx q[3];
rz(-2.5899124) q[3];
sx q[3];
rz(1.967103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42528459) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(2.6728805) q[2];
rz(2.5675755) q[3];
sx q[3];
rz(-1.6711957) q[3];
sx q[3];
rz(0.65042692) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9531517) q[0];
sx q[0];
rz(-1.783239) q[0];
sx q[0];
rz(-0.69001946) q[0];
rz(0.72235876) q[1];
sx q[1];
rz(-1.1411618) q[1];
sx q[1];
rz(-2.3648327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070097797) q[0];
sx q[0];
rz(-0.24452183) q[0];
sx q[0];
rz(0.31425176) q[0];
rz(-2.5525307) q[2];
sx q[2];
rz(-1.8839594) q[2];
sx q[2];
rz(0.89436603) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36464035) q[1];
sx q[1];
rz(-1.9951275) q[1];
sx q[1];
rz(1.2838497) q[1];
rz(0.89069907) q[3];
sx q[3];
rz(-1.5560702) q[3];
sx q[3];
rz(2.9662728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7338099) q[2];
sx q[2];
rz(-1.5668543) q[2];
sx q[2];
rz(0.053827914) q[2];
rz(-0.6221866) q[3];
sx q[3];
rz(-2.6572808) q[3];
sx q[3];
rz(-0.079306451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4902041) q[0];
sx q[0];
rz(-1.1505928) q[0];
sx q[0];
rz(1.401061) q[0];
rz(0.98848629) q[1];
sx q[1];
rz(-1.2683615) q[1];
sx q[1];
rz(2.7338457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7177454) q[0];
sx q[0];
rz(-1.9739986) q[0];
sx q[0];
rz(1.3333596) q[0];
x q[1];
rz(1.7290224) q[2];
sx q[2];
rz(-0.68099772) q[2];
sx q[2];
rz(2.3883377) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.249392) q[1];
sx q[1];
rz(-0.86848753) q[1];
sx q[1];
rz(-2.5125395) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8697168) q[3];
sx q[3];
rz(-1.4375981) q[3];
sx q[3];
rz(-1.3877752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72655788) q[2];
sx q[2];
rz(-2.0099535) q[2];
sx q[2];
rz(-0.69990194) q[2];
rz(2.0578201) q[3];
sx q[3];
rz(-0.86728573) q[3];
sx q[3];
rz(2.9054902) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.938195) q[0];
sx q[0];
rz(-1.6827826) q[0];
sx q[0];
rz(-2.7517125) q[0];
rz(-1.1979206) q[1];
sx q[1];
rz(-0.094466297) q[1];
sx q[1];
rz(0.87497154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1695677) q[0];
sx q[0];
rz(-0.9476757) q[0];
sx q[0];
rz(0.02768691) q[0];
rz(-pi) q[1];
rz(-2.0535371) q[2];
sx q[2];
rz(-0.74626479) q[2];
sx q[2];
rz(-1.3811228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5905206) q[1];
sx q[1];
rz(-1.1391907) q[1];
sx q[1];
rz(1.9533532) q[1];
x q[2];
rz(-1.0040818) q[3];
sx q[3];
rz(-1.9329786) q[3];
sx q[3];
rz(-0.099206533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8209057) q[2];
sx q[2];
rz(-0.69183886) q[2];
sx q[2];
rz(-0.64986491) q[2];
rz(1.1167022) q[3];
sx q[3];
rz(-1.2973123) q[3];
sx q[3];
rz(0.19702774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4138625) q[0];
sx q[0];
rz(-0.05302269) q[0];
sx q[0];
rz(0.19456385) q[0];
rz(2.3692756) q[1];
sx q[1];
rz(-2.0077424) q[1];
sx q[1];
rz(0.17759855) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91084448) q[0];
sx q[0];
rz(-0.0017191942) q[0];
sx q[0];
rz(-1.2079713) q[0];
x q[1];
rz(2.4899876) q[2];
sx q[2];
rz(-2.0984887) q[2];
sx q[2];
rz(1.0678295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.011780364) q[1];
sx q[1];
rz(-0.82158554) q[1];
sx q[1];
rz(-0.66844861) q[1];
x q[2];
rz(-2.4836618) q[3];
sx q[3];
rz(-2.2773491) q[3];
sx q[3];
rz(2.547675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3557768) q[2];
sx q[2];
rz(-0.59712258) q[2];
sx q[2];
rz(-0.43006483) q[2];
rz(1.1108584) q[3];
sx q[3];
rz(-1.4579803) q[3];
sx q[3];
rz(0.41657579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149902) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(-1.1732187) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(2.0963893) q[2];
sx q[2];
rz(-1.3202616) q[2];
sx q[2];
rz(3.1332982) q[2];
rz(-2.2852702) q[3];
sx q[3];
rz(-1.1416937) q[3];
sx q[3];
rz(1.0037224) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
