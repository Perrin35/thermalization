OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(-0.010308417) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(4.5865321) q[1];
sx q[1];
rz(10.001339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86069292) q[0];
sx q[0];
rz(-1.4497978) q[0];
sx q[0];
rz(-0.31009407) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2498807) q[2];
sx q[2];
rz(-2.0713965) q[2];
sx q[2];
rz(-0.5989738) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3763872) q[1];
sx q[1];
rz(-0.55843267) q[1];
sx q[1];
rz(-0.33513481) q[1];
x q[2];
rz(2.9750132) q[3];
sx q[3];
rz(-2.2969679) q[3];
sx q[3];
rz(-2.8960901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.6813797) q[2];
sx q[2];
rz(-1.8560393) q[2];
rz(-1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-1.0514528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-0.24573627) q[0];
rz(2.083678) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(-0.33448514) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933788) q[0];
sx q[0];
rz(-2.0248027) q[0];
sx q[0];
rz(-0.55310849) q[0];
rz(1.9742786) q[2];
sx q[2];
rz(-2.2902787) q[2];
sx q[2];
rz(-0.011842273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2776371) q[1];
sx q[1];
rz(-2.5647128) q[1];
sx q[1];
rz(-0.83630348) q[1];
rz(-2.7354419) q[3];
sx q[3];
rz(-0.63780071) q[3];
sx q[3];
rz(0.33406182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1841396) q[2];
sx q[2];
rz(-2.2440971) q[2];
sx q[2];
rz(-1.755836) q[2];
rz(-2.1615248) q[3];
sx q[3];
rz(-0.46531427) q[3];
sx q[3];
rz(0.050962713) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053213) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(1.7514239) q[0];
rz(2.7888489) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(2.5211451) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0095795) q[0];
sx q[0];
rz(-1.6452098) q[0];
sx q[0];
rz(3.0509909) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98288713) q[2];
sx q[2];
rz(-1.4841784) q[2];
sx q[2];
rz(-0.69332214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3603649) q[1];
sx q[1];
rz(-1.8958127) q[1];
sx q[1];
rz(-0.31045352) q[1];
rz(-pi) q[2];
rz(0.51692112) q[3];
sx q[3];
rz(-2.1955588) q[3];
sx q[3];
rz(0.79510915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0188401) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(-1.8943465) q[2];
rz(2.9774169) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71488798) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(1.2870652) q[0];
rz(2.5413068) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(-0.90369019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0125192) q[0];
sx q[0];
rz(-0.12618574) q[0];
sx q[0];
rz(2.587785) q[0];
x q[1];
rz(-0.37995423) q[2];
sx q[2];
rz(-1.5075052) q[2];
sx q[2];
rz(-1.582887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43634448) q[1];
sx q[1];
rz(-1.3911934) q[1];
sx q[1];
rz(-0.257538) q[1];
rz(-pi) q[2];
rz(-1.8925531) q[3];
sx q[3];
rz(-1.1776973) q[3];
sx q[3];
rz(2.2571795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6220182) q[2];
sx q[2];
rz(-0.833424) q[2];
sx q[2];
rz(0.77345094) q[2];
rz(1.2889688) q[3];
sx q[3];
rz(-2.6123612) q[3];
sx q[3];
rz(-0.73498631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2499823) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(2.6494359) q[0];
rz(-0.53897578) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(1.0999365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019177) q[0];
sx q[0];
rz(-1.2309886) q[0];
sx q[0];
rz(-1.2948582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1824838) q[2];
sx q[2];
rz(-1.5796184) q[2];
sx q[2];
rz(-0.1977284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7485986) q[1];
sx q[1];
rz(-1.4778607) q[1];
sx q[1];
rz(2.3307021) q[1];
rz(-pi) q[2];
rz(-0.57657974) q[3];
sx q[3];
rz(-2.5825273) q[3];
sx q[3];
rz(0.41043974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34053549) q[2];
sx q[2];
rz(-2.8057782) q[2];
sx q[2];
rz(-2.578793) q[2];
rz(2.4417012) q[3];
sx q[3];
rz(-1.8507345) q[3];
sx q[3];
rz(-2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7971147) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(-2.7864454) q[0];
rz(-1.9225325) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(-0.452279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1817976) q[0];
sx q[0];
rz(-0.99645146) q[0];
sx q[0];
rz(-2.6622165) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3269749) q[2];
sx q[2];
rz(-0.53485188) q[2];
sx q[2];
rz(1.8823106) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2792795) q[1];
sx q[1];
rz(-1.8932027) q[1];
sx q[1];
rz(3.1161948) q[1];
x q[2];
rz(-1.7132492) q[3];
sx q[3];
rz(-0.94678426) q[3];
sx q[3];
rz(-2.595866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63048116) q[2];
sx q[2];
rz(-2.6177572) q[2];
sx q[2];
rz(-3.0044978) q[2];
rz(1.5824205) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(0.77663511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95098507) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(1.5964339) q[0];
rz(-2.6243788) q[1];
sx q[1];
rz(-1.9904174) q[1];
sx q[1];
rz(-0.98701611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0407437) q[0];
sx q[0];
rz(-0.13253875) q[0];
sx q[0];
rz(-0.1610456) q[0];
rz(-pi) q[1];
rz(-2.8124078) q[2];
sx q[2];
rz(-2.1527157) q[2];
sx q[2];
rz(2.0792368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2199041) q[1];
sx q[1];
rz(-0.35383407) q[1];
sx q[1];
rz(0.71531957) q[1];
x q[2];
rz(1.889398) q[3];
sx q[3];
rz(-1.650839) q[3];
sx q[3];
rz(-2.6032676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3551657) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-3.1332341) q[2];
rz(-2.9110294) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(-1.6647388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1281328) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(-1.531456) q[0];
rz(-1.5229335) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(1.5628372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3561479) q[0];
sx q[0];
rz(-1.7164937) q[0];
sx q[0];
rz(2.5375318) q[0];
x q[1];
rz(-2.7655536) q[2];
sx q[2];
rz(-1.4362772) q[2];
sx q[2];
rz(-2.0903115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48846515) q[1];
sx q[1];
rz(-1.152895) q[1];
sx q[1];
rz(2.6326219) q[1];
rz(-0.66791299) q[3];
sx q[3];
rz(-1.3156943) q[3];
sx q[3];
rz(0.15984838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83850399) q[2];
sx q[2];
rz(-1.1469301) q[2];
sx q[2];
rz(-3.1040891) q[2];
rz(0.23669067) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(-1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26838747) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(-1.0733806) q[0];
rz(-2.7997596) q[1];
sx q[1];
rz(-0.57917246) q[1];
sx q[1];
rz(2.5198708) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1960012) q[0];
sx q[0];
rz(-2.669793) q[0];
sx q[0];
rz(-0.5990754) q[0];
rz(-0.15839496) q[2];
sx q[2];
rz(-2.4753053) q[2];
sx q[2];
rz(1.0384384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1248884) q[1];
sx q[1];
rz(-1.0811662) q[1];
sx q[1];
rz(-1.917065) q[1];
rz(-pi) q[2];
rz(1.6584431) q[3];
sx q[3];
rz(-1.7109799) q[3];
sx q[3];
rz(-0.46938458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5198034) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(2.5780047) q[2];
rz(-1.311519) q[3];
sx q[3];
rz(-1.9026285) q[3];
sx q[3];
rz(2.3254976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.1821063) q[0];
sx q[0];
rz(-2.2725548) q[0];
sx q[0];
rz(-2.2667789) q[0];
rz(1.733571) q[1];
sx q[1];
rz(-2.4751016) q[1];
sx q[1];
rz(2.5061238) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5498813) q[0];
sx q[0];
rz(-2.3259865) q[0];
sx q[0];
rz(0.41378234) q[0];
rz(-pi) q[1];
rz(-2.8671625) q[2];
sx q[2];
rz(-1.1688902) q[2];
sx q[2];
rz(0.91584554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2848031) q[1];
sx q[1];
rz(-1.1691879) q[1];
sx q[1];
rz(1.0076341) q[1];
x q[2];
rz(-2.3395679) q[3];
sx q[3];
rz(-2.4200508) q[3];
sx q[3];
rz(0.82593838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1982939) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(0.80149209) q[2];
rz(1.21579) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(0.041778684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6112678) q[0];
sx q[0];
rz(-1.7014736) q[0];
sx q[0];
rz(0.3076719) q[0];
rz(1.2904185) q[1];
sx q[1];
rz(-0.94480521) q[1];
sx q[1];
rz(0.31029846) q[1];
rz(1.4353095) q[2];
sx q[2];
rz(-1.7510964) q[2];
sx q[2];
rz(-0.68000678) q[2];
rz(2.967703) q[3];
sx q[3];
rz(-2.4774543) q[3];
sx q[3];
rz(3.1042393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
