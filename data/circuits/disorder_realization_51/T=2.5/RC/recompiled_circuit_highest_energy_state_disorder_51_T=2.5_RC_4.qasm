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
rz(-0.75495523) q[0];
sx q[0];
rz(3.8929953) q[0];
sx q[0];
rz(10.473517) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(-1.7555305) q[1];
sx q[1];
rz(3.0233033) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6615579) q[0];
sx q[0];
rz(-1.2999417) q[0];
sx q[0];
rz(-2.1476905) q[0];
x q[1];
rz(-3.0906244) q[2];
sx q[2];
rz(-1.5870924) q[2];
sx q[2];
rz(2.7560134) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50128776) q[1];
sx q[1];
rz(-0.03420642) q[1];
sx q[1];
rz(0.51987363) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8382073) q[3];
sx q[3];
rz(-1.6441364) q[3];
sx q[3];
rz(2.908542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92363858) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(2.9812319) q[2];
rz(-1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639975) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(0.98597041) q[0];
rz(-1.5522955) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(1.5812965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575657) q[0];
sx q[0];
rz(-1.0959033) q[0];
sx q[0];
rz(2.7090461) q[0];
rz(1.5918757) q[2];
sx q[2];
rz(-2.1651742) q[2];
sx q[2];
rz(3.0836058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4680921) q[1];
sx q[1];
rz(-1.7335658) q[1];
sx q[1];
rz(-2.7767608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2038764) q[3];
sx q[3];
rz(-1.4707397) q[3];
sx q[3];
rz(1.2591187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5220149) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(1.0349234) q[2];
rz(0.50152913) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(-2.1653304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72306776) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.204741) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(-3.0012896) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85976542) q[0];
sx q[0];
rz(-0.66376309) q[0];
sx q[0];
rz(-2.0053932) q[0];
rz(-pi) q[1];
rz(-0.42201747) q[2];
sx q[2];
rz(-2.2520503) q[2];
sx q[2];
rz(1.2764507) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.122815) q[1];
sx q[1];
rz(-0.63624708) q[1];
sx q[1];
rz(0.38738363) q[1];
rz(-pi) q[2];
rz(1.9123593) q[3];
sx q[3];
rz(-1.4384801) q[3];
sx q[3];
rz(-1.7225456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77421618) q[2];
sx q[2];
rz(-2.1420631) q[2];
sx q[2];
rz(2.9191169) q[2];
rz(-3.0761278) q[3];
sx q[3];
rz(-1.8583349) q[3];
sx q[3];
rz(0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881994) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(2.5945493) q[0];
rz(0.25829265) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(2.8000854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40453005) q[0];
sx q[0];
rz(-1.5642691) q[0];
sx q[0];
rz(-0.0013690283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8027868) q[2];
sx q[2];
rz(-1.1735386) q[2];
sx q[2];
rz(0.65639979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8402183) q[1];
sx q[1];
rz(-0.84053381) q[1];
sx q[1];
rz(2.6867742) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6686752) q[3];
sx q[3];
rz(-1.1120136) q[3];
sx q[3];
rz(2.6948351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42939886) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(0.94924259) q[2];
rz(0.74887577) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(-0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6147989) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(1.5507966) q[0];
rz(-1.7855478) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(0.063025085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5274104) q[0];
sx q[0];
rz(-3.0689803) q[0];
sx q[0];
rz(-2.7114948) q[0];
rz(-pi) q[1];
rz(2.1049324) q[2];
sx q[2];
rz(-0.89293081) q[2];
sx q[2];
rz(1.2346083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57471993) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(-1.779056) q[1];
x q[2];
rz(-2.6659417) q[3];
sx q[3];
rz(-2.3010572) q[3];
sx q[3];
rz(2.2851839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4847792) q[2];
sx q[2];
rz(-1.3098837) q[2];
sx q[2];
rz(0.64378929) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463882) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(-0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(2.9933062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7374518) q[0];
sx q[0];
rz(-0.1830398) q[0];
sx q[0];
rz(3.1279081) q[0];
rz(-pi) q[1];
rz(2.9695187) q[2];
sx q[2];
rz(-2.7248757) q[2];
sx q[2];
rz(-2.2234349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4365789) q[1];
sx q[1];
rz(-1.3376029) q[1];
sx q[1];
rz(1.7089273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17354024) q[3];
sx q[3];
rz(-2.1719912) q[3];
sx q[3];
rz(1.3663101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4699012) q[2];
sx q[2];
rz(-2.9000059) q[2];
sx q[2];
rz(3.0174603) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-0.1709443) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588876) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(-0.74147725) q[0];
rz(-2.8575836) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(-2.8264118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18703546) q[0];
sx q[0];
rz(-1.5976359) q[0];
sx q[0];
rz(-1.6363653) q[0];
rz(-pi) q[1];
rz(2.9340247) q[2];
sx q[2];
rz(-1.1313442) q[2];
sx q[2];
rz(-1.0644827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8410346) q[1];
sx q[1];
rz(-1.729064) q[1];
sx q[1];
rz(-2.5421418) q[1];
rz(-pi) q[2];
rz(1.6739082) q[3];
sx q[3];
rz(-2.9167843) q[3];
sx q[3];
rz(1.53035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2726941) q[2];
sx q[2];
rz(-1.0947451) q[2];
sx q[2];
rz(2.4197253) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751936) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(1.5750634) q[0];
rz(-2.9381835) q[1];
sx q[1];
rz(-1.2982439) q[1];
sx q[1];
rz(0.64483109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6015907) q[0];
sx q[0];
rz(-2.9338957) q[0];
sx q[0];
rz(-1.6112441) q[0];
x q[1];
rz(0.60226925) q[2];
sx q[2];
rz(-1.5503484) q[2];
sx q[2];
rz(-1.9589613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7643824) q[1];
sx q[1];
rz(-1.6444211) q[1];
sx q[1];
rz(2.1097357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1172148) q[3];
sx q[3];
rz(-2.0189813) q[3];
sx q[3];
rz(0.77806015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5440392) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(-2.7618347) q[2];
rz(-1.0455658) q[3];
sx q[3];
rz(-1.9087722) q[3];
sx q[3];
rz(1.9426965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7757292) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(-1.7849543) q[0];
rz(-0.44048539) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(0.7007362) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627351) q[0];
sx q[0];
rz(-2.252251) q[0];
sx q[0];
rz(-0.89618857) q[0];
rz(-pi) q[1];
rz(-2.9435025) q[2];
sx q[2];
rz(-0.89025324) q[2];
sx q[2];
rz(-2.1777505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7149188) q[1];
sx q[1];
rz(-0.85619421) q[1];
sx q[1];
rz(1.8780519) q[1];
rz(-2.7089617) q[3];
sx q[3];
rz(-1.5670766) q[3];
sx q[3];
rz(1.5100117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35357722) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(1.8549982) q[2];
rz(0.50518099) q[3];
sx q[3];
rz(-2.6954539) q[3];
sx q[3];
rz(-1.2222458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(1.5995481) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-0.33682987) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3339506) q[0];
sx q[0];
rz(-1.5729289) q[0];
sx q[0];
rz(2.8239408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20614236) q[2];
sx q[2];
rz(-1.2954324) q[2];
sx q[2];
rz(-2.5554267) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1244119) q[1];
sx q[1];
rz(-1.681911) q[1];
sx q[1];
rz(3.0558056) q[1];
x q[2];
rz(-2.7261717) q[3];
sx q[3];
rz(-1.0491228) q[3];
sx q[3];
rz(-3.0565302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51356641) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(0.27964082) q[2];
rz(2.5102992) q[3];
sx q[3];
rz(-0.93835962) q[3];
sx q[3];
rz(2.5173371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414128) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(2.367876) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(2.28188) q[2];
sx q[2];
rz(-1.0621241) q[2];
sx q[2];
rz(-0.56806628) q[2];
rz(-3.0841699) q[3];
sx q[3];
rz(-2.7647655) q[3];
sx q[3];
rz(-0.066508807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
