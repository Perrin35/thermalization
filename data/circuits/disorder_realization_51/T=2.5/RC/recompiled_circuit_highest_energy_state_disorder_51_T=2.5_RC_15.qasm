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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2631419) q[0];
sx q[0];
rz(-2.1241444) q[0];
sx q[0];
rz(-0.31991495) q[0];
rz(-pi) q[1];
rz(3.0906244) q[2];
sx q[2];
rz(-1.5545003) q[2];
sx q[2];
rz(2.7560134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50128776) q[1];
sx q[1];
rz(-0.03420642) q[1];
sx q[1];
rz(2.621719) q[1];
x q[2];
rz(2.9004495) q[3];
sx q[3];
rz(-2.8297348) q[3];
sx q[3];
rz(-1.57392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2179541) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(2.9812319) q[2];
rz(1.977836) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9775951) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(0.98597041) q[0];
rz(1.5522955) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(-1.5812965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47871209) q[0];
sx q[0];
rz(-1.1887738) q[0];
sx q[0];
rz(1.0555313) q[0];
rz(-pi) q[1];
rz(-1.5918757) q[2];
sx q[2];
rz(-2.1651742) q[2];
sx q[2];
rz(0.057986857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5039798) q[1];
sx q[1];
rz(-0.39800522) q[1];
sx q[1];
rz(-0.43136533) q[1];
rz(-pi) q[2];
rz(-0.12388568) q[3];
sx q[3];
rz(-2.200211) q[3];
sx q[3];
rz(-2.7567425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5220149) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(2.1066693) q[2];
rz(2.6400635) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(2.1653304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72306776) q[0];
sx q[0];
rz(-1.1534961) q[0];
sx q[0];
rz(1.204741) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(-0.14030309) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36049309) q[0];
sx q[0];
rz(-1.8331967) q[0];
sx q[0];
rz(-0.95376063) q[0];
rz(-2.7195752) q[2];
sx q[2];
rz(-2.2520503) q[2];
sx q[2];
rz(-1.2764507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.122815) q[1];
sx q[1];
rz(-0.63624708) q[1];
sx q[1];
rz(-2.754209) q[1];
x q[2];
rz(3.0012673) q[3];
sx q[3];
rz(-1.9092536) q[3];
sx q[3];
rz(-2.942977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77421618) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(2.9191169) q[2];
rz(-0.065464822) q[3];
sx q[3];
rz(-1.8583349) q[3];
sx q[3];
rz(2.2953667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533933) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(0.54704332) q[0];
rz(-2.8833) q[1];
sx q[1];
rz(-0.02198418) q[1];
sx q[1];
rz(-2.8000854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753174) q[0];
sx q[0];
rz(-1.5721653) q[0];
sx q[0];
rz(1.5642691) q[0];
rz(1.9893622) q[2];
sx q[2];
rz(-1.2593049) q[2];
sx q[2];
rz(2.3626987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66879067) q[1];
sx q[1];
rz(-2.3040132) q[1];
sx q[1];
rz(1.1146783) q[1];
rz(-pi) q[2];
rz(0.47291748) q[3];
sx q[3];
rz(-2.0295791) q[3];
sx q[3];
rz(2.6948351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7121938) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(0.94924259) q[2];
rz(-2.3927169) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(-3.0794028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267938) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(1.5507966) q[0];
rz(-1.7855478) q[1];
sx q[1];
rz(-3.1372012) q[1];
sx q[1];
rz(-0.063025085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5274104) q[0];
sx q[0];
rz(-0.072612397) q[0];
sx q[0];
rz(-2.7114948) q[0];
rz(-2.5777528) q[2];
sx q[2];
rz(-2.3055674) q[2];
sx q[2];
rz(-2.6631402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57471993) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(1.779056) q[1];
rz(0.78181569) q[3];
sx q[3];
rz(-1.2226579) q[3];
sx q[3];
rz(2.0962417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65681347) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(2.4978034) q[2];
rz(-2.2667609) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952045) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(2.9429759) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(0.14828646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7374518) q[0];
sx q[0];
rz(-2.9585529) q[0];
sx q[0];
rz(3.1279081) q[0];
rz(-2.7303549) q[2];
sx q[2];
rz(-1.6401575) q[2];
sx q[2];
rz(0.81024059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1634341) q[1];
sx q[1];
rz(-0.27039195) q[1];
sx q[1];
rz(-0.52537523) q[1];
x q[2];
rz(0.96252302) q[3];
sx q[3];
rz(-1.4279162) q[3];
sx q[3];
rz(2.8382728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6716914) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(-0.12413231) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-2.997213) q[3];
sx q[3];
rz(-2.9706484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588876) q[0];
sx q[0];
rz(-3.0165065) q[0];
sx q[0];
rz(2.4001154) q[0];
rz(0.28400907) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(0.31518087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3819987) q[0];
sx q[0];
rz(-1.505251) q[0];
sx q[0];
rz(3.1146953) q[0];
x q[1];
rz(1.9839331) q[2];
sx q[2];
rz(-0.48309946) q[2];
sx q[2];
rz(2.536762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6446723) q[1];
sx q[1];
rz(-0.61750353) q[1];
sx q[1];
rz(0.27568494) q[1];
x q[2];
rz(1.3471425) q[3];
sx q[3];
rz(-1.5478494) q[3];
sx q[3];
rz(-3.0815041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(2.4197253) q[2];
rz(0.38247153) q[3];
sx q[3];
rz(-1.1451274) q[3];
sx q[3];
rz(-1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.5751936) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(1.5750634) q[0];
rz(2.9381835) q[1];
sx q[1];
rz(-1.2982439) q[1];
sx q[1];
rz(-0.64483109) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5602556) q[0];
sx q[0];
rz(-1.7783209) q[0];
sx q[0];
rz(3.1330714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5459841) q[2];
sx q[2];
rz(-2.1729219) q[2];
sx q[2];
rz(2.7393722) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1496274) q[1];
sx q[1];
rz(-1.0334762) q[1];
sx q[1];
rz(3.0558636) q[1];
rz(0.73904343) q[3];
sx q[3];
rz(-2.5151281) q[3];
sx q[3];
rz(1.5193957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5440392) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(0.37975797) q[2];
rz(-2.0960268) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(1.9426965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7757292) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(1.3566383) q[0];
rz(-2.7011073) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(-0.7007362) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-2.252251) q[0];
sx q[0];
rz(0.89618857) q[0];
x q[1];
rz(-0.1980902) q[2];
sx q[2];
rz(-2.2513394) q[2];
sx q[2];
rz(-2.1777505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4266738) q[1];
sx q[1];
rz(-0.85619421) q[1];
sx q[1];
rz(1.8780519) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1327206) q[3];
sx q[3];
rz(-2.7089467) q[3];
sx q[3];
rz(3.0888626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35357722) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(1.2865944) q[2];
rz(0.50518099) q[3];
sx q[3];
rz(-2.6954539) q[3];
sx q[3];
rz(-1.2222458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5995481) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-3.1344963) q[1];
sx q[1];
rz(0.33682987) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754691) q[0];
sx q[0];
rz(-1.2531452) q[0];
sx q[0];
rz(-1.5730412) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9438528) q[2];
sx q[2];
rz(-0.34239951) q[2];
sx q[2];
rz(-1.2417718) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5975128) q[1];
sx q[1];
rz(-1.656053) q[1];
sx q[1];
rz(1.6823177) q[1];
rz(0.41542094) q[3];
sx q[3];
rz(-1.0491228) q[3];
sx q[3];
rz(-3.0565302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6280262) q[2];
sx q[2];
rz(-0.94945532) q[2];
sx q[2];
rz(2.8619518) q[2];
rz(2.5102992) q[3];
sx q[3];
rz(-0.93835962) q[3];
sx q[3];
rz(2.5173371) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414128) q[0];
sx q[0];
rz(-1.5914088) q[0];
sx q[0];
rz(1.7803022) q[0];
rz(-0.77371669) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(-2.5071267) q[2];
sx q[2];
rz(-2.1773311) q[2];
sx q[2];
rz(0.60550634) q[2];
rz(-0.057422765) q[3];
sx q[3];
rz(-0.37682711) q[3];
sx q[3];
rz(3.0750838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
