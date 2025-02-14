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
rz(-2.39019) q[0];
sx q[0];
rz(1.0487392) q[0];
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
rz(-1.0174482) q[0];
sx q[0];
rz(0.31991495) q[0];
x q[1];
rz(-0.050968214) q[2];
sx q[2];
rz(-1.5870924) q[2];
sx q[2];
rz(-2.7560134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1227545) q[1];
sx q[1];
rz(-1.5411106) q[1];
sx q[1];
rz(1.5877941) q[1];
x q[2];
rz(-0.24114313) q[3];
sx q[3];
rz(-2.8297348) q[3];
sx q[3];
rz(1.5676726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92363858) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(-2.9812319) q[2];
rz(-1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(-0.81419301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639975) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(-2.1556222) q[0];
rz(-1.5892971) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(-1.5602962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6742822) q[0];
sx q[0];
rz(-2.5105866) q[0];
sx q[0];
rz(-2.254807) q[0];
rz(0.031173988) q[2];
sx q[2];
rz(-2.5468862) q[2];
sx q[2];
rz(0.020356914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1007019) q[1];
sx q[1];
rz(-1.9305848) q[1];
sx q[1];
rz(1.7448107) q[1];
rz(-pi) q[2];
rz(3.017707) q[3];
sx q[3];
rz(-0.94138161) q[3];
sx q[3];
rz(-0.38485011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(-1.0349234) q[2];
rz(-2.6400635) q[3];
sx q[3];
rz(-0.087567121) q[3];
sx q[3];
rz(-0.97626221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72306776) q[0];
sx q[0];
rz(-1.1534961) q[0];
sx q[0];
rz(1.9368517) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(3.0012896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85976542) q[0];
sx q[0];
rz(-2.4778296) q[0];
sx q[0];
rz(1.1361994) q[0];
x q[1];
rz(0.42201747) q[2];
sx q[2];
rz(-2.2520503) q[2];
sx q[2];
rz(-1.2764507) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2724803) q[1];
sx q[1];
rz(-1.3444053) q[1];
sx q[1];
rz(2.5416944) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14032538) q[3];
sx q[3];
rz(-1.232339) q[3];
sx q[3];
rz(-0.1986157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3673765) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(0.2224758) q[2];
rz(-3.0761278) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533933) q[0];
sx q[0];
rz(-0.024217483) q[0];
sx q[0];
rz(2.5945493) q[0];
rz(-2.8833) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(2.8000854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1662752) q[0];
sx q[0];
rz(-1.5721653) q[0];
sx q[0];
rz(-1.5773236) q[0];
x q[1];
rz(0.33880587) q[2];
sx q[2];
rz(-1.1735386) q[2];
sx q[2];
rz(-0.65639979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.472802) q[1];
sx q[1];
rz(-0.83757949) q[1];
sx q[1];
rz(1.1146783) q[1];
rz(0.82587256) q[3];
sx q[3];
rz(-0.64651981) q[3];
sx q[3];
rz(2.7310128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7121938) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(2.1923501) q[2];
rz(2.3927169) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(0.062189814) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6147989) q[0];
sx q[0];
rz(-3.1068046) q[0];
sx q[0];
rz(1.590796) q[0];
rz(-1.7855478) q[1];
sx q[1];
rz(-3.1372012) q[1];
sx q[1];
rz(-0.063025085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7558799) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(3.0755733) q[0];
rz(-0.56383987) q[2];
sx q[2];
rz(-0.83602521) q[2];
sx q[2];
rz(0.47845248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57471993) q[1];
sx q[1];
rz(-1.5974471) q[1];
sx q[1];
rz(1.3625366) q[1];
rz(-pi) q[2];
rz(1.0980561) q[3];
sx q[3];
rz(-2.2946828) q[3];
sx q[3];
rz(2.9427317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4847792) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(0.64378929) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(-0.8005825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463882) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(-2.7860506) q[0];
rz(0.19861673) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(-0.14828646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4041408) q[0];
sx q[0];
rz(-0.1830398) q[0];
sx q[0];
rz(-3.1279081) q[0];
x q[1];
rz(-1.4951493) q[2];
sx q[2];
rz(-1.9809857) q[2];
sx q[2];
rz(0.73034053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0394888) q[1];
sx q[1];
rz(-1.436427) q[1];
sx q[1];
rz(-2.9062382) q[1];
rz(-pi) q[2];
rz(-2.1790696) q[3];
sx q[3];
rz(-1.7136765) q[3];
sx q[3];
rz(-2.8382728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6716914) q[2];
sx q[2];
rz(-2.9000059) q[2];
sx q[2];
rz(-0.12413231) q[2];
rz(2.5668674) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-2.9706484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588876) q[0];
sx q[0];
rz(-3.0165065) q[0];
sx q[0];
rz(-2.4001154) q[0];
rz(-0.28400907) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(-0.31518087) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3819987) q[0];
sx q[0];
rz(-1.6363417) q[0];
sx q[0];
rz(-3.1146953) q[0];
rz(-pi) q[1];
rz(-1.9839331) q[2];
sx q[2];
rz(-0.48309946) q[2];
sx q[2];
rz(-2.536762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9786392) q[1];
sx q[1];
rz(-2.1617315) q[1];
sx q[1];
rz(-1.7617474) q[1];
x q[2];
rz(3.1180598) q[3];
sx q[3];
rz(-1.7943903) q[3];
sx q[3];
rz(1.6361039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2726941) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(0.72186738) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(-2.0731879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751936) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(1.5665293) q[0];
rz(2.9381835) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(0.64483109) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0087850182) q[0];
sx q[0];
rz(-1.5791348) q[0];
sx q[0];
rz(-1.3632644) q[0];
rz(-pi) q[1];
x q[1];
rz(3.105509) q[2];
sx q[2];
rz(-2.5390194) q[2];
sx q[2];
rz(2.7831603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.825612) q[1];
sx q[1];
rz(-0.54345268) q[1];
sx q[1];
rz(1.4280591) q[1];
rz(-2.4025492) q[3];
sx q[3];
rz(-0.62646455) q[3];
sx q[3];
rz(-1.5193957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5975534) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(-0.37975797) q[2];
rz(2.0960268) q[3];
sx q[3];
rz(-1.9087722) q[3];
sx q[3];
rz(-1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3658635) q[0];
sx q[0];
rz(-0.033670306) q[0];
sx q[0];
rz(-1.3566383) q[0];
rz(-0.44048539) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(2.4408565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25862274) q[0];
sx q[0];
rz(-2.077266) q[0];
sx q[0];
rz(-0.80431403) q[0];
rz(1.8092625) q[2];
sx q[2];
rz(-0.70435134) q[2];
sx q[2];
rz(2.4865502) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0807018) q[1];
sx q[1];
rz(-1.3403088) q[1];
sx q[1];
rz(-0.73841612) q[1];
x q[2];
rz(3.1327206) q[3];
sx q[3];
rz(-2.7089467) q[3];
sx q[3];
rz(3.0888626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7880154) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(-1.2865944) q[2];
rz(0.50518099) q[3];
sx q[3];
rz(-2.6954539) q[3];
sx q[3];
rz(-1.2222458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196359) q[0];
sx q[0];
rz(-3.0918047) q[0];
sx q[0];
rz(1.5420445) q[0];
rz(0.75625769) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(2.8047628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23035971) q[0];
sx q[0];
rz(-0.31765881) q[0];
sx q[0];
rz(-0.0068276366) q[0];
x q[1];
rz(-2.9354503) q[2];
sx q[2];
rz(-1.2954324) q[2];
sx q[2];
rz(2.5554267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5440798) q[1];
sx q[1];
rz(-1.656053) q[1];
sx q[1];
rz(1.4592749) q[1];
rz(1.0098861) q[3];
sx q[3];
rz(-1.2133382) q[3];
sx q[3];
rz(-1.7021029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6280262) q[2];
sx q[2];
rz(-0.94945532) q[2];
sx q[2];
rz(2.8619518) q[2];
rz(-0.63129342) q[3];
sx q[3];
rz(-0.93835962) q[3];
sx q[3];
rz(2.5173371) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414128) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(-2.367876) q[1];
sx q[1];
rz(-0.63540375) q[1];
sx q[1];
rz(0.2159963) q[1];
rz(0.86376376) q[2];
sx q[2];
rz(-2.2939199) q[2];
sx q[2];
rz(-1.6242956) q[2];
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
