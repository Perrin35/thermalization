OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6516946) q[0];
sx q[0];
rz(5.4013847) q[0];
sx q[0];
rz(9.2800979) q[0];
rz(0.41809234) q[1];
sx q[1];
rz(-0.10377181) q[1];
sx q[1];
rz(-1.6296847) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967752) q[0];
sx q[0];
rz(-0.51048764) q[0];
sx q[0];
rz(1.8148242) q[0];
rz(-pi) q[1];
rz(0.060188607) q[2];
sx q[2];
rz(-2.4071949) q[2];
sx q[2];
rz(-0.73645624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.712973) q[1];
sx q[1];
rz(-0.52134544) q[1];
sx q[1];
rz(-0.85835427) q[1];
rz(0.51201323) q[3];
sx q[3];
rz(-1.9422934) q[3];
sx q[3];
rz(-3.01513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31643733) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(2.4862508) q[2];
rz(1.7573028) q[3];
sx q[3];
rz(-2.1770848) q[3];
sx q[3];
rz(2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69219387) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(1.7610158) q[1];
sx q[1];
rz(-0.5813798) q[1];
sx q[1];
rz(1.9320206) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0469477) q[0];
sx q[0];
rz(-2.1535465) q[0];
sx q[0];
rz(2.9215982) q[0];
rz(-1.5530994) q[2];
sx q[2];
rz(-1.624384) q[2];
sx q[2];
rz(0.45547661) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3833419) q[1];
sx q[1];
rz(-1.1102601) q[1];
sx q[1];
rz(1.7136445) q[1];
x q[2];
rz(-0.38755667) q[3];
sx q[3];
rz(-2.1394552) q[3];
sx q[3];
rz(-2.89944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9380583) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(-2.5246942) q[2];
rz(1.6054224) q[3];
sx q[3];
rz(-1.0441531) q[3];
sx q[3];
rz(-0.49083403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817552) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(2.4165261) q[0];
rz(1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(2.1824172) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77391439) q[0];
sx q[0];
rz(-1.9579906) q[0];
sx q[0];
rz(-1.0605124) q[0];
x q[1];
rz(-2.3940635) q[2];
sx q[2];
rz(-2.5995214) q[2];
sx q[2];
rz(-2.4666748) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37466418) q[1];
sx q[1];
rz(-0.21682021) q[1];
sx q[1];
rz(-1.0105074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.972288) q[3];
sx q[3];
rz(-1.7506071) q[3];
sx q[3];
rz(-1.7468417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3265257) q[2];
sx q[2];
rz(-1.8748137) q[2];
sx q[2];
rz(2.8845924) q[2];
rz(1.0810931) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(-1.5628373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0839888) q[0];
sx q[0];
rz(-2.8631518) q[0];
sx q[0];
rz(1.9637015) q[0];
rz(-2.9914757) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(-1.6252801) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35203194) q[0];
sx q[0];
rz(-1.2403249) q[0];
sx q[0];
rz(2.8432557) q[0];
rz(-pi) q[1];
rz(-0.84211911) q[2];
sx q[2];
rz(-1.2884022) q[2];
sx q[2];
rz(2.6248464) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5387568) q[1];
sx q[1];
rz(-1.2456166) q[1];
sx q[1];
rz(1.4657337) q[1];
rz(0.64370207) q[3];
sx q[3];
rz(-1.4182404) q[3];
sx q[3];
rz(-1.8495454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84245044) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(2.2184856) q[2];
rz(-0.9592157) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8206896) q[0];
sx q[0];
rz(-0.90383363) q[0];
sx q[0];
rz(-0.88291105) q[0];
rz(1.8461022) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(2.6466218) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6225016) q[0];
sx q[0];
rz(-1.8880196) q[0];
sx q[0];
rz(-2.836801) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60762824) q[2];
sx q[2];
rz(-0.80127412) q[2];
sx q[2];
rz(-0.60333383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40303142) q[1];
sx q[1];
rz(-2.237169) q[1];
sx q[1];
rz(-1.4590461) q[1];
rz(1.048377) q[3];
sx q[3];
rz(-1.256878) q[3];
sx q[3];
rz(-0.81410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8098658) q[2];
sx q[2];
rz(-3.0220384) q[2];
sx q[2];
rz(-1.8734107) q[2];
rz(-0.082503334) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(-0.65703195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4853915) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(1.50151) q[0];
rz(-1.062475) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(-1.1753488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1262681) q[0];
sx q[0];
rz(-1.7909174) q[0];
sx q[0];
rz(-2.8298529) q[0];
rz(-pi) q[1];
rz(-0.69357102) q[2];
sx q[2];
rz(-1.7993002) q[2];
sx q[2];
rz(-0.19738838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3835241) q[1];
sx q[1];
rz(-1.3192156) q[1];
sx q[1];
rz(-1.6918159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31194056) q[3];
sx q[3];
rz(-2.3565203) q[3];
sx q[3];
rz(1.5573481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33092734) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(-2.1539099) q[2];
rz(0.87743131) q[3];
sx q[3];
rz(-0.26724795) q[3];
sx q[3];
rz(-0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-1.1838609) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(3.0027332) q[0];
rz(-1.9271556) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(-0.81659281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3432085) q[0];
sx q[0];
rz(-1.7955762) q[0];
sx q[0];
rz(-3.1220469) q[0];
rz(-pi) q[1];
rz(-1.8310407) q[2];
sx q[2];
rz(-0.81965551) q[2];
sx q[2];
rz(1.2092839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.066799435) q[1];
sx q[1];
rz(-1.6082676) q[1];
sx q[1];
rz(-0.85185677) q[1];
rz(0.84888938) q[3];
sx q[3];
rz(-2.1037641) q[3];
sx q[3];
rz(-0.70895586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3198118) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(0.4846586) q[2];
rz(-2.7759077) q[3];
sx q[3];
rz(-1.3644783) q[3];
sx q[3];
rz(1.9827838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0115688) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(3.044361) q[0];
rz(-3.0367127) q[1];
sx q[1];
rz(-2.361894) q[1];
sx q[1];
rz(-1.8269151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7044367) q[0];
sx q[0];
rz(-1.5644685) q[0];
sx q[0];
rz(3.1392211) q[0];
rz(-1.8123367) q[2];
sx q[2];
rz(-1.353423) q[2];
sx q[2];
rz(-1.5474943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7200249) q[1];
sx q[1];
rz(-1.0526669) q[1];
sx q[1];
rz(-1.6314477) q[1];
rz(-pi) q[2];
rz(1.0339853) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(-1.0190462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9647727) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(2.9883265) q[2];
rz(1.9249453) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(0.6161859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1185054) q[0];
sx q[0];
rz(-0.0054792976) q[0];
sx q[0];
rz(1.6368921) q[0];
rz(2.3432689) q[1];
sx q[1];
rz(-1.5232122) q[1];
sx q[1];
rz(-0.33531478) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.714761) q[0];
sx q[0];
rz(-0.26894595) q[0];
sx q[0];
rz(1.8105276) q[0];
rz(-pi) q[1];
rz(-0.026235418) q[2];
sx q[2];
rz(-2.3910284) q[2];
sx q[2];
rz(0.98294965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7882725) q[1];
sx q[1];
rz(-1.2306552) q[1];
sx q[1];
rz(-1.2659094) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3319546) q[3];
sx q[3];
rz(-2.5510029) q[3];
sx q[3];
rz(-2.1826377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12751427) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(-2.7039995) q[2];
rz(-2.0994999) q[3];
sx q[3];
rz(-1.4053248) q[3];
sx q[3];
rz(2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9893148) q[0];
sx q[0];
rz(-2.2570026) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(0.97688976) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(0.50122112) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64912632) q[0];
sx q[0];
rz(-1.3808131) q[0];
sx q[0];
rz(-1.8590007) q[0];
x q[1];
rz(1.9947204) q[2];
sx q[2];
rz(-1.3739283) q[2];
sx q[2];
rz(-0.19320657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2071562) q[1];
sx q[1];
rz(-2.0324576) q[1];
sx q[1];
rz(-2.8769879) q[1];
rz(2.1887423) q[3];
sx q[3];
rz(-2.5898159) q[3];
sx q[3];
rz(-2.1976041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1180798) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(-0.21772131) q[2];
rz(1.2348385) q[3];
sx q[3];
rz(-0.66399884) q[3];
sx q[3];
rz(-2.2209404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49391838) q[0];
sx q[0];
rz(-1.4287345) q[0];
sx q[0];
rz(-2.7609974) q[0];
rz(0.71161288) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(0.099945036) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(-2.7357581) q[3];
sx q[3];
rz(-1.3530227) q[3];
sx q[3];
rz(1.4128662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
