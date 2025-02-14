OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4898981) q[0];
sx q[0];
rz(-2.2597921) q[0];
sx q[0];
rz(-2.9969126) q[0];
rz(3.559685) q[1];
sx q[1];
rz(3.2453645) q[1];
sx q[1];
rz(11.054463) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4747595) q[0];
sx q[0];
rz(-2.0647683) q[0];
sx q[0];
rz(-3.0071063) q[0];
rz(-0.060188607) q[2];
sx q[2];
rz(-0.73439774) q[2];
sx q[2];
rz(2.4051364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6408893) q[1];
sx q[1];
rz(-1.2391866) q[1];
sx q[1];
rz(-1.9808116) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51201323) q[3];
sx q[3];
rz(-1.9422934) q[3];
sx q[3];
rz(3.01513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31643733) q[2];
sx q[2];
rz(-1.2707767) q[2];
sx q[2];
rz(0.6553418) q[2];
rz(-1.7573028) q[3];
sx q[3];
rz(-0.96450788) q[3];
sx q[3];
rz(2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69219387) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(0.50814116) q[0];
rz(1.3805768) q[1];
sx q[1];
rz(-0.5813798) q[1];
sx q[1];
rz(1.2095721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.094645) q[0];
sx q[0];
rz(-2.1535465) q[0];
sx q[0];
rz(2.9215982) q[0];
x q[1];
rz(-1.5530994) q[2];
sx q[2];
rz(-1.5172086) q[2];
sx q[2];
rz(2.686116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75825071) q[1];
sx q[1];
rz(-2.0313325) q[1];
sx q[1];
rz(1.7136445) q[1];
x q[2];
rz(-0.96662917) q[3];
sx q[3];
rz(-1.2467017) q[3];
sx q[3];
rz(-1.5965727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2035344) q[2];
sx q[2];
rz(-1.396023) q[2];
sx q[2];
rz(2.5246942) q[2];
rz(1.5361702) q[3];
sx q[3];
rz(-1.0441531) q[3];
sx q[3];
rz(-2.6507586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817552) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(-2.4165261) q[0];
rz(-1.2695351) q[1];
sx q[1];
rz(-2.4639362) q[1];
sx q[1];
rz(0.95917541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1364258) q[0];
sx q[0];
rz(-2.040103) q[0];
sx q[0];
rz(-0.43715663) q[0];
x q[1];
rz(-0.41590642) q[2];
sx q[2];
rz(-1.9291483) q[2];
sx q[2];
rz(-2.9168473) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.945248) q[1];
sx q[1];
rz(-1.3875393) q[1];
sx q[1];
rz(-0.11653479) q[1];
x q[2];
rz(-0.82328513) q[3];
sx q[3];
rz(-0.24634493) q[3];
sx q[3];
rz(-0.63194094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3265257) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(0.25700021) q[2];
rz(-1.0810931) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(1.5628373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0576039) q[0];
sx q[0];
rz(-2.8631518) q[0];
sx q[0];
rz(1.1778911) q[0];
rz(-2.9914757) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(-1.6252801) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1105791) q[0];
sx q[0];
rz(-0.4415126) q[0];
sx q[0];
rz(-2.2792321) q[0];
x q[1];
rz(0.84211911) q[2];
sx q[2];
rz(-1.8531905) q[2];
sx q[2];
rz(-0.51674622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8575688) q[1];
sx q[1];
rz(-2.800436) q[1];
sx q[1];
rz(0.30156231) q[1];
x q[2];
rz(0.64370207) q[3];
sx q[3];
rz(-1.7233522) q[3];
sx q[3];
rz(1.8495454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84245044) q[2];
sx q[2];
rz(-0.56988847) q[2];
sx q[2];
rz(0.92310706) q[2];
rz(0.9592157) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(-0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32090309) q[0];
sx q[0];
rz(-2.237759) q[0];
sx q[0];
rz(-2.2586816) q[0];
rz(-1.2954905) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(2.6466218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1495234) q[0];
sx q[0];
rz(-1.2816634) q[0];
sx q[0];
rz(1.2393214) q[0];
x q[1];
rz(-2.4385543) q[2];
sx q[2];
rz(-1.1482757) q[2];
sx q[2];
rz(-2.6247744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58261694) q[1];
sx q[1];
rz(-0.6742653) q[1];
sx q[1];
rz(0.14087454) q[1];
rz(-pi) q[2];
rz(2.1476188) q[3];
sx q[3];
rz(-0.60189542) q[3];
sx q[3];
rz(2.8771597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8098658) q[2];
sx q[2];
rz(-3.0220384) q[2];
sx q[2];
rz(1.8734107) q[2];
rz(-3.0590893) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(-2.4845607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6562011) q[0];
sx q[0];
rz(-0.20683658) q[0];
sx q[0];
rz(-1.50151) q[0];
rz(1.062475) q[1];
sx q[1];
rz(-1.9891918) q[1];
sx q[1];
rz(1.9662439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9906015) q[0];
sx q[0];
rz(-0.37953934) q[0];
sx q[0];
rz(0.63024704) q[0];
rz(-pi) q[1];
rz(-1.2771036) q[2];
sx q[2];
rz(-2.2429464) q[2];
sx q[2];
rz(1.1872684) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8379587) q[1];
sx q[1];
rz(-2.8629747) q[1];
sx q[1];
rz(2.7024804) q[1];
rz(-pi) q[2];
rz(-2.8296521) q[3];
sx q[3];
rz(-2.3565203) q[3];
sx q[3];
rz(1.5842445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8106653) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(-2.1539099) q[2];
rz(2.2641613) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(2.4115244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9577318) q[0];
sx q[0];
rz(-1.5393625) q[0];
sx q[0];
rz(3.0027332) q[0];
rz(1.9271556) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(0.81659281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3735376) q[0];
sx q[0];
rz(-1.5517424) q[0];
sx q[0];
rz(-1.7956177) q[0];
rz(-pi) q[1];
rz(-0.76824378) q[2];
sx q[2];
rz(-1.7599987) q[2];
sx q[2];
rz(-0.18174325) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6703709) q[1];
sx q[1];
rz(-0.85247096) q[1];
sx q[1];
rz(3.0918151) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2927033) q[3];
sx q[3];
rz(-1.0378285) q[3];
sx q[3];
rz(-0.70895586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3198118) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(-0.4846586) q[2];
rz(-2.7759077) q[3];
sx q[3];
rz(-1.3644783) q[3];
sx q[3];
rz(1.9827838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(-0.77969867) q[1];
sx q[1];
rz(1.8269151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7044367) q[0];
sx q[0];
rz(-1.5771241) q[0];
sx q[0];
rz(3.1392211) q[0];
rz(-pi) q[1];
rz(-2.9179321) q[2];
sx q[2];
rz(-1.8065435) q[2];
sx q[2];
rz(-0.076381269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7200249) q[1];
sx q[1];
rz(-2.0889258) q[1];
sx q[1];
rz(1.6314477) q[1];
x q[2];
rz(2.1076074) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(1.0190462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9647727) q[2];
sx q[2];
rz(-1.3975881) q[2];
sx q[2];
rz(-0.1532661) q[2];
rz(-1.9249453) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(2.5254068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0230873) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(-1.5047005) q[0];
rz(2.3432689) q[1];
sx q[1];
rz(-1.5232122) q[1];
sx q[1];
rz(2.8062779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42683168) q[0];
sx q[0];
rz(-2.8726467) q[0];
sx q[0];
rz(1.3310651) q[0];
x q[1];
rz(0.7503926) q[2];
sx q[2];
rz(-1.5886891) q[2];
sx q[2];
rz(0.6070348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.032051) q[1];
sx q[1];
rz(-2.6888043) q[1];
sx q[1];
rz(-2.4381263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3319546) q[3];
sx q[3];
rz(-0.59058978) q[3];
sx q[3];
rz(2.1826377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12751427) q[2];
sx q[2];
rz(-1.1684343) q[2];
sx q[2];
rz(2.7039995) q[2];
rz(-1.0420927) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(-0.54245814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9893148) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(2.1647029) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(2.6403715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2758482) q[0];
sx q[0];
rz(-1.2879218) q[0];
sx q[0];
rz(-0.19794835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1468723) q[2];
sx q[2];
rz(-1.7676643) q[2];
sx q[2];
rz(0.19320657) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3851191) q[1];
sx q[1];
rz(-1.8071399) q[1];
sx q[1];
rz(-2.0467351) q[1];
rz(-2.7990255) q[3];
sx q[3];
rz(-2.0122572) q[3];
sx q[3];
rz(-1.6395237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1180798) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(-0.21772131) q[2];
rz(-1.2348385) q[3];
sx q[3];
rz(-0.66399884) q[3];
sx q[3];
rz(2.2209404) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476743) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(-2.4299798) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(-0.099945036) q[2];
sx q[2];
rz(-0.47158416) q[2];
sx q[2];
rz(-0.33000962) q[2];
rz(1.3344516) q[3];
sx q[3];
rz(-1.1750887) q[3];
sx q[3];
rz(-0.25050161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
