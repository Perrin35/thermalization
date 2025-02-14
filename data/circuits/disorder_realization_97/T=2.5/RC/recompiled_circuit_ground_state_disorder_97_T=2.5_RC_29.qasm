OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41807732) q[0];
sx q[0];
rz(3.1273841) q[0];
sx q[0];
rz(9.4655884) q[0];
rz(-0.040925097) q[1];
sx q[1];
rz(2.0857781) q[1];
sx q[1];
rz(9.981282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9153298) q[0];
sx q[0];
rz(-1.2253083) q[0];
sx q[0];
rz(2.9468263) q[0];
rz(-pi) q[1];
rz(2.7950613) q[2];
sx q[2];
rz(-2.541447) q[2];
sx q[2];
rz(-0.66020667) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79136363) q[1];
sx q[1];
rz(-1.204317) q[1];
sx q[1];
rz(0.54270737) q[1];
rz(1.3338148) q[3];
sx q[3];
rz(-1.4795884) q[3];
sx q[3];
rz(2.3438675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1540264) q[2];
sx q[2];
rz(-1.7511) q[2];
sx q[2];
rz(1.7517368) q[2];
rz(0.050358199) q[3];
sx q[3];
rz(-1.4202838) q[3];
sx q[3];
rz(-1.095298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20767009) q[0];
sx q[0];
rz(-1.4140797) q[0];
sx q[0];
rz(3.0811144) q[0];
rz(2.1531847) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(-0.21313937) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3857065) q[0];
sx q[0];
rz(-1.1637628) q[0];
sx q[0];
rz(-0.72188053) q[0];
rz(-1.7333116) q[2];
sx q[2];
rz(-1.6846416) q[2];
sx q[2];
rz(-1.3448037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2717183) q[1];
sx q[1];
rz(-1.3187746) q[1];
sx q[1];
rz(-1.6099754) q[1];
rz(2.0673236) q[3];
sx q[3];
rz(-1.7317711) q[3];
sx q[3];
rz(2.9825008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1468143) q[2];
sx q[2];
rz(-1.549336) q[2];
sx q[2];
rz(0.0063304338) q[2];
rz(-1.7019042) q[3];
sx q[3];
rz(-0.78841811) q[3];
sx q[3];
rz(-0.44062135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45563662) q[0];
sx q[0];
rz(-0.39304471) q[0];
sx q[0];
rz(-1.3156369) q[0];
rz(1.3872967) q[1];
sx q[1];
rz(-2.5376814) q[1];
sx q[1];
rz(3.0953298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14213698) q[0];
sx q[0];
rz(-1.6007413) q[0];
sx q[0];
rz(1.5718476) q[0];
rz(-2.4607804) q[2];
sx q[2];
rz(-0.57475427) q[2];
sx q[2];
rz(2.0865692) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1686778) q[1];
sx q[1];
rz(-2.4371304) q[1];
sx q[1];
rz(-2.5903914) q[1];
rz(-pi) q[2];
rz(-0.95607676) q[3];
sx q[3];
rz(-1.4341364) q[3];
sx q[3];
rz(-2.7352493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5578654) q[2];
sx q[2];
rz(-1.9670308) q[2];
sx q[2];
rz(-2.2306856) q[2];
rz(1.9117428) q[3];
sx q[3];
rz(-2.3807447) q[3];
sx q[3];
rz(-0.41874829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8987027) q[0];
sx q[0];
rz(-1.3985343) q[0];
sx q[0];
rz(-2.718495) q[0];
rz(-0.93370071) q[1];
sx q[1];
rz(-1.2010778) q[1];
sx q[1];
rz(2.7074714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2796626) q[0];
sx q[0];
rz(-1.8619143) q[0];
sx q[0];
rz(1.5345957) q[0];
rz(-pi) q[1];
rz(1.7291315) q[2];
sx q[2];
rz(-2.2248817) q[2];
sx q[2];
rz(-2.9323062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0435183) q[1];
sx q[1];
rz(-1.6988108) q[1];
sx q[1];
rz(0.56358595) q[1];
x q[2];
rz(-1.0634123) q[3];
sx q[3];
rz(-0.23325167) q[3];
sx q[3];
rz(0.95913974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9551129) q[2];
sx q[2];
rz(-3.0035512) q[2];
sx q[2];
rz(1.0008) q[2];
rz(-0.75438386) q[3];
sx q[3];
rz(-0.84027925) q[3];
sx q[3];
rz(-1.8050571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56483018) q[0];
sx q[0];
rz(-2.5046528) q[0];
sx q[0];
rz(-1.379396) q[0];
rz(2.5881536) q[1];
sx q[1];
rz(-1.484) q[1];
sx q[1];
rz(2.4310506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0481537) q[0];
sx q[0];
rz(-1.7293674) q[0];
sx q[0];
rz(-2.3331654) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24181409) q[2];
sx q[2];
rz(-1.6266357) q[2];
sx q[2];
rz(-2.9602891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.50113) q[1];
sx q[1];
rz(-1.8014596) q[1];
sx q[1];
rz(3.1063558) q[1];
x q[2];
rz(0.54039012) q[3];
sx q[3];
rz(-1.4266664) q[3];
sx q[3];
rz(1.543951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1826627) q[2];
sx q[2];
rz(-1.6561597) q[2];
sx q[2];
rz(0.043969285) q[2];
rz(2.8751657) q[3];
sx q[3];
rz(-0.1259585) q[3];
sx q[3];
rz(-2.9282279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63785315) q[0];
sx q[0];
rz(-1.6197236) q[0];
sx q[0];
rz(2.7730686) q[0];
rz(2.7875994) q[1];
sx q[1];
rz(-2.0123672) q[1];
sx q[1];
rz(1.0634208) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76621687) q[0];
sx q[0];
rz(-1.5581344) q[0];
sx q[0];
rz(-1.4305737) q[0];
rz(-pi) q[1];
rz(-1.9475458) q[2];
sx q[2];
rz(-2.8440335) q[2];
sx q[2];
rz(-2.3524795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8700347) q[1];
sx q[1];
rz(-1.9816625) q[1];
sx q[1];
rz(2.690587) q[1];
x q[2];
rz(0.20078378) q[3];
sx q[3];
rz(-0.75888435) q[3];
sx q[3];
rz(-0.52606612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4949558) q[2];
sx q[2];
rz(-1.2594624) q[2];
sx q[2];
rz(1.4259526) q[2];
rz(-0.87092733) q[3];
sx q[3];
rz(-2.0723074) q[3];
sx q[3];
rz(-0.23437414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17212412) q[0];
sx q[0];
rz(-2.5775462) q[0];
sx q[0];
rz(-0.010757541) q[0];
rz(2.5116008) q[1];
sx q[1];
rz(-1.0456345) q[1];
sx q[1];
rz(-0.90525544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7566959) q[0];
sx q[0];
rz(-1.6357452) q[0];
sx q[0];
rz(2.1048291) q[0];
rz(-pi) q[1];
rz(2.808242) q[2];
sx q[2];
rz(-1.0819737) q[2];
sx q[2];
rz(0.49788633) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7510482) q[1];
sx q[1];
rz(-2.6169951) q[1];
sx q[1];
rz(-0.28769298) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7131374) q[3];
sx q[3];
rz(-0.96142229) q[3];
sx q[3];
rz(0.46655015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4524727) q[2];
sx q[2];
rz(-1.2484739) q[2];
sx q[2];
rz(-3.0242237) q[2];
rz(0.72702414) q[3];
sx q[3];
rz(-0.89594642) q[3];
sx q[3];
rz(1.1400384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3177976) q[0];
sx q[0];
rz(-1.057484) q[0];
sx q[0];
rz(-0.29528883) q[0];
rz(-1.3914039) q[1];
sx q[1];
rz(-2.4877805) q[1];
sx q[1];
rz(-2.6546435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21509296) q[0];
sx q[0];
rz(-1.3736667) q[0];
sx q[0];
rz(-1.6736945) q[0];
x q[1];
rz(-1.0560965) q[2];
sx q[2];
rz(-1.8499415) q[2];
sx q[2];
rz(2.5893958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0978969) q[1];
sx q[1];
rz(-2.0123031) q[1];
sx q[1];
rz(-0.82951464) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5540274) q[3];
sx q[3];
rz(-0.31767818) q[3];
sx q[3];
rz(-0.4099572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9140909) q[2];
sx q[2];
rz(-0.5946244) q[2];
sx q[2];
rz(-1.5458858) q[2];
rz(-0.60901904) q[3];
sx q[3];
rz(-1.1142497) q[3];
sx q[3];
rz(-0.6685037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85536999) q[0];
sx q[0];
rz(-0.90317059) q[0];
sx q[0];
rz(1.5089996) q[0];
rz(1.1434932) q[1];
sx q[1];
rz(-1.9805464) q[1];
sx q[1];
rz(0.11018363) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4297053) q[0];
sx q[0];
rz(-1.674187) q[0];
sx q[0];
rz(0.59075077) q[0];
x q[1];
rz(-1.1082591) q[2];
sx q[2];
rz(-1.8018844) q[2];
sx q[2];
rz(-3.0840735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5563123) q[1];
sx q[1];
rz(-1.9686799) q[1];
sx q[1];
rz(-2.3817987) q[1];
rz(-pi) q[2];
rz(-1.5694782) q[3];
sx q[3];
rz(-2.0681046) q[3];
sx q[3];
rz(0.2494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1630359) q[2];
sx q[2];
rz(-1.0806934) q[2];
sx q[2];
rz(2.2829096) q[2];
rz(-1.5052172) q[3];
sx q[3];
rz(-1.7965763) q[3];
sx q[3];
rz(-1.6987364) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34864146) q[0];
sx q[0];
rz(-2.0703147) q[0];
sx q[0];
rz(2.4559378) q[0];
rz(0.81949702) q[1];
sx q[1];
rz(-1.6623431) q[1];
sx q[1];
rz(2.562838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98812719) q[0];
sx q[0];
rz(-0.99737001) q[0];
sx q[0];
rz(0.98116409) q[0];
rz(-pi) q[1];
rz(-0.24722331) q[2];
sx q[2];
rz(-1.0178538) q[2];
sx q[2];
rz(-1.2485112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6616649) q[1];
sx q[1];
rz(-2.6746422) q[1];
sx q[1];
rz(-0.60162094) q[1];
rz(3.0323735) q[3];
sx q[3];
rz(-1.6192604) q[3];
sx q[3];
rz(-2.8630961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3026352) q[2];
sx q[2];
rz(-2.9557648) q[2];
sx q[2];
rz(0.1798943) q[2];
rz(-0.57080615) q[3];
sx q[3];
rz(-2.4195318) q[3];
sx q[3];
rz(2.4898237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30052341) q[0];
sx q[0];
rz(-2.4916334) q[0];
sx q[0];
rz(-1.4148225) q[0];
rz(-2.6534446) q[1];
sx q[1];
rz(-2.585325) q[1];
sx q[1];
rz(-1.5811031) q[1];
rz(0.93964259) q[2];
sx q[2];
rz(-1.2705391) q[2];
sx q[2];
rz(1.4654797) q[2];
rz(-0.45541422) q[3];
sx q[3];
rz(-0.58697636) q[3];
sx q[3];
rz(-0.1252115) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
