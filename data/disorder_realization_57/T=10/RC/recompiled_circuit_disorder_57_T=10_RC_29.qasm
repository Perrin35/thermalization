OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(-2.3318113) q[0];
sx q[0];
rz(0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(1.3637435) q[1];
sx q[1];
rz(10.627828) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8578313) q[0];
sx q[0];
rz(-0.28781048) q[0];
sx q[0];
rz(-0.83588155) q[0];
rz(-1.3538023) q[2];
sx q[2];
rz(-1.3218244) q[2];
sx q[2];
rz(-0.72460246) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0359218) q[1];
sx q[1];
rz(-2.8997313) q[1];
sx q[1];
rz(-0.24005228) q[1];
rz(-2.8486512) q[3];
sx q[3];
rz(-1.5353234) q[3];
sx q[3];
rz(-2.5778511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8649341) q[2];
sx q[2];
rz(-2.0099137) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(1.6137326) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(0.73408192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.99877015) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(-0.2858513) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430543) q[0];
sx q[0];
rz(-2.1715954) q[0];
sx q[0];
rz(1.1061125) q[0];
rz(-pi) q[1];
rz(2.9026047) q[2];
sx q[2];
rz(-1.1766489) q[2];
sx q[2];
rz(-0.1220526) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8394168) q[1];
sx q[1];
rz(-1.2100879) q[1];
sx q[1];
rz(-1.9990666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.500202) q[3];
sx q[3];
rz(-2.4603599) q[3];
sx q[3];
rz(2.2487946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9849898) q[2];
sx q[2];
rz(-1.9282324) q[2];
sx q[2];
rz(1.9796237) q[2];
rz(-0.087163838) q[3];
sx q[3];
rz(-1.6379084) q[3];
sx q[3];
rz(0.073908977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(-2.4940925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.203813) q[0];
sx q[0];
rz(-0.63596361) q[0];
sx q[0];
rz(-0.075809191) q[0];
rz(-pi) q[1];
rz(2.5325003) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(-0.44037214) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4574617) q[1];
sx q[1];
rz(-1.1161472) q[1];
sx q[1];
rz(2.734114) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87326761) q[3];
sx q[3];
rz(-1.0928565) q[3];
sx q[3];
rz(2.0877116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8422164) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(-1.5807318) q[2];
rz(2.1598699) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(2.2696944) q[0];
rz(-0.38726989) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(0.29104582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8431393) q[0];
sx q[0];
rz(-1.7994587) q[0];
sx q[0];
rz(0.97897162) q[0];
rz(-pi) q[1];
rz(-1.8918858) q[2];
sx q[2];
rz(-1.2920657) q[2];
sx q[2];
rz(2.2433787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2070771) q[1];
sx q[1];
rz(-2.0743437) q[1];
sx q[1];
rz(2.3758923) q[1];
rz(1.8157186) q[3];
sx q[3];
rz(-2.5513253) q[3];
sx q[3];
rz(3.0274689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3691833) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(-0.54405653) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13957025) q[0];
sx q[0];
rz(-2.1152571) q[0];
sx q[0];
rz(0.64506662) q[0];
rz(2.5158665) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7932927) q[0];
sx q[0];
rz(-0.59069809) q[0];
sx q[0];
rz(2.9212055) q[0];
rz(2.7475287) q[2];
sx q[2];
rz(-0.85959083) q[2];
sx q[2];
rz(1.994339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2423357) q[1];
sx q[1];
rz(-1.3471654) q[1];
sx q[1];
rz(2.1171338) q[1];
rz(-pi) q[2];
rz(0.55769063) q[3];
sx q[3];
rz(-1.4413177) q[3];
sx q[3];
rz(-1.492471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3877635) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(1.3641664) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(-0.92145872) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0867778) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(2.4247647) q[0];
rz(-2.7038799) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(0.81370083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1610634) q[0];
sx q[0];
rz(-1.5080394) q[0];
sx q[0];
rz(1.784523) q[0];
x q[1];
rz(-0.72197638) q[2];
sx q[2];
rz(-1.1140031) q[2];
sx q[2];
rz(1.7320088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8061714) q[1];
sx q[1];
rz(-2.2911982) q[1];
sx q[1];
rz(0.45067388) q[1];
rz(-pi) q[2];
rz(-1.2171451) q[3];
sx q[3];
rz(-1.8600703) q[3];
sx q[3];
rz(-1.7464964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8852691) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880144) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(0.6860835) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(0.66326052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8535186) q[0];
sx q[0];
rz(-0.92068866) q[0];
sx q[0];
rz(-3.1115565) q[0];
rz(-pi) q[1];
rz(0.29055039) q[2];
sx q[2];
rz(-0.65813375) q[2];
sx q[2];
rz(-1.5302637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7546332) q[1];
sx q[1];
rz(-2.6135751) q[1];
sx q[1];
rz(-1.9792884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0955986) q[3];
sx q[3];
rz(-2.0810658) q[3];
sx q[3];
rz(1.0218395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-2.3040999) q[2];
sx q[2];
rz(-3.1265124) q[2];
rz(-2.1298501) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052977234) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(0.10738871) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(0.4577786) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7664238) q[0];
sx q[0];
rz(-2.1918115) q[0];
sx q[0];
rz(-0.39874886) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56577487) q[2];
sx q[2];
rz(-0.83116311) q[2];
sx q[2];
rz(-2.4772252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6785504) q[1];
sx q[1];
rz(-1.2150587) q[1];
sx q[1];
rz(-2.8929391) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5910208) q[3];
sx q[3];
rz(-0.62567657) q[3];
sx q[3];
rz(1.613137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7765939) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(1.4252023) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.8553998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872221) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(-2.1881058) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(-2.5568533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027771587) q[0];
sx q[0];
rz(-0.50243938) q[0];
sx q[0];
rz(-0.2343982) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.598001) q[2];
sx q[2];
rz(-1.2541176) q[2];
sx q[2];
rz(-1.8900332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.016420267) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(0.16146407) q[1];
rz(-pi) q[2];
rz(-1.8381848) q[3];
sx q[3];
rz(-1.4328453) q[3];
sx q[3];
rz(-1.4064058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(0.27627036) q[2];
rz(-2.590498) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(2.7579257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0693414) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(-1.1570702) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.6190593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009506) q[0];
sx q[0];
rz(-0.92693936) q[0];
sx q[0];
rz(-1.180091) q[0];
rz(-pi) q[1];
rz(-2.7340639) q[2];
sx q[2];
rz(-0.73020836) q[2];
sx q[2];
rz(-0.66116316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5628221) q[1];
sx q[1];
rz(-0.65014231) q[1];
sx q[1];
rz(-0.27362089) q[1];
rz(-pi) q[2];
rz(0.75928648) q[3];
sx q[3];
rz(-2.5873103) q[3];
sx q[3];
rz(-0.62461531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9860501) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(0.84214169) q[2];
rz(-1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-0.026528927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(-1.7799829) q[2];
sx q[2];
rz(-1.5447415) q[2];
sx q[2];
rz(0.011199657) q[2];
rz(-0.079506569) q[3];
sx q[3];
rz(-2.2600997) q[3];
sx q[3];
rz(-3.1117677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
