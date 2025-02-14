OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(-1.886263) q[0];
sx q[0];
rz(2.6502967) q[0];
rz(2.8331941) q[1];
sx q[1];
rz(-1.3003132) q[1];
sx q[1];
rz(0.058606776) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8347625) q[0];
sx q[0];
rz(-1.0245203) q[0];
sx q[0];
rz(0.52301126) q[0];
rz(-1.116007) q[2];
sx q[2];
rz(-2.3014549) q[2];
sx q[2];
rz(-2.0065713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72926988) q[1];
sx q[1];
rz(-2.7979921) q[1];
sx q[1];
rz(-2.6471828) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0618151) q[3];
sx q[3];
rz(-1.8892131) q[3];
sx q[3];
rz(1.7169881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10597524) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(0.38113511) q[2];
rz(3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(-2.0924675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1586665) q[0];
sx q[0];
rz(-0.52678147) q[0];
sx q[0];
rz(2.7016933) q[0];
rz(-3.0797709) q[1];
sx q[1];
rz(-1.5301751) q[1];
sx q[1];
rz(-0.27438146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44654057) q[0];
sx q[0];
rz(-0.98869222) q[0];
sx q[0];
rz(-2.5869408) q[0];
rz(-pi) q[1];
rz(1.0978866) q[2];
sx q[2];
rz(-0.47240546) q[2];
sx q[2];
rz(0.96443572) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.741863) q[1];
sx q[1];
rz(-1.4951839) q[1];
sx q[1];
rz(-2.0194598) q[1];
x q[2];
rz(0.54631375) q[3];
sx q[3];
rz(-1.7726745) q[3];
sx q[3];
rz(-3.0936732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2939833) q[2];
sx q[2];
rz(-0.01394883) q[2];
sx q[2];
rz(3.0713165) q[2];
rz(-2.516732) q[3];
sx q[3];
rz(-1.3288386) q[3];
sx q[3];
rz(0.074987324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845508) q[0];
sx q[0];
rz(-1.5032737) q[0];
sx q[0];
rz(2.7497838) q[0];
rz(-2.1458972) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(-3.0973184) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366695) q[0];
sx q[0];
rz(-1.2630487) q[0];
sx q[0];
rz(2.1627764) q[0];
x q[1];
rz(2.2703553) q[2];
sx q[2];
rz(-2.2164564) q[2];
sx q[2];
rz(-1.7019367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6926319) q[1];
sx q[1];
rz(-1.3621743) q[1];
sx q[1];
rz(1.7600928) q[1];
rz(-pi) q[2];
rz(-1.6674394) q[3];
sx q[3];
rz(-2.1657888) q[3];
sx q[3];
rz(-2.4293824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6280881) q[2];
sx q[2];
rz(-2.2254483) q[2];
sx q[2];
rz(-1.2543031) q[2];
rz(-3.0301376) q[3];
sx q[3];
rz(-0.97193757) q[3];
sx q[3];
rz(2.286262) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3067538) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(0.87598959) q[0];
rz(0.39729473) q[1];
sx q[1];
rz(-1.5700424) q[1];
sx q[1];
rz(-1.809459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080058424) q[0];
sx q[0];
rz(-1.3803015) q[0];
sx q[0];
rz(-0.47866042) q[0];
rz(-2.6880363) q[2];
sx q[2];
rz(-2.1906914) q[2];
sx q[2];
rz(0.93728055) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6072453) q[1];
sx q[1];
rz(-1.6493454) q[1];
sx q[1];
rz(0.15395366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2355742) q[3];
sx q[3];
rz(-2.0959508) q[3];
sx q[3];
rz(-0.88334879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1490271) q[2];
sx q[2];
rz(-1.488648) q[2];
sx q[2];
rz(1.3943025) q[2];
rz(2.1483138) q[3];
sx q[3];
rz(-1.1634049) q[3];
sx q[3];
rz(1.8786028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(-2.2827523) q[0];
rz(-0.42293388) q[1];
sx q[1];
rz(-0.49086389) q[1];
sx q[1];
rz(-1.4609969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92991023) q[0];
sx q[0];
rz(-2.0693464) q[0];
sx q[0];
rz(-1.7105667) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.00255835) q[2];
sx q[2];
rz(-1.125549) q[2];
sx q[2];
rz(1.3819577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8940058) q[1];
sx q[1];
rz(-1.7628221) q[1];
sx q[1];
rz(-0.2094261) q[1];
x q[2];
rz(-1.8247174) q[3];
sx q[3];
rz(-2.6220206) q[3];
sx q[3];
rz(-2.7412753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6964263) q[2];
sx q[2];
rz(-2.7327765) q[2];
sx q[2];
rz(1.4985296) q[2];
rz(-2.0131352) q[3];
sx q[3];
rz(-2.0339637) q[3];
sx q[3];
rz(0.64277738) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6978825) q[0];
sx q[0];
rz(-1.3297798) q[0];
sx q[0];
rz(2.4238996) q[0];
rz(2.7806661) q[1];
sx q[1];
rz(-0.91054994) q[1];
sx q[1];
rz(0.31401971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7166724) q[0];
sx q[0];
rz(-1.0756399) q[0];
sx q[0];
rz(-1.6177482) q[0];
x q[1];
rz(1.5271565) q[2];
sx q[2];
rz(-0.72481643) q[2];
sx q[2];
rz(-2.5284655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9523588) q[1];
sx q[1];
rz(-1.2554607) q[1];
sx q[1];
rz(-1.0811514) q[1];
rz(-0.76992294) q[3];
sx q[3];
rz(-1.439538) q[3];
sx q[3];
rz(0.24410393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6641984) q[2];
sx q[2];
rz(-1.3496642) q[2];
sx q[2];
rz(-2.9273709) q[2];
rz(-2.2231806) q[3];
sx q[3];
rz(-0.94614202) q[3];
sx q[3];
rz(-1.4001747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65900954) q[0];
sx q[0];
rz(-0.56203401) q[0];
sx q[0];
rz(2.5157978) q[0];
rz(-1.91045) q[1];
sx q[1];
rz(-1.2673255) q[1];
sx q[1];
rz(-2.321718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5441262) q[0];
sx q[0];
rz(-1.0080999) q[0];
sx q[0];
rz(-1.94348) q[0];
rz(2.6447222) q[2];
sx q[2];
rz(-2.3665603) q[2];
sx q[2];
rz(-2.1642223) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.956683) q[1];
sx q[1];
rz(-1.9016097) q[1];
sx q[1];
rz(0.079903101) q[1];
x q[2];
rz(-0.094535703) q[3];
sx q[3];
rz(-1.5526287) q[3];
sx q[3];
rz(-1.2430735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95255533) q[2];
sx q[2];
rz(-1.8841212) q[2];
sx q[2];
rz(0.38743585) q[2];
rz(0.45346013) q[3];
sx q[3];
rz(-2.267434) q[3];
sx q[3];
rz(-2.9816154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74359918) q[0];
sx q[0];
rz(-1.1389808) q[0];
sx q[0];
rz(2.0177662) q[0];
rz(0.75346142) q[1];
sx q[1];
rz(-1.1898142) q[1];
sx q[1];
rz(1.9728647) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334621) q[0];
sx q[0];
rz(-2.4442844) q[0];
sx q[0];
rz(-1.4278317) q[0];
x q[1];
rz(-0.0059466023) q[2];
sx q[2];
rz(-2.0281254) q[2];
sx q[2];
rz(-1.3471239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1064042) q[1];
sx q[1];
rz(-2.4420083) q[1];
sx q[1];
rz(-2.0856218) q[1];
rz(-pi) q[2];
rz(-0.46242373) q[3];
sx q[3];
rz(-0.96721691) q[3];
sx q[3];
rz(1.0360595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0231102) q[2];
sx q[2];
rz(-2.1239943) q[2];
sx q[2];
rz(-2.6712096) q[2];
rz(1.6810301) q[3];
sx q[3];
rz(-1.2602256) q[3];
sx q[3];
rz(2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39101446) q[0];
sx q[0];
rz(-1.3704726) q[0];
sx q[0];
rz(-1.4960666) q[0];
rz(-2.0059313) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(2.6527203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.931536) q[0];
sx q[0];
rz(-1.6029412) q[0];
sx q[0];
rz(1.710762) q[0];
rz(2.082294) q[2];
sx q[2];
rz(-1.6548205) q[2];
sx q[2];
rz(0.78843278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.012163398) q[1];
sx q[1];
rz(-1.5930854) q[1];
sx q[1];
rz(-0.049456923) q[1];
rz(-0.094303207) q[3];
sx q[3];
rz(-0.62178388) q[3];
sx q[3];
rz(2.0981385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44498542) q[2];
sx q[2];
rz(-2.5312436) q[2];
sx q[2];
rz(1.102591) q[2];
rz(1.5251112) q[3];
sx q[3];
rz(-2.3307255) q[3];
sx q[3];
rz(-1.6703687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0082323) q[0];
sx q[0];
rz(-2.6120549) q[0];
sx q[0];
rz(-1.3326921) q[0];
rz(2.7545199) q[1];
sx q[1];
rz(-1.2553071) q[1];
sx q[1];
rz(2.0852087) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1176193) q[0];
sx q[0];
rz(-1.7101544) q[0];
sx q[0];
rz(1.8169212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41802539) q[2];
sx q[2];
rz(-1.4565832) q[2];
sx q[2];
rz(0.25164652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7596408) q[1];
sx q[1];
rz(-2.3084062) q[1];
sx q[1];
rz(-0.52565378) q[1];
rz(-2.4664573) q[3];
sx q[3];
rz(-2.0656485) q[3];
sx q[3];
rz(-2.2371045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1308088) q[2];
sx q[2];
rz(-0.84785145) q[2];
sx q[2];
rz(1.7424142) q[2];
rz(1.0731267) q[3];
sx q[3];
rz(-1.9173887) q[3];
sx q[3];
rz(2.6789902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90049967) q[0];
sx q[0];
rz(-1.6443962) q[0];
sx q[0];
rz(-1.6045438) q[0];
rz(-2.4201139) q[1];
sx q[1];
rz(-0.57054467) q[1];
sx q[1];
rz(1.9346938) q[1];
rz(1.1804562) q[2];
sx q[2];
rz(-1.9705638) q[2];
sx q[2];
rz(0.86314291) q[2];
rz(-0.7714959) q[3];
sx q[3];
rz(-1.6416999) q[3];
sx q[3];
rz(0.083569817) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
