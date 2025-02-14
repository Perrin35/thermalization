OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7795774) q[0];
sx q[0];
rz(-0.22688046) q[0];
sx q[0];
rz(1.7664631) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(0.28092608) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37051113) q[0];
sx q[0];
rz(-2.2265052) q[0];
sx q[0];
rz(-0.20072584) q[0];
rz(-1.8616025) q[2];
sx q[2];
rz(-1.6806597) q[2];
sx q[2];
rz(1.8786877) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.0065767584) q[1];
sx q[1];
rz(-2.6408051) q[1];
sx q[1];
rz(2.1441475) q[1];
x q[2];
rz(-2.6219764) q[3];
sx q[3];
rz(-2.5027788) q[3];
sx q[3];
rz(-2.6878302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0278339) q[2];
sx q[2];
rz(-1.7366624) q[2];
sx q[2];
rz(-0.11715451) q[2];
rz(-2.8095918) q[3];
sx q[3];
rz(-0.74688512) q[3];
sx q[3];
rz(-0.78506708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(-2.7681328) q[0];
rz(-0.39730486) q[1];
sx q[1];
rz(-1.9970857) q[1];
sx q[1];
rz(-0.73748803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1289559) q[0];
sx q[0];
rz(-1.0397433) q[0];
sx q[0];
rz(-2.3319202) q[0];
rz(1.039871) q[2];
sx q[2];
rz(-1.0801407) q[2];
sx q[2];
rz(0.38035989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4026248) q[1];
sx q[1];
rz(-1.0254271) q[1];
sx q[1];
rz(1.7508372) q[1];
rz(0.37009671) q[3];
sx q[3];
rz(-2.3488999) q[3];
sx q[3];
rz(1.1350138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.069933683) q[2];
sx q[2];
rz(-1.6120913) q[2];
sx q[2];
rz(2.7640479) q[2];
rz(2.8318882) q[3];
sx q[3];
rz(-2.0524502) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51652235) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(1.9260433) q[0];
rz(-1.0141605) q[1];
sx q[1];
rz(-1.7182257) q[1];
sx q[1];
rz(-0.97022143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824066) q[0];
sx q[0];
rz(-2.0699887) q[0];
sx q[0];
rz(-0.046120709) q[0];
x q[1];
rz(0.54527905) q[2];
sx q[2];
rz(-0.98731326) q[2];
sx q[2];
rz(2.3346221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41170317) q[1];
sx q[1];
rz(-1.3940934) q[1];
sx q[1];
rz(-1.2248301) q[1];
rz(-0.35200624) q[3];
sx q[3];
rz(-1.4624634) q[3];
sx q[3];
rz(-1.7217404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0936475) q[2];
sx q[2];
rz(-0.84257546) q[2];
sx q[2];
rz(2.688431) q[2];
rz(1.2293182) q[3];
sx q[3];
rz(-1.8639576) q[3];
sx q[3];
rz(0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969084) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(-0.36566439) q[0];
rz(-0.91224313) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(-0.40547392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039818) q[0];
sx q[0];
rz(-1.3657327) q[0];
sx q[0];
rz(0.12139856) q[0];
x q[1];
rz(0.22511668) q[2];
sx q[2];
rz(-0.76864132) q[2];
sx q[2];
rz(-2.3942238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.34205375) q[1];
sx q[1];
rz(-2.5785682) q[1];
sx q[1];
rz(-2.8232226) q[1];
rz(2.502127) q[3];
sx q[3];
rz(-1.829095) q[3];
sx q[3];
rz(-2.4909508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.071986467) q[2];
sx q[2];
rz(-1.4402086) q[2];
sx q[2];
rz(-1.5427422) q[2];
rz(-0.37374464) q[3];
sx q[3];
rz(-1.6662686) q[3];
sx q[3];
rz(1.6811949) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986346) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(2.4454818) q[0];
rz(1.8822582) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(2.7764244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8445963) q[0];
sx q[0];
rz(-1.1410603) q[0];
sx q[0];
rz(0.55083042) q[0];
rz(-2.8234473) q[2];
sx q[2];
rz(-1.1765624) q[2];
sx q[2];
rz(0.94903681) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2303673) q[1];
sx q[1];
rz(-0.54552286) q[1];
sx q[1];
rz(3.1389152) q[1];
x q[2];
rz(-1.9837512) q[3];
sx q[3];
rz(-0.98414183) q[3];
sx q[3];
rz(-2.0362622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68088561) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(0.17670259) q[2];
rz(1.3867311) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(-0.37469125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143167) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(1.5981307) q[0];
rz(1.8371001) q[1];
sx q[1];
rz(-2.0871128) q[1];
sx q[1];
rz(2.3354882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62063187) q[0];
sx q[0];
rz(-1.3343108) q[0];
sx q[0];
rz(-2.3804401) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88509934) q[2];
sx q[2];
rz(-2.4174567) q[2];
sx q[2];
rz(0.76937461) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1404952) q[1];
sx q[1];
rz(-1.6204658) q[1];
sx q[1];
rz(-1.7752663) q[1];
x q[2];
rz(0.57512734) q[3];
sx q[3];
rz(-2.025423) q[3];
sx q[3];
rz(-1.8607339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1690037) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(-0.063610323) q[2];
rz(2.5566067) q[3];
sx q[3];
rz(-1.7413185) q[3];
sx q[3];
rz(-1.6271648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1503898) q[0];
sx q[0];
rz(-1.3911893) q[0];
sx q[0];
rz(0.72917953) q[0];
rz(-1.9623307) q[1];
sx q[1];
rz(-2.3919892) q[1];
sx q[1];
rz(-2.8470305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33074327) q[0];
sx q[0];
rz(-2.0669575) q[0];
sx q[0];
rz(1.8627573) q[0];
rz(0.11004098) q[2];
sx q[2];
rz(-2.0015284) q[2];
sx q[2];
rz(-0.34039341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9443317) q[1];
sx q[1];
rz(-0.75350584) q[1];
sx q[1];
rz(-1.7427518) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90055777) q[3];
sx q[3];
rz(-1.8769185) q[3];
sx q[3];
rz(-0.30508074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4906759) q[2];
sx q[2];
rz(-1.7016405) q[2];
sx q[2];
rz(-2.4093464) q[2];
rz(0.8693153) q[3];
sx q[3];
rz(-1.1704051) q[3];
sx q[3];
rz(-1.7349582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.231584) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(2.3439132) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-0.99635092) q[1];
sx q[1];
rz(1.4083883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080264576) q[0];
sx q[0];
rz(-1.3664075) q[0];
sx q[0];
rz(-0.56046446) q[0];
x q[1];
rz(-1.1154091) q[2];
sx q[2];
rz(-2.6386127) q[2];
sx q[2];
rz(-0.044654559) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8127784) q[1];
sx q[1];
rz(-1.8056889) q[1];
sx q[1];
rz(0.77130227) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3878294) q[3];
sx q[3];
rz(-1.8865693) q[3];
sx q[3];
rz(2.0980199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.306376) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(-0.19732538) q[2];
rz(-2.3399682) q[3];
sx q[3];
rz(-0.97951952) q[3];
sx q[3];
rz(2.7200123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1590969) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(-2.2267447) q[0];
rz(0.21028701) q[1];
sx q[1];
rz(-2.5631914) q[1];
sx q[1];
rz(0.68967825) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1062968) q[0];
sx q[0];
rz(-1.6637422) q[0];
sx q[0];
rz(-2.9151117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014691512) q[2];
sx q[2];
rz(-1.4533051) q[2];
sx q[2];
rz(-1.2047639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27305973) q[1];
sx q[1];
rz(-1.0976205) q[1];
sx q[1];
rz(-0.57762967) q[1];
rz(2.1237064) q[3];
sx q[3];
rz(-0.59846557) q[3];
sx q[3];
rz(-2.726647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0476394) q[2];
sx q[2];
rz(-2.2932105) q[2];
sx q[2];
rz(-0.32996714) q[2];
rz(2.0983569) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030815) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(1.0571085) q[0];
rz(-2.8874176) q[1];
sx q[1];
rz(-1.9994241) q[1];
sx q[1];
rz(-2.4370297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0443665) q[0];
sx q[0];
rz(-1.2597879) q[0];
sx q[0];
rz(-2.0661656) q[0];
rz(-1.754934) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(-2.9764701) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5479991) q[1];
sx q[1];
rz(-1.3915466) q[1];
sx q[1];
rz(1.1186734) q[1];
rz(-1.390144) q[3];
sx q[3];
rz(-1.7757987) q[3];
sx q[3];
rz(-1.6222749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6331943) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(0.50676662) q[2];
rz(-2.9554328) q[3];
sx q[3];
rz(-2.9045744) q[3];
sx q[3];
rz(-2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7645466) q[0];
sx q[0];
rz(-1.6896387) q[0];
sx q[0];
rz(1.1402546) q[0];
rz(-2.2346732) q[1];
sx q[1];
rz(-2.4005371) q[1];
sx q[1];
rz(-1.3225318) q[1];
rz(0.53756164) q[2];
sx q[2];
rz(-2.1765709) q[2];
sx q[2];
rz(-1.4449262) q[2];
rz(-0.68578699) q[3];
sx q[3];
rz(-2.7332173) q[3];
sx q[3];
rz(0.12242534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
