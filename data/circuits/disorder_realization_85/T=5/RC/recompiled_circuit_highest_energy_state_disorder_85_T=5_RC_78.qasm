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
rz(-0.83952251) q[0];
sx q[0];
rz(-1.5877725) q[0];
sx q[0];
rz(-2.177218) q[0];
rz(0.40274611) q[1];
sx q[1];
rz(1.4637113) q[1];
sx q[1];
rz(7.2351815) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453934) q[0];
sx q[0];
rz(-2.0598167) q[0];
sx q[0];
rz(3.0351991) q[0];
rz(-pi) q[1];
rz(1.9929041) q[2];
sx q[2];
rz(-1.6607581) q[2];
sx q[2];
rz(-0.46101704) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0178899) q[1];
sx q[1];
rz(-0.54819878) q[1];
sx q[1];
rz(0.75843673) q[1];
rz(-pi) q[2];
rz(-2.0417781) q[3];
sx q[3];
rz(-0.87207651) q[3];
sx q[3];
rz(-0.93633151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52557785) q[2];
sx q[2];
rz(-1.3495477) q[2];
sx q[2];
rz(2.7363321) q[2];
rz(2.0112093) q[3];
sx q[3];
rz(-2.3519792) q[3];
sx q[3];
rz(-1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.2363215) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(2.3765833) q[0];
rz(-0.24461034) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(-0.6388706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9332704) q[0];
sx q[0];
rz(-0.92678419) q[0];
sx q[0];
rz(1.5390551) q[0];
rz(-0.62221726) q[2];
sx q[2];
rz(-1.8832616) q[2];
sx q[2];
rz(1.086233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32989001) q[1];
sx q[1];
rz(-0.59361688) q[1];
sx q[1];
rz(2.4827186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2876533) q[3];
sx q[3];
rz(-0.75732175) q[3];
sx q[3];
rz(-0.18268302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8581533) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(2.8972674) q[2];
rz(-1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(2.9037156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37794381) q[0];
sx q[0];
rz(-1.0841333) q[0];
sx q[0];
rz(-1.2373244) q[0];
rz(1.7297277) q[1];
sx q[1];
rz(-2.6928163) q[1];
sx q[1];
rz(-1.3317187) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1372105) q[0];
sx q[0];
rz(-1.6586275) q[0];
sx q[0];
rz(0.77465221) q[0];
x q[1];
rz(1.9218512) q[2];
sx q[2];
rz(-1.3655541) q[2];
sx q[2];
rz(-2.308941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6159042) q[1];
sx q[1];
rz(-1.6262904) q[1];
sx q[1];
rz(-1.6041808) q[1];
rz(-pi) q[2];
rz(2.5272433) q[3];
sx q[3];
rz(-2.099937) q[3];
sx q[3];
rz(1.2322616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9776913) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(2.4816371) q[2];
rz(-1.4466977) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25809836) q[0];
sx q[0];
rz(-2.3907008) q[0];
sx q[0];
rz(1.1335491) q[0];
rz(0.51620382) q[1];
sx q[1];
rz(-1.3092923) q[1];
sx q[1];
rz(-2.5925327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3908729) q[0];
sx q[0];
rz(-1.4060806) q[0];
sx q[0];
rz(0.14589845) q[0];
rz(-0.15971036) q[2];
sx q[2];
rz(-2.1545112) q[2];
sx q[2];
rz(0.6681548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.085103206) q[1];
sx q[1];
rz(-1.6360549) q[1];
sx q[1];
rz(0.27537055) q[1];
rz(-2.1766014) q[3];
sx q[3];
rz(-0.47156358) q[3];
sx q[3];
rz(-0.71780598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05222008) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(-2.842556) q[2];
rz(0.61797577) q[3];
sx q[3];
rz(-2.0682014) q[3];
sx q[3];
rz(-1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(2.6755565) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-0.51141089) q[1];
sx q[1];
rz(2.3066511) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73663086) q[0];
sx q[0];
rz(-1.6225272) q[0];
sx q[0];
rz(-1.542227) q[0];
rz(2.4136426) q[2];
sx q[2];
rz(-3.0228428) q[2];
sx q[2];
rz(2.438022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7606733) q[1];
sx q[1];
rz(-1.6242508) q[1];
sx q[1];
rz(1.0905114) q[1];
rz(2.4598875) q[3];
sx q[3];
rz(-0.85452628) q[3];
sx q[3];
rz(-2.0264421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9338108) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(0.031522838) q[2];
rz(-2.3837714) q[3];
sx q[3];
rz(-1.1008681) q[3];
sx q[3];
rz(2.5789564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.23444489) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(0.08983287) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(0.33581844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528362) q[0];
sx q[0];
rz(-1.2169219) q[0];
sx q[0];
rz(1.2214946) q[0];
rz(-pi) q[1];
rz(2.0903942) q[2];
sx q[2];
rz(-2.5246387) q[2];
sx q[2];
rz(-1.3796356) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5427067) q[1];
sx q[1];
rz(-2.5395091) q[1];
sx q[1];
rz(1.1335526) q[1];
x q[2];
rz(-0.29817902) q[3];
sx q[3];
rz(-1.1345053) q[3];
sx q[3];
rz(2.0661708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-2.0270831) q[2];
sx q[2];
rz(-1.8709095) q[2];
rz(3.0830248) q[3];
sx q[3];
rz(-1.6978369) q[3];
sx q[3];
rz(0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(0.2659604) q[0];
rz(2.3174441) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(-1.5305653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8088246) q[0];
sx q[0];
rz(-0.36697373) q[0];
sx q[0];
rz(2.3780608) q[0];
rz(-pi) q[1];
x q[1];
rz(1.62362) q[2];
sx q[2];
rz(-2.3011156) q[2];
sx q[2];
rz(0.91816245) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8083473) q[1];
sx q[1];
rz(-0.63831224) q[1];
sx q[1];
rz(-0.45175938) q[1];
x q[2];
rz(0.77218382) q[3];
sx q[3];
rz(-1.8265542) q[3];
sx q[3];
rz(1.0061177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.019235762) q[2];
sx q[2];
rz(-2.0255721) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(-0.14361778) q[3];
sx q[3];
rz(-2.2031281) q[3];
sx q[3];
rz(-0.6374878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.403991) q[0];
sx q[0];
rz(-0.25594512) q[0];
sx q[0];
rz(-0.10424374) q[0];
rz(2.3999816) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(2.7033477) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377698) q[0];
sx q[0];
rz(-1.4684104) q[0];
sx q[0];
rz(-0.67363221) q[0];
rz(-pi) q[1];
rz(0.23542685) q[2];
sx q[2];
rz(-0.98611438) q[2];
sx q[2];
rz(2.2215121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2230473) q[1];
sx q[1];
rz(-2.1761736) q[1];
sx q[1];
rz(0.041408509) q[1];
rz(-0.53536587) q[3];
sx q[3];
rz(-1.6965877) q[3];
sx q[3];
rz(-0.70044242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7354342) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(0.3212277) q[2];
rz(-0.99212232) q[3];
sx q[3];
rz(-1.0935874) q[3];
sx q[3];
rz(1.7295624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2883437) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(-2.993592) q[0];
rz(-2.0026228) q[1];
sx q[1];
rz(-2.482246) q[1];
sx q[1];
rz(2.151087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73511998) q[0];
sx q[0];
rz(-0.48790259) q[0];
sx q[0];
rz(-2.9814629) q[0];
rz(-pi) q[1];
rz(2.0340781) q[2];
sx q[2];
rz(-2.2474237) q[2];
sx q[2];
rz(-0.42596841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5266853) q[1];
sx q[1];
rz(-0.03163306) q[1];
sx q[1];
rz(-1.7355914) q[1];
x q[2];
rz(2.2049974) q[3];
sx q[3];
rz(-2.5504125) q[3];
sx q[3];
rz(-1.9612966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.972435) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(2.5992744) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(-2.1577468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(0.78508776) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(2.7395172) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063906625) q[0];
sx q[0];
rz(-2.5878667) q[0];
sx q[0];
rz(-0.63204773) q[0];
rz(-pi) q[1];
rz(2.3359012) q[2];
sx q[2];
rz(-0.94599722) q[2];
sx q[2];
rz(-2.2479518) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86568923) q[1];
sx q[1];
rz(-2.4498057) q[1];
sx q[1];
rz(1.2671641) q[1];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.1732297) q[3];
sx q[3];
rz(-1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1606007) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(-2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67113039) q[0];
sx q[0];
rz(-2.0736546) q[0];
sx q[0];
rz(0.14508844) q[0];
rz(-0.078770272) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(-1.5552022) q[2];
sx q[2];
rz(-2.369016) q[2];
sx q[2];
rz(-1.1082802) q[2];
rz(2.0918431) q[3];
sx q[3];
rz(-1.6853942) q[3];
sx q[3];
rz(-2.7484721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
