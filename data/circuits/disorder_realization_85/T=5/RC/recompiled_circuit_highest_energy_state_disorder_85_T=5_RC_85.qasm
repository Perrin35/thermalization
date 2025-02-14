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
rz(-2.7388465) q[1];
sx q[1];
rz(-1.4637113) q[1];
sx q[1];
rz(0.95199624) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453934) q[0];
sx q[0];
rz(-2.0598167) q[0];
sx q[0];
rz(-0.10639356) q[0];
rz(-1.1486886) q[2];
sx q[2];
rz(-1.4808345) q[2];
sx q[2];
rz(0.46101704) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0178899) q[1];
sx q[1];
rz(-2.5933939) q[1];
sx q[1];
rz(-2.3831559) q[1];
x q[2];
rz(-0.75593573) q[3];
sx q[3];
rz(-1.2159676) q[3];
sx q[3];
rz(-0.95099041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6160148) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(2.7363321) q[2];
rz(-2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.214504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2363215) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(-2.3765833) q[0];
rz(-2.8969823) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(0.6388706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8804358) q[0];
sx q[0];
rz(-0.64468274) q[0];
sx q[0];
rz(0.042244367) q[0];
rz(-pi) q[1];
rz(-1.1924001) q[2];
sx q[2];
rz(-0.98289433) q[2];
sx q[2];
rz(0.70158106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42145211) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(-1.1790183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25821553) q[3];
sx q[3];
rz(-0.85047837) q[3];
sx q[3];
rz(0.19816601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2834393) q[2];
sx q[2];
rz(-0.44698134) q[2];
sx q[2];
rz(-0.24432527) q[2];
rz(-1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(-0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7636488) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(-1.7297277) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(1.809874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1372105) q[0];
sx q[0];
rz(-1.6586275) q[0];
sx q[0];
rz(2.3669404) q[0];
rz(-pi) q[1];
rz(-2.1151416) q[2];
sx q[2];
rz(-2.7371001) q[2];
sx q[2];
rz(2.9112688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52568846) q[1];
sx q[1];
rz(-1.6262904) q[1];
sx q[1];
rz(-1.5374119) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5272433) q[3];
sx q[3];
rz(-1.0416556) q[3];
sx q[3];
rz(-1.2322616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1639013) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(0.65995556) q[2];
rz(-1.6948949) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-2.0383294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25809836) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(2.6253888) q[1];
sx q[1];
rz(-1.3092923) q[1];
sx q[1];
rz(-0.54905999) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66019193) q[0];
sx q[0];
rz(-0.21960077) q[0];
sx q[0];
rz(-2.289413) q[0];
x q[1];
rz(-1.334515) q[2];
sx q[2];
rz(-2.5388814) q[2];
sx q[2];
rz(-2.757795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.258865) q[1];
sx q[1];
rz(-2.8587864) q[1];
sx q[1];
rz(0.23587366) q[1];
x q[2];
rz(-2.1766014) q[3];
sx q[3];
rz(-0.47156358) q[3];
sx q[3];
rz(-0.71780598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.05222008) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(2.842556) q[2];
rz(0.61797577) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95737902) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(2.6755565) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(-2.3066511) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4049618) q[0];
sx q[0];
rz(-1.6225272) q[0];
sx q[0];
rz(1.542227) q[0];
rz(-pi) q[1];
rz(-2.4136426) q[2];
sx q[2];
rz(-0.11874983) q[2];
sx q[2];
rz(2.438022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8495167) q[1];
sx q[1];
rz(-2.6585732) q[1];
sx q[1];
rz(1.6860875) q[1];
rz(0.72839177) q[3];
sx q[3];
rz(-2.0660704) q[3];
sx q[3];
rz(-3.107695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20778188) q[2];
sx q[2];
rz(-2.2489579) q[2];
sx q[2];
rz(-0.031522838) q[2];
rz(-2.3837714) q[3];
sx q[3];
rz(-1.1008681) q[3];
sx q[3];
rz(-0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.23444489) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(-3.0517598) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(-2.8057742) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15648297) q[0];
sx q[0];
rz(-1.2439737) q[0];
sx q[0];
rz(-2.7669897) q[0];
rz(-2.1226826) q[2];
sx q[2];
rz(-1.2794211) q[2];
sx q[2];
rz(0.62770972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5988859) q[1];
sx q[1];
rz(-0.60208354) q[1];
sx q[1];
rz(-1.1335526) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0246216) q[3];
sx q[3];
rz(-1.3012851) q[3];
sx q[3];
rz(-2.5170642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-2.0270831) q[2];
sx q[2];
rz(-1.8709095) q[2];
rz(-0.058567889) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(2.7931255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(0.2659604) q[0];
rz(-2.3174441) q[1];
sx q[1];
rz(-1.3105086) q[1];
sx q[1];
rz(1.6110274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3327681) q[0];
sx q[0];
rz(-2.7746189) q[0];
sx q[0];
rz(0.76353188) q[0];
x q[1];
rz(-1.62362) q[2];
sx q[2];
rz(-0.84047707) q[2];
sx q[2];
rz(-2.2234302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87674031) q[1];
sx q[1];
rz(-2.1365668) q[1];
sx q[1];
rz(1.2575722) q[1];
rz(-pi) q[2];
rz(1.2208185) q[3];
sx q[3];
rz(-0.82982291) q[3];
sx q[3];
rz(-2.8184824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1223569) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(0.96987152) q[2];
rz(0.14361778) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73760167) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(-2.3999816) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(2.7033477) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377698) q[0];
sx q[0];
rz(-1.4684104) q[0];
sx q[0];
rz(2.4679604) q[0];
rz(-0.97317071) q[2];
sx q[2];
rz(-1.3750374) q[2];
sx q[2];
rz(2.6225066) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5129152) q[1];
sx q[1];
rz(-1.6048429) q[1];
sx q[1];
rz(2.176575) q[1];
x q[2];
rz(-2.6062268) q[3];
sx q[3];
rz(-1.4450049) q[3];
sx q[3];
rz(-0.70044242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40615842) q[2];
sx q[2];
rz(-2.8696852) q[2];
sx q[2];
rz(-0.3212277) q[2];
rz(2.1494703) q[3];
sx q[3];
rz(-1.0935874) q[3];
sx q[3];
rz(1.7295624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883437) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(-0.14800063) q[0];
rz(2.0026228) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(2.151087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873224) q[0];
sx q[0];
rz(-2.0519216) q[0];
sx q[0];
rz(1.4863798) q[0];
x q[1];
rz(2.6338088) q[2];
sx q[2];
rz(-2.3427026) q[2];
sx q[2];
rz(0.2474167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3504178) q[1];
sx q[1];
rz(-1.5656078) q[1];
sx q[1];
rz(1.6020011) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0664987) q[3];
sx q[3];
rz(-1.9073581) q[3];
sx q[3];
rz(-0.15777212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(-0.54231826) q[2];
rz(1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(-0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9024314) q[0];
sx q[0];
rz(-2.3434174) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(-0.78508776) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(-2.7395172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94984834) q[0];
sx q[0];
rz(-1.254891) q[0];
sx q[0];
rz(2.6788968) q[0];
x q[1];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2520894) q[1];
sx q[1];
rz(-2.2252603) q[1];
sx q[1];
rz(-2.8988125) q[1];
x q[2];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98099199) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(-2.215019) q[2];
rz(2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4704623) q[0];
sx q[0];
rz(-1.067938) q[0];
sx q[0];
rz(-2.9965042) q[0];
rz(-3.0628224) q[1];
sx q[1];
rz(-1.7369743) q[1];
sx q[1];
rz(-1.8440934) q[1];
rz(-0.79828046) q[2];
sx q[2];
rz(-1.5816806) q[2];
sx q[2];
rz(0.47368373) q[2];
rz(1.79803) q[3];
sx q[3];
rz(-2.6092292) q[3];
sx q[3];
rz(1.7672641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
