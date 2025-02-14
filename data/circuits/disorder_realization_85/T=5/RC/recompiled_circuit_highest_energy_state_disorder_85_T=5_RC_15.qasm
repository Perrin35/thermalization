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
rz(0.96437469) q[0];
rz(0.40274611) q[1];
sx q[1];
rz(-1.6778814) q[1];
sx q[1];
rz(2.1895964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6168686) q[0];
sx q[0];
rz(-1.6646806) q[0];
sx q[0];
rz(1.0794207) q[0];
x q[1];
rz(1.3540719) q[2];
sx q[2];
rz(-0.43102362) q[2];
sx q[2];
rz(-2.2292516) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9086985) q[1];
sx q[1];
rz(-1.9373937) q[1];
sx q[1];
rz(-0.41723771) q[1];
x q[2];
rz(2.0417781) q[3];
sx q[3];
rz(-2.2695161) q[3];
sx q[3];
rz(-0.93633151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52557785) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(-0.40526059) q[2];
rz(2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90527117) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(2.3765833) q[0];
rz(2.8969823) q[1];
sx q[1];
rz(-1.8965992) q[1];
sx q[1];
rz(-2.5027221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2611569) q[0];
sx q[0];
rz(-2.4969099) q[0];
sx q[0];
rz(0.042244367) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62221726) q[2];
sx q[2];
rz(-1.2583311) q[2];
sx q[2];
rz(1.086233) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7201405) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(-1.1790183) q[1];
x q[2];
rz(-0.25821553) q[3];
sx q[3];
rz(-0.85047837) q[3];
sx q[3];
rz(-0.19816601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2834393) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-2.8972674) q[2];
rz(1.1089995) q[3];
sx q[3];
rz(-1.8651483) q[3];
sx q[3];
rz(2.9037156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37794381) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(-1.7297277) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(1.809874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6645637) q[0];
sx q[0];
rz(-0.77858276) q[0];
sx q[0];
rz(-3.016359) q[0];
x q[1];
rz(1.2197415) q[2];
sx q[2];
rz(-1.3655541) q[2];
sx q[2];
rz(2.308941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0983372) q[1];
sx q[1];
rz(-1.5374633) q[1];
sx q[1];
rz(0.055524932) q[1];
x q[2];
rz(-0.79257719) q[3];
sx q[3];
rz(-0.78781414) q[3];
sx q[3];
rz(2.8590096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9776913) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(0.65995556) q[2];
rz(1.4466977) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-2.0383294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834943) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(-2.0080436) q[0];
rz(-2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(2.5925327) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66019193) q[0];
sx q[0];
rz(-0.21960077) q[0];
sx q[0];
rz(-0.85217969) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8070776) q[2];
sx q[2];
rz(-2.5388814) q[2];
sx q[2];
rz(0.38379764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0564894) q[1];
sx q[1];
rz(-1.6360549) q[1];
sx q[1];
rz(0.27537055) q[1];
rz(-pi) q[2];
rz(1.9677343) q[3];
sx q[3];
rz(-1.8324495) q[3];
sx q[3];
rz(1.4058827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0893726) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1842136) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(0.46603611) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(0.83494157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83564312) q[0];
sx q[0];
rz(-1.5993274) q[0];
sx q[0];
rz(3.0898407) q[0];
rz(-2.4136426) q[2];
sx q[2];
rz(-3.0228428) q[2];
sx q[2];
rz(0.70357067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38091936) q[1];
sx q[1];
rz(-1.6242508) q[1];
sx q[1];
rz(-2.0510813) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4132009) q[3];
sx q[3];
rz(-1.0755223) q[3];
sx q[3];
rz(0.033897696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9338108) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(-3.1100698) q[2];
rz(0.75782123) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9071478) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(-3.0517598) q[1];
sx q[1];
rz(-0.49982163) q[1];
sx q[1];
rz(2.8057742) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851097) q[0];
sx q[0];
rz(-1.2439737) q[0];
sx q[0];
rz(2.7669897) q[0];
rz(-pi) q[1];
rz(1.0511984) q[2];
sx q[2];
rz(-2.5246387) q[2];
sx q[2];
rz(-1.761957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5988859) q[1];
sx q[1];
rz(-2.5395091) q[1];
sx q[1];
rz(-2.00804) q[1];
rz(2.0246216) q[3];
sx q[3];
rz(-1.8403075) q[3];
sx q[3];
rz(-0.62452841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.971656) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(1.2706832) q[2];
rz(-3.0830248) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(-2.7931255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(-2.8756323) q[0];
rz(-0.82414857) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(1.6110274) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53479311) q[0];
sx q[0];
rz(-1.8329808) q[0];
sx q[0];
rz(1.3110089) q[0];
rz(-pi) q[1];
rz(-0.73101298) q[2];
sx q[2];
rz(-1.5314529) q[2];
sx q[2];
rz(-2.524216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8083473) q[1];
sx q[1];
rz(-0.63831224) q[1];
sx q[1];
rz(-0.45175938) q[1];
rz(-pi) q[2];
rz(-1.2208185) q[3];
sx q[3];
rz(-0.82982291) q[3];
sx q[3];
rz(-0.32311026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.019235762) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(0.96987152) q[2];
rz(-0.14361778) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(0.6374878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73760167) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(0.10424374) q[0];
rz(-2.3999816) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(0.43824497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5271664) q[0];
sx q[0];
rz(-2.2402555) q[0];
sx q[0];
rz(1.4400843) q[0];
rz(-1.2319698) q[2];
sx q[2];
rz(-0.6251339) q[2];
sx q[2];
rz(1.3300611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99122366) q[1];
sx q[1];
rz(-0.60661483) q[1];
sx q[1];
rz(-1.5110509) q[1];
rz(-pi) q[2];
rz(2.8986071) q[3];
sx q[3];
rz(-2.5930517) q[3];
sx q[3];
rz(-2.0627956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7354342) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(2.820365) q[2];
rz(-0.99212232) q[3];
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
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(2.993592) q[0];
rz(-2.0026228) q[1];
sx q[1];
rz(-2.482246) q[1];
sx q[1];
rz(-0.99050561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1642056) q[0];
sx q[0];
rz(-1.6456104) q[0];
sx q[0];
rz(-2.6590024) q[0];
x q[1];
rz(0.50778383) q[2];
sx q[2];
rz(-0.79889005) q[2];
sx q[2];
rz(-2.894176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79117486) q[1];
sx q[1];
rz(-1.5656078) q[1];
sx q[1];
rz(-1.5395916) q[1];
rz(-1.0750939) q[3];
sx q[3];
rz(-1.9073581) q[3];
sx q[3];
rz(-2.9838205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1691576) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(0.54231826) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.9138391) q[3];
sx q[3];
rz(-0.98384583) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9024314) q[0];
sx q[0];
rz(-2.3434174) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(-2.3565049) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(2.7395172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77469414) q[0];
sx q[0];
rz(-2.0089564) q[0];
sx q[0];
rz(1.9209981) q[0];
rz(-pi) q[1];
rz(0.80569141) q[2];
sx q[2];
rz(-2.1955954) q[2];
sx q[2];
rz(-2.2479518) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6732503) q[1];
sx q[1];
rz(-1.3788917) q[1];
sx q[1];
rz(-0.90190355) q[1];
rz(-pi) q[2];
rz(-0.28487408) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(-1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98099199) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(-2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(1.8615104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
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
rz(1.5863905) q[2];
sx q[2];
rz(-2.369016) q[2];
sx q[2];
rz(-1.1082802) q[2];
rz(-0.13194247) q[3];
sx q[3];
rz(-2.0880825) q[3];
sx q[3];
rz(-1.1121398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
