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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096199252) q[0];
sx q[0];
rz(-1.081776) q[0];
sx q[0];
rz(-3.0351991) q[0];
rz(-pi) q[1];
x q[1];
rz(0.098564221) q[2];
sx q[2];
rz(-1.9910887) q[2];
sx q[2];
rz(1.9914876) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9086985) q[1];
sx q[1];
rz(-1.9373937) q[1];
sx q[1];
rz(2.7243549) q[1];
x q[2];
rz(-0.75593573) q[3];
sx q[3];
rz(-1.9256251) q[3];
sx q[3];
rz(-2.1906022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52557785) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(0.40526059) q[2];
rz(-1.1303834) q[3];
sx q[3];
rz(-2.3519792) q[3];
sx q[3];
rz(-1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2363215) q[0];
sx q[0];
rz(-0.032624809) q[0];
sx q[0];
rz(-2.3765833) q[0];
rz(-0.24461034) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(2.5027221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8804358) q[0];
sx q[0];
rz(-0.64468274) q[0];
sx q[0];
rz(3.0993483) q[0];
rz(-pi) q[1];
rz(-2.5193754) q[2];
sx q[2];
rz(-1.2583311) q[2];
sx q[2];
rz(-2.0553596) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8117026) q[1];
sx q[1];
rz(-2.5479758) q[1];
sx q[1];
rz(0.65887405) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83373006) q[3];
sx q[3];
rz(-1.7639188) q[3];
sx q[3];
rz(-1.5964791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2834393) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-2.8972674) q[2];
rz(-2.0325932) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37794381) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(1.2373244) q[0];
rz(-1.411865) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(1.3317187) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.004382172) q[0];
sx q[0];
rz(-1.4829652) q[0];
sx q[0];
rz(-0.77465221) q[0];
rz(-pi) q[1];
rz(1.026451) q[2];
sx q[2];
rz(-2.7371001) q[2];
sx q[2];
rz(-0.23032388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0983372) q[1];
sx q[1];
rz(-1.6041293) q[1];
sx q[1];
rz(0.055524932) q[1];
rz(0.94966763) q[3];
sx q[3];
rz(-2.0916208) q[3];
sx q[3];
rz(0.68062147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9776913) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(0.65995556) q[2];
rz(-1.6948949) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25809836) q[0];
sx q[0];
rz(-2.3907008) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(0.54905999) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66019193) q[0];
sx q[0];
rz(-0.21960077) q[0];
sx q[0];
rz(2.289413) q[0];
rz(-pi) q[1];
rz(-1.8070776) q[2];
sx q[2];
rz(-0.60271128) q[2];
sx q[2];
rz(-2.757795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0564894) q[1];
sx q[1];
rz(-1.6360549) q[1];
sx q[1];
rz(0.27537055) q[1];
rz(-pi) q[2];
rz(0.96499126) q[3];
sx q[3];
rz(-0.47156358) q[3];
sx q[3];
rz(2.4237867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0893726) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(0.29903665) q[2];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-3.1035677) q[0];
sx q[0];
rz(-2.6755565) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(0.83494157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73663086) q[0];
sx q[0];
rz(-1.6225272) q[0];
sx q[0];
rz(1.5993656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72795002) q[2];
sx q[2];
rz(-0.11874983) q[2];
sx q[2];
rz(2.438022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7606733) q[1];
sx q[1];
rz(-1.6242508) q[1];
sx q[1];
rz(1.0905114) q[1];
x q[2];
rz(2.4598875) q[3];
sx q[3];
rz(-0.85452628) q[3];
sx q[3];
rz(1.1151506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9338108) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(-0.031522838) q[2];
rz(0.75782123) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(-2.5789564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9071478) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(-0.031540792) q[0];
rz(-3.0517598) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(-2.8057742) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851097) q[0];
sx q[0];
rz(-1.2439737) q[0];
sx q[0];
rz(-2.7669897) q[0];
rz(-pi) q[1];
rz(-0.33863314) q[2];
sx q[2];
rz(-1.0446608) q[2];
sx q[2];
rz(-2.3735685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5988859) q[1];
sx q[1];
rz(-2.5395091) q[1];
sx q[1];
rz(-1.1335526) q[1];
rz(-2.8434136) q[3];
sx q[3];
rz(-2.0070873) q[3];
sx q[3];
rz(2.0661708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.971656) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(1.8709095) q[2];
rz(-0.058567889) q[3];
sx q[3];
rz(-1.6978369) q[3];
sx q[3];
rz(-2.7931255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.612065) q[0];
sx q[0];
rz(-2.6030774) q[0];
sx q[0];
rz(-2.8756323) q[0];
rz(2.3174441) q[1];
sx q[1];
rz(-1.3105086) q[1];
sx q[1];
rz(1.5305653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96722051) q[0];
sx q[0];
rz(-1.3200813) q[0];
sx q[0];
rz(-2.870737) q[0];
rz(-0.73101298) q[2];
sx q[2];
rz(-1.5314529) q[2];
sx q[2];
rz(-2.524216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2756259) q[1];
sx q[1];
rz(-1.8339364) q[1];
sx q[1];
rz(2.5530173) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3694088) q[3];
sx q[3];
rz(-1.3150384) q[3];
sx q[3];
rz(-1.0061177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.019235762) q[2];
sx q[2];
rz(-2.0255721) q[2];
sx q[2];
rz(2.1717211) q[2];
rz(-2.9979749) q[3];
sx q[3];
rz(-2.2031281) q[3];
sx q[3];
rz(-2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.403991) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(-0.10424374) q[0];
rz(2.3999816) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(0.43824497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40566984) q[0];
sx q[0];
rz(-2.4614262) q[0];
sx q[0];
rz(-0.16323547) q[0];
rz(-pi) q[1];
rz(-2.9061658) q[2];
sx q[2];
rz(-2.1554783) q[2];
sx q[2];
rz(0.92008051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.150369) q[1];
sx q[1];
rz(-0.60661483) q[1];
sx q[1];
rz(1.6305417) q[1];
rz(2.6062268) q[3];
sx q[3];
rz(-1.6965877) q[3];
sx q[3];
rz(2.4411502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40615842) q[2];
sx q[2];
rz(-2.8696852) q[2];
sx q[2];
rz(-2.820365) q[2];
rz(-0.99212232) q[3];
sx q[3];
rz(-1.0935874) q[3];
sx q[3];
rz(1.7295624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2883437) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(0.14800063) q[0];
rz(-1.1389698) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(2.151087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97738701) q[0];
sx q[0];
rz(-1.6456104) q[0];
sx q[0];
rz(2.6590024) q[0];
x q[1];
rz(-1.1075145) q[2];
sx q[2];
rz(-2.2474237) q[2];
sx q[2];
rz(2.7156242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79117486) q[1];
sx q[1];
rz(-1.5656078) q[1];
sx q[1];
rz(-1.5395916) q[1];
rz(-pi) q[2];
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
rz(1.1691576) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(-0.54231826) q[2];
rz(1.4646685) q[3];
sx q[3];
rz(-1.9138391) q[3];
sx q[3];
rz(0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(-0.78508776) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(-0.40207544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94984834) q[0];
sx q[0];
rz(-1.254891) q[0];
sx q[0];
rz(-2.6788968) q[0];
rz(-pi) q[1];
rz(2.3763972) q[2];
sx q[2];
rz(-0.94586654) q[2];
sx q[2];
rz(3.011572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2759034) q[1];
sx q[1];
rz(-0.69178693) q[1];
sx q[1];
rz(-1.2671641) q[1];
rz(-pi) q[2];
rz(1.1583331) q[3];
sx q[3];
rz(-1.8329046) q[3];
sx q[3];
rz(-2.7893808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1606007) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(-0.55046925) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.8615104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.078770272) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(-1.5863905) q[2];
sx q[2];
rz(-0.77257665) q[2];
sx q[2];
rz(2.0333124) q[2];
rz(1.3435626) q[3];
sx q[3];
rz(-0.5323635) q[3];
sx q[3];
rz(-1.3743286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
