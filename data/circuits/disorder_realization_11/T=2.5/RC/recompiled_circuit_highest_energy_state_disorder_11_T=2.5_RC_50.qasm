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
rz(0.053057916) q[0];
sx q[0];
rz(-2.7795656) q[0];
sx q[0];
rz(1.9757353) q[0];
rz(2.2447658) q[1];
sx q[1];
rz(4.5936102) q[1];
sx q[1];
rz(11.138154) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4025164) q[0];
sx q[0];
rz(-1.5105643) q[0];
sx q[0];
rz(-0.51243685) q[0];
x q[1];
rz(-0.25144724) q[2];
sx q[2];
rz(-1.5000952) q[2];
sx q[2];
rz(1.9858907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1686656) q[1];
sx q[1];
rz(-1.4391237) q[1];
sx q[1];
rz(1.3802856) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0479508) q[3];
sx q[3];
rz(-2.6029498) q[3];
sx q[3];
rz(2.8847196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6894138) q[2];
sx q[2];
rz(-0.016760085) q[2];
sx q[2];
rz(2.9407799) q[2];
rz(0.14874841) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(-0.31140056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1397322) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(-1.0396022) q[0];
rz(0.014558583) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(1.5537517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4400875) q[0];
sx q[0];
rz(-2.0962414) q[0];
sx q[0];
rz(-0.82869014) q[0];
rz(-pi) q[1];
rz(1.7660242) q[2];
sx q[2];
rz(-3.0659817) q[2];
sx q[2];
rz(-1.6603254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0224198) q[1];
sx q[1];
rz(-1.3117562) q[1];
sx q[1];
rz(-3.1062192) q[1];
rz(-pi) q[2];
rz(0.33287853) q[3];
sx q[3];
rz(-1.7529704) q[3];
sx q[3];
rz(-1.6265488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62850922) q[2];
sx q[2];
rz(-1.5336978) q[2];
sx q[2];
rz(-1.7536564) q[2];
rz(-1.3758818) q[3];
sx q[3];
rz(-2.0949771) q[3];
sx q[3];
rz(2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.7498216) q[0];
sx q[0];
rz(-2.9048558) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(1.5982184) q[1];
sx q[1];
rz(-0.18074712) q[1];
sx q[1];
rz(-2.1764596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53837073) q[0];
sx q[0];
rz(-3.1206836) q[0];
sx q[0];
rz(-1.7666817) q[0];
rz(-pi) q[1];
rz(-0.94723126) q[2];
sx q[2];
rz(-2.5922225) q[2];
sx q[2];
rz(1.4538764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4016006) q[1];
sx q[1];
rz(-1.6066181) q[1];
sx q[1];
rz(-1.7154681) q[1];
rz(3.0583417) q[3];
sx q[3];
rz(-1.4752582) q[3];
sx q[3];
rz(-0.18886177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33660108) q[2];
sx q[2];
rz(-2.470863) q[2];
sx q[2];
rz(0.86656183) q[2];
rz(2.0349515) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57257819) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(1.6160075) q[0];
rz(-0.010604803) q[1];
sx q[1];
rz(-3.1378101) q[1];
sx q[1];
rz(-2.3912281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.01975) q[0];
sx q[0];
rz(-1.6656309) q[0];
sx q[0];
rz(-2.4767843) q[0];
x q[1];
rz(-0.60549824) q[2];
sx q[2];
rz(-1.7586305) q[2];
sx q[2];
rz(2.2772636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6746241) q[1];
sx q[1];
rz(-1.6858988) q[1];
sx q[1];
rz(0.50927598) q[1];
x q[2];
rz(-1.4441277) q[3];
sx q[3];
rz(-1.2614246) q[3];
sx q[3];
rz(0.29361967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7287207) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.9000165) q[2];
rz(0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-2.2976105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.64549696) q[0];
sx q[0];
rz(-0.050300751) q[0];
sx q[0];
rz(2.2182933) q[0];
rz(-2.3410489) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(0.18005767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31949735) q[0];
sx q[0];
rz(-1.332177) q[0];
sx q[0];
rz(3.0988295) q[0];
x q[1];
rz(-3.0434612) q[2];
sx q[2];
rz(-1.6509102) q[2];
sx q[2];
rz(1.1289885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6517366) q[1];
sx q[1];
rz(-1.8102614) q[1];
sx q[1];
rz(0.36291562) q[1];
rz(-3.0308444) q[3];
sx q[3];
rz(-1.9344653) q[3];
sx q[3];
rz(-1.8520385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8397612) q[2];
sx q[2];
rz(-1.3000891) q[2];
sx q[2];
rz(1.7054455) q[2];
rz(-1.885421) q[3];
sx q[3];
rz(-1.6585645) q[3];
sx q[3];
rz(-0.095631599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65211463) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(-2.6982464) q[0];
rz(-1.1706785) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(-0.43177691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0042808) q[0];
sx q[0];
rz(-2.3385424) q[0];
sx q[0];
rz(-2.7025168) q[0];
rz(-pi) q[1];
x q[1];
rz(1.742892) q[2];
sx q[2];
rz(-1.6531214) q[2];
sx q[2];
rz(2.4351127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0178294) q[1];
sx q[1];
rz(-2.4430954) q[1];
sx q[1];
rz(-2.5269954) q[1];
x q[2];
rz(1.3408991) q[3];
sx q[3];
rz(-1.8935003) q[3];
sx q[3];
rz(1.6663807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7410437) q[2];
sx q[2];
rz(-0.86047518) q[2];
sx q[2];
rz(-1.7575556) q[2];
rz(-2.4436229) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(-2.82011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8292002) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(0.37305748) q[0];
rz(0.27698764) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(-0.75609797) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7519203) q[0];
sx q[0];
rz(-2.396466) q[0];
sx q[0];
rz(-0.57655277) q[0];
rz(-1.937708) q[2];
sx q[2];
rz(-1.9607301) q[2];
sx q[2];
rz(0.45118794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3792808) q[1];
sx q[1];
rz(-1.9729776) q[1];
sx q[1];
rz(2.1625175) q[1];
x q[2];
rz(-0.96090286) q[3];
sx q[3];
rz(-2.641819) q[3];
sx q[3];
rz(-2.5148396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21195352) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(0.92145222) q[2];
rz(-3.0742505) q[3];
sx q[3];
rz(-1.208297) q[3];
sx q[3];
rz(1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9901554) q[0];
sx q[0];
rz(-2.86148) q[0];
sx q[0];
rz(3.0113599) q[0];
rz(-0.79365927) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(0.28009716) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6009499) q[0];
sx q[0];
rz(-1.6495935) q[0];
sx q[0];
rz(-0.053683563) q[0];
rz(-0.21799223) q[2];
sx q[2];
rz(-1.483886) q[2];
sx q[2];
rz(-2.4119968) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6895161) q[1];
sx q[1];
rz(-1.0930287) q[1];
sx q[1];
rz(-0.46124752) q[1];
rz(-pi) q[2];
rz(2.04591) q[3];
sx q[3];
rz(-2.2244448) q[3];
sx q[3];
rz(2.978289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0555931) q[2];
sx q[2];
rz(-1.4693825) q[2];
sx q[2];
rz(0.80297339) q[2];
rz(-1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-2.3101961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(0.10920864) q[0];
rz(-2.7681007) q[1];
sx q[1];
rz(-1.1871352) q[1];
sx q[1];
rz(2.5901897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38574252) q[0];
sx q[0];
rz(-2.0944716) q[0];
sx q[0];
rz(0.68501212) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18420561) q[2];
sx q[2];
rz(-1.4478683) q[2];
sx q[2];
rz(-0.21609989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0950737) q[1];
sx q[1];
rz(-2.1035026) q[1];
sx q[1];
rz(-2.2979876) q[1];
rz(-pi) q[2];
x q[2];
rz(2.738036) q[3];
sx q[3];
rz(-1.3440895) q[3];
sx q[3];
rz(1.8098531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52770829) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(-1.8180397) q[2];
rz(1.8771133) q[3];
sx q[3];
rz(-1.2885965) q[3];
sx q[3];
rz(-3.1356139) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5453813) q[0];
sx q[0];
rz(-0.63544202) q[0];
sx q[0];
rz(-2.4052461) q[0];
rz(-2.9296854) q[1];
sx q[1];
rz(-2.1129463) q[1];
sx q[1];
rz(-1.597499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1101226) q[0];
sx q[0];
rz(-1.2916864) q[0];
sx q[0];
rz(-1.5895248) q[0];
rz(-pi) q[1];
rz(-3.1112291) q[2];
sx q[2];
rz(-1.5734451) q[2];
sx q[2];
rz(0.16636682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36240931) q[1];
sx q[1];
rz(-1.1115555) q[1];
sx q[1];
rz(2.4613441) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61640255) q[3];
sx q[3];
rz(-2.6129301) q[3];
sx q[3];
rz(-1.5746436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.836901) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(-1.8677853) q[2];
rz(-1.4443719) q[3];
sx q[3];
rz(-3.0675409) q[3];
sx q[3];
rz(1.6465638) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.973751) q[0];
sx q[0];
rz(-1.5580307) q[0];
sx q[0];
rz(1.8488484) q[0];
rz(1.6043067) q[1];
sx q[1];
rz(-2.2289386) q[1];
sx q[1];
rz(-2.9569721) q[1];
rz(0.040097728) q[2];
sx q[2];
rz(-1.5798777) q[2];
sx q[2];
rz(-0.29691534) q[2];
rz(1.3956232) q[3];
sx q[3];
rz(-2.5268231) q[3];
sx q[3];
rz(-2.3781621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
