OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(-1.2892105) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754817) q[0];
sx q[0];
rz(-0.55395836) q[0];
sx q[0];
rz(1.0004811) q[0];
rz(1.5775561) q[2];
sx q[2];
rz(-1.5464371) q[2];
sx q[2];
rz(1.7560132) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7154555) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(2.9261158) q[1];
x q[2];
rz(1.8208002) q[3];
sx q[3];
rz(-2.8651926) q[3];
sx q[3];
rz(1.9763415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(3.0060449) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(2.360789) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.3134726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03598729) q[0];
sx q[0];
rz(-2.8337038) q[0];
sx q[0];
rz(-2.3769828) q[0];
rz(-pi) q[1];
rz(1.08554) q[2];
sx q[2];
rz(-1.2962356) q[2];
sx q[2];
rz(-0.95825125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5116509) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(-2.6444142) q[1];
rz(-3.0659862) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(2.9755039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375219) q[0];
sx q[0];
rz(-1.4253758) q[0];
sx q[0];
rz(1.279633) q[0];
rz(-pi) q[1];
rz(0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(2.6454676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7158311) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(-0.75158822) q[1];
x q[2];
rz(1.8649678) q[3];
sx q[3];
rz(-2.4066381) q[3];
sx q[3];
rz(2.8923349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7819536) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.3004998) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0825901) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(0.77484432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9967277) q[0];
sx q[0];
rz(-1.554368) q[0];
sx q[0];
rz(-0.073547151) q[0];
x q[1];
rz(1.61548) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.9862663) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91037726) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(1.3740205) q[1];
x q[2];
rz(-0.1713486) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(-1.6397887) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90785039) q[0];
sx q[0];
rz(-1.478501) q[0];
sx q[0];
rz(-1.924563) q[0];
rz(-pi) q[1];
rz(-2.2568251) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.9212854) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6377423) q[1];
sx q[1];
rz(-0.78617326) q[1];
sx q[1];
rz(2.1450858) q[1];
rz(-2.1343469) q[3];
sx q[3];
rz(-1.2754903) q[3];
sx q[3];
rz(1.2122648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(3.1048408) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058379563) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6769584) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(-0.93924378) q[0];
rz(-pi) q[1];
x q[1];
rz(1.594627) q[2];
sx q[2];
rz(-1.3696635) q[2];
sx q[2];
rz(0.80326524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.444126) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(2.6078348) q[1];
rz(-pi) q[2];
rz(-2.6050657) q[3];
sx q[3];
rz(-2.0201207) q[3];
sx q[3];
rz(0.41538737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(0.84386688) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.2197781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0735002) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(1.5324701) q[0];
x q[1];
rz(-2.2336823) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(-0.31928911) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4525758) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(-0.31022443) q[1];
x q[2];
rz(-1.1620429) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(-0.85588928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.09181) q[0];
sx q[0];
rz(-0.70735332) q[0];
sx q[0];
rz(2.011767) q[0];
x q[1];
rz(-1.8843083) q[2];
sx q[2];
rz(-1.0672617) q[2];
sx q[2];
rz(1.7538278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1815391) q[1];
sx q[1];
rz(-1.0078733) q[1];
sx q[1];
rz(1.7510508) q[1];
x q[2];
rz(0.51729047) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.6329637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976534) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-2.9336714) q[0];
x q[1];
rz(-0.091392322) q[2];
sx q[2];
rz(-2.1835727) q[2];
sx q[2];
rz(-2.6932655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.846261) q[1];
sx q[1];
rz(-1.8896855) q[1];
sx q[1];
rz(-1.0237414) q[1];
rz(-0.056592654) q[3];
sx q[3];
rz(-0.86846272) q[3];
sx q[3];
rz(3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(-1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.3814829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488659) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(-0.086765246) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9276766) q[2];
sx q[2];
rz(-0.48366085) q[2];
sx q[2];
rz(0.96825251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(-2.6020357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8990717) q[3];
sx q[3];
rz(-0.74088135) q[3];
sx q[3];
rz(1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-0.59664574) q[2];
rz(-2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(1.8252774) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(1.4894555) q[3];
sx q[3];
rz(-1.0621536) q[3];
sx q[3];
rz(-2.5928706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];