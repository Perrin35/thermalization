OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(1.8523822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754817) q[0];
sx q[0];
rz(-2.5876343) q[0];
sx q[0];
rz(1.0004811) q[0];
rz(-pi) q[1];
rz(0.024359811) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(2.9565405) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42613712) q[1];
sx q[1];
rz(-1.8278367) q[1];
sx q[1];
rz(-2.9261158) q[1];
rz(-pi) q[2];
rz(-1.3025769) q[3];
sx q[3];
rz(-1.5032288) q[3];
sx q[3];
rz(-0.64642954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(-2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(-0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1056054) q[0];
sx q[0];
rz(-0.30788883) q[0];
sx q[0];
rz(0.76460989) q[0];
x q[1];
rz(1.08554) q[2];
sx q[2];
rz(-1.8453571) q[2];
sx q[2];
rz(-2.1833414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(-2.980742) q[1];
rz(-pi) q[2];
x q[2];
rz(0.075606451) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057782877) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(0.16608873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27667339) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(-2.9898781) q[0];
rz(2.5327352) q[2];
sx q[2];
rz(-1.9945842) q[2];
sx q[2];
rz(0.54178967) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7158311) q[1];
sx q[1];
rz(-0.67486963) q[1];
sx q[1];
rz(2.3900044) q[1];
x q[2];
rz(-0.2563128) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(-3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-0.44149533) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(0.77484432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.144865) q[0];
sx q[0];
rz(-1.554368) q[0];
sx q[0];
rz(-0.073547151) q[0];
x q[1];
rz(1.5261126) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.1553264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6619819) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(1.9740723) q[1];
rz(-pi) q[2];
rz(-0.99153783) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(3.0338874) q[1];
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
rz(0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.9212854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6458466) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-2.2699725) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0884573) q[3];
sx q[3];
rz(-2.5128551) q[3];
sx q[3];
rz(2.351458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(0.010820476) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(-1.0169792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65360083) q[0];
sx q[0];
rz(-1.9032318) q[0];
sx q[0];
rz(1.0792653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20118841) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(-2.3692998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6974666) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(2.6078348) q[1];
rz(-pi) q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(-0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-0.84386688) q[0];
rz(1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780753) q[0];
sx q[0];
rz(-0.099637195) q[0];
sx q[0];
rz(0.3936605) q[0];
rz(-2.8104066) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(-1.0489724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0691094) q[1];
sx q[1];
rz(-1.8169889) q[1];
sx q[1];
rz(-0.90087955) q[1];
x q[2];
rz(1.3656093) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(2.6147571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(0.037242446) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9650139) q[0];
sx q[0];
rz(-1.8518378) q[0];
sx q[0];
rz(2.2289508) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2572844) q[2];
sx q[2];
rz(-1.0672617) q[2];
sx q[2];
rz(-1.7538278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.65539) q[1];
sx q[1];
rz(-1.4185925) q[1];
sx q[1];
rz(2.5712719) q[1];
rz(2.7301634) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(0.065447741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(0.55465737) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2239969) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976534) q[0];
sx q[0];
rz(-0.88730747) q[0];
sx q[0];
rz(-0.20792122) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6998859) q[2];
sx q[2];
rz(-2.522905) q[2];
sx q[2];
rz(2.8512851) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0548045) q[1];
sx q[1];
rz(-1.0541704) q[1];
sx q[1];
rz(0.36887849) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5040594) q[3];
sx q[3];
rz(-0.70422322) q[3];
sx q[3];
rz(0.080554068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(2.7049086) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.3814829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.6185332) q[0];
rz(0.47423116) q[2];
sx q[2];
rz(-1.6696764) q[2];
sx q[2];
rz(2.729051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(-0.53955697) q[1];
rz(-pi) q[2];
rz(-2.284427) q[3];
sx q[3];
rz(-1.3514337) q[3];
sx q[3];
rz(0.21881783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(-2.6560442) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.3960719) q[2];
sx q[2];
rz(-2.8833485) q[2];
sx q[2];
rz(-1.9042518) q[2];
rz(0.14470312) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
