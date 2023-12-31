OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(2.826638) q[1];
sx q[1];
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0955015) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(3.1216937) q[0];
rz(-2.7675682) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(-0.49486578) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.284277) q[1];
sx q[1];
rz(-0.89968649) q[1];
sx q[1];
rz(-1.7637232) q[1];
rz(-pi) q[2];
rz(2.8217444) q[3];
sx q[3];
rz(-1.5187129) q[3];
sx q[3];
rz(-2.6944514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.989495) q[0];
rz(1.2377897) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-2.6706085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7171213) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(0.061901285) q[0];
x q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(1.0242467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86733782) q[1];
sx q[1];
rz(-1.5161637) q[1];
sx q[1];
rz(-2.0591303) q[1];
rz(0.69555517) q[3];
sx q[3];
rz(-1.5091981) q[3];
sx q[3];
rz(-1.0775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0835138) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96734756) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(1.9232737) q[0];
rz(0.36053948) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(0.53256065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5369536) q[1];
sx q[1];
rz(-0.96452689) q[1];
sx q[1];
rz(-2.9301675) q[1];
rz(2.2112591) q[3];
sx q[3];
rz(-2.4822682) q[3];
sx q[3];
rz(1.3673155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(-2.9411194) q[2];
rz(-0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(2.2669852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(1.8391795) q[0];
rz(0.34747296) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(0.70870542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.306327) q[1];
sx q[1];
rz(-1.7976465) q[1];
sx q[1];
rz(-0.38362417) q[1];
rz(1.6455669) q[3];
sx q[3];
rz(-1.4660335) q[3];
sx q[3];
rz(-1.9024224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.4495513) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5071881) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(0.27642823) q[0];
rz(-pi) q[1];
rz(0.50136106) q[2];
sx q[2];
rz(-1.4793581) q[2];
sx q[2];
rz(-0.83133343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44504657) q[1];
sx q[1];
rz(-1.3475218) q[1];
sx q[1];
rz(-2.9162507) q[1];
rz(1.284243) q[3];
sx q[3];
rz(-0.78714579) q[3];
sx q[3];
rz(-0.033165008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(0.26838475) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6049833) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(1.8761294) q[0];
x q[1];
rz(2.2141586) q[2];
sx q[2];
rz(-0.76403996) q[2];
sx q[2];
rz(-1.6709136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.222059) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(-0.49948378) q[1];
x q[2];
rz(2.3985732) q[3];
sx q[3];
rz(-1.7457186) q[3];
sx q[3];
rz(2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053453) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-2.5678182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928497) q[0];
sx q[0];
rz(-0.22073711) q[0];
sx q[0];
rz(1.0462532) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7710118) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(-1.3862762) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80572762) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(-2.2373799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1366828) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(0.79808455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
rz(-2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745236) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(1.2033071) q[0];
rz(-1.2175351) q[2];
sx q[2];
rz(-1.8262987) q[2];
sx q[2];
rz(-0.080163408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1295373) q[1];
sx q[1];
rz(-2.4730198) q[1];
sx q[1];
rz(-1.2052016) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5116974) q[3];
sx q[3];
rz(-0.44145465) q[3];
sx q[3];
rz(-2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4351589) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81460641) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(0.58821046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5327685) q[1];
sx q[1];
rz(-2.3012487) q[1];
sx q[1];
rz(-0.26809147) q[1];
x q[2];
rz(-0.25837274) q[3];
sx q[3];
rz(-1.4900472) q[3];
sx q[3];
rz(-1.0563869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-0.231803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7983539) q[0];
sx q[0];
rz(-0.51734561) q[0];
sx q[0];
rz(-1.377064) q[0];
x q[1];
rz(0.24256369) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(-2.501542) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9634562) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(-1.7174277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9276804) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(-2.5122535) q[2];
sx q[2];
rz(-0.39471252) q[2];
sx q[2];
rz(-0.82804745) q[2];
rz(-2.9536392) q[3];
sx q[3];
rz(-1.5899851) q[3];
sx q[3];
rz(-0.24485484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
