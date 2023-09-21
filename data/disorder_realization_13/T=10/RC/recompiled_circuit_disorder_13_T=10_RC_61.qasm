OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(-0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0955015) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(0.019898947) q[0];
x q[1];
rz(2.1452791) q[2];
sx q[2];
rz(-2.516054) q[2];
sx q[2];
rz(-1.1686981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1617042) q[1];
sx q[1];
rz(-0.6941422) q[1];
sx q[1];
rz(2.90467) q[1];
x q[2];
rz(1.625657) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(-2.0351792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(0.89481568) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.5367616) q[1];
sx q[1];
rz(2.6706085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-2.6563245) q[0];
sx q[0];
rz(1.6886061) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2330301) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(-1.0242467) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5405611) q[1];
sx q[1];
rz(-2.6504576) q[1];
sx q[1];
rz(-1.4547552) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095951565) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-2.5747091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3210179) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(0.029650173) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5191684) q[2];
sx q[2];
rz(-1.9308959) q[2];
sx q[2];
rz(1.0200295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9860552) q[1];
sx q[1];
rz(-1.3974766) q[1];
sx q[1];
rz(-0.9539414) q[1];
rz(2.2112591) q[3];
sx q[3];
rz(-2.4822682) q[3];
sx q[3];
rz(-1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(2.2669852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0229605) q[0];
sx q[0];
rz(-1.7264139) q[0];
sx q[0];
rz(-0.96373425) q[0];
rz(-1.4250408) q[2];
sx q[2];
rz(-1.9521904) q[2];
sx q[2];
rz(0.33304735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.306327) q[1];
sx q[1];
rz(-1.7976465) q[1];
sx q[1];
rz(-0.38362417) q[1];
rz(-pi) q[2];
rz(-1.6455669) q[3];
sx q[3];
rz(-1.4660335) q[3];
sx q[3];
rz(-1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(-0.57058913) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76737228) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(2.2348997) q[0];
rz(2.9530753) q[2];
sx q[2];
rz(-2.6326615) q[2];
sx q[2];
rz(0.9045507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.78391) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(0.79343474) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3371313) q[3];
sx q[3];
rz(-1.3692229) q[3];
sx q[3];
rz(1.3988914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1935929) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282288) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(-0.22426228) q[0];
x q[1];
rz(-2.6199117) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(0.66643836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9195337) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(-0.49948378) q[1];
rz(-2.3985732) q[3];
sx q[3];
rz(-1.7457186) q[3];
sx q[3];
rz(0.35374853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(0.57377446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550802) q[0];
sx q[0];
rz(-1.3801314) q[0];
sx q[0];
rz(-3.029682) q[0];
rz(-1.7710118) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(-1.3862762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7635203) q[1];
sx q[1];
rz(-2.1280648) q[1];
sx q[1];
rz(2.4850363) q[1];
rz(-2.1366828) q[3];
sx q[3];
rz(-0.27540576) q[3];
sx q[3];
rz(-2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2670691) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(-1.9382856) q[0];
rz(0.92408085) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(0.88976394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8746652) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(0.93519559) q[1];
x q[2];
rz(-1.6298953) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(2.303404) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7064338) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(0.81460641) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(2.5533822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5327685) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(0.26809147) q[1];
rz(-pi) q[2];
rz(1.4872876) q[3];
sx q[3];
rz(-1.3132846) q[3];
sx q[3];
rz(2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3964513) q[0];
sx q[0];
rz(-1.6661577) q[0];
sx q[0];
rz(2.0800637) q[0];
rz(-pi) q[1];
rz(-0.78328697) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(2.96539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.787549) q[1];
sx q[1];
rz(-1.7122867) q[1];
sx q[1];
rz(-0.267412) q[1];
rz(-pi) q[2];
rz(0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-2.3378519) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-2.172773) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(0.3248365) q[2];
sx q[2];
rz(-1.3424716) q[2];
sx q[2];
rz(-2.9906103) q[2];
rz(-1.5512636) q[3];
sx q[3];
rz(-1.3828779) q[3];
sx q[3];
rz(-1.8120017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
