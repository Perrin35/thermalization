OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.2831777) q[0];
sx q[0];
rz(-1.5812961) q[0];
sx q[0];
rz(-0.66891447) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73092587) q[0];
sx q[0];
rz(-1.9528292) q[0];
sx q[0];
rz(-2.8305603) q[0];
rz(0.046836179) q[2];
sx q[2];
rz(-1.9745262) q[2];
sx q[2];
rz(2.9138034) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60379825) q[1];
sx q[1];
rz(-1.9532353) q[1];
sx q[1];
rz(2.8549744) q[1];
rz(-pi) q[2];
rz(-1.2143308) q[3];
sx q[3];
rz(-2.2495396) q[3];
sx q[3];
rz(-1.9744622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5657438) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(-3.1080833) q[2];
rz(-1.2014028) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(1.0777773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031438436) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(2.187619) q[0];
rz(1.4454449) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.1711858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4142708) q[0];
sx q[0];
rz(-2.1571299) q[0];
sx q[0];
rz(-0.81887736) q[0];
rz(-2.7523405) q[2];
sx q[2];
rz(-0.67485038) q[2];
sx q[2];
rz(2.9794326) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5038472) q[1];
sx q[1];
rz(-0.23023573) q[1];
sx q[1];
rz(3.0158325) q[1];
x q[2];
rz(2.21255) q[3];
sx q[3];
rz(-0.52601846) q[3];
sx q[3];
rz(-2.0656916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90616068) q[2];
sx q[2];
rz(-1.4305328) q[2];
sx q[2];
rz(-1.99235) q[2];
rz(0.033128459) q[3];
sx q[3];
rz(-1.6413611) q[3];
sx q[3];
rz(-2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843435) q[0];
sx q[0];
rz(-2.3488022) q[0];
sx q[0];
rz(1.0302011) q[0];
rz(1.9016117) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(2.1485567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1112615) q[0];
sx q[0];
rz(-1.6188038) q[0];
sx q[0];
rz(-2.6859849) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4740691) q[2];
sx q[2];
rz(-1.4771531) q[2];
sx q[2];
rz(2.8053631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56089106) q[1];
sx q[1];
rz(-2.1296394) q[1];
sx q[1];
rz(-0.88878312) q[1];
rz(-0.42476023) q[3];
sx q[3];
rz(-2.3094607) q[3];
sx q[3];
rz(-0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21732907) q[2];
sx q[2];
rz(-1.1161048) q[2];
sx q[2];
rz(2.7868311) q[2];
rz(0.99700704) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(-0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.16911258) q[0];
sx q[0];
rz(-1.7153772) q[0];
sx q[0];
rz(-2.0347563) q[0];
rz(2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(-0.69721627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95283629) q[0];
sx q[0];
rz(-1.3719014) q[0];
sx q[0];
rz(-2.5752978) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12733404) q[2];
sx q[2];
rz(-2.4155354) q[2];
sx q[2];
rz(-0.51148326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8839421) q[1];
sx q[1];
rz(-1.9387987) q[1];
sx q[1];
rz(1.1543399) q[1];
x q[2];
rz(1.5955865) q[3];
sx q[3];
rz(-1.5462592) q[3];
sx q[3];
rz(-2.8061574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-0.29402688) q[2];
rz(-2.0781519) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251637) q[0];
sx q[0];
rz(-2.9603781) q[0];
sx q[0];
rz(0.46318769) q[0];
rz(0.40395346) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(2.892866) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28379074) q[0];
sx q[0];
rz(-1.3732855) q[0];
sx q[0];
rz(-0.0017314712) q[0];
rz(-1.8291446) q[2];
sx q[2];
rz(-2.097192) q[2];
sx q[2];
rz(2.0451982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.675093) q[1];
sx q[1];
rz(-2.6394301) q[1];
sx q[1];
rz(0.84169047) q[1];
rz(-pi) q[2];
rz(-2.24108) q[3];
sx q[3];
rz(-2.1888362) q[3];
sx q[3];
rz(-1.4423808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8467466) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(-1.3999636) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(0.2259026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1840709) q[0];
sx q[0];
rz(-0.83860832) q[0];
sx q[0];
rz(0.34570178) q[0];
rz(0.050994571) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(-0.7720224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8851903) q[0];
sx q[0];
rz(-1.2202383) q[0];
sx q[0];
rz(-1.3023443) q[0];
rz(-pi) q[1];
rz(-0.096944158) q[2];
sx q[2];
rz(-2.2205995) q[2];
sx q[2];
rz(-0.61496269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.0061969697) q[1];
sx q[1];
rz(-0.67402285) q[1];
sx q[1];
rz(-1.037302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0431248) q[3];
sx q[3];
rz(-1.3195795) q[3];
sx q[3];
rz(-2.4285897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(-0.15743206) q[2];
rz(2.7031247) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(-2.1854775) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081414374) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(9/(11*pi)) q[0];
rz(2.5632437) q[1];
sx q[1];
rz(-1.7702421) q[1];
sx q[1];
rz(-0.15301212) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082913) q[0];
sx q[0];
rz(-2.387945) q[0];
sx q[0];
rz(0.80259364) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4744841) q[2];
sx q[2];
rz(-1.6553214) q[2];
sx q[2];
rz(2.9241004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12100231) q[1];
sx q[1];
rz(-2.4205225) q[1];
sx q[1];
rz(0.64863689) q[1];
x q[2];
rz(2.8197391) q[3];
sx q[3];
rz(-1.1328837) q[3];
sx q[3];
rz(-1.9066325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1795307) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.3029441) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-2.0343659) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1846979) q[0];
sx q[0];
rz(-0.95320025) q[0];
sx q[0];
rz(1.4554998) q[0];
rz(-0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(0.89967322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036741055) q[0];
sx q[0];
rz(-2.0213958) q[0];
sx q[0];
rz(1.3356871) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3191965) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(1.0902001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65708651) q[1];
sx q[1];
rz(-0.58875798) q[1];
sx q[1];
rz(-1.4935054) q[1];
rz(-3.1042388) q[3];
sx q[3];
rz(-1.3539697) q[3];
sx q[3];
rz(-1.9268394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4453033) q[2];
sx q[2];
rz(-0.60089198) q[2];
sx q[2];
rz(-2.3882833) q[2];
rz(2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(-1.1118579) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488572) q[0];
sx q[0];
rz(-2.1519372) q[0];
sx q[0];
rz(0.85187546) q[0];
rz(-1.4082255) q[1];
sx q[1];
rz(-1.6586761) q[1];
sx q[1];
rz(-0.3826938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9155884) q[0];
sx q[0];
rz(-1.3289641) q[0];
sx q[0];
rz(0.42892021) q[0];
x q[1];
rz(-0.76994728) q[2];
sx q[2];
rz(-1.7367518) q[2];
sx q[2];
rz(2.6024352) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5671317) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(-1.4682795) q[1];
x q[2];
rz(-2.9858573) q[3];
sx q[3];
rz(-2.6504575) q[3];
sx q[3];
rz(2.8409426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3771628) q[2];
sx q[2];
rz(-2.1775776) q[2];
sx q[2];
rz(0.60297472) q[2];
rz(1.0682586) q[3];
sx q[3];
rz(-1.7030741) q[3];
sx q[3];
rz(-0.8409797) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076440796) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(-0.35368791) q[0];
rz(2.0091281) q[1];
sx q[1];
rz(-0.9681038) q[1];
sx q[1];
rz(2.7521334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6923461) q[0];
sx q[0];
rz(-2.5857537) q[0];
sx q[0];
rz(2.6920094) q[0];
x q[1];
rz(-1.9612938) q[2];
sx q[2];
rz(-0.26528851) q[2];
sx q[2];
rz(-0.92084322) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3130541) q[1];
sx q[1];
rz(-1.8450929) q[1];
sx q[1];
rz(2.6649339) q[1];
rz(1.6234446) q[3];
sx q[3];
rz(-1.0431759) q[3];
sx q[3];
rz(1.9260581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13122095) q[2];
sx q[2];
rz(-1.5067357) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(-1.5366588) q[3];
sx q[3];
rz(-0.48303548) q[3];
sx q[3];
rz(-0.35227942) q[3];
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
rz(-0.99240408) q[0];
sx q[0];
rz(-2.438899) q[0];
sx q[0];
rz(-2.6371523) q[0];
rz(-3.0814677) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(0.073645097) q[2];
sx q[2];
rz(-1.6785868) q[2];
sx q[2];
rz(1.6389107) q[2];
rz(-0.89636421) q[3];
sx q[3];
rz(-1.5435565) q[3];
sx q[3];
rz(0.35741318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
