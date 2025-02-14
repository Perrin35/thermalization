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
rz(pi/2) q[2];
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
rz(1.9929041) q[2];
sx q[2];
rz(-1.4808345) q[2];
sx q[2];
rz(-2.6805756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18033548) q[1];
sx q[1];
rz(-1.1828268) q[1];
sx q[1];
rz(1.1731824) q[1];
rz(2.6463616) q[3];
sx q[3];
rz(-2.3217045) q[3];
sx q[3];
rz(2.8748729) q[3];
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
rz(-0.40526059) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.8965992) q[1];
sx q[1];
rz(2.5027221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8804358) q[0];
sx q[0];
rz(-2.4969099) q[0];
sx q[0];
rz(-3.0993483) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62221726) q[2];
sx q[2];
rz(-1.8832616) q[2];
sx q[2];
rz(2.0553596) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32989001) q[1];
sx q[1];
rz(-0.59361688) q[1];
sx q[1];
rz(-0.65887405) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3078626) q[3];
sx q[3];
rz(-1.3776738) q[3];
sx q[3];
rz(1.5964791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8581533) q[2];
sx q[2];
rz(-0.44698134) q[2];
sx q[2];
rz(2.8972674) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7636488) q[0];
sx q[0];
rz(-1.0841333) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(-1.7297277) q[1];
sx q[1];
rz(-2.6928163) q[1];
sx q[1];
rz(-1.809874) q[1];
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
rz(-pi) q[1];
rz(-0.21816605) q[2];
sx q[2];
rz(-1.2274172) q[2];
sx q[2];
rz(-0.81264466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0983372) q[1];
sx q[1];
rz(-1.6041293) q[1];
sx q[1];
rz(-3.0860677) q[1];
rz(-pi) q[2];
rz(-2.3490155) q[3];
sx q[3];
rz(-0.78781414) q[3];
sx q[3];
rz(-2.8590096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1639013) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(0.65995556) q[2];
rz(1.6948949) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834943) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(0.51620382) q[1];
sx q[1];
rz(-1.3092923) q[1];
sx q[1];
rz(-2.5925327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9375795) q[0];
sx q[0];
rz(-1.7147062) q[0];
sx q[0];
rz(1.4043441) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15971036) q[2];
sx q[2];
rz(-0.98708144) q[2];
sx q[2];
rz(-2.4734378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0564894) q[1];
sx q[1];
rz(-1.6360549) q[1];
sx q[1];
rz(2.8662221) q[1];
x q[2];
rz(0.28259773) q[3];
sx q[3];
rz(-1.188084) q[3];
sx q[3];
rz(3.0846918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0893726) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(-0.29903665) q[2];
rz(-2.5236169) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1842136) q[0];
sx q[0];
rz(-0.03802499) q[0];
sx q[0];
rz(0.46603611) q[0];
rz(-2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(2.3066511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3059495) q[0];
sx q[0];
rz(-1.5993274) q[0];
sx q[0];
rz(0.051751922) q[0];
x q[1];
rz(-0.088836121) q[2];
sx q[2];
rz(-1.6497017) q[2];
sx q[2];
rz(0.14278459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1620491) q[1];
sx q[1];
rz(-2.0503373) q[1];
sx q[1];
rz(3.0813345) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4598875) q[3];
sx q[3];
rz(-0.85452628) q[3];
sx q[3];
rz(2.0264421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20778188) q[2];
sx q[2];
rz(-2.2489579) q[2];
sx q[2];
rz(-0.031522838) q[2];
rz(-0.75782123) q[3];
sx q[3];
rz(-2.0407245) q[3];
sx q[3];
rz(-0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.031540792) q[0];
rz(-3.0517598) q[1];
sx q[1];
rz(-0.49982163) q[1];
sx q[1];
rz(-0.33581844) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851097) q[0];
sx q[0];
rz(-1.8976189) q[0];
sx q[0];
rz(-2.7669897) q[0];
rz(-pi) q[1];
rz(2.8029595) q[2];
sx q[2];
rz(-2.0969318) q[2];
sx q[2];
rz(2.3735685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0586186) q[1];
sx q[1];
rz(-1.0320283) q[1];
sx q[1];
rz(0.28317105) q[1];
rz(-1.0085513) q[3];
sx q[3];
rz(-0.52298552) q[3];
sx q[3];
rz(-0.44660911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1699367) q[2];
sx q[2];
rz(-2.0270831) q[2];
sx q[2];
rz(1.8709095) q[2];
rz(-0.058567889) q[3];
sx q[3];
rz(-1.6978369) q[3];
sx q[3];
rz(0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(2.8756323) q[0];
rz(0.82414857) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(-1.6110274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53479311) q[0];
sx q[0];
rz(-1.3086119) q[0];
sx q[0];
rz(1.8305837) q[0];
x q[1];
rz(1.5179727) q[2];
sx q[2];
rz(-2.3011156) q[2];
sx q[2];
rz(2.2234302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33324533) q[1];
sx q[1];
rz(-0.63831224) q[1];
sx q[1];
rz(-0.45175938) q[1];
rz(0.35857753) q[3];
sx q[3];
rz(-2.3365575) q[3];
sx q[3];
rz(-0.81881675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1223569) q[2];
sx q[2];
rz(-2.0255721) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(2.9979749) q[3];
sx q[3];
rz(-2.2031281) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73760167) q[0];
sx q[0];
rz(-0.25594512) q[0];
sx q[0];
rz(0.10424374) q[0];
rz(-0.741611) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(0.43824497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0377698) q[0];
sx q[0];
rz(-1.6731823) q[0];
sx q[0];
rz(-2.4679604) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9096228) q[2];
sx q[2];
rz(-2.5164587) q[2];
sx q[2];
rz(-1.8115316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.150369) q[1];
sx q[1];
rz(-2.5349778) q[1];
sx q[1];
rz(1.6305417) q[1];
rz(-0.24298553) q[3];
sx q[3];
rz(-2.5930517) q[3];
sx q[3];
rz(-2.0627956) q[3];
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
rz(0.3212277) q[2];
rz(2.1494703) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(1.4120302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(0.14800063) q[0];
rz(1.1389698) q[1];
sx q[1];
rz(-0.6593467) q[1];
sx q[1];
rz(0.99050561) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73511998) q[0];
sx q[0];
rz(-0.48790259) q[0];
sx q[0];
rz(2.9814629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1075145) q[2];
sx q[2];
rz(-0.89416891) q[2];
sx q[2];
rz(-0.42596841) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3618092) q[1];
sx q[1];
rz(-1.6020007) q[1];
sx q[1];
rz(-0.0051910944) q[1];
x q[2];
rz(0.3785636) q[3];
sx q[3];
rz(-1.1051911) q[3];
sx q[3];
rz(1.9052802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(-2.5992744) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(2.3565049) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(2.7395172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77469414) q[0];
sx q[0];
rz(-1.1326362) q[0];
sx q[0];
rz(-1.9209981) q[0];
rz(-pi) q[1];
rz(0.78530757) q[2];
sx q[2];
rz(-2.1673137) q[2];
sx q[2];
rz(-1.1889088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86568923) q[1];
sx q[1];
rz(-0.69178693) q[1];
sx q[1];
rz(1.8744286) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98100752) q[3];
sx q[3];
rz(-2.6569603) q[3];
sx q[3];
rz(-2.4576602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98099199) q[2];
sx q[2];
rz(-0.6610142) q[2];
sx q[2];
rz(2.215019) q[2];
rz(0.55046925) q[3];
sx q[3];
rz(-2.2678943) q[3];
sx q[3];
rz(1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(0.015197531) q[2];
sx q[2];
rz(-0.79833818) q[2];
sx q[2];
rz(-1.0865059) q[2];
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
