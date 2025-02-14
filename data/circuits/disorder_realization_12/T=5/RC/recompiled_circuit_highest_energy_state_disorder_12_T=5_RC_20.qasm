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
rz(2.9830018) q[0];
sx q[0];
rz(-2.6245485) q[0];
sx q[0];
rz(-1.3102732) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(0.96460834) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63435455) q[0];
sx q[0];
rz(-2.1956081) q[0];
sx q[0];
rz(1.4039197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3753424) q[2];
sx q[2];
rz(-1.9103622) q[2];
sx q[2];
rz(1.5887141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55935153) q[1];
sx q[1];
rz(-0.808945) q[1];
sx q[1];
rz(1.9870583) q[1];
x q[2];
rz(-1.7157285) q[3];
sx q[3];
rz(-0.67620899) q[3];
sx q[3];
rz(-1.9761666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0059119314) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(-1.4646336) q[2];
rz(1.8320734) q[3];
sx q[3];
rz(-2.194761) q[3];
sx q[3];
rz(0.32002282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86318535) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(-1.7742668) q[0];
rz(-1.4890081) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(-2.8489825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45918834) q[0];
sx q[0];
rz(-1.6688129) q[0];
sx q[0];
rz(-2.8062264) q[0];
rz(-pi) q[1];
rz(-0.66537036) q[2];
sx q[2];
rz(-0.20466802) q[2];
sx q[2];
rz(0.39835793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.039174883) q[1];
sx q[1];
rz(-0.31512773) q[1];
sx q[1];
rz(2.4246895) q[1];
rz(-pi) q[2];
rz(-1.5047362) q[3];
sx q[3];
rz(-1.8365897) q[3];
sx q[3];
rz(3.0817666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9802398) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(0.18388595) q[2];
rz(0.45267496) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9420796) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(2.3364501) q[0];
rz(-1.0792271) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(1.8260746) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.091324) q[0];
sx q[0];
rz(-1.6299695) q[0];
sx q[0];
rz(2.4135804) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82683021) q[2];
sx q[2];
rz(-1.2203802) q[2];
sx q[2];
rz(2.0598799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6650369) q[1];
sx q[1];
rz(-1.8522693) q[1];
sx q[1];
rz(0.095515619) q[1];
x q[2];
rz(0.30267164) q[3];
sx q[3];
rz(-1.4883452) q[3];
sx q[3];
rz(-2.3668805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9341854) q[2];
sx q[2];
rz(-1.6560053) q[2];
sx q[2];
rz(2.5918813) q[2];
rz(0.011119757) q[3];
sx q[3];
rz(-2.34237) q[3];
sx q[3];
rz(1.3219249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(-2.8380561) q[0];
rz(2.983298) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(1.0096445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9859814) q[0];
sx q[0];
rz(-1.8537052) q[0];
sx q[0];
rz(-2.020218) q[0];
rz(-pi) q[1];
rz(-0.60229566) q[2];
sx q[2];
rz(-1.3118852) q[2];
sx q[2];
rz(-2.7342755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6827439) q[1];
sx q[1];
rz(-1.7284596) q[1];
sx q[1];
rz(0.47994061) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49780881) q[3];
sx q[3];
rz(-1.5890997) q[3];
sx q[3];
rz(-0.25204424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2046854) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(2.8724907) q[2];
rz(-1.4540295) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(-1.5788186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8708385) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(0.36218542) q[0];
rz(2.4902792) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(-1.5247033) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0391738) q[0];
sx q[0];
rz(-0.96511594) q[0];
sx q[0];
rz(1.7282979) q[0];
rz(-pi) q[1];
rz(-1.6848504) q[2];
sx q[2];
rz(-2.5156351) q[2];
sx q[2];
rz(-2.2618641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.406639) q[1];
sx q[1];
rz(-1.2124774) q[1];
sx q[1];
rz(-1.0124769) q[1];
rz(2.397419) q[3];
sx q[3];
rz(-1.2726882) q[3];
sx q[3];
rz(1.0013195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0321956) q[2];
sx q[2];
rz(-0.95462644) q[2];
sx q[2];
rz(-1.3738013) q[2];
rz(-0.43921709) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(-0.60905987) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0405025) q[0];
sx q[0];
rz(-0.73223615) q[0];
sx q[0];
rz(0.26443431) q[0];
rz(-2.4624372) q[1];
sx q[1];
rz(-2.2676088) q[1];
sx q[1];
rz(-0.98495475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0593582) q[0];
sx q[0];
rz(-1.2684039) q[0];
sx q[0];
rz(-1.2281657) q[0];
rz(-2.0206775) q[2];
sx q[2];
rz(-1.3123672) q[2];
sx q[2];
rz(-0.33201906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15036525) q[1];
sx q[1];
rz(-2.5179407) q[1];
sx q[1];
rz(0.6536478) q[1];
rz(2.0130751) q[3];
sx q[3];
rz(-1.7566214) q[3];
sx q[3];
rz(-0.020026576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.039845) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(0.43231535) q[2];
rz(-2.1281706) q[3];
sx q[3];
rz(-1.5347967) q[3];
sx q[3];
rz(0.54949808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6559615) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(2.0320758) q[0];
rz(2.3484777) q[1];
sx q[1];
rz(-1.7879281) q[1];
sx q[1];
rz(-1.5951593) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414559) q[0];
sx q[0];
rz(-2.1329892) q[0];
sx q[0];
rz(-0.67397042) q[0];
rz(-2.4199615) q[2];
sx q[2];
rz(-1.7296556) q[2];
sx q[2];
rz(-1.5312486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90518752) q[1];
sx q[1];
rz(-1.7710454) q[1];
sx q[1];
rz(0.88242857) q[1];
rz(-pi) q[2];
rz(-1.6768544) q[3];
sx q[3];
rz(-2.2844188) q[3];
sx q[3];
rz(0.80435637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82211295) q[2];
sx q[2];
rz(-2.959368) q[2];
sx q[2];
rz(1.3694084) q[2];
rz(0.93044126) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(-1.6377009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996138) q[0];
sx q[0];
rz(-1.1811341) q[0];
sx q[0];
rz(0.30366316) q[0];
rz(0.97967255) q[1];
sx q[1];
rz(-2.517608) q[1];
sx q[1];
rz(-2.7365275) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41202085) q[0];
sx q[0];
rz(-1.3776017) q[0];
sx q[0];
rz(-2.7576819) q[0];
rz(-pi) q[1];
rz(1.1978537) q[2];
sx q[2];
rz(-0.44633807) q[2];
sx q[2];
rz(2.9624903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77877264) q[1];
sx q[1];
rz(-1.3937794) q[1];
sx q[1];
rz(-0.57139146) q[1];
x q[2];
rz(0.20196721) q[3];
sx q[3];
rz(-1.9900366) q[3];
sx q[3];
rz(-0.48101411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8651809) q[2];
sx q[2];
rz(-0.96651912) q[2];
sx q[2];
rz(2.442404) q[2];
rz(2.3326323) q[3];
sx q[3];
rz(-1.8374846) q[3];
sx q[3];
rz(-0.32304025) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.1572105) q[0];
sx q[0];
rz(-3.004177) q[0];
rz(0.049086463) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346283) q[0];
sx q[0];
rz(-2.6457743) q[0];
sx q[0];
rz(0.056255503) q[0];
x q[1];
rz(1.5829344) q[2];
sx q[2];
rz(-1.7430787) q[2];
sx q[2];
rz(2.7207295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12700554) q[1];
sx q[1];
rz(-2.0620605) q[1];
sx q[1];
rz(0.36775695) q[1];
rz(-pi) q[2];
rz(0.93799641) q[3];
sx q[3];
rz(-0.15060234) q[3];
sx q[3];
rz(-1.4678528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55234838) q[2];
sx q[2];
rz(-1.2184703) q[2];
sx q[2];
rz(-0.98769665) q[2];
rz(-0.47932953) q[3];
sx q[3];
rz(-1.5440116) q[3];
sx q[3];
rz(-2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.669303) q[0];
sx q[0];
rz(-0.23858128) q[0];
sx q[0];
rz(3.0008089) q[0];
rz(-0.36551481) q[1];
sx q[1];
rz(-0.82217685) q[1];
sx q[1];
rz(-2.4929094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539584) q[0];
sx q[0];
rz(-1.2231069) q[0];
sx q[0];
rz(2.4566922) q[0];
rz(-pi) q[1];
rz(-0.12154433) q[2];
sx q[2];
rz(-2.3732608) q[2];
sx q[2];
rz(-2.2939081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4280871) q[1];
sx q[1];
rz(-2.1660342) q[1];
sx q[1];
rz(-1.5028605) q[1];
rz(1.3289823) q[3];
sx q[3];
rz(-1.8886654) q[3];
sx q[3];
rz(2.1235025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4623798) q[2];
sx q[2];
rz(-1.0054192) q[2];
sx q[2];
rz(-0.072754808) q[2];
rz(-0.66271979) q[3];
sx q[3];
rz(-1.3136256) q[3];
sx q[3];
rz(2.976118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0543906) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(-1.6064593) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(-1.6894111) q[2];
sx q[2];
rz(-1.7977503) q[2];
sx q[2];
rz(1.1888421) q[2];
rz(-2.2519464) q[3];
sx q[3];
rz(-1.2500661) q[3];
sx q[3];
rz(1.4483857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
