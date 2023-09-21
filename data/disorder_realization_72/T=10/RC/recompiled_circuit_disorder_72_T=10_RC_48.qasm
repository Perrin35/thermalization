OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2274269) q[0];
sx q[0];
rz(-1.5072855) q[0];
sx q[0];
rz(0.27557296) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9688641) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(-1.2781065) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0807304) q[1];
sx q[1];
rz(-0.84134358) q[1];
sx q[1];
rz(-1.3688449) q[1];
rz(-pi) q[2];
rz(1.6742168) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(-2.8570989) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53045814) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(-0.64972933) q[0];
rz(1.3833212) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(1.0152917) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36675378) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(1.2971836) q[1];
rz(-1.2172132) q[3];
sx q[3];
rz(-1.0014605) q[3];
sx q[3];
rz(-0.47488892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4602063) q[0];
sx q[0];
rz(-1.9708512) q[0];
sx q[0];
rz(-1.7494739) q[0];
rz(-pi) q[1];
rz(-2.8475464) q[2];
sx q[2];
rz(-1.1883192) q[2];
sx q[2];
rz(1.8665714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(-1.6509389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4530573) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(1.1586787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-0.5501737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.600425) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(3.0279972) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30087154) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(3.0645264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.354047) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(-1.3613308) q[1];
rz(1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(-2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-2.8660529) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(0.91526389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1137431) q[0];
sx q[0];
rz(-2.5115974) q[0];
sx q[0];
rz(-1.4907452) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0643164) q[2];
sx q[2];
rz(-2.6151267) q[2];
sx q[2];
rz(-2.4174945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.535714) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(-1.4153626) q[1];
rz(0.79742934) q[3];
sx q[3];
rz(-2.1505822) q[3];
sx q[3];
rz(-0.42339719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2561589) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(-2.8594349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5846892) q[2];
sx q[2];
rz(-2.5957426) q[2];
sx q[2];
rz(0.57005537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10627667) q[1];
sx q[1];
rz(-0.66674399) q[1];
sx q[1];
rz(-1.7229401) q[1];
rz(-1.8144238) q[3];
sx q[3];
rz(-2.1414087) q[3];
sx q[3];
rz(2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(-1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92857498) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(1.4777044) q[0];
rz(2.0790714) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(-1.6342271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(1.7794442) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0464826) q[3];
sx q[3];
rz(-1.2246119) q[3];
sx q[3];
rz(0.1971052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(-2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93452867) q[0];
sx q[0];
rz(-1.3158187) q[0];
sx q[0];
rz(1.2364788) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92808) q[2];
sx q[2];
rz(-2.4734481) q[2];
sx q[2];
rz(-1.1754787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0667463) q[1];
sx q[1];
rz(-0.42145887) q[1];
sx q[1];
rz(-1.4261817) q[1];
rz(1.950374) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(-0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582726) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
rz(-2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34708729) q[1];
sx q[1];
rz(-0.80650389) q[1];
sx q[1];
rz(-1.1361213) q[1];
x q[2];
rz(1.1205348) q[3];
sx q[3];
rz(-1.0362719) q[3];
sx q[3];
rz(-0.083732097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(-2.5433345) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82829581) q[0];
sx q[0];
rz(-1.7772355) q[0];
sx q[0];
rz(-2.9199827) q[0];
rz(-pi) q[1];
rz(-2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-3.0002909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2559402) q[1];
sx q[1];
rz(-1.2068032) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
rz(-0.15828295) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-0.61959735) q[2];
sx q[2];
rz(-1.4705114) q[2];
sx q[2];
rz(0.32497139) q[2];
rz(2.7491309) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
