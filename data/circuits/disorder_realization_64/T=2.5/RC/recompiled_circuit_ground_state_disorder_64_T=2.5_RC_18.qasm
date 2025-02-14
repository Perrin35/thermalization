OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(-0.10804478) q[0];
sx q[0];
rz(0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(1.4563814) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0385168) q[0];
sx q[0];
rz(-1.7060192) q[0];
sx q[0];
rz(1.9413906) q[0];
rz(-pi) q[1];
rz(-0.13338344) q[2];
sx q[2];
rz(-1.6074153) q[2];
sx q[2];
rz(-0.9549779) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3177569) q[1];
sx q[1];
rz(-1.2878622) q[1];
sx q[1];
rz(-2.1499277) q[1];
rz(-pi) q[2];
rz(-1.637804) q[3];
sx q[3];
rz(-2.5634804) q[3];
sx q[3];
rz(-2.1117977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1504537) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(-1.0452622) q[2];
rz(-2.1446877) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(-1.8256942) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5528706) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(1.7887822) q[0];
rz(-3.1006587) q[1];
sx q[1];
rz(-1.2177688) q[1];
sx q[1];
rz(-1.5997684) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7307161) q[0];
sx q[0];
rz(-1.7109032) q[0];
sx q[0];
rz(-0.1430169) q[0];
rz(-0.016288494) q[2];
sx q[2];
rz(-1.5920581) q[2];
sx q[2];
rz(2.4580815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2165228) q[1];
sx q[1];
rz(-2.1669186) q[1];
sx q[1];
rz(-3.104167) q[1];
rz(-pi) q[2];
rz(-2.1349234) q[3];
sx q[3];
rz(-2.4897277) q[3];
sx q[3];
rz(1.8396371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.505595) q[2];
sx q[2];
rz(-0.032568585) q[2];
sx q[2];
rz(2.6231498) q[2];
rz(-1.1238267) q[3];
sx q[3];
rz(-0.80945194) q[3];
sx q[3];
rz(2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951303) q[0];
sx q[0];
rz(-3.0914682) q[0];
sx q[0];
rz(-1.5396402) q[0];
rz(-2.4165972) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(2.4620788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5711229) q[0];
sx q[0];
rz(-0.75665604) q[0];
sx q[0];
rz(-2.4837912) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.021596507) q[2];
sx q[2];
rz(-0.21449247) q[2];
sx q[2];
rz(-0.15211764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.589349) q[1];
sx q[1];
rz(-1.2584983) q[1];
sx q[1];
rz(2.8683788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98285003) q[3];
sx q[3];
rz(-1.790739) q[3];
sx q[3];
rz(-1.4877121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9069549) q[2];
sx q[2];
rz(-2.8997771) q[2];
sx q[2];
rz(2.6406636) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2723715) q[0];
sx q[0];
rz(-0.048373241) q[0];
sx q[0];
rz(2.3424171) q[0];
rz(2.8436106) q[1];
sx q[1];
rz(-0.25845343) q[1];
sx q[1];
rz(-0.90202773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.018343) q[0];
sx q[0];
rz(-0.76597528) q[0];
sx q[0];
rz(2.9428456) q[0];
rz(-2.600014) q[2];
sx q[2];
rz(-1.4409833) q[2];
sx q[2];
rz(2.1101348) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9060256) q[1];
sx q[1];
rz(-2.9936571) q[1];
sx q[1];
rz(-1.6355913) q[1];
rz(-1.6482192) q[3];
sx q[3];
rz(-1.6257878) q[3];
sx q[3];
rz(-2.1907326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(-1.4704977) q[2];
rz(-2.9521613) q[3];
sx q[3];
rz(-3.093284) q[3];
sx q[3];
rz(2.3725177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.9162132) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(2.7498229) q[0];
rz(2.0562992) q[1];
sx q[1];
rz(-3.1317874) q[1];
sx q[1];
rz(0.36878961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15426668) q[0];
sx q[0];
rz(-1.0079103) q[0];
sx q[0];
rz(-1.994654) q[0];
rz(-0.26924765) q[2];
sx q[2];
rz(-1.9740065) q[2];
sx q[2];
rz(-1.1336114) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5318067) q[1];
sx q[1];
rz(-1.5641677) q[1];
sx q[1];
rz(3.1408491) q[1];
rz(-0.7097575) q[3];
sx q[3];
rz(-1.3922979) q[3];
sx q[3];
rz(-0.75542007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5607249) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(-1.3690534) q[2];
rz(3.0154058) q[3];
sx q[3];
rz(-2.7044665) q[3];
sx q[3];
rz(2.0807467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5996025) q[0];
sx q[0];
rz(-2.3805711) q[0];
sx q[0];
rz(-1.5916995) q[0];
rz(0.97269336) q[1];
sx q[1];
rz(-2.9284271) q[1];
sx q[1];
rz(2.1125643) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3297076) q[0];
sx q[0];
rz(-0.68674478) q[0];
sx q[0];
rz(-0.48573692) q[0];
rz(-pi) q[1];
rz(-1.712079) q[2];
sx q[2];
rz(-1.8064692) q[2];
sx q[2];
rz(2.6038246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29518471) q[1];
sx q[1];
rz(-1.4746397) q[1];
sx q[1];
rz(-0.00043649159) q[1];
x q[2];
rz(-0.036110445) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(0.13425628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.782393) q[2];
sx q[2];
rz(-1.7155557) q[2];
sx q[2];
rz(-2.7101809) q[2];
rz(-0.99724489) q[3];
sx q[3];
rz(-0.037171818) q[3];
sx q[3];
rz(-0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.7998841) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(-2.0836015) q[0];
rz(-2.3175088) q[1];
sx q[1];
rz(-3.1415756) q[1];
sx q[1];
rz(-0.81738671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1202755) q[0];
sx q[0];
rz(-0.45545721) q[0];
sx q[0];
rz(-2.0648455) q[0];
x q[1];
rz(2.779473) q[2];
sx q[2];
rz(-0.019500168) q[2];
sx q[2];
rz(-1.8060613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89725607) q[1];
sx q[1];
rz(-0.24459363) q[1];
sx q[1];
rz(-1.7331428) q[1];
rz(-pi) q[2];
rz(0.55900662) q[3];
sx q[3];
rz(-0.83448258) q[3];
sx q[3];
rz(-2.112039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.151513) q[2];
sx q[2];
rz(-1.238287) q[2];
sx q[2];
rz(1.6286758) q[2];
rz(-0.40142909) q[3];
sx q[3];
rz(-3.1135058) q[3];
sx q[3];
rz(-2.1872971) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7734739) q[0];
sx q[0];
rz(-3.0768657) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(0.041944567) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(-2.1379474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40239375) q[0];
sx q[0];
rz(-0.93951997) q[0];
sx q[0];
rz(1.4556134) q[0];
rz(1.233439) q[2];
sx q[2];
rz(-1.5784831) q[2];
sx q[2];
rz(2.6111662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.072362547) q[1];
sx q[1];
rz(-1.3389999) q[1];
sx q[1];
rz(0.20834558) q[1];
rz(-pi) q[2];
rz(0.65559994) q[3];
sx q[3];
rz(-1.3274945) q[3];
sx q[3];
rz(0.52025627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4369138) q[2];
sx q[2];
rz(-3.0824326) q[2];
sx q[2];
rz(1.7436279) q[2];
rz(-0.49020234) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(-1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040084664) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-0.21139938) q[0];
rz(2.0231817) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(-1.6308019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64563066) q[0];
sx q[0];
rz(-1.2617153) q[0];
sx q[0];
rz(-1.8739971) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5671179) q[2];
sx q[2];
rz(-0.7099289) q[2];
sx q[2];
rz(1.0021051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76912266) q[1];
sx q[1];
rz(-1.5536397) q[1];
sx q[1];
rz(2.9869798) q[1];
rz(-pi) q[2];
rz(-0.16638593) q[3];
sx q[3];
rz(-1.659044) q[3];
sx q[3];
rz(-1.7235173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3596892) q[2];
sx q[2];
rz(-3.1236881) q[2];
sx q[2];
rz(2.8622799) q[2];
rz(-0.8638047) q[3];
sx q[3];
rz(-0.0044718663) q[3];
sx q[3];
rz(-2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(1.3155235) q[0];
rz(2.3663991) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(-1.7528037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1895386) q[0];
sx q[0];
rz(-1.7197968) q[0];
sx q[0];
rz(2.8556264) q[0];
rz(-pi) q[1];
rz(2.7371762) q[2];
sx q[2];
rz(-1.7787654) q[2];
sx q[2];
rz(-0.31851092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1422279) q[1];
sx q[1];
rz(-1.5714297) q[1];
sx q[1];
rz(-3.1406979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1517066) q[3];
sx q[3];
rz(-1.1858318) q[3];
sx q[3];
rz(-0.46483913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3751601) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(-1.4520221) q[2];
rz(0.25894138) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(1.1567206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5158841) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(-1.6547849) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(-2.941359) q[2];
sx q[2];
rz(-1.5641134) q[2];
sx q[2];
rz(-1.3431298) q[2];
rz(-0.026080118) q[3];
sx q[3];
rz(-1.2240388) q[3];
sx q[3];
rz(0.025561533) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
