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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(-2.562685) q[1];
sx q[1];
rz(0.97715598) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8492964) q[0];
sx q[0];
rz(-1.9417282) q[0];
sx q[0];
rz(-1.7975259) q[0];
x q[1];
rz(-1.5571655) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(1.4135828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.4358913) q[1];
sx q[1];
rz(-0.6026661) q[1];
x q[2];
rz(2.7746088) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(0.56753502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5620293) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(2.5207632) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(2.5625693) q[0];
rz(2.31965) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-0.45375219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762208) q[0];
sx q[0];
rz(-2.1848618) q[0];
sx q[0];
rz(1.8472415) q[0];
x q[1];
rz(1.683382) q[2];
sx q[2];
rz(-2.0346918) q[2];
sx q[2];
rz(0.91700208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1316741) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(-2.591955) q[1];
rz(-pi) q[2];
rz(-2.1165127) q[3];
sx q[3];
rz(-1.1835872) q[3];
sx q[3];
rz(2.6301165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.73432505) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(0.56488758) q[2];
rz(-0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(0.55150223) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(-2.636748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0612978) q[0];
sx q[0];
rz(-0.43879959) q[0];
sx q[0];
rz(-0.3831692) q[0];
rz(-pi) q[1];
rz(-1.6644434) q[2];
sx q[2];
rz(-1.159707) q[2];
sx q[2];
rz(1.9699485) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4703731) q[1];
sx q[1];
rz(-2.3335636) q[1];
sx q[1];
rz(0.68617448) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48074333) q[3];
sx q[3];
rz(-1.9078443) q[3];
sx q[3];
rz(-0.8686921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(3.0943387) q[2];
rz(-0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(-1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4243917) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(1.3264054) q[0];
rz(-1.6560417) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43760532) q[0];
sx q[0];
rz(-3.1235031) q[0];
sx q[0];
rz(1.6419069) q[0];
rz(-pi) q[1];
rz(1.0718173) q[2];
sx q[2];
rz(-1.8868251) q[2];
sx q[2];
rz(1.8763148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.078770854) q[1];
sx q[1];
rz(-2.7981542) q[1];
sx q[1];
rz(1.7867286) q[1];
x q[2];
rz(0.62795378) q[3];
sx q[3];
rz(-1.1276922) q[3];
sx q[3];
rz(-1.7114044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.7690313) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(0.27455583) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127301) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(0.10511705) q[0];
rz(-2.8016727) q[1];
sx q[1];
rz(-0.9143908) q[1];
sx q[1];
rz(-1.0629268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2539627) q[0];
sx q[0];
rz(-0.85940532) q[0];
sx q[0];
rz(-1.5894058) q[0];
x q[1];
rz(1.3219484) q[2];
sx q[2];
rz(-2.3169059) q[2];
sx q[2];
rz(0.4363554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8430994) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(1.5850409) q[1];
rz(-pi) q[2];
rz(1.5664728) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(-2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.8049392) q[2];
rz(3.0873599) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(1.0460151) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(-0.85375839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7542242) q[0];
sx q[0];
rz(-1.4367391) q[0];
sx q[0];
rz(-0.06019528) q[0];
rz(-2.1519203) q[2];
sx q[2];
rz(-2.2377491) q[2];
sx q[2];
rz(-3.0158531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5484838) q[1];
sx q[1];
rz(-2.1442772) q[1];
sx q[1];
rz(-0.67326633) q[1];
rz(-pi) q[2];
rz(2.7925909) q[3];
sx q[3];
rz(-1.4274538) q[3];
sx q[3];
rz(-2.3239612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(-0.74756527) q[2];
rz(-0.93303624) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-0.30195495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(2.5001496) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(2.0279121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9819558) q[0];
sx q[0];
rz(-0.41712077) q[0];
sx q[0];
rz(1.812716) q[0];
rz(-pi) q[1];
rz(0.895787) q[2];
sx q[2];
rz(-0.75715827) q[2];
sx q[2];
rz(-0.96162187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8889474) q[1];
sx q[1];
rz(-1.3919437) q[1];
sx q[1];
rz(2.0986544) q[1];
rz(2.6212531) q[3];
sx q[3];
rz(-1.6005777) q[3];
sx q[3];
rz(-2.7958718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-0.75801545) q[2];
rz(2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4161943) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(0.92195177) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(3.0070378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050722402) q[0];
sx q[0];
rz(-2.7742552) q[0];
sx q[0];
rz(1.814117) q[0];
rz(2.806611) q[2];
sx q[2];
rz(-0.6689531) q[2];
sx q[2];
rz(1.9425336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1052195) q[1];
sx q[1];
rz(-1.4810307) q[1];
sx q[1];
rz(-1.5516993) q[1];
x q[2];
rz(-2.4681925) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(-1.0750225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(-1.1456683) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(1.0660508) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.3409748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0308834) q[0];
sx q[0];
rz(-1.481505) q[0];
sx q[0];
rz(-0.16020328) q[0];
x q[1];
rz(2.3843896) q[2];
sx q[2];
rz(-1.5067889) q[2];
sx q[2];
rz(-0.369095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.080440259) q[1];
sx q[1];
rz(-2.3576479) q[1];
sx q[1];
rz(0.32439167) q[1];
rz(-pi) q[2];
rz(-1.789647) q[3];
sx q[3];
rz(-2.7776981) q[3];
sx q[3];
rz(-0.6536676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0880903) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(3.0255393) q[2];
rz(-1.8807489) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(-1.457823) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55105227) q[0];
sx q[0];
rz(-0.11369471) q[0];
sx q[0];
rz(0.1846479) q[0];
rz(0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1383863) q[0];
sx q[0];
rz(-1.0995563) q[0];
sx q[0];
rz(-0.22881656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72214224) q[2];
sx q[2];
rz(-1.5173679) q[2];
sx q[2];
rz(-2.4268399) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25333693) q[1];
sx q[1];
rz(-1.3789021) q[1];
sx q[1];
rz(-0.95805053) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2380794) q[3];
sx q[3];
rz(-2.2196688) q[3];
sx q[3];
rz(1.0024662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47716466) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(-2.7208929) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(1.8506321) q[2];
sx q[2];
rz(-1.6879514) q[2];
sx q[2];
rz(0.25456706) q[2];
rz(1.609997) q[3];
sx q[3];
rz(-1.3269674) q[3];
sx q[3];
rz(1.5404601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
