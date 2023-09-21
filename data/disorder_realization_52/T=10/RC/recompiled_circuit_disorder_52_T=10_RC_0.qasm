OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9031154) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(1.205501) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56514481) q[2];
sx q[2];
rz(-1.5752715) q[2];
sx q[2];
rz(-0.43484136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86815392) q[1];
sx q[1];
rz(-1.8167129) q[1];
sx q[1];
rz(1.6383365) q[1];
x q[2];
rz(-2.2803335) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(2.0303126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(0.65594977) q[2];
rz(1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(2.7080652) q[0];
rz(-2.9128089) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4097071) q[0];
sx q[0];
rz(-1.9623555) q[0];
sx q[0];
rz(0.10086985) q[0];
rz(-2.3184001) q[2];
sx q[2];
rz(-2.7232812) q[2];
sx q[2];
rz(0.38564607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65456395) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(-1.8064206) q[1];
x q[2];
rz(-0.98874436) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(0.83272979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(0.69765222) q[2];
rz(0.12156045) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.7720222) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-0.26161584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944377) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(-1.5887512) q[0];
rz(1.4357655) q[2];
sx q[2];
rz(-2.3927852) q[2];
sx q[2];
rz(0.43040648) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1640472) q[1];
sx q[1];
rz(-2.3512212) q[1];
sx q[1];
rz(0.17069823) q[1];
rz(2.6286078) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(2.5740636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-0.27954277) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(-1.2483695) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6682537) q[0];
sx q[0];
rz(-3.0248397) q[0];
sx q[0];
rz(-0.97572414) q[0];
x q[1];
rz(-1.7522881) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(-0.69603053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.215938) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(-1.4309806) q[1];
rz(-0.72327153) q[3];
sx q[3];
rz(-1.372882) q[3];
sx q[3];
rz(-2.7151782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(-1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(-2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-2.2391589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31061253) q[0];
sx q[0];
rz(-1.6164301) q[0];
sx q[0];
rz(2.5020585) q[0];
rz(-pi) q[1];
rz(-0.47307737) q[2];
sx q[2];
rz(-2.6396751) q[2];
sx q[2];
rz(2.9830473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24689281) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(-0.65319368) q[1];
rz(-3.1387572) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(2.1383274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3301795) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70226442) q[0];
sx q[0];
rz(-1.4685255) q[0];
sx q[0];
rz(-0.066145397) q[0];
x q[1];
rz(2.3927275) q[2];
sx q[2];
rz(-0.61486926) q[2];
sx q[2];
rz(2.2392337) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0338248) q[1];
sx q[1];
rz(-1.0948085) q[1];
sx q[1];
rz(1.0955217) q[1];
x q[2];
rz(0.73268907) q[3];
sx q[3];
rz(-1.2142039) q[3];
sx q[3];
rz(-0.69571146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(-1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-2.1988791) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38521117) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(-0.60473196) q[0];
rz(-pi) q[1];
rz(-0.96037453) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(-0.20197091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9224291) q[1];
sx q[1];
rz(-1.9814098) q[1];
sx q[1];
rz(-0.34904042) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0122635) q[3];
sx q[3];
rz(-2.5848476) q[3];
sx q[3];
rz(2.7260775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(2.5121636) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-0.28392917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2487508) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(-2.6177004) q[0];
rz(-pi) q[1];
rz(-2.5416652) q[2];
sx q[2];
rz(-1.4495965) q[2];
sx q[2];
rz(-1.8184513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64360196) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(-2.1010146) q[1];
rz(-pi) q[2];
rz(2.4817804) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(-0.82994474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0652086) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47855908) q[0];
sx q[0];
rz(-1.8543108) q[0];
sx q[0];
rz(-1.5525596) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(-1.1073081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8252392) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(-0.52492001) q[1];
x q[2];
rz(-0.38735729) q[3];
sx q[3];
rz(-0.82327561) q[3];
sx q[3];
rz(-2.2203317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4298657) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(2.4263884) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.573695) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(2.7609603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56363737) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(-1.7726583) q[0];
x q[1];
rz(-2.5194174) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(-0.13571339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3818647) q[1];
sx q[1];
rz(-1.3761531) q[1];
sx q[1];
rz(2.0774283) q[1];
x q[2];
rz(-0.59488876) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(-0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0011065817) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(-1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052335) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(-1.3409875) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(2.9410578) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];