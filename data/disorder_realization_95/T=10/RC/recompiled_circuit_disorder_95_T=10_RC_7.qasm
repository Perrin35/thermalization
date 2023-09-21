OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(2.907213) q[1];
sx q[1];
rz(-2.8657764) q[1];
sx q[1];
rz(-2.0770567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2250741) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(0.9401456) q[0];
rz(2.7877596) q[2];
sx q[2];
rz(-0.96828038) q[2];
sx q[2];
rz(-1.8941855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(-2.4982846) q[1];
x q[2];
rz(-1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-2.2251341) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40936138) q[0];
sx q[0];
rz(-1.1890113) q[0];
sx q[0];
rz(-0.9071) q[0];
rz(-pi) q[1];
rz(-2.7777113) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(1.770307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0757054) q[1];
sx q[1];
rz(-1.1527219) q[1];
sx q[1];
rz(1.651152) q[1];
x q[2];
rz(-1.8804178) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(-2.0203967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(-2.0545161) q[0];
x q[1];
rz(-1.9636743) q[2];
sx q[2];
rz(-1.5117206) q[2];
sx q[2];
rz(-1.3438366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.725863) q[1];
sx q[1];
rz(-1.2408857) q[1];
sx q[1];
rz(1.2582474) q[1];
rz(-pi) q[2];
rz(0.49384533) q[3];
sx q[3];
rz(-1.4673127) q[3];
sx q[3];
rz(2.8007357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(2.1901954) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(0.11638164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560293) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(1.5390736) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(-1.7184005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0105646) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(2.4478711) q[1];
rz(-2.651752) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.9929569) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0557077) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(-1.4447312) q[0];
x q[1];
rz(-1.0257452) q[2];
sx q[2];
rz(-0.75136853) q[2];
sx q[2];
rz(0.92912208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7733113) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.7721304) q[1];
rz(-0.43907459) q[3];
sx q[3];
rz(-0.77507654) q[3];
sx q[3];
rz(0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8800031) q[0];
sx q[0];
rz(-2.3149009) q[0];
sx q[0];
rz(1.4522626) q[0];
x q[1];
rz(-0.88271898) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(-2.8628652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8099765) q[1];
sx q[1];
rz(-2.9990701) q[1];
sx q[1];
rz(-2.317821) q[1];
rz(0.19595512) q[3];
sx q[3];
rz(-1.8926419) q[3];
sx q[3];
rz(-0.67924196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241087) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(0.80362513) q[0];
rz(1.4326101) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(1.1952343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8163029) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(-1.7273278) q[1];
x q[2];
rz(-2.1738449) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655501) q[0];
sx q[0];
rz(-1.3001469) q[0];
sx q[0];
rz(-2.3483777) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4079354) q[2];
sx q[2];
rz(-1.7454299) q[2];
sx q[2];
rz(2.3101431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(0.68316858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82151316) q[3];
sx q[3];
rz(-1.320208) q[3];
sx q[3];
rz(-0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.6220629) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1544979) q[0];
sx q[0];
rz(-2.8821766) q[0];
sx q[0];
rz(0.82018606) q[0];
rz(-2.8081886) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(2.8508027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77482624) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(1.1714539) q[1];
rz(-pi) q[2];
rz(1.3689234) q[3];
sx q[3];
rz(-1.8608421) q[3];
sx q[3];
rz(-1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-0.15670776) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(-2.8740846) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88975554) q[0];
sx q[0];
rz(-1.4714843) q[0];
sx q[0];
rz(1.8343783) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.033038) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(-3.0548981) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1117489) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(2.7010121) q[1];
rz(2.3481028) q[3];
sx q[3];
rz(-1.3440545) q[3];
sx q[3];
rz(-1.9179118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7662979) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-2.5557774) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];