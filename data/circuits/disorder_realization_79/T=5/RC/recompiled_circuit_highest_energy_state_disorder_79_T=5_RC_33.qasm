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
rz(1.5647178) q[0];
sx q[0];
rz(-1.325542) q[0];
sx q[0];
rz(2.5817885) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(-1.3275194) q[1];
sx q[1];
rz(-1.3865857) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35980269) q[0];
sx q[0];
rz(-1.823678) q[0];
sx q[0];
rz(-1.0183615) q[0];
rz(-2.785478) q[2];
sx q[2];
rz(-0.31415161) q[2];
sx q[2];
rz(-0.34700307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69405729) q[1];
sx q[1];
rz(-2.9172998) q[1];
sx q[1];
rz(-2.7560225) q[1];
rz(1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(2.5293722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5379415) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(2.3407827) q[3];
sx q[3];
rz(-1.2468612) q[3];
sx q[3];
rz(-2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2342728) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(0.32670879) q[0];
rz(-1.7960583) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(2.2950642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162982) q[0];
sx q[0];
rz(-2.7413053) q[0];
sx q[0];
rz(-1.8497224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88794215) q[2];
sx q[2];
rz(-1.789621) q[2];
sx q[2];
rz(-1.8162948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7134756) q[1];
sx q[1];
rz(-1.5813811) q[1];
sx q[1];
rz(-1.0416563) q[1];
x q[2];
rz(1.8655103) q[3];
sx q[3];
rz(-0.8407774) q[3];
sx q[3];
rz(-1.4286014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9091866) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(2.6573913) q[2];
rz(1.3062612) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(-1.9561214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4803798) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(-2.6015306) q[0];
rz(1.1506608) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(2.5475492) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10162486) q[0];
sx q[0];
rz(-1.4007845) q[0];
sx q[0];
rz(0.11084686) q[0];
x q[1];
rz(2.743529) q[2];
sx q[2];
rz(-2.4467565) q[2];
sx q[2];
rz(-0.24206012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6642528) q[1];
sx q[1];
rz(-2.3788733) q[1];
sx q[1];
rz(-1.699317) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.698231) q[3];
sx q[3];
rz(-0.48526796) q[3];
sx q[3];
rz(-1.3290153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8192886) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(0.6907531) q[2];
rz(-2.6639719) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(2.7121108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531042) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(2.389287) q[0];
rz(2.2350156) q[1];
sx q[1];
rz(-2.7598359) q[1];
sx q[1];
rz(-0.22901542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9045712) q[0];
sx q[0];
rz(-1.2777877) q[0];
sx q[0];
rz(-1.3364393) q[0];
rz(-1.3439517) q[2];
sx q[2];
rz(-1.0431492) q[2];
sx q[2];
rz(-0.91560748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.467207) q[1];
sx q[1];
rz(-1.8388646) q[1];
sx q[1];
rz(-1.8169136) q[1];
rz(-pi) q[2];
rz(-0.58521478) q[3];
sx q[3];
rz(-1.8793595) q[3];
sx q[3];
rz(1.2206961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(-1.0329186) q[2];
rz(-1.0659263) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(1.8377931) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8835939) q[0];
sx q[0];
rz(-3.0553387) q[0];
sx q[0];
rz(-1.8121383) q[0];
rz(-0.52779245) q[1];
sx q[1];
rz(-1.8254435) q[1];
sx q[1];
rz(-2.7425308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2891727) q[0];
sx q[0];
rz(-1.6483083) q[0];
sx q[0];
rz(-1.2184675) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9827323) q[2];
sx q[2];
rz(-2.2592681) q[2];
sx q[2];
rz(2.7044123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4265824) q[1];
sx q[1];
rz(-1.2599328) q[1];
sx q[1];
rz(2.3764627) q[1];
rz(-0.19029096) q[3];
sx q[3];
rz(-2.2566593) q[3];
sx q[3];
rz(-3.0113285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1479147) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(-2.1368775) q[2];
rz(-2.8435454) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(2.0731488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3598969) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(1.1915278) q[0];
rz(-2.6021992) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(-2.3753812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.51821) q[0];
sx q[0];
rz(-2.1124387) q[0];
sx q[0];
rz(0.40671273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4922392) q[2];
sx q[2];
rz(-2.1198065) q[2];
sx q[2];
rz(1.8446814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2431406) q[1];
sx q[1];
rz(-0.86075538) q[1];
sx q[1];
rz(0.30177994) q[1];
x q[2];
rz(-1.3260854) q[3];
sx q[3];
rz(-0.38077106) q[3];
sx q[3];
rz(-2.0920193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3300276) q[2];
sx q[2];
rz(-1.6946239) q[2];
sx q[2];
rz(-0.88328973) q[2];
rz(-2.6106994) q[3];
sx q[3];
rz(-0.58404946) q[3];
sx q[3];
rz(2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4009878) q[0];
sx q[0];
rz(-0.29380909) q[0];
sx q[0];
rz(-0.64939943) q[0];
rz(-0.3736639) q[1];
sx q[1];
rz(-2.3039736) q[1];
sx q[1];
rz(1.1423473) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4209265) q[0];
sx q[0];
rz(-2.775179) q[0];
sx q[0];
rz(-1.1144251) q[0];
rz(-pi) q[1];
rz(0.20988237) q[2];
sx q[2];
rz(-1.58605) q[2];
sx q[2];
rz(-1.5535056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2616407) q[1];
sx q[1];
rz(-1.8790207) q[1];
sx q[1];
rz(2.37642) q[1];
rz(-pi) q[2];
rz(-3.0965831) q[3];
sx q[3];
rz(-1.112072) q[3];
sx q[3];
rz(0.59797244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1335527) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(-0.35888654) q[2];
rz(0.60838962) q[3];
sx q[3];
rz(-1.5549992) q[3];
sx q[3];
rz(-2.2801094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95978874) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(0.96610075) q[0];
rz(2.6995662) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(0.40837049) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7743384) q[0];
sx q[0];
rz(-1.5975195) q[0];
sx q[0];
rz(1.4382118) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28200217) q[2];
sx q[2];
rz(-0.15371727) q[2];
sx q[2];
rz(1.745818) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30701559) q[1];
sx q[1];
rz(-2.53741) q[1];
sx q[1];
rz(1.3175563) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73159742) q[3];
sx q[3];
rz(-1.7642191) q[3];
sx q[3];
rz(3.132427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8926706) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(1.9425617) q[2];
rz(1.8438953) q[3];
sx q[3];
rz(-2.0632931) q[3];
sx q[3];
rz(-0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2099828) q[0];
sx q[0];
rz(-1.7117806) q[0];
sx q[0];
rz(1.1071052) q[0];
rz(2.9529086) q[1];
sx q[1];
rz(-0.37982267) q[1];
sx q[1];
rz(1.3245378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640776) q[0];
sx q[0];
rz(-1.0594089) q[0];
sx q[0];
rz(-1.2426503) q[0];
x q[1];
rz(0.97365427) q[2];
sx q[2];
rz(-1.5882379) q[2];
sx q[2];
rz(-0.23250015) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81428567) q[1];
sx q[1];
rz(-1.5734356) q[1];
sx q[1];
rz(1.0304818) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4408165) q[3];
sx q[3];
rz(-1.1808529) q[3];
sx q[3];
rz(-1.5696456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3482427) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(-0.14774518) q[2];
rz(2.31367) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.70656908) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(1.6178004) q[0];
rz(-1.891547) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(-2.7203383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68116248) q[0];
sx q[0];
rz(-1.958485) q[0];
sx q[0];
rz(1.1739385) q[0];
rz(-pi) q[1];
rz(-1.4827316) q[2];
sx q[2];
rz(-0.95092623) q[2];
sx q[2];
rz(-1.6234835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7819216) q[1];
sx q[1];
rz(-3.0542448) q[1];
sx q[1];
rz(-0.10876258) q[1];
rz(-pi) q[2];
rz(-0.45205595) q[3];
sx q[3];
rz(-1.789942) q[3];
sx q[3];
rz(-0.62458405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0201575) q[2];
sx q[2];
rz(-2.406106) q[2];
sx q[2];
rz(0.45041931) q[2];
rz(1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12238518) q[0];
sx q[0];
rz(-1.5337802) q[0];
sx q[0];
rz(1.5564729) q[0];
rz(-2.4284594) q[1];
sx q[1];
rz(-1.5596371) q[1];
sx q[1];
rz(0.023991931) q[1];
rz(2.0959185) q[2];
sx q[2];
rz(-1.1633506) q[2];
sx q[2];
rz(1.9539471) q[2];
rz(-0.080056277) q[3];
sx q[3];
rz(-1.6430942) q[3];
sx q[3];
rz(-2.5157635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
