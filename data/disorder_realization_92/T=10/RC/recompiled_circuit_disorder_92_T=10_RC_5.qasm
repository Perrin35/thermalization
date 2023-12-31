OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(-3.1080973) q[0];
sx q[0];
rz(-1.7749696) q[0];
rz(2.0904436) q[1];
sx q[1];
rz(-1.6519974) q[1];
sx q[1];
rz(-1.1319914) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7819408) q[0];
sx q[0];
rz(-1.2255166) q[0];
sx q[0];
rz(-2.3495673) q[0];
rz(-pi) q[1];
rz(2.3425383) q[2];
sx q[2];
rz(-0.22944268) q[2];
sx q[2];
rz(-0.45861751) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72343091) q[1];
sx q[1];
rz(-1.9730113) q[1];
sx q[1];
rz(2.1869786) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2855929) q[3];
sx q[3];
rz(-2.1072901) q[3];
sx q[3];
rz(0.59629089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(-1.1532016) q[2];
rz(0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(-0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.3020878) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(1.5867651) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(-0.15393004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4707697) q[0];
sx q[0];
rz(-1.5851067) q[0];
sx q[0];
rz(-1.0679246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2203127) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(-2.9608179) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31175266) q[1];
sx q[1];
rz(-2.0318188) q[1];
sx q[1];
rz(0.47383576) q[1];
x q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(1.9382167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2543891) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(-1.3228234) q[2];
rz(1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-0.71050182) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674295) q[0];
sx q[0];
rz(-2.0236334) q[0];
sx q[0];
rz(2.2316566) q[0];
rz(0.77727708) q[1];
sx q[1];
rz(-2.3060019) q[1];
sx q[1];
rz(-2.1562703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9604208) q[0];
sx q[0];
rz(-0.67502484) q[0];
sx q[0];
rz(1.3135364) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2276332) q[2];
sx q[2];
rz(-1.3105536) q[2];
sx q[2];
rz(0.39011207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4640376) q[1];
sx q[1];
rz(-2.8863393) q[1];
sx q[1];
rz(-0.21932253) q[1];
rz(2.3722234) q[3];
sx q[3];
rz(-0.83988512) q[3];
sx q[3];
rz(1.6826671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2788006) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(1.8939691) q[2];
rz(-0.4425846) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(-2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2414395) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(-1.1234269) q[0];
rz(0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.3935864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508586) q[0];
sx q[0];
rz(-0.77461857) q[0];
sx q[0];
rz(-0.61268341) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4030928) q[2];
sx q[2];
rz(-1.1412914) q[2];
sx q[2];
rz(-1.6657366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0949888) q[1];
sx q[1];
rz(-0.91218439) q[1];
sx q[1];
rz(0.98595001) q[1];
x q[2];
rz(-2.9080503) q[3];
sx q[3];
rz(-2.2925348) q[3];
sx q[3];
rz(-2.8535709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5559876) q[2];
sx q[2];
rz(3.139479) q[2];
rz(0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(-0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.026022) q[0];
sx q[0];
rz(-1.1239115) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(-1.4670124) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(-1.7747169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05058771) q[0];
sx q[0];
rz(-2.707621) q[0];
sx q[0];
rz(2.6758053) q[0];
rz(-1.5405802) q[2];
sx q[2];
rz(-2.317791) q[2];
sx q[2];
rz(0.012416427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83134507) q[1];
sx q[1];
rz(-0.86569769) q[1];
sx q[1];
rz(0.10465937) q[1];
rz(-pi) q[2];
rz(-2.6392691) q[3];
sx q[3];
rz(-1.886743) q[3];
sx q[3];
rz(-1.5497006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(-2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(-0.50338411) q[0];
rz(0.21884306) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(2.4898081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3394649) q[0];
sx q[0];
rz(-1.3314684) q[0];
sx q[0];
rz(1.7665187) q[0];
rz(-pi) q[1];
rz(-1.6794372) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(-1.7679364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5946878) q[1];
sx q[1];
rz(-1.5190131) q[1];
sx q[1];
rz(-1.4361708) q[1];
x q[2];
rz(-2.8586646) q[3];
sx q[3];
rz(-1.471721) q[3];
sx q[3];
rz(-2.2968452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.3971389) q[2];
sx q[2];
rz(-0.39166489) q[2];
rz(-0.20646778) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713292) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(2.1719334) q[0];
rz(0.58352739) q[1];
sx q[1];
rz(-1.1279761) q[1];
sx q[1];
rz(3.0113509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7182065) q[0];
sx q[0];
rz(-1.5209746) q[0];
sx q[0];
rz(-2.7888265) q[0];
x q[1];
rz(2.4415605) q[2];
sx q[2];
rz(-0.80169741) q[2];
sx q[2];
rz(1.9966372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8091734) q[1];
sx q[1];
rz(-1.0671339) q[1];
sx q[1];
rz(-1.3219576) q[1];
rz(-0.98826615) q[3];
sx q[3];
rz(-1.8997972) q[3];
sx q[3];
rz(0.44579166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5600081) q[2];
sx q[2];
rz(-1.0858034) q[2];
rz(0.052224934) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(-2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645585) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(-1.9751256) q[0];
rz(-2.4160066) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(-2.3988147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22314534) q[0];
sx q[0];
rz(-2.3728275) q[0];
sx q[0];
rz(-2.8630775) q[0];
rz(2.1002662) q[2];
sx q[2];
rz(-0.96457446) q[2];
sx q[2];
rz(0.60339061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0689773) q[1];
sx q[1];
rz(-1.9779357) q[1];
sx q[1];
rz(-0.8168656) q[1];
x q[2];
rz(-2.2320896) q[3];
sx q[3];
rz(-1.0854183) q[3];
sx q[3];
rz(-0.73736008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(-2.880704) q[2];
rz(-2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(-2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7850007) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(0.93604952) q[0];
rz(-3.127457) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(0.65151185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148841) q[0];
sx q[0];
rz(-1.5944905) q[0];
sx q[0];
rz(-0.004304927) q[0];
rz(0.75479836) q[2];
sx q[2];
rz(-2.3279466) q[2];
sx q[2];
rz(2.1458643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41496719) q[1];
sx q[1];
rz(-1.8888998) q[1];
sx q[1];
rz(3.0148562) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43990073) q[3];
sx q[3];
rz(-1.9045881) q[3];
sx q[3];
rz(1.4378289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8462048) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(2.8335617) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(-1.3380922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4073407) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(2.667528) q[1];
sx q[1];
rz(-1.4539376) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5546075) q[0];
sx q[0];
rz(-0.73903144) q[0];
sx q[0];
rz(0.87297312) q[0];
x q[1];
rz(0.7147185) q[2];
sx q[2];
rz(-0.93313365) q[2];
sx q[2];
rz(1.4710466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.6369032) q[1];
sx q[1];
rz(-0.46823374) q[1];
sx q[1];
rz(0.56028985) q[1];
x q[2];
rz(-2.8836807) q[3];
sx q[3];
rz(-2.5690418) q[3];
sx q[3];
rz(0.27975988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-2.0861574) q[3];
sx q[3];
rz(2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67125852) q[0];
sx q[0];
rz(-2.8699734) q[0];
sx q[0];
rz(1.097453) q[0];
rz(0.63256565) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(-0.40614265) q[2];
sx q[2];
rz(-2.4293025) q[2];
sx q[2];
rz(0.14355125) q[2];
rz(-2.408151) q[3];
sx q[3];
rz(-2.1248795) q[3];
sx q[3];
rz(-2.636573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
