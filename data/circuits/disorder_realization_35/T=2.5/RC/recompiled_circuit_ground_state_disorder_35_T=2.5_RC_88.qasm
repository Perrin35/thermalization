OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(2.3647302) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(1.3219272) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8213733) q[0];
sx q[0];
rz(-0.8614653) q[0];
sx q[0];
rz(0.18289645) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58682449) q[2];
sx q[2];
rz(-0.57673645) q[2];
sx q[2];
rz(1.5962034) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.894134) q[1];
sx q[1];
rz(-1.1620635) q[1];
sx q[1];
rz(-0.031886727) q[1];
rz(2.9453927) q[3];
sx q[3];
rz(-1.8561188) q[3];
sx q[3];
rz(-1.5232435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(0.35749164) q[2];
rz(-2.5743971) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95328632) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(-1.8081007) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(-0.83998799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6463722) q[0];
sx q[0];
rz(-0.53421181) q[0];
sx q[0];
rz(1.2665126) q[0];
x q[1];
rz(0.2510906) q[2];
sx q[2];
rz(-2.5320964) q[2];
sx q[2];
rz(1.239038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31038302) q[1];
sx q[1];
rz(-2.134517) q[1];
sx q[1];
rz(-2.4604164) q[1];
x q[2];
rz(1.5072823) q[3];
sx q[3];
rz(-2.5908785) q[3];
sx q[3];
rz(-0.39001071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(1.4060414) q[2];
rz(2.9794335) q[3];
sx q[3];
rz(-1.5608965) q[3];
sx q[3];
rz(-2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(0.42174569) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.755736) q[1];
sx q[1];
rz(2.8657894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0244747) q[0];
sx q[0];
rz(-1.8016651) q[0];
sx q[0];
rz(-0.33322115) q[0];
rz(-1.3490729) q[2];
sx q[2];
rz(-1.1100551) q[2];
sx q[2];
rz(2.0784476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2308064) q[1];
sx q[1];
rz(-1.2503997) q[1];
sx q[1];
rz(-1.5159392) q[1];
x q[2];
rz(1.4546605) q[3];
sx q[3];
rz(-0.7051055) q[3];
sx q[3];
rz(-3.1381315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3383823) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(-2.4513643) q[2];
rz(-0.35987443) q[3];
sx q[3];
rz(-1.3709603) q[3];
sx q[3];
rz(2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91367078) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-2.2699455) q[0];
rz(-0.40920416) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86996468) q[0];
sx q[0];
rz(-1.6883435) q[0];
sx q[0];
rz(2.6948476) q[0];
rz(2.6957507) q[2];
sx q[2];
rz(-2.9135357) q[2];
sx q[2];
rz(2.7331309) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.732199) q[1];
sx q[1];
rz(-2.7229561) q[1];
sx q[1];
rz(2.6217209) q[1];
rz(-2.8866037) q[3];
sx q[3];
rz(-0.6643799) q[3];
sx q[3];
rz(2.5685316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(2.3166166) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(2.2671674) q[0];
rz(2.8127316) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(-2.0430476) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335137) q[0];
sx q[0];
rz(-2.7493434) q[0];
sx q[0];
rz(0.29229887) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6107164) q[2];
sx q[2];
rz(-2.4062059) q[2];
sx q[2];
rz(0.11271726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5034628) q[1];
sx q[1];
rz(-2.0582804) q[1];
sx q[1];
rz(-2.8007228) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8642578) q[3];
sx q[3];
rz(-2.1430052) q[3];
sx q[3];
rz(-0.25857833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99981368) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-2.8580247) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11033002) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(1.7562235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2792086) q[0];
sx q[0];
rz(-0.68717903) q[0];
sx q[0];
rz(1.1670668) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4848861) q[2];
sx q[2];
rz(-2.067778) q[2];
sx q[2];
rz(0.53705207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54810235) q[1];
sx q[1];
rz(-1.2498651) q[1];
sx q[1];
rz(-1.3671419) q[1];
x q[2];
rz(-2.9870944) q[3];
sx q[3];
rz(-1.6067926) q[3];
sx q[3];
rz(0.29496962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.063227) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(0.37609491) q[2];
rz(-2.0070576) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.1217839) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(-0.58492297) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(1.9580511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4485156) q[0];
sx q[0];
rz(-1.7127556) q[0];
sx q[0];
rz(-1.4407115) q[0];
x q[1];
rz(0.41924119) q[2];
sx q[2];
rz(-2.4346874) q[2];
sx q[2];
rz(-2.8173994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21321024) q[1];
sx q[1];
rz(-2.2021535) q[1];
sx q[1];
rz(-0.31698314) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46681301) q[3];
sx q[3];
rz(-1.3315505) q[3];
sx q[3];
rz(1.6166663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(2.9417876) q[2];
rz(2.0810769) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(-0.04960355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878433) q[0];
sx q[0];
rz(-1.6596154) q[0];
sx q[0];
rz(0.31295452) q[0];
rz(0.9494268) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.688028) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87438238) q[0];
sx q[0];
rz(-1.6613164) q[0];
sx q[0];
rz(0.1874013) q[0];
rz(3.0309148) q[2];
sx q[2];
rz(-2.5191436) q[2];
sx q[2];
rz(1.770293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1968699) q[1];
sx q[1];
rz(-0.88366449) q[1];
sx q[1];
rz(-1.8623167) q[1];
x q[2];
rz(-0.62409921) q[3];
sx q[3];
rz(-0.54275504) q[3];
sx q[3];
rz(1.8763157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7877385) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(-1.3346416) q[2];
rz(-1.3245964) q[3];
sx q[3];
rz(-2.4748804) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036309328) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(2.4435254) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-1.061729) q[1];
sx q[1];
rz(0.36044136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1669641) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(-1.0247158) q[0];
rz(-pi) q[1];
rz(-2.7588604) q[2];
sx q[2];
rz(-1.6305411) q[2];
sx q[2];
rz(-2.003423) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1680582) q[1];
sx q[1];
rz(-2.413001) q[1];
sx q[1];
rz(0.97961564) q[1];
x q[2];
rz(0.27020755) q[3];
sx q[3];
rz(-0.5654618) q[3];
sx q[3];
rz(-2.5266441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2299819) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(2.3243375) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(-0.219492) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5068186) q[0];
sx q[0];
rz(-0.84832484) q[0];
sx q[0];
rz(-3.1095374) q[0];
rz(1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(-0.65418902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830156) q[0];
sx q[0];
rz(-2.9111283) q[0];
sx q[0];
rz(2.5925267) q[0];
x q[1];
rz(-1.1172574) q[2];
sx q[2];
rz(-2.3058866) q[2];
sx q[2];
rz(-1.0356366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6950049) q[1];
sx q[1];
rz(-2.4420629) q[1];
sx q[1];
rz(2.3746731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93864949) q[3];
sx q[3];
rz(-1.3397927) q[3];
sx q[3];
rz(1.7836026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0111982) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(2.1350258) q[2];
rz(2.448163) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7158159) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(-0.88059942) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(2.2702552) q[2];
sx q[2];
rz(-1.9819145) q[2];
sx q[2];
rz(-2.869538) q[2];
rz(1.1893336) q[3];
sx q[3];
rz(-0.91100024) q[3];
sx q[3];
rz(-0.80445214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
