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
rz(0.053057916) q[0];
sx q[0];
rz(0.36202708) q[0];
sx q[0];
rz(10.590635) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(1.4282164) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7390763) q[0];
sx q[0];
rz(-1.5105643) q[0];
sx q[0];
rz(-0.51243685) q[0];
x q[1];
rz(-0.25144724) q[2];
sx q[2];
rz(-1.6414974) q[2];
sx q[2];
rz(1.155702) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5184091) q[1];
sx q[1];
rz(-1.759638) q[1];
sx q[1];
rz(3.0075226) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8737439) q[3];
sx q[3];
rz(-1.0976657) q[3];
sx q[3];
rz(-0.79896636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6894138) q[2];
sx q[2];
rz(-3.1248326) q[2];
sx q[2];
rz(2.9407799) q[2];
rz(-2.9928442) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1397322) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(1.0396022) q[0];
rz(-0.014558583) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(-1.5537517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7717465) q[0];
sx q[0];
rz(-0.87939191) q[0];
sx q[0];
rz(0.86172523) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3755685) q[2];
sx q[2];
rz(-3.0659817) q[2];
sx q[2];
rz(-1.6603254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0224198) q[1];
sx q[1];
rz(-1.3117562) q[1];
sx q[1];
rz(-3.1062192) q[1];
rz(2.6282477) q[3];
sx q[3];
rz(-0.37783315) q[3];
sx q[3];
rz(-0.53841089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5130834) q[2];
sx q[2];
rz(-1.5336978) q[2];
sx q[2];
rz(-1.7536564) q[2];
rz(1.7657109) q[3];
sx q[3];
rz(-2.0949771) q[3];
sx q[3];
rz(2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3917711) q[0];
sx q[0];
rz(-0.23673683) q[0];
sx q[0];
rz(2.5340875) q[0];
rz(1.5982184) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(-0.9651331) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133234) q[0];
sx q[0];
rz(-1.566727) q[0];
sx q[0];
rz(-1.550287) q[0];
rz(0.34334646) q[2];
sx q[2];
rz(-2.008524) q[2];
sx q[2];
rz(0.9870607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7399921) q[1];
sx q[1];
rz(-1.5349746) q[1];
sx q[1];
rz(1.7154681) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2854961) q[3];
sx q[3];
rz(-3.0149547) q[3];
sx q[3];
rz(0.52982012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33660108) q[2];
sx q[2];
rz(-2.470863) q[2];
sx q[2];
rz(-2.2750308) q[2];
rz(2.0349515) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5690145) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(1.5255852) q[0];
rz(-3.1309879) q[1];
sx q[1];
rz(-3.1378101) q[1];
sx q[1];
rz(2.3912281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5724154) q[0];
sx q[0];
rz(-2.4710725) q[0];
sx q[0];
rz(-2.9886093) q[0];
rz(-0.60549824) q[2];
sx q[2];
rz(-1.3829621) q[2];
sx q[2];
rz(-2.2772636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1018095) q[1];
sx q[1];
rz(-1.0652115) q[1];
sx q[1];
rz(-1.702448) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8298805) q[3];
sx q[3];
rz(-1.4501713) q[3];
sx q[3];
rz(-1.3159304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7287207) q[2];
sx q[2];
rz(-1.0888638) q[2];
sx q[2];
rz(-1.2415761) q[2];
rz(0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-2.2976105) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64549696) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(-2.2182933) q[0];
rz(-0.80054379) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(2.961535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6430017) q[0];
sx q[0];
rz(-0.24234903) q[0];
sx q[0];
rz(1.7447628) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6864941) q[2];
sx q[2];
rz(-0.1265993) q[2];
sx q[2];
rz(1.1243658) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99111588) q[1];
sx q[1];
rz(-1.2186945) q[1];
sx q[1];
rz(1.8262509) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0308444) q[3];
sx q[3];
rz(-1.2071273) q[3];
sx q[3];
rz(1.2895541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8397612) q[2];
sx q[2];
rz(-1.3000891) q[2];
sx q[2];
rz(-1.4361471) q[2];
rz(-1.885421) q[3];
sx q[3];
rz(-1.6585645) q[3];
sx q[3];
rz(-0.095631599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489478) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(2.6982464) q[0];
rz(-1.1706785) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(-0.43177691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0042808) q[0];
sx q[0];
rz(-2.3385424) q[0];
sx q[0];
rz(-0.43907586) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3987006) q[2];
sx q[2];
rz(-1.6531214) q[2];
sx q[2];
rz(-0.70647994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1990423) q[1];
sx q[1];
rz(-1.1909134) q[1];
sx q[1];
rz(2.540285) q[1];
rz(-0.33080868) q[3];
sx q[3];
rz(-1.7886297) q[3];
sx q[3];
rz(3.1200925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7410437) q[2];
sx q[2];
rz(-0.86047518) q[2];
sx q[2];
rz(-1.7575556) q[2];
rz(2.4436229) q[3];
sx q[3];
rz(-0.75708404) q[3];
sx q[3];
rz(-2.82011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8292002) q[0];
sx q[0];
rz(-1.7990524) q[0];
sx q[0];
rz(-0.37305748) q[0];
rz(-0.27698764) q[1];
sx q[1];
rz(-0.00033683446) q[1];
sx q[1];
rz(2.3854947) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6269659) q[0];
sx q[0];
rz(-1.9494162) q[0];
sx q[0];
rz(0.65831229) q[0];
x q[1];
rz(-2.4239642) q[2];
sx q[2];
rz(-2.6127079) q[2];
sx q[2];
rz(2.8019049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0758772) q[1];
sx q[1];
rz(-1.0317894) q[1];
sx q[1];
rz(2.6679941) q[1];
x q[2];
rz(2.1806898) q[3];
sx q[3];
rz(-2.641819) q[3];
sx q[3];
rz(0.62675301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9296391) q[2];
sx q[2];
rz(-2.5888207) q[2];
sx q[2];
rz(2.2201404) q[2];
rz(-3.0742505) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(-1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15143722) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(0.13023278) q[0];
rz(-2.3479334) q[1];
sx q[1];
rz(-3.1399813) q[1];
sx q[1];
rz(-2.8614955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.139629) q[0];
sx q[0];
rz(-3.0462777) q[0];
sx q[0];
rz(-2.167666) q[0];
rz(-2.75861) q[2];
sx q[2];
rz(-0.23442253) q[2];
sx q[2];
rz(1.2146666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4520766) q[1];
sx q[1];
rz(-2.048564) q[1];
sx q[1];
rz(2.6803451) q[1];
rz(1.0956826) q[3];
sx q[3];
rz(-0.91714782) q[3];
sx q[3];
rz(-0.16330367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0555931) q[2];
sx q[2];
rz(-1.6722101) q[2];
sx q[2];
rz(-2.3386193) q[2];
rz(-1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-2.3101961) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(0.10920864) q[0];
rz(-0.373492) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(2.5901897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5057988) q[0];
sx q[0];
rz(-2.3058545) q[0];
sx q[0];
rz(-2.4018025) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5481651) q[2];
sx q[2];
rz(-2.9205236) q[2];
sx q[2];
rz(1.9367557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0903328) q[1];
sx q[1];
rz(-0.96104927) q[1];
sx q[1];
rz(2.4734797) q[1];
x q[2];
rz(0.53111303) q[3];
sx q[3];
rz(-0.45980849) q[3];
sx q[3];
rz(-2.8961033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52770829) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(-1.323553) q[2];
rz(-1.8771133) q[3];
sx q[3];
rz(-1.2885965) q[3];
sx q[3];
rz(3.1356139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962113) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(0.7363466) q[0];
rz(2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(-1.597499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1101226) q[0];
sx q[0];
rz(-1.2916864) q[0];
sx q[0];
rz(-1.5895248) q[0];
rz(-pi) q[1];
rz(-1.5734463) q[2];
sx q[2];
rz(-1.6011597) q[2];
sx q[2];
rz(-1.4043491) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5527344) q[1];
sx q[1];
rz(-2.169796) q[1];
sx q[1];
rz(1.0042648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2451363) q[3];
sx q[3];
rz(-1.146637) q[3];
sx q[3];
rz(-2.2617634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.836901) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(1.2738073) q[2];
rz(-1.4443719) q[3];
sx q[3];
rz(-3.0675409) q[3];
sx q[3];
rz(-1.4950289) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16784167) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(1.6043067) q[1];
sx q[1];
rz(-2.2289386) q[1];
sx q[1];
rz(-2.9569721) q[1];
rz(0.040097728) q[2];
sx q[2];
rz(-1.5798777) q[2];
sx q[2];
rz(-0.29691534) q[2];
rz(0.96327412) q[3];
sx q[3];
rz(-1.6714851) q[3];
sx q[3];
rz(2.1906399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
