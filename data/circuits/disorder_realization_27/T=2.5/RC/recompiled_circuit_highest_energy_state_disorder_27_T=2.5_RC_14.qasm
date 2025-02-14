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
rz(-1.6838411) q[0];
sx q[0];
rz(-0.10310752) q[0];
sx q[0];
rz(0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65771253) q[0];
sx q[0];
rz(-1.9498511) q[0];
sx q[0];
rz(-2.0188278) q[0];
rz(-1.3288055) q[2];
sx q[2];
rz(-2.274161) q[2];
sx q[2];
rz(2.6405328) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9903914) q[1];
sx q[1];
rz(-1.5631274) q[1];
sx q[1];
rz(2.5911314) q[1];
x q[2];
rz(-1.9363251) q[3];
sx q[3];
rz(-0.46896471) q[3];
sx q[3];
rz(-0.381857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.890581) q[2];
sx q[2];
rz(-0.0958395) q[2];
sx q[2];
rz(1.7354234) q[2];
rz(0.74875325) q[3];
sx q[3];
rz(-2.483832) q[3];
sx q[3];
rz(0.21193084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9894079) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(0.33921355) q[0];
rz(-2.5665414) q[1];
sx q[1];
rz(-1.9564068) q[1];
sx q[1];
rz(-0.38160479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4313916) q[0];
sx q[0];
rz(-1.1553191) q[0];
sx q[0];
rz(0.25449591) q[0];
rz(-1.9100045) q[2];
sx q[2];
rz(-1.3217862) q[2];
sx q[2];
rz(-2.9033962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72605342) q[1];
sx q[1];
rz(-1.9697297) q[1];
sx q[1];
rz(0.010460214) q[1];
x q[2];
rz(-1.0142542) q[3];
sx q[3];
rz(-2.294971) q[3];
sx q[3];
rz(-2.4374645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40994689) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(-2.6551969) q[2];
rz(0.54388034) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(0.16415183) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.630702) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(-2.1422332) q[0];
rz(0.73127812) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(-0.17459248) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0431049) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(1.5763603) q[0];
rz(1.393494) q[2];
sx q[2];
rz(-1.7288704) q[2];
sx q[2];
rz(2.6690935) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0831956) q[1];
sx q[1];
rz(-3.0388012) q[1];
sx q[1];
rz(1.0344857) q[1];
rz(1.1483795) q[3];
sx q[3];
rz(-0.49634051) q[3];
sx q[3];
rz(-1.6272735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38873765) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(-1.8641776) q[2];
rz(0.82469213) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81553066) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(1.1482358) q[0];
rz(0.3321906) q[1];
sx q[1];
rz(-2.5416608) q[1];
sx q[1];
rz(1.0848328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116219) q[0];
sx q[0];
rz(-2.8154439) q[0];
sx q[0];
rz(1.1205313) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81697322) q[2];
sx q[2];
rz(-1.3187871) q[2];
sx q[2];
rz(-1.9793881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12382896) q[1];
sx q[1];
rz(-2.4014086) q[1];
sx q[1];
rz(-0.077880903) q[1];
rz(-2.3938136) q[3];
sx q[3];
rz(-2.5430395) q[3];
sx q[3];
rz(-2.7876496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79160249) q[2];
sx q[2];
rz(-0.76452667) q[2];
sx q[2];
rz(0.0066268607) q[2];
rz(2.918112) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(-2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70581907) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(-1.9594132) q[0];
rz(1.301282) q[1];
sx q[1];
rz(-2.2186406) q[1];
sx q[1];
rz(-0.15929793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2625339) q[0];
sx q[0];
rz(-0.28994432) q[0];
sx q[0];
rz(3.0109809) q[0];
x q[1];
rz(-0.18929357) q[2];
sx q[2];
rz(-1.2603501) q[2];
sx q[2];
rz(0.60223641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0634126) q[1];
sx q[1];
rz(-1.0920807) q[1];
sx q[1];
rz(-0.53005752) q[1];
x q[2];
rz(-0.42493762) q[3];
sx q[3];
rz(-2.1930088) q[3];
sx q[3];
rz(3.034301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84676877) q[2];
sx q[2];
rz(-2.939665) q[2];
sx q[2];
rz(-1.798604) q[2];
rz(2.790847) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(-0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89966929) q[0];
sx q[0];
rz(-0.82829183) q[0];
sx q[0];
rz(-0.1668461) q[0];
rz(-0.22653656) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(0.47795263) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.307823) q[0];
sx q[0];
rz(-1.2306899) q[0];
sx q[0];
rz(0.97619636) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9327546) q[2];
sx q[2];
rz(-0.88506341) q[2];
sx q[2];
rz(-0.63206965) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30580995) q[1];
sx q[1];
rz(-1.0210345) q[1];
sx q[1];
rz(2.1533986) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80676956) q[3];
sx q[3];
rz(-0.58611996) q[3];
sx q[3];
rz(2.7703843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6106674) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(-1.8492071) q[2];
rz(2.2409706) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(-1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6677299) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(1.5245755) q[0];
rz(-1.698311) q[1];
sx q[1];
rz(-0.89569211) q[1];
sx q[1];
rz(2.6406094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9568142) q[0];
sx q[0];
rz(-1.9405455) q[0];
sx q[0];
rz(0.92147227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8800182) q[2];
sx q[2];
rz(-2.1078601) q[2];
sx q[2];
rz(1.352965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0308548) q[1];
sx q[1];
rz(-2.1424865) q[1];
sx q[1];
rz(2.5377889) q[1];
rz(-pi) q[2];
rz(0.37338169) q[3];
sx q[3];
rz(-0.33003673) q[3];
sx q[3];
rz(-0.16597834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3524126) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(2.6782356) q[2];
rz(-3.0739259) q[3];
sx q[3];
rz(-1.3031518) q[3];
sx q[3];
rz(2.9786003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891334) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(0.64396089) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(-2.8535829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2525507) q[0];
sx q[0];
rz(-1.5324549) q[0];
sx q[0];
rz(-2.949723) q[0];
x q[1];
rz(-2.5019849) q[2];
sx q[2];
rz(-2.158463) q[2];
sx q[2];
rz(-1.7220875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5254613) q[1];
sx q[1];
rz(-1.0122293) q[1];
sx q[1];
rz(1.2458879) q[1];
x q[2];
rz(0.27585101) q[3];
sx q[3];
rz(-0.39182651) q[3];
sx q[3];
rz(1.0759169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94706941) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(-0.21491773) q[2];
rz(-1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(-0.18068331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27338481) q[0];
sx q[0];
rz(-0.93623638) q[0];
sx q[0];
rz(2.0599763) q[0];
rz(0.99010211) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(2.8066011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43083999) q[0];
sx q[0];
rz(-2.8169605) q[0];
sx q[0];
rz(-1.4650833) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0063517687) q[2];
sx q[2];
rz(-1.9926903) q[2];
sx q[2];
rz(-0.61074257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27410848) q[1];
sx q[1];
rz(-1.6699526) q[1];
sx q[1];
rz(1.656437) q[1];
rz(3.1331691) q[3];
sx q[3];
rz(-1.971774) q[3];
sx q[3];
rz(0.13051912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.058502402) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(-2.7527909) q[2];
rz(2.7723516) q[3];
sx q[3];
rz(-0.25596127) q[3];
sx q[3];
rz(-2.2970439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19287547) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(2.4396851) q[0];
rz(-1.5319872) q[1];
sx q[1];
rz(-0.67370266) q[1];
sx q[1];
rz(0.38175499) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2546291) q[0];
sx q[0];
rz(-1.5077295) q[0];
sx q[0];
rz(3.0596759) q[0];
rz(1.4456621) q[2];
sx q[2];
rz(-1.9788673) q[2];
sx q[2];
rz(3.0681075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.020479105) q[1];
sx q[1];
rz(-1.8217297) q[1];
sx q[1];
rz(-0.24848715) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66017229) q[3];
sx q[3];
rz(-2.2260465) q[3];
sx q[3];
rz(1.2869664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4664885) q[2];
sx q[2];
rz(-1.0305104) q[2];
sx q[2];
rz(-0.64811903) q[2];
rz(1.0068007) q[3];
sx q[3];
rz(-1.1454134) q[3];
sx q[3];
rz(-3.0692611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.61160144) q[0];
sx q[0];
rz(-1.7399104) q[0];
sx q[0];
rz(-2.4763784) q[0];
rz(2.8499659) q[1];
sx q[1];
rz(-1.7886152) q[1];
sx q[1];
rz(2.0033966) q[1];
rz(-1.2164581) q[2];
sx q[2];
rz(-1.0783429) q[2];
sx q[2];
rz(1.1781296) q[2];
rz(-1.9129561) q[3];
sx q[3];
rz(-0.4426601) q[3];
sx q[3];
rz(1.9667251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
