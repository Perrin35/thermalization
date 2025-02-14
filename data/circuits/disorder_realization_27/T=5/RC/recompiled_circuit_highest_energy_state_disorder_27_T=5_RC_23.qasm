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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(0.21188307) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(2.7636757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48403063) q[0];
sx q[0];
rz(-3.1320509) q[0];
sx q[0];
rz(-1.6173167) q[0];
x q[1];
rz(-2.7013999) q[2];
sx q[2];
rz(-2.4131916) q[2];
sx q[2];
rz(0.0022526646) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7841872) q[1];
sx q[1];
rz(-2.2012246) q[1];
sx q[1];
rz(1.9131843) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7573646) q[3];
sx q[3];
rz(-0.93776449) q[3];
sx q[3];
rz(0.6642793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8523031) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(-3.0526414) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(-2.037183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1386344) q[0];
sx q[0];
rz(-1.3413652) q[0];
sx q[0];
rz(-0.63284875) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3588793) q[2];
sx q[2];
rz(-1.209895) q[2];
sx q[2];
rz(0.95907839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30518499) q[1];
sx q[1];
rz(-2.6713604) q[1];
sx q[1];
rz(-1.2364619) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16544754) q[3];
sx q[3];
rz(-2.5962127) q[3];
sx q[3];
rz(0.65217962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(3.1125715) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-0.40036449) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5059197) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(-2.1257504) q[0];
x q[1];
rz(-1.3157) q[2];
sx q[2];
rz(-2.2410903) q[2];
sx q[2];
rz(0.75302659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4246042) q[1];
sx q[1];
rz(-1.4115872) q[1];
sx q[1];
rz(-2.6648494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34137643) q[3];
sx q[3];
rz(-1.2814643) q[3];
sx q[3];
rz(-3.0738664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6843188) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(2.686783) q[2];
rz(-1.2416035) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(-0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860745) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-0.00051001471) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(2.9972163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7847608) q[0];
sx q[0];
rz(-2.5481173) q[0];
sx q[0];
rz(-1.9463825) q[0];
rz(-2.8881489) q[2];
sx q[2];
rz(-1.4254693) q[2];
sx q[2];
rz(-3.1325454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6295602) q[1];
sx q[1];
rz(-1.343439) q[1];
sx q[1];
rz(2.4870706) q[1];
rz(-pi) q[2];
rz(1.6779283) q[3];
sx q[3];
rz(-1.0262353) q[3];
sx q[3];
rz(-0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(-0.96251881) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(2.8008154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(-1.0386946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13030355) q[0];
sx q[0];
rz(-1.4892007) q[0];
sx q[0];
rz(2.5528583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4034385) q[2];
sx q[2];
rz(-2.0655491) q[2];
sx q[2];
rz(3.11657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9925527) q[1];
sx q[1];
rz(-0.91616183) q[1];
sx q[1];
rz(2.4813985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.065537621) q[3];
sx q[3];
rz(-1.1315232) q[3];
sx q[3];
rz(-0.36079839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.3795615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370711) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(1.1035408) q[0];
rz(-0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(-0.59741098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99973122) q[0];
sx q[0];
rz(-1.3628565) q[0];
sx q[0];
rz(-0.19241649) q[0];
rz(2.8812203) q[2];
sx q[2];
rz(-0.7625167) q[2];
sx q[2];
rz(2.0338361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.589041) q[1];
sx q[1];
rz(-1.8218166) q[1];
sx q[1];
rz(1.4000721) q[1];
rz(1.9416945) q[3];
sx q[3];
rz(-1.4850332) q[3];
sx q[3];
rz(-0.24468064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(-0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48653212) q[0];
sx q[0];
rz(-1.9831816) q[0];
sx q[0];
rz(-0.038473815) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9486987) q[0];
sx q[0];
rz(-1.7561551) q[0];
sx q[0];
rz(-2.9167487) q[0];
rz(-pi) q[1];
rz(-2.269417) q[2];
sx q[2];
rz(-2.6769014) q[2];
sx q[2];
rz(-1.8605491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0044627) q[1];
sx q[1];
rz(-2.4871768) q[1];
sx q[1];
rz(0.14738247) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41681186) q[3];
sx q[3];
rz(-0.14188611) q[3];
sx q[3];
rz(-2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(2.5433507) q[2];
rz(0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365731) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-0.2463499) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(-2.6447703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(-2.4606649) q[0];
rz(-pi) q[1];
rz(0.17077568) q[2];
sx q[2];
rz(-0.71238067) q[2];
sx q[2];
rz(0.26584372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0223479) q[1];
sx q[1];
rz(-1.8059219) q[1];
sx q[1];
rz(-2.0658595) q[1];
rz(1.0786177) q[3];
sx q[3];
rz(-1.4340622) q[3];
sx q[3];
rz(-2.0355952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42314998) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(1.1540958) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(-0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(-3.0754454) q[0];
rz(-1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(-2.537421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19949958) q[0];
sx q[0];
rz(-1.8460994) q[0];
sx q[0];
rz(1.3574187) q[0];
rz(-0.2653052) q[2];
sx q[2];
rz(-2.5992706) q[2];
sx q[2];
rz(3.0408183) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0688404) q[1];
sx q[1];
rz(-1.3786331) q[1];
sx q[1];
rz(0.56984624) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.484453) q[3];
sx q[3];
rz(-2.1517188) q[3];
sx q[3];
rz(2.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8784647) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(0.43295941) q[2];
rz(-1.912502) q[3];
sx q[3];
rz(-1.2058328) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(-1.2731592) q[0];
rz(0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(-2.926631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5962385) q[0];
sx q[0];
rz(-2.4430877) q[0];
sx q[0];
rz(2.7424988) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0029055) q[2];
sx q[2];
rz(-2.0615163) q[2];
sx q[2];
rz(0.034772074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8394154) q[1];
sx q[1];
rz(-1.068699) q[1];
sx q[1];
rz(0.79151293) q[1];
rz(-pi) q[2];
rz(0.56193476) q[3];
sx q[3];
rz(-1.3519545) q[3];
sx q[3];
rz(-2.9666025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(-0.95883933) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6017629) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(-0.70855793) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(2.5790527) q[2];
sx q[2];
rz(-1.2934791) q[2];
sx q[2];
rz(0.51359609) q[2];
rz(0.72193969) q[3];
sx q[3];
rz(-1.4293213) q[3];
sx q[3];
rz(-2.8534129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
