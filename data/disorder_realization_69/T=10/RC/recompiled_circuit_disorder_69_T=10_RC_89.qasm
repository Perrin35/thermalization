OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(2.6740958) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43359783) q[0];
sx q[0];
rz(-0.85277075) q[0];
sx q[0];
rz(0.17311592) q[0];
rz(2.1490287) q[2];
sx q[2];
rz(-2.7928068) q[2];
sx q[2];
rz(2.0860096) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.8036098) q[1];
sx q[1];
rz(2.0884872) q[1];
x q[2];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(-2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(-3.0541259) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-2.6175245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0511765) q[0];
sx q[0];
rz(-0.91958445) q[0];
sx q[0];
rz(2.766846) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2696813) q[2];
sx q[2];
rz(-1.2362091) q[2];
sx q[2];
rz(-1.7988234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8088278) q[1];
sx q[1];
rz(-1.1953224) q[1];
sx q[1];
rz(1.3198225) q[1];
rz(-1.4161795) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(-0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-0.18951167) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42650578) q[0];
sx q[0];
rz(-0.77930342) q[0];
sx q[0];
rz(-0.39003178) q[0];
x q[1];
rz(0.26280975) q[2];
sx q[2];
rz(-0.38729471) q[2];
sx q[2];
rz(-0.58236052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7889001) q[1];
sx q[1];
rz(-0.2460203) q[1];
sx q[1];
rz(0.23060631) q[1];
rz(-pi) q[2];
rz(0.17574163) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1069964) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-1.0881902) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(2.0514354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50463146) q[0];
sx q[0];
rz(-1.9177027) q[0];
sx q[0];
rz(0.42731254) q[0];
rz(-2.1074739) q[2];
sx q[2];
rz(-2.3148429) q[2];
sx q[2];
rz(0.85988753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62413299) q[1];
sx q[1];
rz(-2.1412666) q[1];
sx q[1];
rz(-1.073451) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.777321) q[3];
sx q[3];
rz(-1.9897451) q[3];
sx q[3];
rz(-1.9358313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-0.85582716) q[2];
rz(1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677508) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(2.2348316) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1779551) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(2.9658085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7829195) q[1];
sx q[1];
rz(-1.9649319) q[1];
sx q[1];
rz(-0.47788099) q[1];
rz(-pi) q[2];
rz(-1.3768457) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(0.84996163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.871792) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.2924682) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(2.4687016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56129365) q[0];
sx q[0];
rz(-1.6742047) q[0];
sx q[0];
rz(-2.1456477) q[0];
rz(-pi) q[1];
rz(0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(1.6753472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8012961) q[1];
sx q[1];
rz(-2.166689) q[1];
sx q[1];
rz(1.8227541) q[1];
rz(-pi) q[2];
rz(2.7298922) q[3];
sx q[3];
rz(-1.6414335) q[3];
sx q[3];
rz(0.27765805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.792753) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(-1.0621747) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(1.7395082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3081449) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(-1.863198) q[0];
x q[1];
rz(-1.8066605) q[2];
sx q[2];
rz(-0.81996041) q[2];
sx q[2];
rz(-1.5136528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2981616) q[1];
sx q[1];
rz(-1.9318252) q[1];
sx q[1];
rz(-1.7346738) q[1];
rz(1.7756895) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(-0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(0.68022234) q[2];
rz(2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(2.6002398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66051018) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(-0.020945992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(2.9609749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.036877) q[1];
sx q[1];
rz(-2.2762205) q[1];
sx q[1];
rz(-2.5887262) q[1];
rz(0.3803216) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(-2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-0.0099649075) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(1.6962956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1709258) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(-2.4404581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2575349) q[2];
sx q[2];
rz(-2.3381655) q[2];
sx q[2];
rz(2.6510309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6362793) q[1];
sx q[1];
rz(-1.5863824) q[1];
sx q[1];
rz(-0.63025766) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5310998) q[3];
sx q[3];
rz(-0.99184147) q[3];
sx q[3];
rz(1.0510774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69892591) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(0.60738579) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(-1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-2.7493431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0454355) q[0];
sx q[0];
rz(-0.7924315) q[0];
sx q[0];
rz(-0.57604726) q[0];
rz(-pi) q[1];
rz(-0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(-3.0160883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2892462) q[1];
sx q[1];
rz(-2.1834071) q[1];
sx q[1];
rz(-2.5958519) q[1];
rz(-pi) q[2];
rz(-1.5997821) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(2.8857185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-1.2148576) q[2];
sx q[2];
rz(-2.3050953) q[2];
sx q[2];
rz(-0.54511025) q[2];
rz(-2.5946887) q[3];
sx q[3];
rz(-0.96596598) q[3];
sx q[3];
rz(2.1706497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
