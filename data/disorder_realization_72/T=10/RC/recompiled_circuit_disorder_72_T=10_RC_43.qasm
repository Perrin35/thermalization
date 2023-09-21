OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2274269) q[0];
sx q[0];
rz(-1.5072855) q[0];
sx q[0];
rz(0.27557296) q[0];
x q[1];
rz(2.9688641) q[2];
sx q[2];
rz(-0.85731259) q[2];
sx q[2];
rz(-1.2781065) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5159113) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(-2.4019269) q[1];
rz(-pi) q[2];
rz(-1.6742168) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(-0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(-2.1388334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8589477) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(-2.100201) q[0];
rz(2.9687256) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(2.6181521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7111295) q[1];
sx q[1];
rz(-1.4158447) q[1];
sx q[1];
rz(-2.1610545) q[1];
rz(0.59870436) q[3];
sx q[3];
rz(-1.2748534) q[3];
sx q[3];
rz(2.242089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(0.88358203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(-1.3921188) q[0];
x q[1];
rz(2.8475464) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(-1.2750212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.8557465) q[1];
rz(1.6885353) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54116762) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(0.11359544) q[0];
rz(2.8470464) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(1.3603269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.354047) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(1.3613308) q[1];
rz(-pi) q[2];
rz(-2.9123141) q[3];
sx q[3];
rz(-1.6212654) q[3];
sx q[3];
rz(2.5073187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-0.91526389) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.6179248) q[0];
sx q[0];
rz(-0.94232725) q[0];
x q[1];
rz(0.5251685) q[2];
sx q[2];
rz(-1.6095973) q[2];
sx q[2];
rz(0.9135439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.535714) q[1];
sx q[1];
rz(-1.6026346) q[1];
sx q[1];
rz(1.4153626) q[1];
x q[2];
rz(-0.79742934) q[3];
sx q[3];
rz(-2.1505822) q[3];
sx q[3];
rz(0.42339719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(-2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3298033) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(-1.1008218) q[0];
x q[1];
rz(0.47607143) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(1.4897886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5571594) q[1];
sx q[1];
rz(-1.6646619) q[1];
sx q[1];
rz(-0.90969109) q[1];
x q[2];
rz(0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-1.126359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61486812) q[0];
sx q[0];
rz(-1.65979) q[0];
sx q[0];
rz(0.29863775) q[0];
x q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4451051) q[1];
sx q[1];
rz(-1.7587874) q[1];
sx q[1];
rz(2.6871215) q[1];
rz(-pi) q[2];
rz(-0.095110006) q[3];
sx q[3];
rz(-1.9169807) q[3];
sx q[3];
rz(0.1971052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-3.1088366) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93452867) q[0];
sx q[0];
rz(-1.8257739) q[0];
sx q[0];
rz(-1.2364788) q[0];
x q[1];
rz(1.007349) q[2];
sx q[2];
rz(-1.1903561) q[2];
sx q[2];
rz(0.13605875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(1.9883518) q[1];
rz(-pi) q[2];
rz(1.950374) q[3];
sx q[3];
rz(-1.4879585) q[3];
sx q[3];
rz(0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(0.73293066) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982614) q[0];
sx q[0];
rz(-1.4364463) q[0];
sx q[0];
rz(-1.6249379) q[0];
x q[1];
rz(-0.23156667) q[2];
sx q[2];
rz(-0.61083691) q[2];
sx q[2];
rz(-1.4996741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-0.81307462) q[1];
x q[2];
rz(-2.5599307) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(1.7074701) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(-0.5982582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.3539679) q[0];
sx q[0];
rz(-1.3593332) q[0];
rz(-pi) q[1];
rz(2.5355849) q[2];
sx q[2];
rz(-0.25642828) q[2];
sx q[2];
rz(-2.3026349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88565247) q[1];
sx q[1];
rz(-1.2068032) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
rz(0.15828295) q[3];
sx q[3];
rz(-1.4083107) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(1.6937704) q[2];
sx q[2];
rz(-0.95477827) q[2];
sx q[2];
rz(1.9670602) q[2];
rz(-2.7491309) q[3];
sx q[3];
rz(-2.2069208) q[3];
sx q[3];
rz(1.8451286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];