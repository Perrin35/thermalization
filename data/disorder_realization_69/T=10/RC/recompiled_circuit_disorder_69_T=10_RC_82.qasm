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
rz(-0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2517393) q[0];
sx q[0];
rz(-1.4407053) q[0];
sx q[0];
rz(-0.84530172) q[0];
rz(-pi) q[1];
rz(1.8664076) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(-1.0653898) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(1.0531055) q[1];
rz(-pi) q[2];
rz(1.7476139) q[3];
sx q[3];
rz(-1.8999945) q[3];
sx q[3];
rz(0.65229177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(-2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4755718) q[0];
sx q[0];
rz(-0.7374987) q[0];
sx q[0];
rz(-2.0185508) q[0];
rz(1.0753724) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(-2.5410595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8088278) q[1];
sx q[1];
rz(-1.9462703) q[1];
sx q[1];
rz(-1.3198225) q[1];
x q[2];
rz(1.4161795) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(-2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(-2.5168193) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-0.18951167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097475826) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(-1.9301027) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37522845) q[2];
sx q[2];
rz(-1.4725176) q[2];
sx q[2];
rz(0.74429846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5901958) q[1];
sx q[1];
rz(-1.8101748) q[1];
sx q[1];
rz(1.6281284) q[1];
rz(-pi) q[2];
rz(-0.96954815) q[3];
sx q[3];
rz(-1.7161955) q[3];
sx q[3];
rz(1.6325412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1069964) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(0.43131367) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(-2.003147) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(1.0901573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074993) q[0];
sx q[0];
rz(-2.5980237) q[0];
sx q[0];
rz(0.71732934) q[0];
rz(0.50699373) q[2];
sx q[2];
rz(-2.255313) q[2];
sx q[2];
rz(0.13912858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.412147) q[1];
sx q[1];
rz(-0.73819654) q[1];
sx q[1];
rz(-0.63936887) q[1];
rz(-pi) q[2];
rz(-2.7146043) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-0.85582716) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.1700464) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677508) q[0];
sx q[0];
rz(-1.2497328) q[0];
sx q[0];
rz(-2.2348316) q[0];
rz(-pi) q[1];
rz(-2.8180426) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(-2.8138585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7331446) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(-2.0088197) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31801362) q[3];
sx q[3];
rz(-0.56088305) q[3];
sx q[3];
rz(-2.6624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.871792) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(1.8491245) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56129365) q[0];
sx q[0];
rz(-1.467388) q[0];
sx q[0];
rz(-2.1456477) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1168675) q[2];
sx q[2];
rz(-2.284986) q[2];
sx q[2];
rz(-1.6753472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7676047) q[1];
sx q[1];
rz(-1.3629706) q[1];
sx q[1];
rz(0.61088224) q[1];
x q[2];
rz(0.41170044) q[3];
sx q[3];
rz(-1.5001591) q[3];
sx q[3];
rz(-2.8639346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.792753) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(1.0621747) q[0];
rz(1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7715493) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(2.7605961) q[0];
x q[1];
rz(2.3767396) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(-2.9219251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7359888) q[1];
sx q[1];
rz(-0.39499184) q[1];
sx q[1];
rz(-2.7337381) q[1];
rz(-pi) q[2];
rz(3.0038463) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4222251) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(0.24958615) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-2.6002398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4562948) q[0];
sx q[0];
rz(-2.1349499) q[0];
sx q[0];
rz(-1.5575404) q[0];
x q[1];
rz(-1.2301684) q[2];
sx q[2];
rz(-1.1414141) q[2];
sx q[2];
rz(-0.18061772) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.036877) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(-2.5887262) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5726611) q[3];
sx q[3];
rz(-1.1904753) q[3];
sx q[3];
rz(2.5170381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-0.0099649075) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24628595) q[0];
sx q[0];
rz(-0.9406957) q[0];
sx q[0];
rz(-0.90755983) q[0];
x q[1];
rz(2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(-2.5788139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.054113764) q[1];
sx q[1];
rz(-2.2009654) q[1];
sx q[1];
rz(-1.5900882) q[1];
x q[2];
rz(-0.61049283) q[3];
sx q[3];
rz(-0.99184147) q[3];
sx q[3];
rz(-2.0905153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(-2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.8764953) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(-2.4368068) q[0];
rz(-pi) q[1];
rz(0.99598155) q[2];
sx q[2];
rz(-2.1157584) q[2];
sx q[2];
rz(-1.121322) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0872005) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(-0.88263504) q[1];
rz(-3.0030389) q[3];
sx q[3];
rz(-0.20692736) q[3];
sx q[3];
rz(3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(2.0742119) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(2.773182) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(0.89041238) q[3];
sx q[3];
rz(-1.1287516) q[3];
sx q[3];
rz(0.26656084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
