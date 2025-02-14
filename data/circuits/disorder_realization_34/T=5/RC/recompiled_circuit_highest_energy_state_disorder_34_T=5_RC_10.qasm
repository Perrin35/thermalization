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
rz(-1.7582769) q[0];
sx q[0];
rz(-1.7051237) q[0];
sx q[0];
rz(-0.96487784) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(5.8536018) q[1];
sx q[1];
rz(11.516973) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222264) q[0];
sx q[0];
rz(-0.99577409) q[0];
sx q[0];
rz(2.7336043) q[0];
rz(-pi) q[1];
rz(2.9994316) q[2];
sx q[2];
rz(-1.249525) q[2];
sx q[2];
rz(-0.69883332) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2856372) q[1];
sx q[1];
rz(-0.95889839) q[1];
sx q[1];
rz(2.6735071) q[1];
rz(1.4780294) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(0.64914671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6875978) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(-0.018608658) q[2];
rz(2.5971557) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740771) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(-2.6065705) q[0];
rz(2.6016443) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.7417057) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2311418) q[0];
sx q[0];
rz(-0.60908356) q[0];
sx q[0];
rz(2.2876491) q[0];
x q[1];
rz(-2.5783456) q[2];
sx q[2];
rz(-2.0958825) q[2];
sx q[2];
rz(-1.6890749) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3185841) q[1];
sx q[1];
rz(-2.6806147) q[1];
sx q[1];
rz(2.5557842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0008282) q[3];
sx q[3];
rz(-2.1080351) q[3];
sx q[3];
rz(-1.9443823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(2.2155217) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(-2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7922908) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(1.6287623) q[0];
rz(2.8577562) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(-2.0294752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130198) q[0];
sx q[0];
rz(-1.5637102) q[0];
sx q[0];
rz(-1.5779693) q[0];
x q[1];
rz(-1.6792504) q[2];
sx q[2];
rz(-1.4949833) q[2];
sx q[2];
rz(2.6159942) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49431153) q[1];
sx q[1];
rz(-1.0930645) q[1];
sx q[1];
rz(-1.7812438) q[1];
rz(-pi) q[2];
rz(1.2073969) q[3];
sx q[3];
rz(-1.8135241) q[3];
sx q[3];
rz(-3.0128765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(0.88895041) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(2.8269732) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(1.7291732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88453) q[0];
sx q[0];
rz(-1.8842959) q[0];
sx q[0];
rz(0.78961163) q[0];
rz(-1.6666404) q[2];
sx q[2];
rz(-1.7828701) q[2];
sx q[2];
rz(-0.98835683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38217332) q[1];
sx q[1];
rz(-1.628211) q[1];
sx q[1];
rz(-2.7613954) q[1];
rz(-1.0869157) q[3];
sx q[3];
rz(-1.3984507) q[3];
sx q[3];
rz(0.92007557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(0.55595428) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(-0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81412643) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(0.59447527) q[0];
rz(-2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-0.97602239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3626427) q[0];
sx q[0];
rz(-0.57641006) q[0];
sx q[0];
rz(-2.2767785) q[0];
rz(-pi) q[1];
rz(0.75404928) q[2];
sx q[2];
rz(-1.591914) q[2];
sx q[2];
rz(1.4625975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0343209) q[1];
sx q[1];
rz(-1.8684505) q[1];
sx q[1];
rz(-2.4536282) q[1];
x q[2];
rz(2.5247252) q[3];
sx q[3];
rz(-0.36412334) q[3];
sx q[3];
rz(2.0193651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(2.5308934) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(-1.8122199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139451) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(2.3620391) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(2.4028042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29863623) q[0];
sx q[0];
rz(-1.2901352) q[0];
sx q[0];
rz(2.929351) q[0];
x q[1];
rz(0.03348695) q[2];
sx q[2];
rz(-2.8346363) q[2];
sx q[2];
rz(-2.592776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-2.457239) q[1];
sx q[1];
rz(-2.1696287) q[1];
rz(-pi) q[2];
rz(-2.7802864) q[3];
sx q[3];
rz(-1.2194467) q[3];
sx q[3];
rz(-1.6811937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5256727) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(-2.2789148) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(1.1427243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(1.020485) q[0];
rz(3.0842969) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(0.78757706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6660992) q[0];
sx q[0];
rz(-1.8084053) q[0];
sx q[0];
rz(-0.65639021) q[0];
rz(-2.8767005) q[2];
sx q[2];
rz(-2.8745626) q[2];
sx q[2];
rz(-1.9445813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5503834) q[1];
sx q[1];
rz(-2.0616643) q[1];
sx q[1];
rz(-0.20259133) q[1];
x q[2];
rz(0.034463783) q[3];
sx q[3];
rz(-0.98317819) q[3];
sx q[3];
rz(-1.8922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3202177) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(2.5569432) q[2];
rz(0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24340165) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(0.32824326) q[0];
rz(0.39168656) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(1.4788871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2770689) q[0];
sx q[0];
rz(-1.3448937) q[0];
sx q[0];
rz(-0.27033131) q[0];
x q[1];
rz(-0.52430341) q[2];
sx q[2];
rz(-1.9881762) q[2];
sx q[2];
rz(0.5988754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6819671) q[1];
sx q[1];
rz(-0.29469583) q[1];
sx q[1];
rz(1.9035643) q[1];
rz(-pi) q[2];
rz(-0.21027474) q[3];
sx q[3];
rz(-1.9326903) q[3];
sx q[3];
rz(2.0023605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(-2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(-2.1797144) q[0];
rz(-0.23421639) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(-2.403517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9316406) q[0];
sx q[0];
rz(-2.6745213) q[0];
sx q[0];
rz(0.54801236) q[0];
x q[1];
rz(2.0162651) q[2];
sx q[2];
rz(-0.98410749) q[2];
sx q[2];
rz(2.4054804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9752561) q[1];
sx q[1];
rz(-1.8615906) q[1];
sx q[1];
rz(1.1329805) q[1];
x q[2];
rz(1.290019) q[3];
sx q[3];
rz(-2.2178136) q[3];
sx q[3];
rz(3.0386277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6568079) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(-2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.242908) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(-0.86724487) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(1.0702466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8673422) q[0];
sx q[0];
rz(-2.946642) q[0];
sx q[0];
rz(-2.250953) q[0];
rz(-pi) q[1];
rz(0.82735586) q[2];
sx q[2];
rz(-2.1340886) q[2];
sx q[2];
rz(2.2652596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5336736) q[1];
sx q[1];
rz(-1.1864788) q[1];
sx q[1];
rz(-0.55748765) q[1];
x q[2];
rz(-0.66859122) q[3];
sx q[3];
rz(-2.9644659) q[3];
sx q[3];
rz(1.9884381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.9099859) q[2];
rz(-1.6407137) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-2.7267743) q[1];
sx q[1];
rz(-1.7638313) q[1];
sx q[1];
rz(1.5324963) q[1];
rz(3.0267486) q[2];
sx q[2];
rz(-0.44841246) q[2];
sx q[2];
rz(2.8186225) q[2];
rz(0.98607705) q[3];
sx q[3];
rz(-0.90860962) q[3];
sx q[3];
rz(-1.2398401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
