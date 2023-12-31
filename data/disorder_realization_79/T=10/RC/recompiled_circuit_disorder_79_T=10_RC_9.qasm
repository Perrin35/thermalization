OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642759) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(0.067406128) q[0];
x q[1];
rz(-2.9002951) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(1.8482006) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.623917) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(-0.6448402) q[1];
x q[2];
rz(-1.0055553) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(-0.934787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(-2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-2.0400955) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0365021) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(-0.0037395517) q[0];
rz(-pi) q[1];
rz(2.0287839) q[2];
sx q[2];
rz(-0.53456842) q[2];
sx q[2];
rz(2.2528258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.022545594) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(2.4382298) q[1];
rz(-pi) q[2];
x q[2];
rz(2.39605) q[3];
sx q[3];
rz(-1.3127483) q[3];
sx q[3];
rz(-0.52449709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24414177) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(-1.1600526) q[0];
rz(-pi) q[1];
rz(0.86513743) q[2];
sx q[2];
rz(-2.1102998) q[2];
sx q[2];
rz(1.7973961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.264818) q[1];
sx q[1];
rz(-0.41300981) q[1];
sx q[1];
rz(-3.1289711) q[1];
x q[2];
rz(-0.86977264) q[3];
sx q[3];
rz(-1.8043892) q[3];
sx q[3];
rz(-1.9529395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67733726) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1151428) q[0];
sx q[0];
rz(-2.0949674) q[0];
sx q[0];
rz(1.3036222) q[0];
x q[1];
rz(-1.9374574) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(-1.0166849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81740582) q[1];
sx q[1];
rz(-1.4293912) q[1];
sx q[1];
rz(-1.6721339) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5303909) q[3];
sx q[3];
rz(-1.2308321) q[3];
sx q[3];
rz(-1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(0.51182169) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.408668) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(2.1599105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3571346) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(-2.5255192) q[0];
rz(-pi) q[1];
rz(0.9972516) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(-2.2821102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86451929) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(-2.8237052) q[1];
rz(-pi) q[2];
rz(-1.6793628) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(-1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.813252) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-2.7691675) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344791) q[0];
sx q[0];
rz(-1.665859) q[0];
sx q[0];
rz(1.6001742) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0566453) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(-1.2598318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2416934) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(2.4844869) q[1];
x q[2];
rz(-1.1214439) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.526028) q[0];
sx q[0];
rz(-3.0897339) q[0];
sx q[0];
rz(-1.9066914) q[0];
rz(-pi) q[1];
rz(2.8924106) q[2];
sx q[2];
rz(-1.8998002) q[2];
sx q[2];
rz(2.5788384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6474825) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(2.8857735) q[1];
rz(-0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(-2.2198026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9471654) q[0];
sx q[0];
rz(-0.88473407) q[0];
sx q[0];
rz(-1.3961193) q[0];
x q[1];
rz(1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(-1.2459754) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7505155) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-0.68117546) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17955762) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(-1.8732656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.4484891) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(-2.6152339) q[0];
rz(-pi) q[1];
rz(-0.87551261) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(-1.8347486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.2320227) q[1];
sx q[1];
rz(-2.4376274) q[1];
rz(0.67919517) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991853) q[0];
sx q[0];
rz(-2.9710794) q[0];
sx q[0];
rz(1.2810345) q[0];
x q[1];
rz(0.54729692) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(1.0123569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9709819) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(1.7112205) q[1];
x q[2];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(-2.836543) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(-0.35076326) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
