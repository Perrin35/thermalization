OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7895176) q[0];
sx q[0];
rz(-1.2791469) q[0];
sx q[0];
rz(-1.350499) q[0];
x q[1];
rz(2.9169693) q[2];
sx q[2];
rz(-0.42806872) q[2];
sx q[2];
rz(-0.12878865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6341056) q[1];
sx q[1];
rz(-1.7934985) q[1];
sx q[1];
rz(2.1102064) q[1];
x q[2];
rz(0.033693245) q[3];
sx q[3];
rz(-1.9851306) q[3];
sx q[3];
rz(1.3106443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(0.12456482) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(1.7547866) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215866) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(0.63175732) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60123721) q[0];
sx q[0];
rz(-1.3155126) q[0];
sx q[0];
rz(-0.028098696) q[0];
rz(-pi) q[1];
rz(2.3089356) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(-0.51072272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1214952) q[1];
sx q[1];
rz(-2.0165682) q[1];
sx q[1];
rz(-0.14264588) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5242819) q[3];
sx q[3];
rz(-0.97913137) q[3];
sx q[3];
rz(-3.068963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.867884) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231409) q[0];
sx q[0];
rz(-1.8225192) q[0];
sx q[0];
rz(1.203042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4726228) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(1.836118) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2883403) q[1];
sx q[1];
rz(-0.77008343) q[1];
sx q[1];
rz(-2.6283162) q[1];
x q[2];
rz(2.418163) q[3];
sx q[3];
rz(-1.9577141) q[3];
sx q[3];
rz(-0.5667516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0597824) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(2.9147193) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(3.0674556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4081057) q[0];
sx q[0];
rz(-1.3164078) q[0];
sx q[0];
rz(-1.0803726) q[0];
rz(-1.9937236) q[2];
sx q[2];
rz(-2.4175156) q[2];
sx q[2];
rz(-0.99087447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6756235) q[1];
sx q[1];
rz(-0.33756653) q[1];
sx q[1];
rz(-1.2077431) q[1];
rz(0.62319237) q[3];
sx q[3];
rz(-0.29286256) q[3];
sx q[3];
rz(-2.5588536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6073109) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.779153) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(1.978925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8487932) q[0];
sx q[0];
rz(-1.4593562) q[0];
sx q[0];
rz(3.1116027) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0983724) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(0.098509468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3001544) q[1];
sx q[1];
rz(-1.9342124) q[1];
sx q[1];
rz(-0.20519786) q[1];
x q[2];
rz(-0.88697042) q[3];
sx q[3];
rz(-2.5081222) q[3];
sx q[3];
rz(-2.5630643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.4676771) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(2.8335559) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.98826) q[0];
sx q[0];
rz(-1.4446265) q[0];
sx q[0];
rz(1.6660965) q[0];
rz(-pi) q[1];
rz(-0.78139379) q[2];
sx q[2];
rz(-0.77971824) q[2];
sx q[2];
rz(1.7674854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8923556) q[1];
sx q[1];
rz(-0.71214572) q[1];
sx q[1];
rz(-1.874079) q[1];
rz(-pi) q[2];
rz(0.34835784) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(-2.0406046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(0.1299783) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(2.470509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13043159) q[0];
sx q[0];
rz(-1.5357619) q[0];
sx q[0];
rz(-0.053671562) q[0];
rz(1.7550049) q[2];
sx q[2];
rz(-1.5407908) q[2];
sx q[2];
rz(-1.9629994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56747251) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(1.1979539) q[1];
rz(-pi) q[2];
rz(-2.2460849) q[3];
sx q[3];
rz(-1.5338147) q[3];
sx q[3];
rz(-1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(-3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1639444) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(-2.7291765) q[0];
rz(-1.6917797) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(-1.1669881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3827688) q[0];
sx q[0];
rz(-2.1817657) q[0];
sx q[0];
rz(-2.8805034) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8896066) q[2];
sx q[2];
rz(-2.406771) q[2];
sx q[2];
rz(-0.97359818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.941277) q[1];
sx q[1];
rz(-1.3822767) q[1];
sx q[1];
rz(1.539906) q[1];
rz(-pi) q[2];
rz(-2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-0.92700672) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.4896726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324069) q[0];
sx q[0];
rz(-2.272185) q[0];
sx q[0];
rz(-2.5138445) q[0];
x q[1];
rz(-1.0987368) q[2];
sx q[2];
rz(-1.5257116) q[2];
sx q[2];
rz(-0.91143196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5786963) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(-1.3955411) q[1];
x q[2];
rz(-2.1065815) q[3];
sx q[3];
rz(-1.1318558) q[3];
sx q[3];
rz(2.5715695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.4302953) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(-1.3841217) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5883314) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(0.28840315) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518236) q[0];
sx q[0];
rz(-2.4252486) q[0];
sx q[0];
rz(-2.1728974) q[0];
rz(-2.8617919) q[2];
sx q[2];
rz(-0.24926148) q[2];
sx q[2];
rz(1.3485497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3210996) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(2.189373) q[1];
rz(0.23707323) q[3];
sx q[3];
rz(-2.5519538) q[3];
sx q[3];
rz(1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(-0.4521162) q[2];
sx q[2];
rz(-1.5099667) q[2];
sx q[2];
rz(-2.920334) q[2];
rz(1.0088624) q[3];
sx q[3];
rz(-1.6812134) q[3];
sx q[3];
rz(-2.5440661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
