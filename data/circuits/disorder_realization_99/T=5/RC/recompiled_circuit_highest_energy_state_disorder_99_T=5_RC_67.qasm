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
rz(-1.7412269) q[0];
sx q[0];
rz(-2.9018612) q[0];
sx q[0];
rz(2.8170407) q[0];
rz(2.1895154) q[1];
sx q[1];
rz(-2.9157186) q[1];
sx q[1];
rz(1.8332551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75057331) q[0];
sx q[0];
rz(-0.68882525) q[0];
sx q[0];
rz(-0.21521615) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.657349) q[2];
sx q[2];
rz(-1.6855557) q[2];
sx q[2];
rz(1.6387303) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0066915) q[1];
sx q[1];
rz(-0.54062245) q[1];
sx q[1];
rz(-1.0512133) q[1];
rz(-pi) q[2];
rz(-1.4300554) q[3];
sx q[3];
rz(-2.2651197) q[3];
sx q[3];
rz(-1.3545881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7294881) q[2];
sx q[2];
rz(-1.2594014) q[2];
sx q[2];
rz(0.35120249) q[2];
rz(1.8515733) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(1.2235519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(2.2112041) q[0];
rz(1.8366086) q[1];
sx q[1];
rz(-2.3001859) q[1];
sx q[1];
rz(-2.5321541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6822426) q[0];
sx q[0];
rz(-1.1455904) q[0];
sx q[0];
rz(0.11328477) q[0];
rz(-pi) q[1];
rz(3.0514293) q[2];
sx q[2];
rz(-1.2740286) q[2];
sx q[2];
rz(-2.0903843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6516782) q[1];
sx q[1];
rz(-2.6082572) q[1];
sx q[1];
rz(2.3467605) q[1];
x q[2];
rz(1.3136765) q[3];
sx q[3];
rz(-0.58565307) q[3];
sx q[3];
rz(2.3975092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.096752) q[2];
sx q[2];
rz(-2.1614306) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(0.52037248) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(0.61409942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-2.4245872) q[0];
sx q[0];
rz(0.30211788) q[0];
rz(2.865454) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(-1.709323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3870704) q[0];
sx q[0];
rz(-1.1898479) q[0];
sx q[0];
rz(0.98636287) q[0];
rz(0.37636855) q[2];
sx q[2];
rz(-2.512815) q[2];
sx q[2];
rz(2.9719458) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36784962) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(1.2817205) q[1];
x q[2];
rz(0.071050771) q[3];
sx q[3];
rz(-1.1059575) q[3];
sx q[3];
rz(2.0185061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3064731) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.845537) q[2];
rz(2.6895788) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(-0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52962676) q[0];
sx q[0];
rz(-1.9744248) q[0];
sx q[0];
rz(-1.3324598) q[0];
rz(-1.0491071) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(1.5886935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4047239) q[0];
sx q[0];
rz(-1.8717614) q[0];
sx q[0];
rz(-2.7851581) q[0];
rz(-pi) q[1];
x q[1];
rz(0.013601549) q[2];
sx q[2];
rz(-2.0715393) q[2];
sx q[2];
rz(0.73464962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9178333) q[1];
sx q[1];
rz(-0.67549113) q[1];
sx q[1];
rz(1.2509327) q[1];
rz(0.53557204) q[3];
sx q[3];
rz(-2.2102997) q[3];
sx q[3];
rz(-2.4909693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4730452) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(2.891053) q[2];
rz(2.1446832) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(2.7610049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32960358) q[0];
sx q[0];
rz(-0.488509) q[0];
sx q[0];
rz(-1.3478152) q[0];
rz(2.7373121) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(-1.1291198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4816138) q[0];
sx q[0];
rz(-0.82757512) q[0];
sx q[0];
rz(-1.0218234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2329699) q[2];
sx q[2];
rz(-0.69402908) q[2];
sx q[2];
rz(0.90210669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2426274) q[1];
sx q[1];
rz(-2.4614442) q[1];
sx q[1];
rz(1.4037474) q[1];
rz(2.8267838) q[3];
sx q[3];
rz(-1.6757351) q[3];
sx q[3];
rz(-0.85280692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79444844) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(1.244119) q[2];
rz(1.0602903) q[3];
sx q[3];
rz(-0.62842193) q[3];
sx q[3];
rz(-0.30317831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4853972) q[0];
sx q[0];
rz(-0.10529101) q[0];
sx q[0];
rz(0.42718497) q[0];
rz(0.034491388) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(-3.131386) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611864) q[0];
sx q[0];
rz(-1.5690593) q[0];
sx q[0];
rz(-0.0023567452) q[0];
x q[1];
rz(1.7882657) q[2];
sx q[2];
rz(-0.64369338) q[2];
sx q[2];
rz(-2.3967607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.806655) q[1];
sx q[1];
rz(-0.68938556) q[1];
sx q[1];
rz(1.4695808) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0766076) q[3];
sx q[3];
rz(-2.0772604) q[3];
sx q[3];
rz(-2.9714366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5222142) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(0.79356066) q[2];
rz(-0.86269745) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4392387) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(-2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(2.8394707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7493499) q[0];
sx q[0];
rz(-1.2686522) q[0];
sx q[0];
rz(0.35373117) q[0];
x q[1];
rz(-2.5100915) q[2];
sx q[2];
rz(-0.69984791) q[2];
sx q[2];
rz(-1.8818784) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4430284) q[1];
sx q[1];
rz(-2.5686567) q[1];
sx q[1];
rz(-2.7354476) q[1];
rz(-pi) q[2];
rz(-2.4551943) q[3];
sx q[3];
rz(-1.0960311) q[3];
sx q[3];
rz(-1.7976185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2332396) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(1.3471777) q[2];
rz(1.2498648) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(-0.10716042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(0.10840848) q[0];
rz(-2.0938865) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7664095) q[0];
sx q[0];
rz(-1.0489223) q[0];
sx q[0];
rz(-0.96065582) q[0];
x q[1];
rz(-3.1347796) q[2];
sx q[2];
rz(-0.76131135) q[2];
sx q[2];
rz(-2.6369264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26886235) q[1];
sx q[1];
rz(-0.62727562) q[1];
sx q[1];
rz(2.4226047) q[1];
rz(-pi) q[2];
rz(0.15411994) q[3];
sx q[3];
rz(-0.74609038) q[3];
sx q[3];
rz(2.8733159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5041647) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(-0.49312433) q[3];
sx q[3];
rz(-0.92032856) q[3];
sx q[3];
rz(1.5639308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2757932) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(0.42617282) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0539249) q[0];
sx q[0];
rz(-3.1172981) q[0];
sx q[0];
rz(0.75456516) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2913646) q[2];
sx q[2];
rz(-2.2022708) q[2];
sx q[2];
rz(2.9852418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81282367) q[1];
sx q[1];
rz(-1.6881094) q[1];
sx q[1];
rz(-1.3552865) q[1];
rz(-pi) q[2];
rz(2.4091899) q[3];
sx q[3];
rz(-2.1653145) q[3];
sx q[3];
rz(-1.5791095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49491945) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.6152629) q[3];
sx q[3];
rz(1.0927965) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786355) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(0.92078513) q[0];
rz(0.50865632) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(-2.7210534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79179278) q[0];
sx q[0];
rz(-2.0689488) q[0];
sx q[0];
rz(2.34116) q[0];
x q[1];
rz(-0.24193544) q[2];
sx q[2];
rz(-0.95429776) q[2];
sx q[2];
rz(-0.37504196) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.072445446) q[1];
sx q[1];
rz(-2.0305567) q[1];
sx q[1];
rz(-0.88847499) q[1];
x q[2];
rz(1.693142) q[3];
sx q[3];
rz(-0.69465717) q[3];
sx q[3];
rz(2.8784424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.73021182) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(0.31663695) q[2];
rz(0.5591048) q[3];
sx q[3];
rz(-0.63566256) q[3];
sx q[3];
rz(-2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6560787) q[0];
sx q[0];
rz(-1.0777127) q[0];
sx q[0];
rz(1.0832473) q[0];
rz(-2.6651233) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(1.6638577) q[2];
sx q[2];
rz(-0.75329594) q[2];
sx q[2];
rz(2.6483173) q[2];
rz(-0.093802916) q[3];
sx q[3];
rz(-1.9524912) q[3];
sx q[3];
rz(-2.5759202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
