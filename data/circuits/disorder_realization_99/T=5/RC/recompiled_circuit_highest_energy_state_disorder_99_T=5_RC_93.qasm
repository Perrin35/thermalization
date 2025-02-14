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
rz(-1.3083375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65303923) q[0];
sx q[0];
rz(-1.7069611) q[0];
sx q[0];
rz(0.67739886) q[0];
x q[1];
rz(2.657349) q[2];
sx q[2];
rz(-1.6855557) q[2];
sx q[2];
rz(1.5028624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89186016) q[1];
sx q[1];
rz(-1.3123871) q[1];
sx q[1];
rz(-2.051146) q[1];
rz(-0.69921391) q[3];
sx q[3];
rz(-1.4627856) q[3];
sx q[3];
rz(2.8349769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41210458) q[2];
sx q[2];
rz(-1.2594014) q[2];
sx q[2];
rz(0.35120249) q[2];
rz(-1.8515733) q[3];
sx q[3];
rz(-2.0100287) q[3];
sx q[3];
rz(-1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38607645) q[0];
sx q[0];
rz(-1.9367243) q[0];
sx q[0];
rz(-0.93038857) q[0];
rz(-1.8366086) q[1];
sx q[1];
rz(-2.3001859) q[1];
sx q[1];
rz(-0.60943857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.45935) q[0];
sx q[0];
rz(-1.1455904) q[0];
sx q[0];
rz(3.0283079) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.090163354) q[2];
sx q[2];
rz(-1.8675641) q[2];
sx q[2];
rz(2.0903843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7623065) q[1];
sx q[1];
rz(-1.9348782) q[1];
sx q[1];
rz(-1.9696139) q[1];
rz(-pi) q[2];
rz(2.1411544) q[3];
sx q[3];
rz(-1.7118239) q[3];
sx q[3];
rz(1.0424249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(-0.074782221) q[2];
rz(0.52037248) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(2.5274932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60270131) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-2.8394748) q[0];
rz(-2.865454) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(-1.4322697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3870704) q[0];
sx q[0];
rz(-1.9517448) q[0];
sx q[0];
rz(-2.1552298) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37636855) q[2];
sx q[2];
rz(-2.512815) q[2];
sx q[2];
rz(-0.16964682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70967445) q[1];
sx q[1];
rz(-2.1248105) q[1];
sx q[1];
rz(-2.9565503) q[1];
rz(-1.4301705) q[3];
sx q[3];
rz(-0.46984497) q[3];
sx q[3];
rz(-1.8610561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8351195) q[2];
sx q[2];
rz(-2.6758631) q[2];
sx q[2];
rz(1.845537) q[2];
rz(0.45201388) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.1671678) q[0];
sx q[0];
rz(-1.8091328) q[0];
rz(-1.0491071) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5528991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73686872) q[0];
sx q[0];
rz(-1.2698313) q[0];
sx q[0];
rz(-2.7851581) q[0];
rz(1.5459486) q[2];
sx q[2];
rz(-0.50091195) q[2];
sx q[2];
rz(2.3786169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59999787) q[1];
sx q[1];
rz(-1.3728956) q[1];
sx q[1];
rz(0.92055121) q[1];
rz(0.96943198) q[3];
sx q[3];
rz(-0.80925377) q[3];
sx q[3];
rz(-1.4331897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4730452) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(0.2505396) q[2];
rz(0.9969095) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(2.7610049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
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
rz(-0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.7937775) q[0];
rz(-2.7373121) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(1.1291198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6599789) q[0];
sx q[0];
rz(-0.82757512) q[0];
sx q[0];
rz(-2.1197693) q[0];
x q[1];
rz(0.90520827) q[2];
sx q[2];
rz(-1.357175) q[2];
sx q[2];
rz(-0.40494949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89896527) q[1];
sx q[1];
rz(-2.4614442) q[1];
sx q[1];
rz(-1.4037474) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8137098) q[3];
sx q[3];
rz(-2.8103069) q[3];
sx q[3];
rz(1.0292018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3471442) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(1.8974737) q[2];
rz(1.0602903) q[3];
sx q[3];
rz(-0.62842193) q[3];
sx q[3];
rz(-0.30317831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(-2.7144077) q[0];
rz(0.034491388) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(-0.010206612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804062) q[0];
sx q[0];
rz(-1.5690593) q[0];
sx q[0];
rz(0.0023567452) q[0];
x q[1];
rz(-0.16049196) q[2];
sx q[2];
rz(-0.94466034) q[2];
sx q[2];
rz(-2.6662835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.806655) q[1];
sx q[1];
rz(-0.68938556) q[1];
sx q[1];
rz(1.6720119) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7179435) q[3];
sx q[3];
rz(-2.4417851) q[3];
sx q[3];
rz(1.0216624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61937845) q[2];
sx q[2];
rz(-2.7284315) q[2];
sx q[2];
rz(2.348032) q[2];
rz(-0.86269745) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4392387) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(-2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(-0.30212197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0724807) q[0];
sx q[0];
rz(-1.2337419) q[0];
sx q[0];
rz(-1.891579) q[0];
rz(-pi) q[1];
rz(0.63150117) q[2];
sx q[2];
rz(-2.4417447) q[2];
sx q[2];
rz(1.8818784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4430284) q[1];
sx q[1];
rz(-0.57293597) q[1];
sx q[1];
rz(-0.40614508) q[1];
rz(-pi) q[2];
rz(-0.68639836) q[3];
sx q[3];
rz(-1.0960311) q[3];
sx q[3];
rz(1.7976185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2332396) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(1.7944149) q[2];
rz(1.8917278) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(-3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5938479) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(0.10840848) q[0];
rz(-2.0938865) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.326556) q[0];
sx q[0];
rz(-2.360965) q[0];
sx q[0];
rz(-2.3579979) q[0];
x q[1];
rz(-0.0068130612) q[2];
sx q[2];
rz(-2.3802813) q[2];
sx q[2];
rz(0.50466621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.093204) q[1];
sx q[1];
rz(-1.1133514) q[1];
sx q[1];
rz(2.0162575) q[1];
rz(2.9874727) q[3];
sx q[3];
rz(-0.74609038) q[3];
sx q[3];
rz(-2.8733159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5041647) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(-2.6484683) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(-1.5776618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86579943) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(0.42617282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0539249) q[0];
sx q[0];
rz(-0.024294596) q[0];
sx q[0];
rz(2.3870275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4075899) q[2];
sx q[2];
rz(-2.2224769) q[2];
sx q[2];
rz(2.3190448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.328769) q[1];
sx q[1];
rz(-1.4534833) q[1];
sx q[1];
rz(1.7863062) q[1];
rz(2.3087193) q[3];
sx q[3];
rz(-2.1578721) q[3];
sx q[3];
rz(-0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49491945) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(1.0529998) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(-1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786355) q[0];
sx q[0];
rz(-2.9243922) q[0];
sx q[0];
rz(2.2208075) q[0];
rz(2.6329363) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(-2.7210534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9284436) q[0];
sx q[0];
rz(-0.9125114) q[0];
sx q[0];
rz(2.4930605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24193544) q[2];
sx q[2];
rz(-0.95429776) q[2];
sx q[2];
rz(0.37504196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0691472) q[1];
sx q[1];
rz(-1.1110359) q[1];
sx q[1];
rz(-0.88847499) q[1];
rz(3.0402571) q[3];
sx q[3];
rz(-2.2592415) q[3];
sx q[3];
rz(-2.7197414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(2.8249557) q[2];
rz(2.5824879) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(0.81812286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.485514) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(0.47646933) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(1.477735) q[2];
sx q[2];
rz(-2.3882967) q[2];
sx q[2];
rz(-0.49327539) q[2];
rz(1.3415402) q[3];
sx q[3];
rz(-0.39250249) q[3];
sx q[3];
rz(-2.8233118) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
