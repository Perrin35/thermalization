OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(-2.9380517) q[0];
sx q[0];
rz(-1.8939053) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(-2.4434659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0188698) q[0];
sx q[0];
rz(-1.9779045) q[0];
sx q[0];
rz(-2.5886445) q[0];
rz(2.6883672) q[2];
sx q[2];
rz(-1.5648603) q[2];
sx q[2];
rz(0.56528948) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9351077) q[1];
sx q[1];
rz(-0.84151959) q[1];
sx q[1];
rz(0.40928276) q[1];
rz(-pi) q[2];
rz(-1.3157482) q[3];
sx q[3];
rz(-1.5849893) q[3];
sx q[3];
rz(-3.1007781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(0.87948925) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5010928) q[0];
sx q[0];
rz(-1.1859897) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(-2.1462006) q[1];
sx q[1];
rz(-1.6973015) q[1];
sx q[1];
rz(0.46517864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8841726) q[0];
sx q[0];
rz(-1.3883988) q[0];
sx q[0];
rz(2.718363) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7654573) q[2];
sx q[2];
rz(-0.97797624) q[2];
sx q[2];
rz(1.3786292) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77762802) q[1];
sx q[1];
rz(-2.2644357) q[1];
sx q[1];
rz(-0.65366807) q[1];
rz(1.3690794) q[3];
sx q[3];
rz(-1.8276102) q[3];
sx q[3];
rz(1.6981924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62749949) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(0.56525362) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(-2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8239215) q[0];
sx q[0];
rz(-1.7963821) q[0];
sx q[0];
rz(0.70096651) q[0];
rz(-0.40755454) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(1.0754546) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0757119) q[0];
sx q[0];
rz(-1.2383108) q[0];
sx q[0];
rz(-0.57384579) q[0];
rz(-pi) q[1];
rz(2.9161152) q[2];
sx q[2];
rz(-0.92712159) q[2];
sx q[2];
rz(0.89356092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6523931) q[1];
sx q[1];
rz(-1.4402188) q[1];
sx q[1];
rz(1.9741535) q[1];
rz(1.8355708) q[3];
sx q[3];
rz(-1.522101) q[3];
sx q[3];
rz(2.8935561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0740697) q[2];
sx q[2];
rz(-1.3347551) q[2];
sx q[2];
rz(1.734181) q[2];
rz(-2.6635026) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(2.5673229) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(2.0950914) q[0];
rz(2.8872755) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(-0.3124803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0346261) q[0];
sx q[0];
rz(-0.43145942) q[0];
sx q[0];
rz(-1.8851243) q[0];
rz(-0.66019571) q[2];
sx q[2];
rz(-1.1822299) q[2];
sx q[2];
rz(-1.6903433) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4832125) q[1];
sx q[1];
rz(-1.8129895) q[1];
sx q[1];
rz(-0.32434742) q[1];
rz(-1.0279827) q[3];
sx q[3];
rz(-0.58251721) q[3];
sx q[3];
rz(1.5912671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8326524) q[2];
sx q[2];
rz(-2.4231484) q[2];
sx q[2];
rz(2.9529875) q[2];
rz(-0.23379937) q[3];
sx q[3];
rz(-0.89582396) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69498649) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(-0.0083228668) q[0];
rz(-2.7024929) q[1];
sx q[1];
rz(-1.7726026) q[1];
sx q[1];
rz(-1.8956553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2089691) q[0];
sx q[0];
rz(-1.1198988) q[0];
sx q[0];
rz(-1.3572925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9925084) q[2];
sx q[2];
rz(-1.8051762) q[2];
sx q[2];
rz(2.5917376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.627927) q[1];
sx q[1];
rz(-2.6509319) q[1];
sx q[1];
rz(2.2918244) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2034791) q[3];
sx q[3];
rz(-2.4532336) q[3];
sx q[3];
rz(1.3753152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7370854) q[2];
sx q[2];
rz(-3.0948907) q[2];
sx q[2];
rz(1.5076293) q[2];
rz(-0.59219939) q[3];
sx q[3];
rz(-1.7630354) q[3];
sx q[3];
rz(2.4018905) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208039) q[0];
sx q[0];
rz(-1.692166) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(-1.9256915) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(-1.1908092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95071973) q[0];
sx q[0];
rz(-0.83774532) q[0];
sx q[0];
rz(-0.0035973408) q[0];
rz(-1.5120248) q[2];
sx q[2];
rz(-0.55355367) q[2];
sx q[2];
rz(0.89828459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3336181) q[1];
sx q[1];
rz(-2.3876973) q[1];
sx q[1];
rz(-3.1303166) q[1];
rz(-pi) q[2];
rz(-1.6281288) q[3];
sx q[3];
rz(-1.2964037) q[3];
sx q[3];
rz(-2.969034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33209458) q[2];
sx q[2];
rz(-1.3846493) q[2];
sx q[2];
rz(0.14796999) q[2];
rz(0.45251265) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(-2.7361659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1389403) q[0];
sx q[0];
rz(-1.2661221) q[0];
sx q[0];
rz(1.3462322) q[0];
rz(0.0070618709) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(0.5074358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.00822) q[0];
sx q[0];
rz(-0.79303654) q[0];
sx q[0];
rz(-0.80312463) q[0];
rz(-pi) q[1];
rz(1.1603786) q[2];
sx q[2];
rz(-0.37330356) q[2];
sx q[2];
rz(-2.6609535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9964357) q[1];
sx q[1];
rz(-1.5801589) q[1];
sx q[1];
rz(-0.3939751) q[1];
x q[2];
rz(-1.8431013) q[3];
sx q[3];
rz(-1.5617126) q[3];
sx q[3];
rz(-2.7177725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4355882) q[2];
sx q[2];
rz(-1.0554577) q[2];
sx q[2];
rz(-2.9662507) q[2];
rz(-0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(-2.6800938) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984556) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(1.3125516) q[0];
rz(0.10874272) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(-2.463602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4512166) q[0];
sx q[0];
rz(-1.548598) q[0];
sx q[0];
rz(-1.0772086) q[0];
rz(3.1197879) q[2];
sx q[2];
rz(-0.74912375) q[2];
sx q[2];
rz(-2.5072012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40174088) q[1];
sx q[1];
rz(-2.0900407) q[1];
sx q[1];
rz(0.3335764) q[1];
rz(-0.71700134) q[3];
sx q[3];
rz(-1.9878519) q[3];
sx q[3];
rz(-1.4456229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(-1.9847974) q[2];
rz(-0.50446883) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2061763) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(-0.61844283) q[0];
rz(2.2930324) q[1];
sx q[1];
rz(-0.51785523) q[1];
sx q[1];
rz(2.3458164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37695161) q[0];
sx q[0];
rz(-1.5466154) q[0];
sx q[0];
rz(2.0410246) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3168542) q[2];
sx q[2];
rz(-1.6055664) q[2];
sx q[2];
rz(0.29037133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2256945) q[1];
sx q[1];
rz(-1.8437598) q[1];
sx q[1];
rz(1.2429487) q[1];
x q[2];
rz(0.14988092) q[3];
sx q[3];
rz(-1.9871885) q[3];
sx q[3];
rz(-0.021377953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15647469) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(1.1953243) q[2];
rz(2.7965503) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(-0.8663469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7095551) q[0];
sx q[0];
rz(-0.23831743) q[0];
sx q[0];
rz(3.0639783) q[0];
rz(1.9084825) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(-0.79201039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7562078) q[0];
sx q[0];
rz(-1.0288219) q[0];
sx q[0];
rz(0.78800027) q[0];
rz(-pi) q[1];
rz(0.15969212) q[2];
sx q[2];
rz(-1.7241242) q[2];
sx q[2];
rz(3.0042269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92160901) q[1];
sx q[1];
rz(-1.7243544) q[1];
sx q[1];
rz(-1.1034225) q[1];
rz(-2.3378793) q[3];
sx q[3];
rz(-1.6559542) q[3];
sx q[3];
rz(-2.5309895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21976694) q[2];
sx q[2];
rz(-0.91138387) q[2];
sx q[2];
rz(2.0683973) q[2];
rz(2.8940767) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(0.016935067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76219227) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(2.2629867) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(0.20181561) q[2];
sx q[2];
rz(-2.0463662) q[2];
sx q[2];
rz(-3.0007888) q[2];
rz(0.84390072) q[3];
sx q[3];
rz(-1.6132728) q[3];
sx q[3];
rz(-2.1584395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
