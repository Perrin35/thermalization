OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(-2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7919851) q[0];
sx q[0];
rz(-1.926593) q[0];
sx q[0];
rz(-0.7138568) q[0];
x q[1];
rz(-0.77180441) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(0.31847218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3633903) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.1671288) q[1];
rz(-2.8291679) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-0.47505501) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99012016) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-2.5734076) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3677164) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(-2.3351923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2270826) q[1];
sx q[1];
rz(-0.40725476) q[1];
sx q[1];
rz(1.6362908) q[1];
x q[2];
rz(-0.81289566) q[3];
sx q[3];
rz(-2.025369) q[3];
sx q[3];
rz(-1.2113435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.976982) q[0];
sx q[0];
rz(-1.6342388) q[0];
sx q[0];
rz(-2.001686) q[0];
rz(1.8884044) q[2];
sx q[2];
rz(-2.3177958) q[2];
sx q[2];
rz(-0.93712805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14904505) q[1];
sx q[1];
rz(-0.48679513) q[1];
sx q[1];
rz(0.64736127) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40804789) q[3];
sx q[3];
rz(-1.9983665) q[3];
sx q[3];
rz(-1.6392631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1283778) q[0];
sx q[0];
rz(-0.88071874) q[0];
sx q[0];
rz(0.31353686) q[0];
rz(-pi) q[1];
rz(1.5393125) q[2];
sx q[2];
rz(-1.5891979) q[2];
sx q[2];
rz(-3.1061663) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4094761) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(-0.17318053) q[1];
rz(-pi) q[2];
rz(-3.0691931) q[3];
sx q[3];
rz(-0.94745938) q[3];
sx q[3];
rz(0.49923957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-0.0030227946) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057102324) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(-1.3296207) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51909165) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(1.5965243) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0121213) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(0.71838897) q[1];
x q[2];
rz(-0.32096433) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370699) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(-1.5674595) q[0];
rz(-pi) q[1];
rz(1.6476829) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(-2.4424057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2885292) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(-2.2437375) q[1];
rz(0.45913978) q[3];
sx q[3];
rz(-0.79677478) q[3];
sx q[3];
rz(-1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60550624) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-0.58553186) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(0.62430635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3793959) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(0.72797738) q[0];
x q[1];
rz(-2.675266) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(-1.8898659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.045947) q[1];
sx q[1];
rz(-1.026517) q[1];
sx q[1];
rz(0.73927684) q[1];
x q[2];
rz(-0.47847139) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(-2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(3.0016622) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77057225) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(1.8380941) q[0];
x q[1];
rz(-0.89739563) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(-2.430254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.8433246) q[1];
sx q[1];
rz(-2.6998991) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(-3.0144514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.3508266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5577561) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(-1.3895967) q[0];
rz(1.9753014) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(1.5750969) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.13547922) q[1];
sx q[1];
rz(-2.0524426) q[1];
sx q[1];
rz(2.3856132) q[1];
rz(-pi) q[2];
rz(-0.75307122) q[3];
sx q[3];
rz(-1.3967447) q[3];
sx q[3];
rz(-2.2032602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7231862) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(1.5002804) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17391275) q[1];
sx q[1];
rz(-2.850107) q[1];
sx q[1];
rz(0.12853865) q[1];
x q[2];
rz(0.54639001) q[3];
sx q[3];
rz(-2.8791109) q[3];
sx q[3];
rz(-0.7149834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(0.004301087) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.2658723) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(0.98106445) q[3];
sx q[3];
rz(-0.55064252) q[3];
sx q[3];
rz(-2.8513089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
