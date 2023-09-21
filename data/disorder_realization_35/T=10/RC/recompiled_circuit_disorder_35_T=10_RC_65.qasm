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
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071872358) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(-1.1138492) q[0];
rz(2.3697882) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(0.31847218) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3633903) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.1671288) q[1];
x q[2];
rz(-1.6463046) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.422721) q[0];
sx q[0];
rz(-2.123721) q[0];
sx q[0];
rz(-1.3119112) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77387626) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-2.3351923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2270826) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(1.5053019) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.188835) q[3];
sx q[3];
rz(-0.85988322) q[3];
sx q[3];
rz(-0.074912138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(2.5684165) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1646106) q[0];
sx q[0];
rz(-1.6342388) q[0];
sx q[0];
rz(-2.001686) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77261749) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(0.85341838) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5830071) q[1];
sx q[1];
rz(-1.9531986) q[1];
sx q[1];
rz(1.8797727) q[1];
rz(-pi) q[2];
rz(-2.0315941) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(-1.213221) q[0];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.5359495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4421694) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(2.5073754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(-2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(0.80250424) q[3];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0844903) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(1.3296207) q[0];
rz(0.95894496) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(0.81629717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1369049) q[1];
sx q[1];
rz(-1.3740174) q[1];
sx q[1];
rz(-0.23006567) q[1];
rz(-pi) q[2];
rz(2.8206283) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(-1.9110514) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753871) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-3.1211981) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6476829) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(-2.4424057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8530635) q[1];
sx q[1];
rz(-1.001734) q[1];
sx q[1];
rz(-0.89785518) q[1];
rz(-pi) q[2];
rz(-0.45913978) q[3];
sx q[3];
rz(-2.3448179) q[3];
sx q[3];
rz(-1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(0.62430635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(1.8486345) q[0];
rz(0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(2.6401273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.034060409) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(0.88453102) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1950486) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873868) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(-1.9681853) q[0];
rz(-2.244197) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(-0.71133864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32313777) q[1];
sx q[1];
rz(-1.8433246) q[1];
sx q[1];
rz(2.6998991) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4017879) q[3];
sx q[3];
rz(-1.5825669) q[3];
sx q[3];
rz(-1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(-1.3895967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9753014) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(-1.5664958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8920039) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(0.65111098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3885214) q[3];
sx q[3];
rz(-1.744848) q[3];
sx q[3];
rz(2.2032602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0944259) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4184065) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(-1.5002804) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5200561) q[1];
sx q[1];
rz(-1.6076419) q[1];
sx q[1];
rz(0.28921339) q[1];
rz(2.91594) q[3];
sx q[3];
rz(-1.435558) q[3];
sx q[3];
rz(1.3868563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.5236241) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
rz(-0.32904939) q[3];
sx q[3];
rz(-2.0206738) q[3];
sx q[3];
rz(-2.1856367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];