OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95028967) q[0];
sx q[0];
rz(-0.27009717) q[0];
sx q[0];
rz(0.88859963) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5262573) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(-2.944817) q[0];
x q[1];
rz(-1.0384473) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(2.7388245) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7270679) q[1];
sx q[1];
rz(-1.5804075) q[1];
sx q[1];
rz(-2.8060786) q[1];
rz(-2.9167261) q[3];
sx q[3];
rz(-1.2149842) q[3];
sx q[3];
rz(1.0575305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8218653) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-2.3925609) q[2];
rz(-2.9876515) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(0.099893959) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(-2.7144576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7593133) q[0];
sx q[0];
rz(-2.1370892) q[0];
sx q[0];
rz(1.1575538) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92556503) q[2];
sx q[2];
rz(-1.7665909) q[2];
sx q[2];
rz(-2.2044971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1041873) q[1];
sx q[1];
rz(-1.7113285) q[1];
sx q[1];
rz(1.2804081) q[1];
rz(1.8187567) q[3];
sx q[3];
rz(-0.21020791) q[3];
sx q[3];
rz(1.5402286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(1.3986577) q[2];
rz(-2.9178197) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021521213) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(1.9728164) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(0.57560903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110014) q[0];
sx q[0];
rz(-1.5819893) q[0];
sx q[0];
rz(-0.63909792) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.017150684) q[2];
sx q[2];
rz(-2.2313925) q[2];
sx q[2];
rz(-1.9018381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.528426) q[1];
sx q[1];
rz(-1.3537145) q[1];
sx q[1];
rz(-1.0124799) q[1];
rz(0.024954114) q[3];
sx q[3];
rz(-2.2234699) q[3];
sx q[3];
rz(-0.11634532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(2.1843145) q[2];
rz(0.57473985) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.15605536) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(-1.2611457) q[0];
rz(0.082322923) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(-1.7877158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533723) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(-1.3277256) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8026695) q[2];
sx q[2];
rz(-2.6312575) q[2];
sx q[2];
rz(-1.5477382) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8396047) q[1];
sx q[1];
rz(-2.0652373) q[1];
sx q[1];
rz(1.5093263) q[1];
rz(-pi) q[2];
rz(-2.2500751) q[3];
sx q[3];
rz(-1.1783021) q[3];
sx q[3];
rz(-0.52549997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42133078) q[2];
sx q[2];
rz(-0.91840863) q[2];
sx q[2];
rz(-1.8850231) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(1.4168581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9067136) q[0];
sx q[0];
rz(-0.32787927) q[0];
sx q[0];
rz(1.0176786) q[0];
x q[1];
rz(-1.896865) q[2];
sx q[2];
rz(-0.89674258) q[2];
sx q[2];
rz(-1.3671966) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38450634) q[1];
sx q[1];
rz(-1.0913186) q[1];
sx q[1];
rz(-2.9887385) q[1];
rz(-pi) q[2];
rz(1.9958901) q[3];
sx q[3];
rz(-1.7853961) q[3];
sx q[3];
rz(0.44073013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(-2.055638) q[2];
rz(0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74349657) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(-0.79291517) q[0];
rz(-0.74053699) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(-0.83121306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3073458) q[0];
sx q[0];
rz(-1.1957268) q[0];
sx q[0];
rz(-2.9177279) q[0];
rz(0.0074499091) q[2];
sx q[2];
rz(-1.8478571) q[2];
sx q[2];
rz(1.7817093) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86262396) q[1];
sx q[1];
rz(-1.875307) q[1];
sx q[1];
rz(-2.3079951) q[1];
x q[2];
rz(1.100086) q[3];
sx q[3];
rz(-1.5522955) q[3];
sx q[3];
rz(0.61533606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0703766) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(1.1191204) q[2];
rz(1.8708771) q[3];
sx q[3];
rz(-1.1071353) q[3];
sx q[3];
rz(-1.3814111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041466) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(-1.5420089) q[0];
rz(2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-2.9687845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240066) q[0];
sx q[0];
rz(-2.3476971) q[0];
sx q[0];
rz(-0.85063299) q[0];
rz(1.2020993) q[2];
sx q[2];
rz(-0.92075821) q[2];
sx q[2];
rz(1.4024613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.648293) q[1];
sx q[1];
rz(-1.9956335) q[1];
sx q[1];
rz(-2.8684435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5739848) q[3];
sx q[3];
rz(-2.2877573) q[3];
sx q[3];
rz(-0.53138083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66199866) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(-0.57502037) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32790312) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(1.188259) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(-1.9421008) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48391438) q[0];
sx q[0];
rz(-1.3356613) q[0];
sx q[0];
rz(1.1170618) q[0];
x q[1];
rz(-0.052414465) q[2];
sx q[2];
rz(-1.0804861) q[2];
sx q[2];
rz(2.3971967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6024692) q[1];
sx q[1];
rz(-1.1014465) q[1];
sx q[1];
rz(2.5825809) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0301129) q[3];
sx q[3];
rz(-0.62255854) q[3];
sx q[3];
rz(-0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1477995) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(2.5177453) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(2.0417716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(2.4483335) q[0];
rz(-0.26501003) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.6185435) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361621) q[0];
sx q[0];
rz(-0.79728617) q[0];
sx q[0];
rz(-1.5862341) q[0];
rz(0.9829282) q[2];
sx q[2];
rz(-1.0318021) q[2];
sx q[2];
rz(1.3608152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6553626) q[1];
sx q[1];
rz(-1.2312447) q[1];
sx q[1];
rz(1.5454253) q[1];
rz(0.98966815) q[3];
sx q[3];
rz(-1.5115941) q[3];
sx q[3];
rz(-2.7855685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8388464) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(-1.4576853) q[2];
rz(-2.1863106) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(0.48267522) q[0];
rz(0.92974281) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(-2.7453056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38448725) q[0];
sx q[0];
rz(-1.1570017) q[0];
sx q[0];
rz(-1.1767469) q[0];
rz(-pi) q[1];
rz(-2.9626289) q[2];
sx q[2];
rz(-1.4216627) q[2];
sx q[2];
rz(-0.75527945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63782802) q[1];
sx q[1];
rz(-2.2866268) q[1];
sx q[1];
rz(2.9010335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72209218) q[3];
sx q[3];
rz(-2.6240337) q[3];
sx q[3];
rz(-0.90921569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3489909) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(-1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1031716) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(1.5255899) q[2];
sx q[2];
rz(-2.2635985) q[2];
sx q[2];
rz(0.23512693) q[2];
rz(-0.38402186) q[3];
sx q[3];
rz(-2.1441318) q[3];
sx q[3];
rz(-1.287788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
