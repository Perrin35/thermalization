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
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(2.10973) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2159696) q[0];
sx q[0];
rz(-1.5711938) q[0];
sx q[0];
rz(0.0041603869) q[0];
rz(-pi) q[1];
rz(-2.5450942) q[2];
sx q[2];
rz(-2.7257127) q[2];
sx q[2];
rz(-2.2592827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84656395) q[1];
sx q[1];
rz(-1.16638) q[1];
sx q[1];
rz(1.1961028) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9239465) q[3];
sx q[3];
rz(-0.34322327) q[3];
sx q[3];
rz(-1.9453334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2304113) q[2];
sx q[2];
rz(-1.958467) q[2];
sx q[2];
rz(0.46768701) q[2];
rz(2.6905401) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(-2.8635645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69538799) q[0];
sx q[0];
rz(-0.22286335) q[0];
sx q[0];
rz(-0.15431246) q[0];
rz(0.34126869) q[1];
sx q[1];
rz(-0.35366615) q[1];
sx q[1];
rz(-0.42010677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9720917) q[0];
sx q[0];
rz(-1.7636714) q[0];
sx q[0];
rz(-0.64298301) q[0];
rz(-pi) q[1];
rz(-2.2809199) q[2];
sx q[2];
rz(-1.6063074) q[2];
sx q[2];
rz(-1.4461609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5771835) q[1];
sx q[1];
rz(-1.100005) q[1];
sx q[1];
rz(-1.5692577) q[1];
x q[2];
rz(-0.16815925) q[3];
sx q[3];
rz(-0.75440591) q[3];
sx q[3];
rz(1.2446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5306065) q[2];
sx q[2];
rz(-2.2425118) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(0.84345877) q[3];
sx q[3];
rz(-1.5777028) q[3];
sx q[3];
rz(2.5696866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17324363) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-0.20208836) q[0];
rz(2.659722) q[1];
sx q[1];
rz(-2.9402969) q[1];
sx q[1];
rz(-0.906382) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4545091) q[0];
sx q[0];
rz(-2.2453899) q[0];
sx q[0];
rz(2.1620511) q[0];
rz(-1.1163122) q[2];
sx q[2];
rz(-1.2589728) q[2];
sx q[2];
rz(0.21909595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6828595) q[1];
sx q[1];
rz(-2.5762229) q[1];
sx q[1];
rz(-0.302178) q[1];
rz(-pi) q[2];
rz(1.866147) q[3];
sx q[3];
rz(-1.4100299) q[3];
sx q[3];
rz(-2.8637342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0704982) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(-0.58007288) q[2];
rz(1.5813367) q[3];
sx q[3];
rz(-1.3908849) q[3];
sx q[3];
rz(1.4238547) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(3.0938003) q[0];
rz(1.8675249) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(0.045225708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3078318) q[0];
sx q[0];
rz(-2.2068496) q[0];
sx q[0];
rz(0.18524203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2061646) q[2];
sx q[2];
rz(-1.3029939) q[2];
sx q[2];
rz(-1.8647461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3056335) q[1];
sx q[1];
rz(-0.0063414185) q[1];
sx q[1];
rz(-2.0672227) q[1];
rz(2.2269656) q[3];
sx q[3];
rz(-0.65319971) q[3];
sx q[3];
rz(0.40756215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4101326) q[2];
sx q[2];
rz(-1.8803909) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(2.5951923) q[3];
sx q[3];
rz(-2.140464) q[3];
sx q[3];
rz(-0.31118292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6379717) q[0];
sx q[0];
rz(-1.4542955) q[0];
sx q[0];
rz(-1.7244435) q[0];
rz(-2.0218938) q[1];
sx q[1];
rz(-1.7599301) q[1];
sx q[1];
rz(-1.959257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45641023) q[0];
sx q[0];
rz(-1.5014002) q[0];
sx q[0];
rz(-2.7539537) q[0];
rz(-2.0242906) q[2];
sx q[2];
rz(-1.9879447) q[2];
sx q[2];
rz(-2.3309649) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2675543) q[1];
sx q[1];
rz(-0.8505162) q[1];
sx q[1];
rz(-1.6362067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47207802) q[3];
sx q[3];
rz(-2.3340477) q[3];
sx q[3];
rz(-0.20281747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22415367) q[2];
sx q[2];
rz(-0.76442337) q[2];
sx q[2];
rz(0.42382851) q[2];
rz(1.1118927) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(-2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424585) q[0];
sx q[0];
rz(-0.0722216) q[0];
sx q[0];
rz(2.1891731) q[0];
rz(0.77596387) q[1];
sx q[1];
rz(-0.9351848) q[1];
sx q[1];
rz(-1.5752569) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0512038) q[0];
sx q[0];
rz(-0.04310933) q[0];
sx q[0];
rz(-2.6531522) q[0];
x q[1];
rz(2.1774954) q[2];
sx q[2];
rz(-1.6025873) q[2];
sx q[2];
rz(-2.3111985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.414622) q[1];
sx q[1];
rz(-1.3515673) q[1];
sx q[1];
rz(2.5269233) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9569472) q[3];
sx q[3];
rz(-1.0870516) q[3];
sx q[3];
rz(2.8205706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39885193) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(2.3952132) q[2];
rz(2.0505203) q[3];
sx q[3];
rz(-1.7079791) q[3];
sx q[3];
rz(0.55238849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2090476) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(0.5433425) q[0];
rz(0.53687334) q[1];
sx q[1];
rz(-2.8165292) q[1];
sx q[1];
rz(3.0184025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7167599) q[0];
sx q[0];
rz(-0.93532978) q[0];
sx q[0];
rz(-0.038095564) q[0];
rz(-pi) q[1];
rz(-1.1511717) q[2];
sx q[2];
rz(-2.0438355) q[2];
sx q[2];
rz(1.9422188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.035288485) q[1];
sx q[1];
rz(-2.114944) q[1];
sx q[1];
rz(-0.72334163) q[1];
rz(-2.3513138) q[3];
sx q[3];
rz(-2.1690282) q[3];
sx q[3];
rz(1.3682531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3049551) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(0.48509625) q[2];
rz(1.9890316) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(-3.0936354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3845859) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(0.49569976) q[0];
rz(-0.063830201) q[1];
sx q[1];
rz(-2.3921236) q[1];
sx q[1];
rz(3.0705423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97489965) q[0];
sx q[0];
rz(-1.3314795) q[0];
sx q[0];
rz(-2.6684851) q[0];
rz(-1.7456648) q[2];
sx q[2];
rz(-2.1687458) q[2];
sx q[2];
rz(-1.9964465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20375511) q[1];
sx q[1];
rz(-1.3295638) q[1];
sx q[1];
rz(-0.15934847) q[1];
x q[2];
rz(-1.6556173) q[3];
sx q[3];
rz(-2.0165947) q[3];
sx q[3];
rz(-0.63157192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0157328) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(2.2802172) q[2];
rz(1.3043978) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(-1.5484352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9651589) q[0];
sx q[0];
rz(-0.49089828) q[0];
sx q[0];
rz(1.7023671) q[0];
rz(-3.0610415) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(1.390994) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4347067) q[0];
sx q[0];
rz(-1.9932184) q[0];
sx q[0];
rz(-2.1377853) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.020419196) q[2];
sx q[2];
rz(-1.2660742) q[2];
sx q[2];
rz(-2.8430276) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1999368) q[1];
sx q[1];
rz(-0.98683954) q[1];
sx q[1];
rz(-0.55364174) q[1];
x q[2];
rz(-2.3426314) q[3];
sx q[3];
rz(-0.50833251) q[3];
sx q[3];
rz(1.5768743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3535658) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(0.12727748) q[2];
rz(0.92653972) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(1.658879) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38630286) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(-1.5007098) q[0];
rz(-0.81037784) q[1];
sx q[1];
rz(-0.85226285) q[1];
sx q[1];
rz(-2.3694029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92896748) q[0];
sx q[0];
rz(-1.0662931) q[0];
sx q[0];
rz(0.8302159) q[0];
rz(-pi) q[1];
rz(0.29532642) q[2];
sx q[2];
rz(-1.775303) q[2];
sx q[2];
rz(-0.4819862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82622426) q[1];
sx q[1];
rz(-1.5824707) q[1];
sx q[1];
rz(-1.0213901) q[1];
x q[2];
rz(0.89086074) q[3];
sx q[3];
rz(-1.5308342) q[3];
sx q[3];
rz(-0.077629493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3810252) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(-0.62925657) q[2];
rz(-2.4106846) q[3];
sx q[3];
rz(-0.84354246) q[3];
sx q[3];
rz(-1.6985016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6998941) q[0];
sx q[0];
rz(-2.6581673) q[0];
sx q[0];
rz(2.7375258) q[0];
rz(-0.40456698) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(2.8597833) q[2];
sx q[2];
rz(-1.5055613) q[2];
sx q[2];
rz(0.67508634) q[2];
rz(-1.2734781) q[3];
sx q[3];
rz(-1.0549637) q[3];
sx q[3];
rz(-2.0217287) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
