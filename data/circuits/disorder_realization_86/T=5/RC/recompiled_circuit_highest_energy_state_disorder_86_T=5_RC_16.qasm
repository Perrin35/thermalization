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
rz(-1.4671833) q[0];
sx q[0];
rz(-1.6011342) q[0];
sx q[0];
rz(1.849548) q[0];
rz(2.0436824) q[1];
sx q[1];
rz(-0.75610375) q[1];
sx q[1];
rz(2.4296711) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24589892) q[0];
sx q[0];
rz(-2.993949) q[0];
sx q[0];
rz(-2.9013322) q[0];
x q[1];
rz(2.3825157) q[2];
sx q[2];
rz(-1.402194) q[2];
sx q[2];
rz(-1.027866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56150964) q[1];
sx q[1];
rz(-2.4697008) q[1];
sx q[1];
rz(-0.36999934) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51070945) q[3];
sx q[3];
rz(-1.0282955) q[3];
sx q[3];
rz(1.6245019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11382515) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(-0.85565957) q[2];
rz(1.2565819) q[3];
sx q[3];
rz(-0.62658566) q[3];
sx q[3];
rz(-1.1945061) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41334316) q[0];
sx q[0];
rz(-3.0276868) q[0];
sx q[0];
rz(2.2746427) q[0];
rz(1.4504245) q[1];
sx q[1];
rz(-1.9536641) q[1];
sx q[1];
rz(3.0303722) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635083) q[0];
sx q[0];
rz(-2.1136381) q[0];
sx q[0];
rz(-1.1397902) q[0];
x q[1];
rz(-0.84663518) q[2];
sx q[2];
rz(-0.66613448) q[2];
sx q[2];
rz(0.29633477) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73248374) q[1];
sx q[1];
rz(-0.95327158) q[1];
sx q[1];
rz(-2.3557622) q[1];
x q[2];
rz(1.7584959) q[3];
sx q[3];
rz(-0.71569136) q[3];
sx q[3];
rz(0.3970428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0337246) q[2];
sx q[2];
rz(-2.3794231) q[2];
sx q[2];
rz(2.6503358) q[2];
rz(-1.8852425) q[3];
sx q[3];
rz(-1.7137073) q[3];
sx q[3];
rz(0.45891941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1678109) q[0];
sx q[0];
rz(-1.7901006) q[0];
sx q[0];
rz(-0.0034045086) q[0];
rz(-0.34524521) q[1];
sx q[1];
rz(-0.41989851) q[1];
sx q[1];
rz(2.2664864) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21247877) q[0];
sx q[0];
rz(-1.0153595) q[0];
sx q[0];
rz(-1.4856797) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11981182) q[2];
sx q[2];
rz(-1.3917599) q[2];
sx q[2];
rz(-1.8869635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6236804) q[1];
sx q[1];
rz(-2.0225683) q[1];
sx q[1];
rz(-0.47924115) q[1];
rz(-pi) q[2];
x q[2];
rz(1.544041) q[3];
sx q[3];
rz(-2.3418276) q[3];
sx q[3];
rz(1.2497219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.379999) q[2];
sx q[2];
rz(-2.9271409) q[2];
sx q[2];
rz(-0.3271412) q[2];
rz(1.3867779) q[3];
sx q[3];
rz(-1.1969457) q[3];
sx q[3];
rz(3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4057994) q[0];
sx q[0];
rz(-1.0013591) q[0];
sx q[0];
rz(-1.557198) q[0];
rz(2.7994075) q[1];
sx q[1];
rz(-1.0289959) q[1];
sx q[1];
rz(0.4272961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2262242) q[0];
sx q[0];
rz(-1.6037723) q[0];
sx q[0];
rz(3.0874662) q[0];
rz(0.9924381) q[2];
sx q[2];
rz(-0.61884862) q[2];
sx q[2];
rz(0.98471314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72361833) q[1];
sx q[1];
rz(-0.86043859) q[1];
sx q[1];
rz(2.3478226) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8713636) q[3];
sx q[3];
rz(-1.0247314) q[3];
sx q[3];
rz(-0.48496839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7067318) q[2];
sx q[2];
rz(-0.92919436) q[2];
sx q[2];
rz(0.84563196) q[2];
rz(3.1137858) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(-1.4771247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.528275) q[0];
sx q[0];
rz(-1.0266101) q[0];
sx q[0];
rz(3.0897621) q[0];
rz(1.9822281) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(0.61028284) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51175115) q[0];
sx q[0];
rz(-1.399222) q[0];
sx q[0];
rz(-0.20769329) q[0];
rz(-pi) q[1];
rz(0.29783037) q[2];
sx q[2];
rz(-1.9541085) q[2];
sx q[2];
rz(-0.9155067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38182872) q[1];
sx q[1];
rz(-2.0501158) q[1];
sx q[1];
rz(-1.3717531) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60812833) q[3];
sx q[3];
rz(-1.5481037) q[3];
sx q[3];
rz(0.80858999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(-2.5547011) q[2];
rz(-1.759257) q[3];
sx q[3];
rz(-1.047225) q[3];
sx q[3];
rz(0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1047644) q[0];
sx q[0];
rz(-0.47127518) q[0];
sx q[0];
rz(0.18071827) q[0];
rz(1.7432927) q[1];
sx q[1];
rz(-1.2340052) q[1];
sx q[1];
rz(1.7416551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64547951) q[0];
sx q[0];
rz(-0.81532691) q[0];
sx q[0];
rz(-0.3078356) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79480457) q[2];
sx q[2];
rz(-0.53404885) q[2];
sx q[2];
rz(2.7164823) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3524947) q[1];
sx q[1];
rz(-2.2667655) q[1];
sx q[1];
rz(-0.11274479) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67541452) q[3];
sx q[3];
rz(-2.3590909) q[3];
sx q[3];
rz(0.43849643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7515298) q[2];
sx q[2];
rz(-0.84948245) q[2];
sx q[2];
rz(-2.0709822) q[2];
rz(1.6603598) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(1.7695561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2299131) q[0];
sx q[0];
rz(-0.80753082) q[0];
sx q[0];
rz(-2.5813778) q[0];
rz(2.920976) q[1];
sx q[1];
rz(-2.1897774) q[1];
sx q[1];
rz(-2.2645948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2219243) q[0];
sx q[0];
rz(-3.1238365) q[0];
sx q[0];
rz(-1.71868) q[0];
rz(-0.58135245) q[2];
sx q[2];
rz(-0.46581163) q[2];
sx q[2];
rz(-1.9663481) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24755281) q[1];
sx q[1];
rz(-0.28619372) q[1];
sx q[1];
rz(-1.9412089) q[1];
x q[2];
rz(-3.0362066) q[3];
sx q[3];
rz(-2.2544754) q[3];
sx q[3];
rz(-0.87382853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.666854) q[2];
sx q[2];
rz(-0.89484221) q[2];
sx q[2];
rz(1.3372927) q[2];
rz(2.9232591) q[3];
sx q[3];
rz(-1.0206914) q[3];
sx q[3];
rz(1.7201503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.6364994) q[0];
sx q[0];
rz(-1.8890843) q[0];
sx q[0];
rz(-3.002758) q[0];
rz(-0.88169634) q[1];
sx q[1];
rz(-1.1602743) q[1];
sx q[1];
rz(-2.3169611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52378319) q[0];
sx q[0];
rz(-1.6629761) q[0];
sx q[0];
rz(3.1193798) q[0];
rz(-1.6753372) q[2];
sx q[2];
rz(-1.4741159) q[2];
sx q[2];
rz(0.51890512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0141586) q[1];
sx q[1];
rz(-0.52406543) q[1];
sx q[1];
rz(2.8253943) q[1];
rz(-0.91757284) q[3];
sx q[3];
rz(-1.6503449) q[3];
sx q[3];
rz(-0.099640007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(2.3080431) q[2];
rz(0.27093568) q[3];
sx q[3];
rz(-1.1214333) q[3];
sx q[3];
rz(-0.84735316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9529097) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(-2.6771925) q[0];
rz(2.4480827) q[1];
sx q[1];
rz(-2.214407) q[1];
sx q[1];
rz(-2.4448591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90307921) q[0];
sx q[0];
rz(-3.0173324) q[0];
sx q[0];
rz(1.1073698) q[0];
rz(-1.3515662) q[2];
sx q[2];
rz(-2.9570509) q[2];
sx q[2];
rz(0.74583714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33524698) q[1];
sx q[1];
rz(-1.5434524) q[1];
sx q[1];
rz(-1.8198479) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64033033) q[3];
sx q[3];
rz(-1.5818095) q[3];
sx q[3];
rz(2.2785694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7742179) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(0.83887678) q[2];
rz(1.4539666) q[3];
sx q[3];
rz(-1.5435012) q[3];
sx q[3];
rz(1.669223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(-1.5244315) q[0];
rz(-0.35628191) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(-1.3391395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55521783) q[0];
sx q[0];
rz(-1.5065743) q[0];
sx q[0];
rz(2.1703266) q[0];
x q[1];
rz(-2.8216178) q[2];
sx q[2];
rz(-1.6386445) q[2];
sx q[2];
rz(-2.0463129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7308867) q[1];
sx q[1];
rz(-2.2458898) q[1];
sx q[1];
rz(1.7275586) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4285092) q[3];
sx q[3];
rz(-1.5350966) q[3];
sx q[3];
rz(1.7762828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7431006) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.2644281) q[2];
rz(0.33665952) q[3];
sx q[3];
rz(-2.2008937) q[3];
sx q[3];
rz(1.3338859) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7548512) q[0];
sx q[0];
rz(-1.6254397) q[0];
sx q[0];
rz(2.0425015) q[0];
rz(-0.59973888) q[1];
sx q[1];
rz(-2.6987684) q[1];
sx q[1];
rz(-0.59785688) q[1];
rz(1.9827531) q[2];
sx q[2];
rz(-1.2745274) q[2];
sx q[2];
rz(1.6714255) q[2];
rz(-1.5883432) q[3];
sx q[3];
rz(-2.1033136) q[3];
sx q[3];
rz(3.1284214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
