OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44242087) q[0];
sx q[0];
rz(-2.3306263) q[0];
sx q[0];
rz(-0.45642689) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(-0.18263291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73494324) q[0];
sx q[0];
rz(-1.8514086) q[0];
sx q[0];
rz(-2.2191802) q[0];
rz(-0.47477291) q[2];
sx q[2];
rz(-2.2150196) q[2];
sx q[2];
rz(2.7930773) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4306972) q[1];
sx q[1];
rz(-0.17696807) q[1];
sx q[1];
rz(0.21459337) q[1];
x q[2];
rz(-1.1785281) q[3];
sx q[3];
rz(-2.0571821) q[3];
sx q[3];
rz(-1.9928336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7543588) q[2];
sx q[2];
rz(-2.0626455) q[2];
sx q[2];
rz(-0.31164247) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-0.25185549) q[3];
sx q[3];
rz(0.18865147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99609128) q[0];
sx q[0];
rz(-1.4459193) q[0];
sx q[0];
rz(-1.2392932) q[0];
rz(-0.39341012) q[1];
sx q[1];
rz(-2.0937803) q[1];
sx q[1];
rz(-2.1107103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1686954) q[0];
sx q[0];
rz(-0.6873695) q[0];
sx q[0];
rz(-1.0052135) q[0];
rz(-pi) q[1];
rz(-1.3412881) q[2];
sx q[2];
rz(-1.3131427) q[2];
sx q[2];
rz(0.15783707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0075025) q[1];
sx q[1];
rz(-2.1975027) q[1];
sx q[1];
rz(-2.7078663) q[1];
rz(-pi) q[2];
rz(-0.88406422) q[3];
sx q[3];
rz(-2.2746448) q[3];
sx q[3];
rz(-2.2313909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82765141) q[2];
sx q[2];
rz(-2.2590019) q[2];
sx q[2];
rz(-2.7871056) q[2];
rz(-0.83141023) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(-1.6262511) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671234) q[0];
sx q[0];
rz(-2.9614083) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(-0.88227415) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(0.54214111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71996337) q[0];
sx q[0];
rz(-2.6609592) q[0];
sx q[0];
rz(-2.292657) q[0];
rz(-pi) q[1];
rz(0.19511055) q[2];
sx q[2];
rz(-2.8144022) q[2];
sx q[2];
rz(-1.3766152) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24607813) q[1];
sx q[1];
rz(-2.3195004) q[1];
sx q[1];
rz(-2.8991634) q[1];
x q[2];
rz(-1.8291953) q[3];
sx q[3];
rz(-0.2451788) q[3];
sx q[3];
rz(0.49743957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3279646) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(-2.1072809) q[2];
rz(1.7799001) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(2.5907607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4807602) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(-2.0939636) q[0];
rz(1.818559) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(-0.38527647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40257257) q[0];
sx q[0];
rz(-2.077335) q[0];
sx q[0];
rz(0.49474025) q[0];
rz(-pi) q[1];
rz(-2.4249737) q[2];
sx q[2];
rz(-1.6277378) q[2];
sx q[2];
rz(-0.54887923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4874607) q[1];
sx q[1];
rz(-1.7791554) q[1];
sx q[1];
rz(-1.6435433) q[1];
rz(3.0317467) q[3];
sx q[3];
rz(-2.2513736) q[3];
sx q[3];
rz(2.3504013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72383991) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(-0.27457944) q[2];
rz(0.56139055) q[3];
sx q[3];
rz(-1.0654819) q[3];
sx q[3];
rz(1.6339462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0059589) q[0];
sx q[0];
rz(-0.57666403) q[0];
sx q[0];
rz(-0.86719257) q[0];
rz(0.37569702) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(1.9120749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73342434) q[0];
sx q[0];
rz(-2.6717792) q[0];
sx q[0];
rz(2.7759659) q[0];
x q[1];
rz(-1.0631353) q[2];
sx q[2];
rz(-0.57786059) q[2];
sx q[2];
rz(2.7323728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83513488) q[1];
sx q[1];
rz(-1.1145488) q[1];
sx q[1];
rz(2.2133676) q[1];
rz(-pi) q[2];
rz(1.983316) q[3];
sx q[3];
rz(-1.2591921) q[3];
sx q[3];
rz(1.7330979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1162954) q[2];
sx q[2];
rz(-1.9090434) q[2];
sx q[2];
rz(0.06289014) q[2];
rz(2.7316015) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(-2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9922239) q[0];
sx q[0];
rz(-3.0719482) q[0];
sx q[0];
rz(-2.3058291) q[0];
rz(1.958485) q[1];
sx q[1];
rz(-1.8712021) q[1];
sx q[1];
rz(2.3979208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2398674) q[0];
sx q[0];
rz(-2.5081303) q[0];
sx q[0];
rz(2.6153436) q[0];
rz(-1.2116777) q[2];
sx q[2];
rz(-0.15378498) q[2];
sx q[2];
rz(1.6395456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53568038) q[1];
sx q[1];
rz(-2.4125068) q[1];
sx q[1];
rz(-1.2950743) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9414976) q[3];
sx q[3];
rz(-2.8809014) q[3];
sx q[3];
rz(-2.835915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(-2.3480603) q[2];
rz(2.7225336) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(-2.6194173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94988743) q[0];
sx q[0];
rz(-2.1623623) q[0];
sx q[0];
rz(1.1631843) q[0];
rz(0.78041068) q[1];
sx q[1];
rz(-2.8051832) q[1];
sx q[1];
rz(1.5283248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6097117) q[0];
sx q[0];
rz(-1.6356633) q[0];
sx q[0];
rz(2.2676629) q[0];
rz(-pi) q[1];
rz(-1.9204936) q[2];
sx q[2];
rz(-1.7426859) q[2];
sx q[2];
rz(0.34687172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1133729) q[1];
sx q[1];
rz(-1.372678) q[1];
sx q[1];
rz(-3.1199725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9933543) q[3];
sx q[3];
rz(-1.8236287) q[3];
sx q[3];
rz(-1.8971083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10716001) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(3.072928) q[2];
rz(-0.56985235) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(-2.7710052) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625075) q[0];
sx q[0];
rz(-2.469049) q[0];
sx q[0];
rz(1.8261209) q[0];
rz(-1.8585809) q[1];
sx q[1];
rz(-2.7048769) q[1];
sx q[1];
rz(0.049093094) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4955208) q[0];
sx q[0];
rz(-1.6611413) q[0];
sx q[0];
rz(0.077383872) q[0];
rz(0.95048381) q[2];
sx q[2];
rz(-1.1178218) q[2];
sx q[2];
rz(0.17826232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.093542592) q[1];
sx q[1];
rz(-0.88359264) q[1];
sx q[1];
rz(0.50205135) q[1];
rz(0.92691874) q[3];
sx q[3];
rz(-1.4689494) q[3];
sx q[3];
rz(0.97130126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7580238) q[2];
sx q[2];
rz(-2.263676) q[2];
sx q[2];
rz(-0.54086584) q[2];
rz(1.0848378) q[3];
sx q[3];
rz(-0.70298755) q[3];
sx q[3];
rz(1.7613523) q[3];
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
rz(-2.2329907) q[0];
sx q[0];
rz(-2.1451696) q[0];
sx q[0];
rz(0.10502271) q[0];
rz(2.5999293) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(-0.35596102) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.134577) q[0];
sx q[0];
rz(-1.7099713) q[0];
sx q[0];
rz(0.037875847) q[0];
rz(-pi) q[1];
rz(-2.2860724) q[2];
sx q[2];
rz(-0.86453712) q[2];
sx q[2];
rz(2.962213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40416607) q[1];
sx q[1];
rz(-1.9475137) q[1];
sx q[1];
rz(0.95312037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2139236) q[3];
sx q[3];
rz(-2.8017453) q[3];
sx q[3];
rz(1.353689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9291222) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(-1.3193725) q[2];
rz(1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(0.96778473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.20973715) q[0];
sx q[0];
rz(-0.034448817) q[0];
sx q[0];
rz(-1.6784278) q[0];
rz(-1.4596918) q[1];
sx q[1];
rz(-1.9821143) q[1];
sx q[1];
rz(-0.34585888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5423842) q[0];
sx q[0];
rz(-0.092723474) q[0];
sx q[0];
rz(-2.1044162) q[0];
x q[1];
rz(-2.8246872) q[2];
sx q[2];
rz(-1.5090885) q[2];
sx q[2];
rz(1.2807434) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0557424) q[1];
sx q[1];
rz(-0.80778236) q[1];
sx q[1];
rz(-2.4959688) q[1];
rz(-pi) q[2];
rz(-2.8000051) q[3];
sx q[3];
rz(-2.1862967) q[3];
sx q[3];
rz(0.08344354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.8701376) q[2];
sx q[2];
rz(-0.26091179) q[2];
rz(-1.6711309) q[3];
sx q[3];
rz(-1.7390395) q[3];
sx q[3];
rz(1.4298593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83229257) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(0.55083864) q[1];
sx q[1];
rz(-1.6549587) q[1];
sx q[1];
rz(0.5400198) q[1];
rz(-1.3030686) q[2];
sx q[2];
rz(-0.74543759) q[2];
sx q[2];
rz(-1.8882295) q[2];
rz(-1.8455918) q[3];
sx q[3];
rz(-0.85920371) q[3];
sx q[3];
rz(0.68026713) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
