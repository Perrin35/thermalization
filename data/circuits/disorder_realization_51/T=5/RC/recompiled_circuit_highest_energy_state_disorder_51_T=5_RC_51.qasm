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
rz(0.60740745) q[0];
sx q[0];
rz(-0.4439126) q[0];
sx q[0];
rz(-2.8191415) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(-2.9236887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3716952) q[0];
sx q[0];
rz(-1.2697233) q[0];
sx q[0];
rz(-2.6825344) q[0];
rz(-pi) q[1];
rz(2.9548151) q[2];
sx q[2];
rz(-1.7295656) q[2];
sx q[2];
rz(-1.1940317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7806381) q[1];
sx q[1];
rz(-1.8285969) q[1];
sx q[1];
rz(1.4514489) q[1];
rz(-2.8508223) q[3];
sx q[3];
rz(-1.8760081) q[3];
sx q[3];
rz(-1.206516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-2.71038) q[2];
sx q[2];
rz(-2.9475589) q[2];
rz(-2.7583097) q[3];
sx q[3];
rz(-2.2563939) q[3];
sx q[3];
rz(1.5017728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5782769) q[0];
sx q[0];
rz(-2.0864154) q[0];
sx q[0];
rz(-2.03736) q[0];
rz(-0.68179321) q[1];
sx q[1];
rz(-1.7742523) q[1];
sx q[1];
rz(1.404748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8786312) q[0];
sx q[0];
rz(-1.4436973) q[0];
sx q[0];
rz(-0.2574347) q[0];
rz(-pi) q[1];
rz(-2.4787729) q[2];
sx q[2];
rz(-0.79403764) q[2];
sx q[2];
rz(-0.3074078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7182285) q[1];
sx q[1];
rz(-1.8830788) q[1];
sx q[1];
rz(-0.52257089) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7252411) q[3];
sx q[3];
rz(-0.27746323) q[3];
sx q[3];
rz(0.94547887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76689395) q[2];
sx q[2];
rz(-0.99504343) q[2];
sx q[2];
rz(1.9083171) q[2];
rz(2.342566) q[3];
sx q[3];
rz(-2.7965751) q[3];
sx q[3];
rz(-3.0123805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98053539) q[0];
sx q[0];
rz(-1.072847) q[0];
sx q[0];
rz(-2.5824353) q[0];
rz(-2.0237538) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(3.0303755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2905544) q[0];
sx q[0];
rz(-1.1946329) q[0];
sx q[0];
rz(1.4567503) q[0];
rz(-pi) q[1];
rz(0.62456496) q[2];
sx q[2];
rz(-1.7918158) q[2];
sx q[2];
rz(0.92568892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4065985) q[1];
sx q[1];
rz(-1.3822228) q[1];
sx q[1];
rz(-0.85822924) q[1];
x q[2];
rz(-2.7989332) q[3];
sx q[3];
rz(-2.7268134) q[3];
sx q[3];
rz(0.77769731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8289651) q[2];
sx q[2];
rz(-1.4451005) q[2];
sx q[2];
rz(1.1980537) q[2];
rz(1.7184006) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(-2.6660582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6296185) q[0];
sx q[0];
rz(-2.8422575) q[0];
sx q[0];
rz(2.3053115) q[0];
rz(1.5654927) q[1];
sx q[1];
rz(-2.0471768) q[1];
sx q[1];
rz(0.85087585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353514) q[0];
sx q[0];
rz(-1.6041821) q[0];
sx q[0];
rz(2.1436611) q[0];
rz(-pi) q[1];
rz(0.81340547) q[2];
sx q[2];
rz(-2.4508173) q[2];
sx q[2];
rz(-1.9936313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5961413) q[1];
sx q[1];
rz(-2.343698) q[1];
sx q[1];
rz(-1.376838) q[1];
rz(1.9662764) q[3];
sx q[3];
rz(-1.3787601) q[3];
sx q[3];
rz(-0.78498299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1026844) q[2];
sx q[2];
rz(-0.85997283) q[2];
sx q[2];
rz(-0.50722185) q[2];
rz(0.66633362) q[3];
sx q[3];
rz(-0.75080502) q[3];
sx q[3];
rz(0.91608086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0136593) q[0];
sx q[0];
rz(-1.4480696) q[0];
sx q[0];
rz(1.6130945) q[0];
rz(-1.5616034) q[1];
sx q[1];
rz(-1.0630307) q[1];
sx q[1];
rz(0.16435057) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30475475) q[0];
sx q[0];
rz(-1.4753818) q[0];
sx q[0];
rz(1.4897904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5098624) q[2];
sx q[2];
rz(-0.99020489) q[2];
sx q[2];
rz(2.3931062) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.937035) q[1];
sx q[1];
rz(-1.9144844) q[1];
sx q[1];
rz(0.84852201) q[1];
rz(-0.89495792) q[3];
sx q[3];
rz(-2.7719422) q[3];
sx q[3];
rz(1.0357451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5932172) q[2];
sx q[2];
rz(-2.2943353) q[2];
sx q[2];
rz(0.66620052) q[2];
rz(-3.0001452) q[3];
sx q[3];
rz(-1.5902767) q[3];
sx q[3];
rz(2.1470054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2675466) q[0];
sx q[0];
rz(-0.62804896) q[0];
sx q[0];
rz(-0.48126599) q[0];
rz(-0.73356837) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(0.092863277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32321445) q[0];
sx q[0];
rz(-2.0131454) q[0];
sx q[0];
rz(-2.0967006) q[0];
x q[1];
rz(3.1029557) q[2];
sx q[2];
rz(-1.5687571) q[2];
sx q[2];
rz(-2.3389727) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67750061) q[1];
sx q[1];
rz(-0.94225303) q[1];
sx q[1];
rz(-1.0586865) q[1];
rz(-pi) q[2];
rz(-0.26236292) q[3];
sx q[3];
rz(-1.922058) q[3];
sx q[3];
rz(-2.2134804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7218472) q[2];
sx q[2];
rz(-1.6215308) q[2];
sx q[2];
rz(1.9160371) q[2];
rz(3.0626152) q[3];
sx q[3];
rz(-1.1342528) q[3];
sx q[3];
rz(1.2823766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6135898) q[0];
sx q[0];
rz(-0.71857518) q[0];
sx q[0];
rz(2.3336616) q[0];
rz(1.6174053) q[1];
sx q[1];
rz(-2.9831191) q[1];
sx q[1];
rz(-0.70404109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28939279) q[0];
sx q[0];
rz(-1.5796229) q[0];
sx q[0];
rz(1.5905963) q[0];
x q[1];
rz(-3.1374409) q[2];
sx q[2];
rz(-1.5791681) q[2];
sx q[2];
rz(0.75941759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6345811) q[1];
sx q[1];
rz(-1.5935362) q[1];
sx q[1];
rz(-0.94971595) q[1];
rz(2.9765997) q[3];
sx q[3];
rz(-0.96226245) q[3];
sx q[3];
rz(-0.42770619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.070270553) q[2];
sx q[2];
rz(-0.68838781) q[2];
sx q[2];
rz(0.88453156) q[2];
rz(0.5136579) q[3];
sx q[3];
rz(-1.9653886) q[3];
sx q[3];
rz(0.48532143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.1270444) q[0];
sx q[0];
rz(-0.7651279) q[0];
sx q[0];
rz(0.86084086) q[0];
rz(-0.035482081) q[1];
sx q[1];
rz(-1.9396962) q[1];
sx q[1];
rz(1.4883581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5527715) q[0];
sx q[0];
rz(-2.5591922) q[0];
sx q[0];
rz(-2.3584473) q[0];
rz(1.3272918) q[2];
sx q[2];
rz(-1.3483682) q[2];
sx q[2];
rz(-0.13130638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1775629) q[1];
sx q[1];
rz(-2.7164251) q[1];
sx q[1];
rz(1.1336826) q[1];
rz(-pi) q[2];
rz(-0.27542563) q[3];
sx q[3];
rz(-1.6622433) q[3];
sx q[3];
rz(-2.8851654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52582467) q[2];
sx q[2];
rz(-1.1349698) q[2];
sx q[2];
rz(-1.6305249) q[2];
rz(2.3903971) q[3];
sx q[3];
rz(-2.5978751) q[3];
sx q[3];
rz(-2.3403919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649331) q[0];
sx q[0];
rz(-1.3929921) q[0];
sx q[0];
rz(0.15889731) q[0];
rz(1.502602) q[1];
sx q[1];
rz(-1.8105806) q[1];
sx q[1];
rz(-1.9749036) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.027307) q[0];
sx q[0];
rz(-0.6356568) q[0];
sx q[0];
rz(-0.42814769) q[0];
x q[1];
rz(2.2750365) q[2];
sx q[2];
rz(-2.7242887) q[2];
sx q[2];
rz(-2.6587554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.769628) q[1];
sx q[1];
rz(-1.5329016) q[1];
sx q[1];
rz(-0.016788646) q[1];
rz(-2.1871879) q[3];
sx q[3];
rz(-3.1052178) q[3];
sx q[3];
rz(-2.006435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4022973) q[2];
sx q[2];
rz(-2.4810956) q[2];
sx q[2];
rz(-0.26563773) q[2];
rz(1.9499251) q[3];
sx q[3];
rz(-1.8040413) q[3];
sx q[3];
rz(0.55408293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3091076) q[0];
sx q[0];
rz(-0.6159679) q[0];
sx q[0];
rz(2.4564504) q[0];
rz(2.5849672) q[1];
sx q[1];
rz(-1.6770505) q[1];
sx q[1];
rz(0.83651465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13930411) q[0];
sx q[0];
rz(-1.2320077) q[0];
sx q[0];
rz(2.0071507) q[0];
x q[1];
rz(-1.0511974) q[2];
sx q[2];
rz(-1.300989) q[2];
sx q[2];
rz(-0.74120159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8879509) q[1];
sx q[1];
rz(-0.79756036) q[1];
sx q[1];
rz(-2.5152999) q[1];
rz(-pi) q[2];
rz(1.8440866) q[3];
sx q[3];
rz(-1.4200084) q[3];
sx q[3];
rz(2.085872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4461925) q[2];
sx q[2];
rz(-1.4245028) q[2];
sx q[2];
rz(2.2851473) q[2];
rz(2.942318) q[3];
sx q[3];
rz(-0.98237413) q[3];
sx q[3];
rz(0.91163951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0675426) q[0];
sx q[0];
rz(-1.0150801) q[0];
sx q[0];
rz(-1.6652921) q[0];
rz(-2.603727) q[1];
sx q[1];
rz(-1.270351) q[1];
sx q[1];
rz(-1.7958633) q[1];
rz(0.79277586) q[2];
sx q[2];
rz(-2.1332827) q[2];
sx q[2];
rz(2.573043) q[2];
rz(2.8665681) q[3];
sx q[3];
rz(-2.2755819) q[3];
sx q[3];
rz(-1.6439846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
