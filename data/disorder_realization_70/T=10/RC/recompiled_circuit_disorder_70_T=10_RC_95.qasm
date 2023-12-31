OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(-0.78599077) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91416042) q[0];
sx q[0];
rz(-0.81613805) q[0];
sx q[0];
rz(-3.1022275) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91648957) q[2];
sx q[2];
rz(-1.5208897) q[2];
sx q[2];
rz(0.33545845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.12373911) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(0.53637335) q[1];
rz(1.5473458) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.263164) q[2];
rz(1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(0.29299709) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9843922) q[0];
sx q[0];
rz(-1.7565787) q[0];
sx q[0];
rz(-2.0322582) q[0];
rz(-pi) q[1];
rz(-1.8684623) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(3.1306981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-1.1138798) q[1];
sx q[1];
rz(1.4355852) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-2.0935521) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(2.3838938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0962778) q[0];
sx q[0];
rz(-2.0925539) q[0];
sx q[0];
rz(-1.9406712) q[0];
x q[1];
rz(-0.65501113) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(-0.26817817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.12049) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(2.2282269) q[1];
rz(-0.63852255) q[3];
sx q[3];
rz(-1.6604074) q[3];
sx q[3];
rz(-1.5935957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.3285332) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.7865208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6287597) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(-1.2733449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9424142) q[1];
sx q[1];
rz(-1.5713099) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(-pi) q[2];
rz(-3.0321211) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(-2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.550068) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(-1.0250807) q[0];
x q[1];
rz(-1.1930824) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(-1.2010241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5711172) q[1];
sx q[1];
rz(-1.332453) q[1];
sx q[1];
rz(-1.2586602) q[1];
rz(0.65977804) q[3];
sx q[3];
rz(-1.0228844) q[3];
sx q[3];
rz(1.9632615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(-0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81297368) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(-2.5286753) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8421474) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(3.0999822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4420538) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-2.8738408) q[1];
rz(-pi) q[2];
rz(2.8998428) q[3];
sx q[3];
rz(-1.3931735) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.798511) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(-0.084246158) q[0];
x q[1];
rz(1.9786644) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(2.5025764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6143601) q[1];
sx q[1];
rz(-0.34253201) q[1];
sx q[1];
rz(0.50369461) q[1];
rz(-pi) q[2];
rz(-2.2961388) q[3];
sx q[3];
rz(-2.3618556) q[3];
sx q[3];
rz(2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-0.11809764) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(1.8483298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30771502) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(1.1307554) q[0];
x q[1];
rz(1.6663315) q[2];
sx q[2];
rz(-0.97180688) q[2];
sx q[2];
rz(0.73087382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5871208) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(-0.15092571) q[1];
x q[2];
rz(1.771365) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2244789) q[0];
sx q[0];
rz(-1.0815485) q[0];
sx q[0];
rz(-3.1130303) q[0];
rz(-pi) q[1];
rz(1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(-2.929504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3661256) q[1];
sx q[1];
rz(-2.2062416) q[1];
sx q[1];
rz(2.8833564) q[1];
x q[2];
rz(0.82448126) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(1.256475) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-0.36718711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-0.62717168) q[2];
sx q[2];
rz(-1.7998724) q[2];
sx q[2];
rz(-1.0692182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90875188) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(-2.8244551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.9566386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(2.8425343) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-0.085061442) q[2];
sx q[2];
rz(-2.3903575) q[2];
sx q[2];
rz(0.059546197) q[2];
rz(-0.84898938) q[3];
sx q[3];
rz(-1.6860387) q[3];
sx q[3];
rz(0.54308346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
