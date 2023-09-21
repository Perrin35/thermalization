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
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68361002) q[0];
sx q[0];
rz(-1.542122) q[0];
sx q[0];
rz(0.81575127) q[0];
x q[1];
rz(-0.91648957) q[2];
sx q[2];
rz(-1.5208897) q[2];
sx q[2];
rz(2.8061342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4418728) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-2.2295879) q[1];
rz(-pi) q[2];
rz(-1.5942469) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.263164) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-0.29299709) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-0.22110573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1572004) q[0];
sx q[0];
rz(-1.385014) q[0];
sx q[0];
rz(-1.1093344) q[0];
rz(1.8684623) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(0.01089451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(-1.4355852) q[1];
rz(2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(2.5533) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(1.0480405) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(2.1411238) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0962778) q[0];
sx q[0];
rz(-2.0925539) q[0];
sx q[0];
rz(-1.9406712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7311086) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(2.2018873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0211027) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(-2.2282269) q[1];
rz(-pi) q[2];
rz(-2.5030701) q[3];
sx q[3];
rz(-1.6604074) q[3];
sx q[3];
rz(-1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(-1.8130594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(1.3550718) q[0];
x q[1];
rz(0.51283299) q[2];
sx q[2];
rz(-0.83216681) q[2];
sx q[2];
rz(-1.2733449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36973876) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(1.5688483) q[1];
rz(-0.10947157) q[3];
sx q[3];
rz(-1.7913892) q[3];
sx q[3];
rz(-2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.093734309) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11855928) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(-2.8708007) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(1.9405685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5096163) q[1];
sx q[1];
rz(-0.39034931) q[1];
sx q[1];
rz(2.2401287) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(2.1577948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.323558) q[0];
sx q[0];
rz(-1.8206882) q[0];
sx q[0];
rz(-0.37139335) q[0];
rz(-pi) q[1];
rz(-1.2994453) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(-0.041610418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97555679) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(-1.5029328) q[1];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.3791929) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(-1.8968614) q[0];
rz(-1.9786644) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-0.63901627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(-1.741239) q[1];
x q[2];
rz(0.58052766) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(2.011727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.8483298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338776) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(-2.0108372) q[0];
rz(-pi) q[1];
rz(-0.60111945) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(-2.2476946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5871208) q[1];
sx q[1];
rz(-1.1677824) q[1];
sx q[1];
rz(2.9906669) q[1];
rz(-pi) q[2];
rz(1.771365) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(-1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.8539799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33289136) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-1.0813792) q[0];
rz(0.31918819) q[2];
sx q[2];
rz(-0.50779283) q[2];
sx q[2];
rz(-2.9920981) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3661256) q[1];
sx q[1];
rz(-2.2062416) q[1];
sx q[1];
rz(0.25823621) q[1];
rz(0.82448126) q[3];
sx q[3];
rz(-1.609625) q[3];
sx q[3];
rz(-0.83394921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(0.36718711) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(-3.0935862) q[0];
rz(0.37819241) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(2.3364002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50670934) q[1];
sx q[1];
rz(-0.92223972) q[1];
sx q[1];
rz(-1.3193921) q[1];
rz(-pi) q[2];
rz(-0.96079798) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-3.0565312) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(1.3973665) q[3];
sx q[3];
rz(-0.72931029) q[3];
sx q[3];
rz(1.9839877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
