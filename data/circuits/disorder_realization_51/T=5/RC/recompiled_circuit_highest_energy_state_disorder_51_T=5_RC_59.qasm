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
rz(0.32245114) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(0.21790394) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94643773) q[0];
sx q[0];
rz(-1.1338521) q[0];
sx q[0];
rz(-1.9042356) q[0];
x q[1];
rz(-1.409265) q[2];
sx q[2];
rz(-1.7551975) q[2];
sx q[2];
rz(-2.7946975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.340992) q[1];
sx q[1];
rz(-0.28352724) q[1];
sx q[1];
rz(0.42415828) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8508223) q[3];
sx q[3];
rz(-1.2655846) q[3];
sx q[3];
rz(-1.9350767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-2.71038) q[2];
sx q[2];
rz(-2.9475589) q[2];
rz(-0.38328299) q[3];
sx q[3];
rz(-0.88519874) q[3];
sx q[3];
rz(1.5017728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5782769) q[0];
sx q[0];
rz(-2.0864154) q[0];
sx q[0];
rz(-2.03736) q[0];
rz(0.68179321) q[1];
sx q[1];
rz(-1.3673404) q[1];
sx q[1];
rz(-1.7368447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0007219) q[0];
sx q[0];
rz(-2.8551176) q[0];
sx q[0];
rz(0.46517828) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4656336) q[2];
sx q[2];
rz(-2.0251207) q[2];
sx q[2];
rz(2.3788962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65741828) q[1];
sx q[1];
rz(-2.5403416) q[1];
sx q[1];
rz(2.5674393) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41635155) q[3];
sx q[3];
rz(-0.27746323) q[3];
sx q[3];
rz(2.1961138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76689395) q[2];
sx q[2];
rz(-0.99504343) q[2];
sx q[2];
rz(-1.2332756) q[2];
rz(-2.342566) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(-3.0123805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610573) q[0];
sx q[0];
rz(-2.0687456) q[0];
sx q[0];
rz(2.5824353) q[0];
rz(-1.1178389) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(-3.0303755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32229508) q[0];
sx q[0];
rz(-1.4647554) q[0];
sx q[0];
rz(-2.7631971) q[0];
rz(-pi) q[1];
rz(0.62456496) q[2];
sx q[2];
rz(-1.7918158) q[2];
sx q[2];
rz(0.92568892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7349941) q[1];
sx q[1];
rz(-1.7593699) q[1];
sx q[1];
rz(-0.85822924) q[1];
rz(-pi) q[2];
rz(2.7989332) q[3];
sx q[3];
rz(-0.4147793) q[3];
sx q[3];
rz(0.77769731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3126276) q[2];
sx q[2];
rz(-1.4451005) q[2];
sx q[2];
rz(-1.1980537) q[2];
rz(-1.4231921) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(0.47553441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51197416) q[0];
sx q[0];
rz(-2.8422575) q[0];
sx q[0];
rz(-0.83628118) q[0];
rz(-1.5760999) q[1];
sx q[1];
rz(-1.0944159) q[1];
sx q[1];
rz(-0.85087585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880575) q[0];
sx q[0];
rz(-1.5374105) q[0];
sx q[0];
rz(-0.99793156) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1117029) q[2];
sx q[2];
rz(-2.0238681) q[2];
sx q[2];
rz(1.0525931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2712114) q[1];
sx q[1];
rz(-2.3496501) q[1];
sx q[1];
rz(0.19510896) q[1];
x q[2];
rz(0.20765813) q[3];
sx q[3];
rz(-1.9586143) q[3];
sx q[3];
rz(0.70630276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1026844) q[2];
sx q[2];
rz(-2.2816198) q[2];
sx q[2];
rz(0.50722185) q[2];
rz(0.66633362) q[3];
sx q[3];
rz(-2.3907876) q[3];
sx q[3];
rz(2.2255118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136593) q[0];
sx q[0];
rz(-1.4480696) q[0];
sx q[0];
rz(-1.6130945) q[0];
rz(-1.5616034) q[1];
sx q[1];
rz(-1.0630307) q[1];
sx q[1];
rz(-2.9772421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30475475) q[0];
sx q[0];
rz(-1.4753818) q[0];
sx q[0];
rz(-1.4897904) q[0];
rz(-pi) q[1];
rz(-0.63173024) q[2];
sx q[2];
rz(-2.1513878) q[2];
sx q[2];
rz(2.3931062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4867465) q[1];
sx q[1];
rz(-0.89895144) q[1];
sx q[1];
rz(0.44507546) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2772395) q[3];
sx q[3];
rz(-1.7987712) q[3];
sx q[3];
rz(-3.0346532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5932172) q[2];
sx q[2];
rz(-0.84725738) q[2];
sx q[2];
rz(-2.4753921) q[2];
rz(0.14144746) q[3];
sx q[3];
rz(-1.5513159) q[3];
sx q[3];
rz(0.9945873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.87404609) q[0];
sx q[0];
rz(-2.5135437) q[0];
sx q[0];
rz(0.48126599) q[0];
rz(-0.73356837) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(-3.0487294) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61213148) q[0];
sx q[0];
rz(-2.4680637) q[0];
sx q[0];
rz(-0.81443925) q[0];
rz(-pi) q[1];
rz(-0.052744432) q[2];
sx q[2];
rz(-3.102902) q[2];
sx q[2];
rz(0.82088137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.085437955) q[1];
sx q[1];
rz(-2.3533818) q[1];
sx q[1];
rz(0.5931535) q[1];
rz(-pi) q[2];
rz(2.8792297) q[3];
sx q[3];
rz(-1.922058) q[3];
sx q[3];
rz(-2.2134804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4197454) q[2];
sx q[2];
rz(-1.5200619) q[2];
sx q[2];
rz(-1.2255555) q[2];
rz(-0.078977481) q[3];
sx q[3];
rz(-1.1342528) q[3];
sx q[3];
rz(1.2823766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52800286) q[0];
sx q[0];
rz(-2.4230175) q[0];
sx q[0];
rz(0.80793107) q[0];
rz(1.5241874) q[1];
sx q[1];
rz(-2.9831191) q[1];
sx q[1];
rz(-2.4375516) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86211246) q[0];
sx q[0];
rz(-3.1199146) q[0];
sx q[0];
rz(1.9901748) q[0];
x q[1];
rz(0.0041517082) q[2];
sx q[2];
rz(-1.5624245) q[2];
sx q[2];
rz(2.3821751) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0615394) q[1];
sx q[1];
rz(-0.94990094) q[1];
sx q[1];
rz(-3.1136334) q[1];
x q[2];
rz(0.16499293) q[3];
sx q[3];
rz(-2.1793302) q[3];
sx q[3];
rz(-0.42770619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0713221) q[2];
sx q[2];
rz(-2.4532048) q[2];
sx q[2];
rz(-0.88453156) q[2];
rz(-2.6279348) q[3];
sx q[3];
rz(-1.9653886) q[3];
sx q[3];
rz(-2.6562712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.1270444) q[0];
sx q[0];
rz(-0.7651279) q[0];
sx q[0];
rz(0.86084086) q[0];
rz(-3.1061106) q[1];
sx q[1];
rz(-1.2018964) q[1];
sx q[1];
rz(-1.6532345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5888212) q[0];
sx q[0];
rz(-0.5824005) q[0];
sx q[0];
rz(2.3584473) q[0];
rz(2.9126337) q[2];
sx q[2];
rz(-1.3334074) q[2];
sx q[2];
rz(1.6473532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9640298) q[1];
sx q[1];
rz(-2.7164251) q[1];
sx q[1];
rz(-1.1336826) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6658038) q[3];
sx q[3];
rz(-1.2965513) q[3];
sx q[3];
rz(1.340171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52582467) q[2];
sx q[2];
rz(-2.0066228) q[2];
sx q[2];
rz(-1.6305249) q[2];
rz(0.75119558) q[3];
sx q[3];
rz(-0.54371756) q[3];
sx q[3];
rz(-2.3403919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5649331) q[0];
sx q[0];
rz(-1.7486005) q[0];
sx q[0];
rz(-2.9826953) q[0];
rz(-1.6389906) q[1];
sx q[1];
rz(-1.331012) q[1];
sx q[1];
rz(-1.166689) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5983655) q[0];
sx q[0];
rz(-2.1413714) q[0];
sx q[0];
rz(1.8680509) q[0];
rz(-pi) q[1];
rz(-2.2750365) q[2];
sx q[2];
rz(-0.41730395) q[2];
sx q[2];
rz(0.48283726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9433971) q[1];
sx q[1];
rz(-1.5875729) q[1];
sx q[1];
rz(-1.6086964) q[1];
rz(2.1871879) q[3];
sx q[3];
rz(-3.1052178) q[3];
sx q[3];
rz(-1.1351577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4022973) q[2];
sx q[2];
rz(-0.6604971) q[2];
sx q[2];
rz(2.8759549) q[2];
rz(-1.9499251) q[3];
sx q[3];
rz(-1.3375514) q[3];
sx q[3];
rz(-2.5875097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83248508) q[0];
sx q[0];
rz(-2.5256248) q[0];
sx q[0];
rz(2.4564504) q[0];
rz(0.55662545) q[1];
sx q[1];
rz(-1.6770505) q[1];
sx q[1];
rz(2.305078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5563399) q[0];
sx q[0];
rz(-1.9808021) q[0];
sx q[0];
rz(2.7707731) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0903953) q[2];
sx q[2];
rz(-1.300989) q[2];
sx q[2];
rz(-2.4003911) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54933465) q[1];
sx q[1];
rz(-0.95227949) q[1];
sx q[1];
rz(1.0299512) q[1];
rz(2.9850858) q[3];
sx q[3];
rz(-1.3006852) q[3];
sx q[3];
rz(0.47299415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4461925) q[2];
sx q[2];
rz(-1.4245028) q[2];
sx q[2];
rz(0.85644537) q[2];
rz(-0.19927464) q[3];
sx q[3];
rz(-0.98237413) q[3];
sx q[3];
rz(-2.2299531) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0740501) q[0];
sx q[0];
rz(-1.0150801) q[0];
sx q[0];
rz(-1.6652921) q[0];
rz(-2.603727) q[1];
sx q[1];
rz(-1.270351) q[1];
sx q[1];
rz(-1.7958633) q[1];
rz(0.83897787) q[2];
sx q[2];
rz(-2.2175023) q[2];
sx q[2];
rz(-1.6432696) q[2];
rz(1.8798697) q[3];
sx q[3];
rz(-2.3936987) q[3];
sx q[3];
rz(1.0868418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
