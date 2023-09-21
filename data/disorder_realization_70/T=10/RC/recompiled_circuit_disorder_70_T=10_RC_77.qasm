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
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68361002) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(-2.3258414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2251031) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-0.33545845) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7635599) q[1];
sx q[1];
rz(-2.0114007) q[1];
sx q[1];
rz(-0.9159169) q[1];
x q[2];
rz(1.5942469) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.263164) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-0.22110573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195933) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(-2.9346912) q[0];
x q[1];
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
rz(3.0529777) q[1];
sx q[1];
rz(-2.6664475) q[1];
sx q[1];
rz(-2.8739724) q[1];
x q[2];
rz(-2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-2.3838938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0453148) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(-1.2009215) q[0];
rz(-pi) q[1];
rz(-0.20747848) q[2];
sx q[2];
rz(-1.4138655) q[2];
sx q[2];
rz(-0.66397882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0231087) q[1];
sx q[1];
rz(-0.99220905) q[1];
sx q[1];
rz(-2.5793377) q[1];
x q[2];
rz(1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-0.089126822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-1.3805566) q[0];
sx q[0];
rz(-1.7839) q[0];
sx q[0];
rz(-0.15844945) q[0];
x q[1];
rz(-0.51283299) q[2];
sx q[2];
rz(-0.83216681) q[2];
sx q[2];
rz(1.2733449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37175825) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(3.1410602) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0321211) q[3];
sx q[3];
rz(-1.7913892) q[3];
sx q[3];
rz(-2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(0.85737491) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11855928) q[0];
sx q[0];
rz(-1.0414062) q[0];
sx q[0];
rz(2.8708007) q[0];
x q[1];
rz(0.46576969) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(-0.53616947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5711172) q[1];
sx q[1];
rz(-1.332453) q[1];
sx q[1];
rz(1.2586602) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(-2.1577948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990373) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(1.3034526) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2115057) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(1.2668244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4420538) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(0.26775189) q[1];
rz(-2.8998428) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.3791929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3430816) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(-0.084246158) q[0];
rz(-1.9786644) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(0.63901627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(1.4003537) q[1];
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
rz(-0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.8483298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30771502) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(-2.0108372) q[0];
x q[1];
rz(1.4752611) q[2];
sx q[2];
rz(-0.97180688) q[2];
sx q[2];
rz(2.4107188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1845491) q[1];
sx q[1];
rz(-2.7126985) q[1];
sx q[1];
rz(-1.2317608) q[1];
rz(-1.3702277) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.8539799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8563961) q[0];
sx q[0];
rz(-2.6515793) q[0];
sx q[0];
rz(1.6243837) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48607562) q[2];
sx q[2];
rz(-1.4176148) q[2];
sx q[2];
rz(1.4391522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1907562) q[1];
sx q[1];
rz(-1.7777998) q[1];
sx q[1];
rz(0.91916577) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3171114) q[3];
sx q[3];
rz(-1.609625) q[3];
sx q[3];
rz(-0.83394921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(2.7744055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(0.048006417) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2904097) q[2];
sx q[2];
rz(-2.1791611) q[2];
sx q[2];
rz(0.3384564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2328408) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(2.8244551) q[1];
x q[2];
rz(-2.6606584) q[3];
sx q[3];
rz(-0.98610611) q[3];
sx q[3];
rz(1.9422873) q[3];
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
rz(0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(0.74942855) q[2];
sx q[2];
rz(-1.5127758) q[2];
sx q[2];
rz(1.5681058) q[2];
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