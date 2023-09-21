OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85672985) q[0];
sx q[0];
rz(-0.75548178) q[0];
sx q[0];
rz(1.5289686) q[0];
x q[1];
rz(-0.062866048) q[2];
sx q[2];
rz(-0.91744423) q[2];
sx q[2];
rz(1.1970929) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7635599) q[1];
sx q[1];
rz(-1.130192) q[1];
sx q[1];
rz(-0.9159169) q[1];
rz(-pi) q[2];
rz(0.05539031) q[3];
sx q[3];
rz(-0.4007383) q[3];
sx q[3];
rz(1.8000359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195933) q[0];
sx q[0];
rz(-1.1178734) q[0];
sx q[0];
rz(-0.20690147) q[0];
rz(-pi) q[1];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.9176033) q[2];
sx q[2];
rz(0.32804104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21048966) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(1.4355852) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(-2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.706447) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-0.56133095) q[0];
x q[1];
rz(-1.4104841) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(-2.2018873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11848395) q[1];
sx q[1];
rz(-2.1493836) q[1];
sx q[1];
rz(-2.5793377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63852255) q[3];
sx q[3];
rz(-1.6604074) q[3];
sx q[3];
rz(1.5935957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(2.6304723) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.3285332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.7865208) q[0];
rz(-0.51283299) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(1.8682478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9424142) q[1];
sx q[1];
rz(-1.5713099) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(-pi) q[2];
rz(-1.3489181) q[3];
sx q[3];
rz(-1.4639877) q[3];
sx q[3];
rz(-2.2282003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5915247) q[0];
sx q[0];
rz(-1.3378157) q[0];
sx q[0];
rz(-1.0250807) q[0];
rz(-pi) q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(-1.9405685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.075715847) q[1];
sx q[1];
rz(-1.8738235) q[1];
sx q[1];
rz(2.8916343) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4818146) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(-1.9632615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81297368) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(2.5286753) q[0];
rz(-pi) q[1];
rz(-1.2115057) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(1.2668244) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1660359) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(1.5029328) q[1];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-0.32399696) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.3791929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0066854) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(-1.2447312) q[0];
rz(-1.1629282) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-2.5025764) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5768347) q[1];
sx q[1];
rz(-1.4079637) q[1];
sx q[1];
rz(0.30270438) q[1];
rz(-pi) q[2];
rz(-2.561065) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(-1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30771502) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(1.1307554) q[0];
rz(-pi) q[1];
rz(1.4752611) q[2];
sx q[2];
rz(-2.1697858) q[2];
sx q[2];
rz(0.73087382) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9570436) q[1];
sx q[1];
rz(-2.7126985) q[1];
sx q[1];
rz(1.9098319) q[1];
x q[2];
rz(-1.3702277) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(2.5602706) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(-2.1447694) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9171137) q[0];
sx q[0];
rz(-2.0600442) q[0];
sx q[0];
rz(-3.1130303) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31918819) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(0.14949456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-2.8833564) q[1];
x q[2];
rz(-0.052863315) q[3];
sx q[3];
rz(-2.3164146) q[3];
sx q[3];
rz(2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8866667) q[0];
sx q[0];
rz(-1.5228132) q[0];
sx q[0];
rz(1.6019437) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62717168) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(-1.0692182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50670934) q[1];
sx q[1];
rz(-2.2193529) q[1];
sx q[1];
rz(1.8222005) q[1];
rz(-pi) q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.9566386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.665303) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-0.74942855) q[2];
sx q[2];
rz(-1.6288169) q[2];
sx q[2];
rz(-1.5734869) q[2];
rz(-1.3973665) q[3];
sx q[3];
rz(-2.4122824) q[3];
sx q[3];
rz(-1.1576049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
