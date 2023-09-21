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
rz(2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91416042) q[0];
sx q[0];
rz(-0.81613805) q[0];
sx q[0];
rz(3.1022275) q[0];
rz(1.6526821) q[2];
sx q[2];
rz(-0.6559283) q[2];
sx q[2];
rz(1.3002849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4418728) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(2.2295879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5942469) q[3];
sx q[3];
rz(-1.970885) q[3];
sx q[3];
rz(1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9156076) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(2.9204869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76925812) q[0];
sx q[0];
rz(-2.6466469) q[0];
sx q[0];
rz(1.9702205) q[0];
rz(-pi) q[1];
rz(-1.2731304) q[2];
sx q[2];
rz(-0.3634828) q[2];
sx q[2];
rz(-0.01089451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7213388) q[1];
sx q[1];
rz(-1.4495279) q[1];
sx q[1];
rz(2.6810357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6684181) q[3];
sx q[3];
rz(-1.0704346) q[3];
sx q[3];
rz(0.2122768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0962778) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.2009215) q[0];
x q[1];
rz(0.65501113) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(-0.26817817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.12049) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(0.91336577) q[1];
rz(2.9919639) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(-0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.8130594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.7865208) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3782016) q[2];
sx q[2];
rz(-1.1995458) q[2];
sx q[2];
rz(-0.064918092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37175825) q[1];
sx q[1];
rz(-1.8375988) q[1];
sx q[1];
rz(0.00053243551) q[1];
rz(1.7926746) q[3];
sx q[3];
rz(-1.6776049) q[3];
sx q[3];
rz(-0.91339236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-0.50450528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38406661) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(1.9996044) q[0];
rz(-pi) q[1];
rz(2.675823) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(-2.6054232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0658768) q[1];
sx q[1];
rz(-1.2677691) q[1];
sx q[1];
rz(-2.8916343) q[1];
rz(-pi) q[2];
rz(0.65977804) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15121962) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(1.8381401) q[0];
x q[1];
rz(-1.2994453) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(3.0999822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4420538) q[1];
sx q[1];
rz(-2.8901254) q[1];
sx q[1];
rz(2.8738408) q[1];
rz(-2.8998428) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.7623998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.798511) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(0.084246158) q[0];
x q[1];
rz(-1.9786644) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(0.63901627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56475793) q[1];
sx q[1];
rz(-1.4079637) q[1];
sx q[1];
rz(0.30270438) q[1];
rz(-pi) q[2];
rz(-0.84545387) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-0.29512063) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30771502) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(2.0108372) q[0];
rz(-pi) q[1];
rz(0.60111945) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(-0.89389801) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0656933) q[1];
sx q[1];
rz(-1.7095487) q[1];
sx q[1];
rz(-1.9779512) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6929018) q[3];
sx q[3];
rz(-0.43827) q[3];
sx q[3];
rz(1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.8539799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8087013) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-2.0602134) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(-0.21208866) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95083642) q[1];
sx q[1];
rz(-1.7777998) q[1];
sx q[1];
rz(-0.91916577) q[1];
rz(-pi) q[2];
rz(-1.5136396) q[3];
sx q[3];
rz(-2.3944629) q[3];
sx q[3];
rz(-0.77880083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-2.7744055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254926) q[0];
sx q[0];
rz(-1.5228132) q[0];
sx q[0];
rz(1.6019437) q[0];
x q[1];
rz(1.2904097) q[2];
sx q[2];
rz(-2.1791611) q[2];
sx q[2];
rz(2.8031363) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2314184) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(-2.477596) q[1];
rz(-2.6606584) q[3];
sx q[3];
rz(-0.98610611) q[3];
sx q[3];
rz(-1.1993053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(3.0565312) q[2];
sx q[2];
rz(-2.3903575) q[2];
sx q[2];
rz(0.059546197) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
