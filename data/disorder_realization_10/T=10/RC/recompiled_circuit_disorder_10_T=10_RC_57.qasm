OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(-2.8432863) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059271185) q[0];
sx q[0];
rz(-2.4743609) q[0];
sx q[0];
rz(-0.088952347) q[0];
rz(-pi) q[1];
rz(-1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(-0.36117902) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5391985) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(2.5799275) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0575849) q[3];
sx q[3];
rz(-1.6735895) q[3];
sx q[3];
rz(1.8889129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(-0.47168628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52662151) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(-1.660166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6366028) q[2];
sx q[2];
rz(-1.3244751) q[2];
sx q[2];
rz(-1.9976975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7447409) q[1];
sx q[1];
rz(-2.5137915) q[1];
sx q[1];
rz(0.46698924) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7945812) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(-2.1686045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(2.1726051) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(2.95978) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91711125) q[0];
sx q[0];
rz(-2.3492976) q[0];
sx q[0];
rz(0.86865058) q[0];
rz(0.74080148) q[2];
sx q[2];
rz(-1.891279) q[2];
sx q[2];
rz(-1.6522811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.599217) q[1];
sx q[1];
rz(-1.5812751) q[1];
sx q[1];
rz(0.12565617) q[1];
x q[2];
rz(2.5163243) q[3];
sx q[3];
rz(-1.4843974) q[3];
sx q[3];
rz(0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9905332) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(2.997056) q[0];
rz(-pi) q[1];
rz(-2.5523283) q[2];
sx q[2];
rz(-1.1009842) q[2];
sx q[2];
rz(1.9643009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9150881) q[1];
sx q[1];
rz(-0.18208948) q[1];
sx q[1];
rz(-1.7712797) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2992371) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(-1.3815051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(2.3275862) q[2];
rz(1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-2.2498851) q[0];
rz(-1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-2.9290501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034886995) q[0];
sx q[0];
rz(-1.5505152) q[0];
sx q[0];
rz(3.1253392) q[0];
rz(-pi) q[1];
rz(1.5151305) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(-2.2948613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37413874) q[1];
sx q[1];
rz(-2.1530188) q[1];
sx q[1];
rz(-2.2351082) q[1];
rz(-pi) q[2];
rz(-0.84380031) q[3];
sx q[3];
rz(-1.2146597) q[3];
sx q[3];
rz(-0.86405495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(0.22512063) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(-2.7640142) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8416653) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(-1.4873234) q[0];
rz(-pi) q[1];
rz(-2.8541366) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(3.1099144) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7403455) q[1];
sx q[1];
rz(-2.0013072) q[1];
sx q[1];
rz(3.0262448) q[1];
x q[2];
rz(-0.029929786) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(-2.6605576) q[2];
rz(0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(0.19454923) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(2.887168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5146778) q[0];
sx q[0];
rz(-1.8906381) q[0];
sx q[0];
rz(3.0266648) q[0];
x q[1];
rz(2.9872586) q[2];
sx q[2];
rz(-0.96417226) q[2];
sx q[2];
rz(-0.91147214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19941482) q[1];
sx q[1];
rz(-2.605038) q[1];
sx q[1];
rz(-0.55120991) q[1];
rz(-pi) q[2];
rz(-2.6894327) q[3];
sx q[3];
rz(-1.1434571) q[3];
sx q[3];
rz(-2.6219581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(-2.1379437) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5664068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9109089) q[0];
sx q[0];
rz(-0.77373234) q[0];
sx q[0];
rz(0.45549972) q[0];
rz(-pi) q[1];
rz(-1.8640395) q[2];
sx q[2];
rz(-2.5817421) q[2];
sx q[2];
rz(-0.44550371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34711429) q[1];
sx q[1];
rz(-1.9273259) q[1];
sx q[1];
rz(1.9630678) q[1];
rz(-pi) q[2];
rz(-1.1493707) q[3];
sx q[3];
rz(-1.1024269) q[3];
sx q[3];
rz(-2.3295662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87551293) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(-2.5121571) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4515848) q[0];
sx q[0];
rz(-2.0802054) q[0];
sx q[0];
rz(-0.33126979) q[0];
x q[1];
rz(2.7669737) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6362308) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(1.2760217) q[1];
x q[2];
rz(0.40739079) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(2.9152169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(-0.0020290931) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(-0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(2.1955406) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34110585) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(-2.1675046) q[0];
rz(2.7962748) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(-0.40030865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(0.65995364) q[1];
rz(-1.1935812) q[3];
sx q[3];
rz(-2.0743138) q[3];
sx q[3];
rz(-0.52770381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0845906) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(0.65336147) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54031298) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(2.1379708) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.3238974) q[3];
sx q[3];
rz(-2.8084451) q[3];
sx q[3];
rz(0.13767903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];