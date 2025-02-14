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
rz(1.2020741) q[0];
sx q[0];
rz(-2.095686) q[0];
sx q[0];
rz(2.3910971) q[0];
rz(0.22076386) q[1];
sx q[1];
rz(-2.0808487) q[1];
sx q[1];
rz(-1.1836675) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96776828) q[0];
sx q[0];
rz(-0.90961752) q[0];
sx q[0];
rz(-1.038616) q[0];
x q[1];
rz(1.4040825) q[2];
sx q[2];
rz(-1.445747) q[2];
sx q[2];
rz(0.84593433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69931251) q[1];
sx q[1];
rz(-0.89095014) q[1];
sx q[1];
rz(-1.555843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1424155) q[3];
sx q[3];
rz(-2.7795254) q[3];
sx q[3];
rz(1.3130275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3262647) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(-3.0495138) q[2];
rz(-1.1312671) q[3];
sx q[3];
rz(-1.0695894) q[3];
sx q[3];
rz(0.85589516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6446514) q[0];
sx q[0];
rz(-2.6583789) q[0];
sx q[0];
rz(2.605873) q[0];
rz(0.016853111) q[1];
sx q[1];
rz(-2.2622175) q[1];
sx q[1];
rz(-3.1244315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760289) q[0];
sx q[0];
rz(-2.5390609) q[0];
sx q[0];
rz(2.1621611) q[0];
x q[1];
rz(0.13249915) q[2];
sx q[2];
rz(-2.6713058) q[2];
sx q[2];
rz(2.7072226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2338154) q[1];
sx q[1];
rz(-1.8026979) q[1];
sx q[1];
rz(1.4091989) q[1];
rz(2.5480367) q[3];
sx q[3];
rz(-2.9800804) q[3];
sx q[3];
rz(0.42481977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0993593) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(0.24822203) q[2];
rz(0.92784268) q[3];
sx q[3];
rz(-2.35858) q[3];
sx q[3];
rz(0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8645653) q[0];
sx q[0];
rz(-1.9643354) q[0];
sx q[0];
rz(2.9789341) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-2.9167852) q[1];
sx q[1];
rz(-0.83388296) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637109) q[0];
sx q[0];
rz(-2.1378008) q[0];
sx q[0];
rz(-1.7162697) q[0];
rz(-2.8188848) q[2];
sx q[2];
rz(-1.6889179) q[2];
sx q[2];
rz(0.28216991) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0030886) q[1];
sx q[1];
rz(-0.89996594) q[1];
sx q[1];
rz(-0.53536949) q[1];
rz(-pi) q[2];
rz(1.6271851) q[3];
sx q[3];
rz(-1.1306376) q[3];
sx q[3];
rz(1.4408011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0791846) q[2];
sx q[2];
rz(-1.1391613) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(2.0832113) q[3];
sx q[3];
rz(-0.94401413) q[3];
sx q[3];
rz(3.1112166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.77085483) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.4282861) q[0];
rz(2.3772073) q[1];
sx q[1];
rz(-1.1420219) q[1];
sx q[1];
rz(1.3216602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.747731) q[0];
sx q[0];
rz(-1.9698304) q[0];
sx q[0];
rz(1.1927496) q[0];
rz(-pi) q[1];
rz(1.2390529) q[2];
sx q[2];
rz(-1.52196) q[2];
sx q[2];
rz(1.2570753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.053406) q[1];
sx q[1];
rz(-1.2866396) q[1];
sx q[1];
rz(-2.8233145) q[1];
rz(-0.85614397) q[3];
sx q[3];
rz(-1.2309181) q[3];
sx q[3];
rz(1.4902756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6617714) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(2.856355) q[2];
rz(1.8789004) q[3];
sx q[3];
rz(-1.217239) q[3];
sx q[3];
rz(1.4234022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406463) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(0.88229156) q[0];
rz(-0.67059416) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(-0.98463279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66259164) q[0];
sx q[0];
rz(-0.07018319) q[0];
sx q[0];
rz(2.6778738) q[0];
rz(-pi) q[1];
rz(1.4996455) q[2];
sx q[2];
rz(-0.54279581) q[2];
sx q[2];
rz(-1.1384447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24422503) q[1];
sx q[1];
rz(-0.80522081) q[1];
sx q[1];
rz(2.3697788) q[1];
rz(0.71159848) q[3];
sx q[3];
rz(-2.4715373) q[3];
sx q[3];
rz(-1.3942476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.275445) q[2];
sx q[2];
rz(-2.6068942) q[2];
sx q[2];
rz(2.5246485) q[2];
rz(-0.75211891) q[3];
sx q[3];
rz(-1.813846) q[3];
sx q[3];
rz(0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98944703) q[0];
sx q[0];
rz(-0.716827) q[0];
sx q[0];
rz(-0.77051198) q[0];
rz(-0.72692263) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(-0.70473421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73913888) q[0];
sx q[0];
rz(-2.1069035) q[0];
sx q[0];
rz(-2.2741063) q[0];
x q[1];
rz(-2.0915085) q[2];
sx q[2];
rz(-1.2479223) q[2];
sx q[2];
rz(2.8276557) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26624666) q[1];
sx q[1];
rz(-0.35142144) q[1];
sx q[1];
rz(-0.09697341) q[1];
x q[2];
rz(1.1868531) q[3];
sx q[3];
rz(-0.60572165) q[3];
sx q[3];
rz(1.5662409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3720234) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(0.43323576) q[2];
rz(-2.9406934) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(0.5350298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83605003) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-0.38247633) q[0];
rz(2.1005232) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(-2.6952851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5322558) q[0];
sx q[0];
rz(-0.26383886) q[0];
sx q[0];
rz(-0.20045073) q[0];
x q[1];
rz(0.38241295) q[2];
sx q[2];
rz(-1.3578041) q[2];
sx q[2];
rz(1.2455924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12536167) q[1];
sx q[1];
rz(-0.70763677) q[1];
sx q[1];
rz(-2.1317922) q[1];
x q[2];
rz(-0.84844037) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(-0.27001303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3216766) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(0.72810158) q[2];
rz(0.32650945) q[3];
sx q[3];
rz(-1.510334) q[3];
sx q[3];
rz(0.85730332) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9051055) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(1.2976728) q[0];
rz(-1.0009276) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(0.40850684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9169642) q[0];
sx q[0];
rz(-3.1273807) q[0];
sx q[0];
rz(1.5111708) q[0];
x q[1];
rz(0.82168545) q[2];
sx q[2];
rz(-0.61640451) q[2];
sx q[2];
rz(-0.44611713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0169819) q[1];
sx q[1];
rz(-1.4370185) q[1];
sx q[1];
rz(0.11511187) q[1];
rz(-pi) q[2];
rz(-0.37398963) q[3];
sx q[3];
rz(-1.9714103) q[3];
sx q[3];
rz(0.11696045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11477509) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(0.82359037) q[2];
rz(-2.1987727) q[3];
sx q[3];
rz(-1.7190245) q[3];
sx q[3];
rz(-1.5028809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5470062) q[0];
sx q[0];
rz(-0.44438812) q[0];
sx q[0];
rz(1.2521) q[0];
rz(-1.3638672) q[1];
sx q[1];
rz(-1.5003279) q[1];
sx q[1];
rz(-1.3346671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20864381) q[0];
sx q[0];
rz(-1.874039) q[0];
sx q[0];
rz(-1.6088845) q[0];
rz(-pi) q[1];
rz(-0.23128839) q[2];
sx q[2];
rz(-2.4866001) q[2];
sx q[2];
rz(-0.97823036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6872944) q[1];
sx q[1];
rz(-2.0210938) q[1];
sx q[1];
rz(-1.6520834) q[1];
rz(-1.077721) q[3];
sx q[3];
rz(-2.6209313) q[3];
sx q[3];
rz(0.23948224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.021726457) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(0.29652706) q[2];
rz(1.6622539) q[3];
sx q[3];
rz(-1.6300423) q[3];
sx q[3];
rz(-2.9843946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088365) q[0];
sx q[0];
rz(-0.65082508) q[0];
sx q[0];
rz(-2.346709) q[0];
rz(2.5859313) q[1];
sx q[1];
rz(-1.8652922) q[1];
sx q[1];
rz(1.4478252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705161) q[0];
sx q[0];
rz(-0.85502842) q[0];
sx q[0];
rz(2.0809196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6332288) q[2];
sx q[2];
rz(-2.0678021) q[2];
sx q[2];
rz(1.7071498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8965079) q[1];
sx q[1];
rz(-1.0963529) q[1];
sx q[1];
rz(-0.22405476) q[1];
rz(-pi) q[2];
rz(-2.1673739) q[3];
sx q[3];
rz(-2.0250626) q[3];
sx q[3];
rz(-2.5429986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29622233) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(1.7108062) q[2];
rz(-0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13567781) q[0];
sx q[0];
rz(-2.4527241) q[0];
sx q[0];
rz(1.2149568) q[0];
rz(-1.6750492) q[1];
sx q[1];
rz(-0.92862447) q[1];
sx q[1];
rz(-2.5359572) q[1];
rz(-2.5765924) q[2];
sx q[2];
rz(-1.0641268) q[2];
sx q[2];
rz(0.63882154) q[2];
rz(-1.7411655) q[3];
sx q[3];
rz(-1.6765249) q[3];
sx q[3];
rz(-0.75679814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
