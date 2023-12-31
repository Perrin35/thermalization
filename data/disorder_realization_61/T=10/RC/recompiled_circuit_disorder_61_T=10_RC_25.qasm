OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(1.6917317) q[1];
sx q[1];
rz(5.6258968) q[1];
sx q[1];
rz(13.110553) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13574164) q[0];
sx q[0];
rz(-3.0648181) q[0];
sx q[0];
rz(2.6430921) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7317023) q[2];
sx q[2];
rz(-3.0311243) q[2];
sx q[2];
rz(0.2875178) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6737353) q[1];
sx q[1];
rz(-2.2506672) q[1];
sx q[1];
rz(0.52635898) q[1];
x q[2];
rz(3.1129684) q[3];
sx q[3];
rz(-1.588436) q[3];
sx q[3];
rz(-2.6580435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0698174) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(-0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(-0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4966999) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(1.9702966) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(-1.7374977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728714) q[0];
sx q[0];
rz(-0.88604468) q[0];
sx q[0];
rz(-2.1268334) q[0];
rz(0.48268433) q[2];
sx q[2];
rz(-1.7182554) q[2];
sx q[2];
rz(2.5978616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.055981936) q[1];
sx q[1];
rz(-1.2877052) q[1];
sx q[1];
rz(-1.594747) q[1];
rz(-2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(1.6625422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-0.81958333) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(0.56366411) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6113341) q[0];
sx q[0];
rz(-1.5000492) q[0];
sx q[0];
rz(0.57460873) q[0];
rz(-pi) q[1];
rz(-1.3450422) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(3.0622481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0020395) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(1.0242277) q[1];
rz(0.54791252) q[3];
sx q[3];
rz(-2.4439545) q[3];
sx q[3];
rz(1.1988977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(1.0602661) q[2];
rz(3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.32325) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(-2.9484205) q[0];
rz(-1.5974143) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(0.65778041) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13947978) q[0];
sx q[0];
rz(-1.5407028) q[0];
sx q[0];
rz(0.17480236) q[0];
x q[1];
rz(1.1409608) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(0.62361275) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0187877) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(2.8469574) q[1];
x q[2];
rz(-1.2781906) q[3];
sx q[3];
rz(-1.8183823) q[3];
sx q[3];
rz(0.86865559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0458935) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7970153) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(-1.0239333) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-2.5147298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9526564) q[0];
sx q[0];
rz(-3.1075826) q[0];
sx q[0];
rz(-1.0556428) q[0];
rz(-pi) q[1];
rz(-1.8711898) q[2];
sx q[2];
rz(-1.2671748) q[2];
sx q[2];
rz(2.6189569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31506854) q[1];
sx q[1];
rz(-1.8693923) q[1];
sx q[1];
rz(1.8930757) q[1];
rz(-3.0179126) q[3];
sx q[3];
rz(-1.9607753) q[3];
sx q[3];
rz(-2.9434413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(-0.61156887) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.942379) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(2.8170524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8400612) q[0];
sx q[0];
rz(-1.1467629) q[0];
sx q[0];
rz(1.016894) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3338823) q[2];
sx q[2];
rz(-2.0732023) q[2];
sx q[2];
rz(-0.051256996) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4536344) q[1];
sx q[1];
rz(-2.1046241) q[1];
sx q[1];
rz(0.028793528) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59431521) q[3];
sx q[3];
rz(-0.50110498) q[3];
sx q[3];
rz(-1.8584115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(1.2021525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81284886) q[0];
sx q[0];
rz(-0.7757196) q[0];
sx q[0];
rz(-0.98529718) q[0];
rz(-0.75041109) q[2];
sx q[2];
rz(-1.2806891) q[2];
sx q[2];
rz(1.4587221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89029965) q[1];
sx q[1];
rz(-1.6541462) q[1];
sx q[1];
rz(-3.0492196) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7552745) q[3];
sx q[3];
rz(-1.8309621) q[3];
sx q[3];
rz(-2.9329712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1414286) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(2.2953575) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7164562) q[0];
sx q[0];
rz(-1.3279337) q[0];
sx q[0];
rz(0.48432414) q[0];
rz(2.0267018) q[2];
sx q[2];
rz(-2.0156983) q[2];
sx q[2];
rz(1.0783431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.487405) q[1];
sx q[1];
rz(-2.020917) q[1];
sx q[1];
rz(0.022189157) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9701482) q[3];
sx q[3];
rz(-0.60895863) q[3];
sx q[3];
rz(0.82933784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1476851) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-0.85419401) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-0.6643995) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(-1.2845576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8878471) q[0];
sx q[0];
rz(-0.77639025) q[0];
sx q[0];
rz(-1.1573769) q[0];
x q[1];
rz(0.9719073) q[2];
sx q[2];
rz(-0.45058695) q[2];
sx q[2];
rz(-0.50859261) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5018651) q[1];
sx q[1];
rz(-2.0149391) q[1];
sx q[1];
rz(-0.72413866) q[1];
rz(-2.7515718) q[3];
sx q[3];
rz(-1.5454925) q[3];
sx q[3];
rz(-0.43061531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(0.70358706) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(1.528953) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3846489) q[0];
sx q[0];
rz(-2.1063519) q[0];
sx q[0];
rz(-1.7352362) q[0];
rz(-pi) q[1];
rz(-1.0829955) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(-1.4873193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70982198) q[1];
sx q[1];
rz(-0.55574544) q[1];
sx q[1];
rz(-2.7906448) q[1];
x q[2];
rz(-0.46882792) q[3];
sx q[3];
rz(-0.83823293) q[3];
sx q[3];
rz(-1.1956786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(2.4370082) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4298532) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(-1.7977057) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(0.96991878) q[2];
sx q[2];
rz(-0.51545943) q[2];
sx q[2];
rz(-3.1209844) q[2];
rz(-2.5034954) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
