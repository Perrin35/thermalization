OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(0.71917978) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0702695) q[0];
sx q[0];
rz(-1.9914314) q[0];
sx q[0];
rz(3.0778528) q[0];
rz(-2.1626986) q[2];
sx q[2];
rz(-1.4375293) q[2];
sx q[2];
rz(-2.8507581) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4157384) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(1.0456677) q[1];
x q[2];
rz(0.037976102) q[3];
sx q[3];
rz(-1.3278604) q[3];
sx q[3];
rz(1.2167041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44895479) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(-2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4028567) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-2.7785595) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(1.8889069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4607836) q[0];
sx q[0];
rz(-1.4810116) q[0];
sx q[0];
rz(2.9814274) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86979903) q[2];
sx q[2];
rz(-0.79248488) q[2];
sx q[2];
rz(-3.062641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0879841) q[1];
sx q[1];
rz(-1.0454185) q[1];
sx q[1];
rz(-0.56478516) q[1];
rz(-pi) q[2];
rz(0.18685347) q[3];
sx q[3];
rz(-1.1396176) q[3];
sx q[3];
rz(0.29160515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0210555) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-2.8072642) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-0.25046644) q[0];
sx q[0];
rz(-0.0070455889) q[0];
rz(2.7650611) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(-2.4287756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3526488) q[0];
sx q[0];
rz(-1.7907356) q[0];
sx q[0];
rz(-2.9262732) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17194925) q[2];
sx q[2];
rz(-1.1977344) q[2];
sx q[2];
rz(-0.821515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3908927) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(-0.23551029) q[1];
rz(-2.0531822) q[3];
sx q[3];
rz(-1.4995121) q[3];
sx q[3];
rz(-0.34382581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5793844) q[2];
sx q[2];
rz(0.66398579) q[2];
rz(-1.9021696) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5749213) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(2.5675039) q[0];
rz(-2.8494049) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(-0.92837292) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1855477) q[0];
sx q[0];
rz(-1.5717686) q[0];
sx q[0];
rz(2.2768343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2296154) q[2];
sx q[2];
rz(-2.3301635) q[2];
sx q[2];
rz(-2.2974643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54138994) q[1];
sx q[1];
rz(-1.1536414) q[1];
sx q[1];
rz(0.27020176) q[1];
rz(-pi) q[2];
rz(2.6392322) q[3];
sx q[3];
rz(-0.97548786) q[3];
sx q[3];
rz(-0.82752284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(1.1553923) q[2];
rz(-0.14285764) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(-2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0578385) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-3.0673448) q[0];
rz(-1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-0.9517076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1002016) q[0];
sx q[0];
rz(-0.75267422) q[0];
sx q[0];
rz(-1.4647096) q[0];
rz(-pi) q[1];
rz(-2.3299626) q[2];
sx q[2];
rz(-1.9981761) q[2];
sx q[2];
rz(0.4610093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61102066) q[1];
sx q[1];
rz(-1.7360577) q[1];
sx q[1];
rz(0.9475488) q[1];
x q[2];
rz(0.40162556) q[3];
sx q[3];
rz(-1.0394319) q[3];
sx q[3];
rz(-1.7925966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(1.1523694) q[2];
rz(-2.0955829) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9933269) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(0.48962012) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(-1.5354059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.673296) q[0];
sx q[0];
rz(-1.6185986) q[0];
sx q[0];
rz(0.099138069) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(0.51598179) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8363758) q[1];
sx q[1];
rz(-1.4123385) q[1];
sx q[1];
rz(1.2837019) q[1];
rz(-pi) q[2];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(2.0819506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-0.17573389) q[2];
rz(-1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(1.3659182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7570092) q[0];
sx q[0];
rz(-0.88338822) q[0];
sx q[0];
rz(1.13152) q[0];
rz(0.62909158) q[2];
sx q[2];
rz(-0.39789879) q[2];
sx q[2];
rz(1.2049904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3733702) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(2.8819487) q[1];
x q[2];
rz(-0.52929438) q[3];
sx q[3];
rz(-1.6779473) q[3];
sx q[3];
rz(1.2308434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0697249) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(2.9979624) q[3];
sx q[3];
rz(-1.8840021) q[3];
sx q[3];
rz(0.34415054) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-2.8198077) q[0];
rz(0.92542648) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(2.5245573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2935534) q[0];
sx q[0];
rz(-0.32395054) q[0];
sx q[0];
rz(-2.241961) q[0];
x q[1];
rz(-1.355079) q[2];
sx q[2];
rz(-1.6032013) q[2];
sx q[2];
rz(1.0499133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1083793) q[1];
sx q[1];
rz(-1.9404267) q[1];
sx q[1];
rz(1.0934248) q[1];
x q[2];
rz(1.6156322) q[3];
sx q[3];
rz(-2.3196844) q[3];
sx q[3];
rz(-2.3007948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-3.0275596) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(1.6581416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720855) q[0];
sx q[0];
rz(-1.5811231) q[0];
sx q[0];
rz(-1.6540065) q[0];
rz(-2.9786417) q[2];
sx q[2];
rz(-2.5451247) q[2];
sx q[2];
rz(-0.60566723) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9121453) q[1];
sx q[1];
rz(-1.427269) q[1];
sx q[1];
rz(-1.9043282) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5914707) q[3];
sx q[3];
rz(-1.9024444) q[3];
sx q[3];
rz(1.5105997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.2443939) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72220951) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(-0.37049946) q[0];
rz(2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(-1.3409021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6751624) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(1.1943597) q[0];
x q[1];
rz(0.15721639) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(-2.9027028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95883137) q[1];
sx q[1];
rz(-1.3211831) q[1];
sx q[1];
rz(-1.1968489) q[1];
rz(1.5503928) q[3];
sx q[3];
rz(-1.4495747) q[3];
sx q[3];
rz(-2.5308479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6293388) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(-1.6188999) q[2];
rz(-2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2387977) q[0];
sx q[0];
rz(-1.4807985) q[0];
sx q[0];
rz(2.280622) q[0];
rz(0.32342708) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.6265097) q[2];
sx q[2];
rz(-2.734676) q[2];
sx q[2];
rz(2.717072) q[2];
rz(2.2446185) q[3];
sx q[3];
rz(-0.57861181) q[3];
sx q[3];
rz(1.2275916) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
