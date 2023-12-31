OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(-0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7615258) q[0];
sx q[0];
rz(-2.4702284) q[0];
sx q[0];
rz(1.6094366) q[0];
rz(2.5976074) q[2];
sx q[2];
rz(-1.3731125) q[2];
sx q[2];
rz(2.939784) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3687009) q[1];
sx q[1];
rz(-1.4979616) q[1];
sx q[1];
rz(1.1198977) q[1];
rz(2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651543) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(2.6692078) q[0];
rz(-pi) q[1];
rz(0.86958779) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(-0.7676917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6234399) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(0.63697302) q[1];
x q[2];
rz(-0.20110735) q[3];
sx q[3];
rz(-1.4521986) q[3];
sx q[3];
rz(-0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-0.45062137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0656933) q[0];
sx q[0];
rz(-1.6487299) q[0];
sx q[0];
rz(3.1083641) q[0];
rz(2.2497555) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(2.4613949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33774099) q[1];
sx q[1];
rz(-2.1954814) q[1];
sx q[1];
rz(1.9546263) q[1];
rz(1.9335453) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(-2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(2.8362714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-2.5532789) q[0];
x q[1];
rz(0.46378739) q[2];
sx q[2];
rz(-1.9518491) q[2];
sx q[2];
rz(2.2708937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1034643) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(-0.14726463) q[1];
rz(-pi) q[2];
rz(2.0852854) q[3];
sx q[3];
rz(-1.1122276) q[3];
sx q[3];
rz(1.4460627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(2.0736407) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021279871) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(0.41723199) q[0];
x q[1];
rz(-2.876725) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(2.7153646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0438784) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(-1.9240727) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.898786) q[3];
sx q[3];
rz(-0.77951509) q[3];
sx q[3];
rz(-1.3047583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.4542788) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(-1.3735166) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5645063) q[1];
sx q[1];
rz(1.3173332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23478) q[0];
sx q[0];
rz(-2.5981075) q[0];
sx q[0];
rz(1.5275005) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66770422) q[2];
sx q[2];
rz(-2.0470847) q[2];
sx q[2];
rz(-0.27031937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.091192186) q[1];
sx q[1];
rz(-0.7351774) q[1];
sx q[1];
rz(1.5749732) q[1];
rz(-2.7139211) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0417827) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(-0.89522925) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(-1.1279001) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.2316661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8679778) q[0];
sx q[0];
rz(-2.5036739) q[0];
sx q[0];
rz(-0.62423737) q[0];
rz(-2.5008051) q[2];
sx q[2];
rz(-1.0377585) q[2];
sx q[2];
rz(-1.8719045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98037887) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(-1.4646261) q[1];
rz(-pi) q[2];
rz(-1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(2.5404239) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(-0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6153529) q[0];
sx q[0];
rz(-0.59229367) q[0];
sx q[0];
rz(-3.0583529) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7166491) q[2];
sx q[2];
rz(-2.4307494) q[2];
sx q[2];
rz(-2.5568642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2295099) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(-1.0850701) q[1];
x q[2];
rz(2.8665364) q[3];
sx q[3];
rz(-1.1797136) q[3];
sx q[3];
rz(1.0712136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.3046718) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(2.5794199) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(-1.1086585) q[0];
rz(-2.1015342) q[2];
sx q[2];
rz(-1.4897523) q[2];
sx q[2];
rz(-0.20866742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0253898) q[1];
sx q[1];
rz(-2.2491978) q[1];
sx q[1];
rz(-0.62979001) q[1];
rz(-0.4060679) q[3];
sx q[3];
rz(-1.9417524) q[3];
sx q[3];
rz(-0.0036811034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(-2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(-0.047853619) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536516) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(-2.2257462) q[0];
x q[1];
rz(2.0132952) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(1.0525345) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.88356) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(2.2116823) q[1];
x q[2];
rz(-3.0129274) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-0.50280747) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(-3.1059601) q[2];
sx q[2];
rz(-2.590066) q[2];
sx q[2];
rz(-1.8077015) q[2];
rz(1.1669284) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
