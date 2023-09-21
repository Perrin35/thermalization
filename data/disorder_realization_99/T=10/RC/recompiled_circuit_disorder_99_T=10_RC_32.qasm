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
rz(2.1652048) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3800669) q[0];
sx q[0];
rz(-2.4702284) q[0];
sx q[0];
rz(-1.5321561) q[0];
rz(-0.54398529) q[2];
sx q[2];
rz(-1.7684801) q[2];
sx q[2];
rz(-2.939784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(3.0607037) q[1];
rz(-pi) q[2];
rz(0.84150746) q[3];
sx q[3];
rz(-0.91472799) q[3];
sx q[3];
rz(-2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67316002) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(3.1288778) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(2.1495842) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734946) q[0];
sx q[0];
rz(-2.0272278) q[0];
sx q[0];
rz(-1.2903851) q[0];
rz(-pi) q[1];
rz(-2.2426376) q[2];
sx q[2];
rz(-2.0546753) q[2];
sx q[2];
rz(1.8091786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18842998) q[1];
sx q[1];
rz(-0.97524446) q[1];
sx q[1];
rz(1.9838536) q[1];
rz(-pi) q[2];
rz(-2.6037381) q[3];
sx q[3];
rz(-0.23306498) q[3];
sx q[3];
rz(2.2190998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(1.0773405) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(0.45062137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492836) q[0];
sx q[0];
rz(-1.5376687) q[0];
sx q[0];
rz(-1.4928198) q[0];
rz(0.89183715) q[2];
sx q[2];
rz(-2.6747181) q[2];
sx q[2];
rz(-0.6801978) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1994985) q[1];
sx q[1];
rz(-0.71951413) q[1];
sx q[1];
rz(-0.47902963) q[1];
rz(-pi) q[2];
rz(1.2080473) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(-0.7286275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(1.191656) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-2.0003831) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-2.8362714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13320623) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-2.5532789) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1496454) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(-2.2575833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0381283) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(0.14726463) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78379811) q[3];
sx q[3];
rz(-2.4664306) q[3];
sx q[3];
rz(0.53962196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(-0.40571037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.273542) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(2.4540841) q[0];
x q[1];
rz(1.1666127) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(-3.0886138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.09771422) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(-1.21752) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3769576) q[3];
sx q[3];
rz(-1.7406165) q[3];
sx q[3];
rz(-3.0499383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.4542788) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(-2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23478) q[0];
sx q[0];
rz(-2.5981075) q[0];
sx q[0];
rz(-1.6140922) q[0];
x q[1];
rz(0.98951927) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(1.6473824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.091192186) q[1];
sx q[1];
rz(-0.7351774) q[1];
sx q[1];
rz(1.5666195) q[1];
rz(-pi) q[2];
rz(-0.96885724) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(-2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.09981) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(0.89522925) q[2];
rz(1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-1.1279001) q[0];
rz(2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.9099265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989482) q[0];
sx q[0];
rz(-2.0751187) q[0];
sx q[0];
rz(-1.1619316) q[0];
x q[1];
rz(-2.3627794) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(-0.29701172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5994508) q[1];
sx q[1];
rz(-0.25240024) q[1];
sx q[1];
rz(2.7155994) q[1];
rz(-pi) q[2];
x q[2];
rz(0.021846847) q[3];
sx q[3];
rz(-2.1751746) q[3];
sx q[3];
rz(-1.6250087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-0.59246078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62646482) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(-1.5149087) q[0];
rz(-pi) q[1];
rz(-1.0560889) q[2];
sx q[2];
rz(-2.0851496) q[2];
sx q[2];
rz(-0.27013847) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8471624) q[1];
sx q[1];
rz(-1.8995598) q[1];
sx q[1];
rz(2.2734103) q[1];
rz(-pi) q[2];
rz(1.1660277) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(-0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(-2.2413065) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(-3.1282848) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(2.5794199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87555128) q[0];
sx q[0];
rz(-0.79715675) q[0];
sx q[0];
rz(2.6334727) q[0];
x q[1];
rz(-1.0400585) q[2];
sx q[2];
rz(-1.6518403) q[2];
sx q[2];
rz(-0.20866742) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1162029) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(2.5118026) q[1];
rz(1.1702873) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(-1.7290982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(3.093739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809966) q[0];
sx q[0];
rz(-0.91616917) q[0];
sx q[0];
rz(3.1050443) q[0];
rz(-pi) q[1];
rz(-3.0907862) q[2];
sx q[2];
rz(-1.1287969) q[2];
sx q[2];
rz(-0.49651422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0117482) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(0.5188491) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1286653) q[3];
sx q[3];
rz(-2.7358958) q[3];
sx q[3];
rz(2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(2.6387852) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(0.55124333) q[2];
sx q[2];
rz(-1.5894645) q[2];
sx q[2];
rz(-0.20655256) q[2];
rz(2.7754178) q[3];
sx q[3];
rz(-1.9504642) q[3];
sx q[3];
rz(-1.4049243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];