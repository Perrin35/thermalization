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
rz(1.4662161) q[0];
sx q[0];
rz(-0.94453064) q[0];
sx q[0];
rz(0.20139995) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(-1.720517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22682193) q[0];
sx q[0];
rz(-2.1807007) q[0];
sx q[0];
rz(-0.066909153) q[0];
rz(-pi) q[1];
rz(-0.95008738) q[2];
sx q[2];
rz(-0.7751152) q[2];
sx q[2];
rz(-2.895854) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3609429) q[1];
sx q[1];
rz(-2.2113178) q[1];
sx q[1];
rz(2.065413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1789315) q[3];
sx q[3];
rz(-1.8800991) q[3];
sx q[3];
rz(-2.7397857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308373) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-0.24638677) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.6349207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55913299) q[0];
sx q[0];
rz(-1.8163067) q[0];
sx q[0];
rz(1.1051635) q[0];
rz(0.042889281) q[2];
sx q[2];
rz(-1.873462) q[2];
sx q[2];
rz(-1.3009225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7333201) q[1];
sx q[1];
rz(-2.773914) q[1];
sx q[1];
rz(1.4094844) q[1];
rz(-pi) q[2];
rz(-0.89119567) q[3];
sx q[3];
rz(-2.5053484) q[3];
sx q[3];
rz(-1.6980069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.9096036) q[2];
rz(-2.039382) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(-1.3365041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44322893) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(0.90782905) q[0];
rz(-2.5746131) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(0.10279113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1452652) q[0];
sx q[0];
rz(-1.3609556) q[0];
sx q[0];
rz(2.5681433) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2305713) q[2];
sx q[2];
rz(-3.0098923) q[2];
sx q[2];
rz(-2.1610799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72702209) q[1];
sx q[1];
rz(-2.1466781) q[1];
sx q[1];
rz(-2.058279) q[1];
rz(-1.0166753) q[3];
sx q[3];
rz(-1.301487) q[3];
sx q[3];
rz(2.179972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(-1.7041448) q[2];
rz(-2.7490859) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(-2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.0933541) q[0];
sx q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-1.9527973) q[0];
rz(0.28611723) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(-1.6292705) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.286519) q[0];
sx q[0];
rz(-2.5382156) q[0];
sx q[0];
rz(1.0666749) q[0];
rz(-0.34660201) q[2];
sx q[2];
rz(-0.79814974) q[2];
sx q[2];
rz(2.5426898) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0864073) q[1];
sx q[1];
rz(-2.1595104) q[1];
sx q[1];
rz(-2.9179395) q[1];
rz(-pi) q[2];
rz(2.9305601) q[3];
sx q[3];
rz(-1.0730626) q[3];
sx q[3];
rz(-0.56609234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-2.4646087) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(1.5164794) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077654) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(3.1166792) q[0];
rz(-2.086153) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(2.8581462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.442207) q[0];
sx q[0];
rz(-1.3106579) q[0];
sx q[0];
rz(-2.8752863) q[0];
rz(-2.1589958) q[2];
sx q[2];
rz(-2.9074906) q[2];
sx q[2];
rz(-2.7603619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1920212) q[1];
sx q[1];
rz(-1.5402435) q[1];
sx q[1];
rz(1.0699238) q[1];
rz(-2.2437566) q[3];
sx q[3];
rz(-1.8404318) q[3];
sx q[3];
rz(-1.2854888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8289566) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(0.32583315) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(-1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7376937) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(-0.37495908) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(0.37014827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86287303) q[0];
sx q[0];
rz(-1.5939043) q[0];
sx q[0];
rz(1.5591168) q[0];
rz(-pi) q[1];
rz(-0.13932087) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(0.56688165) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1620336) q[1];
sx q[1];
rz(-2.0311072) q[1];
sx q[1];
rz(-0.76785894) q[1];
x q[2];
rz(2.7137854) q[3];
sx q[3];
rz(-1.5662346) q[3];
sx q[3];
rz(-0.13739861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39773539) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(-2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(2.2174368) q[0];
rz(0.91668516) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7766905) q[0];
sx q[0];
rz(-2.6441537) q[0];
sx q[0];
rz(0.12125347) q[0];
rz(1.3989053) q[2];
sx q[2];
rz(-2.2187833) q[2];
sx q[2];
rz(1.6164219) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3360916) q[1];
sx q[1];
rz(-0.84118836) q[1];
sx q[1];
rz(1.355624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0823233) q[3];
sx q[3];
rz(-1.7632177) q[3];
sx q[3];
rz(-2.6317962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.224996) q[2];
rz(-0.0082958881) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55815721) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(2.2772148) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9821765) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(1.4070562) q[0];
rz(0.094674663) q[2];
sx q[2];
rz(-1.5370768) q[2];
sx q[2];
rz(1.124749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56039366) q[1];
sx q[1];
rz(-2.542109) q[1];
sx q[1];
rz(-1.1349212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2914574) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(-0.13388982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9353443) q[2];
sx q[2];
rz(-1.0103005) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-2.3681971) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(2.7763413) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(-1.431538) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5968558) q[0];
sx q[0];
rz(-1.7628306) q[0];
sx q[0];
rz(-1.3903244) q[0];
x q[1];
rz(1.0042436) q[2];
sx q[2];
rz(-2.1449617) q[2];
sx q[2];
rz(-0.023921704) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.86399779) q[1];
sx q[1];
rz(-0.91121948) q[1];
sx q[1];
rz(-1.9602174) q[1];
rz(-pi) q[2];
rz(-2.6142653) q[3];
sx q[3];
rz(-1.1865215) q[3];
sx q[3];
rz(-1.7155855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(-1.9688152) q[2];
rz(-1.734599) q[3];
sx q[3];
rz(-1.3318136) q[3];
sx q[3];
rz(1.2146568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(0.59463516) q[0];
rz(2.1356964) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(0.88919052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.517006) q[0];
sx q[0];
rz(-1.3505624) q[0];
sx q[0];
rz(-3.0896679) q[0];
rz(-0.070637704) q[2];
sx q[2];
rz(-1.1114235) q[2];
sx q[2];
rz(-2.2769711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31510776) q[1];
sx q[1];
rz(-1.4764305) q[1];
sx q[1];
rz(0.06006518) q[1];
rz(-pi) q[2];
rz(2.6252006) q[3];
sx q[3];
rz(-1.5377196) q[3];
sx q[3];
rz(0.75504485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1903926) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(2.5962043) q[2];
rz(0.96327463) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(-1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6496898) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(2.2979965) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(-0.19914535) q[2];
sx q[2];
rz(-2.0024588) q[2];
sx q[2];
rz(-1.3695516) q[2];
rz(2.5582377) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
