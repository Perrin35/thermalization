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
rz(-0.7926994) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108609) q[0];
sx q[0];
rz(-2.2415677) q[0];
sx q[0];
rz(-0.030681507) q[0];
x q[1];
rz(-1.8007457) q[2];
sx q[2];
rz(-2.1030428) q[2];
sx q[2];
rz(1.4872273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1945038) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(-1.7366921) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.760261) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(-2.965062) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(2.7761249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680981) q[0];
sx q[0];
rz(-2.0272278) q[0];
sx q[0];
rz(-1.8512076) q[0];
x q[1];
rz(-0.59132691) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(-3.025625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47479113) q[1];
sx q[1];
rz(-0.7102237) q[1];
sx q[1];
rz(0.53479654) q[1];
rz(-pi) q[2];
rz(2.9404853) q[3];
sx q[3];
rz(-1.6893941) q[3];
sx q[3];
rz(0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(2.5307632) q[0];
rz(1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(2.6909713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6621646) q[0];
sx q[0];
rz(-3.056884) q[0];
sx q[0];
rz(-1.9730294) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2497555) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(-2.4613949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.140481) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(0.66092296) q[1];
x q[2];
rz(-0.30626014) q[3];
sx q[3];
rz(-0.89950409) q[3];
sx q[3];
rz(-1.1989532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(-0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-2.0003831) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(0.30532125) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353377) q[0];
sx q[0];
rz(-2.4843785) q[0];
sx q[0];
rz(2.6150136) q[0];
x q[1];
rz(2.4112169) q[2];
sx q[2];
rz(-2.5502898) q[2];
sx q[2];
rz(-0.060746047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.893559) q[1];
sx q[1];
rz(-2.3602924) q[1];
sx q[1];
rz(1.7209189) q[1];
rz(-pi) q[2];
rz(-0.78379811) q[3];
sx q[3];
rz(-2.4664306) q[3];
sx q[3];
rz(2.6019707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(-0.40571037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021279871) q[0];
sx q[0];
rz(-2.4097674) q[0];
sx q[0];
rz(-0.41723199) q[0];
rz(-pi) q[1];
rz(2.876725) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(-2.7153646) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.09771422) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(1.21752) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3769576) q[3];
sx q[3];
rz(-1.4009762) q[3];
sx q[3];
rz(0.09165435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562397) q[0];
sx q[0];
rz(-2.1137153) q[0];
sx q[0];
rz(0.026144233) q[0];
rz(-2.4738884) q[2];
sx q[2];
rz(-2.0470847) q[2];
sx q[2];
rz(2.8712733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4827022) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(-0.83562327) q[1];
rz(-0.96885724) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(-2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0417827) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(-2.2463634) q[2];
rz(1.2081395) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(-2.0136925) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3198277) q[0];
sx q[0];
rz(-1.9263096) q[0];
sx q[0];
rz(2.6000644) q[0];
rz(-pi) q[1];
rz(-2.3627794) q[2];
sx q[2];
rz(-0.80873571) q[2];
sx q[2];
rz(0.29701172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5994508) q[1];
sx q[1];
rz(-0.25240024) q[1];
sx q[1];
rz(0.42599328) q[1];
rz(-pi) q[2];
rz(-1.6024186) q[3];
sx q[3];
rz(-0.60472371) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(0.59246078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1661467) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(0.59068824) q[0];
x q[1];
rz(-0.57581298) q[2];
sx q[2];
rz(-2.0137219) q[2];
sx q[2];
rz(1.5720313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2295099) q[1];
sx q[1];
rz(-0.76369897) q[1];
sx q[1];
rz(1.0850701) q[1];
x q[2];
rz(-0.27505625) q[3];
sx q[3];
rz(-1.1797136) q[3];
sx q[3];
rz(-2.070379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.028037926) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.3046718) q[0];
rz(1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-0.56217271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8174443) q[0];
sx q[0];
rz(-1.9262909) q[0];
sx q[0];
rz(-2.4118511) q[0];
rz(-1.7299037) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(-1.2250587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8852946) q[1];
sx q[1];
rz(-2.2513299) q[1];
sx q[1];
rz(-2.2018432) q[1];
x q[2];
rz(0.77769827) q[3];
sx q[3];
rz(-2.5986528) q[3];
sx q[3];
rz(0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(-2.1600058) q[2];
rz(3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809966) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(0.036548338) q[0];
x q[1];
rz(-1.1282975) q[2];
sx q[2];
rz(-1.5248761) q[2];
sx q[2];
rz(-1.0525345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0117482) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(2.6227436) q[1];
rz(-pi) q[2];
rz(1.6258532) q[3];
sx q[3];
rz(-1.9729455) q[3];
sx q[3];
rz(-1.0226585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8235648) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(0.55124333) q[2];
sx q[2];
rz(-1.5894645) q[2];
sx q[2];
rz(-0.20655256) q[2];
rz(0.36617483) q[3];
sx q[3];
rz(-1.1911285) q[3];
sx q[3];
rz(1.7366684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
