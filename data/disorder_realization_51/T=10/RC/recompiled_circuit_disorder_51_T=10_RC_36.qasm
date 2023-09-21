OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(2.060086) q[1];
sx q[1];
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232789) q[0];
sx q[0];
rz(-1.3594472) q[0];
sx q[0];
rz(-1.1002512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62280957) q[2];
sx q[2];
rz(-1.3277167) q[2];
sx q[2];
rz(-0.18261766) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3095113) q[1];
sx q[1];
rz(-0.81334269) q[1];
sx q[1];
rz(1.3002212) q[1];
x q[2];
rz(0.16206046) q[3];
sx q[3];
rz(-1.9063623) q[3];
sx q[3];
rz(0.1048564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(3.0701385) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(-1.6275303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058115) q[0];
sx q[0];
rz(-1.0732871) q[0];
sx q[0];
rz(-2.6742629) q[0];
rz(-pi) q[1];
rz(0.02818429) q[2];
sx q[2];
rz(-2.3079254) q[2];
sx q[2];
rz(0.29836269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4501614) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(1.3732234) q[1];
rz(-pi) q[2];
rz(1.1176923) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(0.97186175) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(0.17175737) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(3.0774975) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(1.6832738) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38760936) q[0];
sx q[0];
rz(-1.7180654) q[0];
sx q[0];
rz(3.0940042) q[0];
rz(-2.1737353) q[2];
sx q[2];
rz(-2.1328916) q[2];
sx q[2];
rz(0.92852718) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7160733) q[1];
sx q[1];
rz(-2.7287373) q[1];
sx q[1];
rz(-3.0070261) q[1];
rz(-2.8267616) q[3];
sx q[3];
rz(-1.587084) q[3];
sx q[3];
rz(-0.44486886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(2.0171719) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.4845928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(-3.1111858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(-0.53211624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1160108) q[1];
sx q[1];
rz(-2.3911747) q[1];
sx q[1];
rz(-2.2376899) q[1];
rz(-pi) q[2];
rz(-2.4184166) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(-0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-0.39988363) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(0.17366017) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325718) q[0];
sx q[0];
rz(-3.0637494) q[0];
sx q[0];
rz(-2.7709333) q[0];
x q[1];
rz(2.6661751) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.6944483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34827161) q[1];
sx q[1];
rz(-1.6760567) q[1];
sx q[1];
rz(-3.1131016) q[1];
x q[2];
rz(-2.5912656) q[3];
sx q[3];
rz(-2.7573857) q[3];
sx q[3];
rz(1.3889544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(2.1693726) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(2.629783) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(1.9981729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91771942) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6799455) q[2];
sx q[2];
rz(-2.4857773) q[2];
sx q[2];
rz(1.1193502) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1868362) q[1];
sx q[1];
rz(-1.7446767) q[1];
sx q[1];
rz(-1.1980921) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1120758) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(-1.908386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.5303622) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0798577) q[0];
sx q[0];
rz(-1.2485866) q[0];
sx q[0];
rz(1.6155924) q[0];
rz(-pi) q[1];
rz(0.71224833) q[2];
sx q[2];
rz(-1.5603666) q[2];
sx q[2];
rz(-2.1005837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28753528) q[1];
sx q[1];
rz(-0.71811986) q[1];
sx q[1];
rz(-1.1487886) q[1];
rz(-pi) q[2];
rz(-0.26182884) q[3];
sx q[3];
rz(-1.7228408) q[3];
sx q[3];
rz(2.9933628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-0.64129889) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(-3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(-1.2088998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6541518) q[0];
sx q[0];
rz(-1.3432661) q[0];
sx q[0];
rz(-0.80765101) q[0];
x q[1];
rz(2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.9343455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4603685) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(-0.016101109) q[1];
rz(-1.204793) q[3];
sx q[3];
rz(-0.98099785) q[3];
sx q[3];
rz(2.3596711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-0.87336826) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-0.31162509) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.5100381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85743839) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(0.5666825) q[0];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(2.0575112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9540625) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(1.3872836) q[3];
sx q[3];
rz(-1.7997051) q[3];
sx q[3];
rz(2.6218417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-2.4831916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347682) q[0];
sx q[0];
rz(-2.0567729) q[0];
sx q[0];
rz(-2.4776239) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.621643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2645996) q[1];
sx q[1];
rz(-1.0942642) q[1];
sx q[1];
rz(2.9546253) q[1];
rz(-1.9202616) q[3];
sx q[3];
rz(-2.0798827) q[3];
sx q[3];
rz(1.5290608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(-2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.9256928) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6170549) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(-1.9032003) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-2.1970489) q[2];
sx q[2];
rz(-1.4251475) q[2];
sx q[2];
rz(0.057694358) q[2];
rz(3.0425439) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];