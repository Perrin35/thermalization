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
rz(0.44904798) q[0];
sx q[0];
rz(-2.0023876) q[0];
sx q[0];
rz(0.44266242) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(3.9332665) q[1];
sx q[1];
rz(10.153833) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75890827) q[0];
sx q[0];
rz(-2.2331862) q[0];
sx q[0];
rz(2.0871967) q[0];
rz(0.0021931113) q[2];
sx q[2];
rz(-2.7518715) q[2];
sx q[2];
rz(-2.268923) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52965611) q[1];
sx q[1];
rz(-1.8126329) q[1];
sx q[1];
rz(-0.87301689) q[1];
rz(-pi) q[2];
rz(0.9425052) q[3];
sx q[3];
rz(-1.9071336) q[3];
sx q[3];
rz(-2.9021341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2810104) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(0.78544593) q[2];
rz(2.1752518) q[3];
sx q[3];
rz(-2.0735425) q[3];
sx q[3];
rz(-0.6399703) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6099089) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(-2.5058643) q[0];
rz(0.6334148) q[1];
sx q[1];
rz(-0.76786357) q[1];
sx q[1];
rz(2.8107218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71459717) q[0];
sx q[0];
rz(-1.9007555) q[0];
sx q[0];
rz(0.64100463) q[0];
x q[1];
rz(-1.517968) q[2];
sx q[2];
rz(-0.46762782) q[2];
sx q[2];
rz(0.29948452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2229022) q[1];
sx q[1];
rz(-1.5371593) q[1];
sx q[1];
rz(1.5123537) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14563609) q[3];
sx q[3];
rz(-0.085509954) q[3];
sx q[3];
rz(2.189429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.50515455) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(0.46112296) q[2];
rz(2.4165706) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(-2.4729474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(1.4918936) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(0.53832501) q[0];
rz(-3.1321943) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(-1.9422772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56012404) q[0];
sx q[0];
rz(-0.5096137) q[0];
sx q[0];
rz(1.4790003) q[0];
rz(-1.7117768) q[2];
sx q[2];
rz(-2.2491697) q[2];
sx q[2];
rz(-1.8110665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4863094) q[1];
sx q[1];
rz(-1.9260672) q[1];
sx q[1];
rz(2.1706697) q[1];
x q[2];
rz(0.84405293) q[3];
sx q[3];
rz(-1.7622602) q[3];
sx q[3];
rz(1.307631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7678396) q[2];
sx q[2];
rz(-1.6497296) q[2];
sx q[2];
rz(2.4191432) q[2];
rz(0.59822285) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(1.9237062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5825321) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(-0.019056888) q[0];
rz(-1.4590774) q[1];
sx q[1];
rz(-2.9190813) q[1];
sx q[1];
rz(2.6313475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8718349) q[0];
sx q[0];
rz(-1.3063916) q[0];
sx q[0];
rz(-0.94839903) q[0];
rz(-0.72204496) q[2];
sx q[2];
rz(-1.3347484) q[2];
sx q[2];
rz(-1.8473193) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31381346) q[1];
sx q[1];
rz(-0.76906068) q[1];
sx q[1];
rz(0.60375795) q[1];
rz(-pi) q[2];
rz(0.41401569) q[3];
sx q[3];
rz(-1.7018205) q[3];
sx q[3];
rz(-1.4690635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9973008) q[2];
sx q[2];
rz(-1.4520626) q[2];
sx q[2];
rz(0.22225456) q[2];
rz(3.0620767) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(0.63681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6154489) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(1.8002321) q[0];
rz(2.6306131) q[1];
sx q[1];
rz(-1.363441) q[1];
sx q[1];
rz(0.43295369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2596672) q[0];
sx q[0];
rz(-2.8779463) q[0];
sx q[0];
rz(2.9211098) q[0];
rz(1.1002925) q[2];
sx q[2];
rz(-2.6950201) q[2];
sx q[2];
rz(2.4960325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49952415) q[1];
sx q[1];
rz(-2.7244901) q[1];
sx q[1];
rz(2.2184371) q[1];
x q[2];
rz(-1.8515685) q[3];
sx q[3];
rz(-2.3237356) q[3];
sx q[3];
rz(0.75033891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13010919) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(-0.9978655) q[2];
rz(2.4326371) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3157432) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(-2.24776) q[0];
rz(-1.029344) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(3.0267874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.478677) q[0];
sx q[0];
rz(-1.3548242) q[0];
sx q[0];
rz(-1.8195282) q[0];
rz(0.64115142) q[2];
sx q[2];
rz(-1.8617407) q[2];
sx q[2];
rz(3.0464622) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8821174) q[1];
sx q[1];
rz(-0.39813706) q[1];
sx q[1];
rz(-2.3601301) q[1];
x q[2];
rz(-2.8607058) q[3];
sx q[3];
rz(-1.6320737) q[3];
sx q[3];
rz(2.3063502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43742314) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(0.21370055) q[2];
rz(-2.5050488) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32695025) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(0.07668177) q[1];
sx q[1];
rz(-2.2814467) q[1];
sx q[1];
rz(2.0032047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58041328) q[0];
sx q[0];
rz(-1.3512159) q[0];
sx q[0];
rz(-0.35803573) q[0];
rz(-1.5745671) q[2];
sx q[2];
rz(-2.3585883) q[2];
sx q[2];
rz(-1.4418999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79087154) q[1];
sx q[1];
rz(-2.6199503) q[1];
sx q[1];
rz(1.9463368) q[1];
rz(-2.4289598) q[3];
sx q[3];
rz(-2.8205835) q[3];
sx q[3];
rz(0.4450624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8818605) q[2];
sx q[2];
rz(-1.0819819) q[2];
sx q[2];
rz(0.38121769) q[2];
rz(-0.11738736) q[3];
sx q[3];
rz(-0.50506794) q[3];
sx q[3];
rz(3.05486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748902) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(2.3122834) q[1];
sx q[1];
rz(-0.81559759) q[1];
sx q[1];
rz(0.097749762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691188) q[0];
sx q[0];
rz(-1.1398672) q[0];
sx q[0];
rz(2.8043967) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4011863) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(0.40415472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48673442) q[1];
sx q[1];
rz(-2.2796101) q[1];
sx q[1];
rz(0.72388078) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5427715) q[3];
sx q[3];
rz(-2.8198979) q[3];
sx q[3];
rz(3.0672586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30663124) q[2];
sx q[2];
rz(-1.1502879) q[2];
sx q[2];
rz(-2.5508896) q[2];
rz(-3.0513884) q[3];
sx q[3];
rz(-2.7019751) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.10839323) q[0];
sx q[0];
rz(-0.84302253) q[0];
sx q[0];
rz(-2.3868308) q[0];
rz(0.33590487) q[1];
sx q[1];
rz(-0.55481189) q[1];
sx q[1];
rz(-1.0364484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514318) q[0];
sx q[0];
rz(-2.2519172) q[0];
sx q[0];
rz(1.7220884) q[0];
rz(-pi) q[1];
rz(-1.5503564) q[2];
sx q[2];
rz(-1.0896434) q[2];
sx q[2];
rz(2.1060973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3295262) q[1];
sx q[1];
rz(-2.6234851) q[1];
sx q[1];
rz(-0.078703367) q[1];
rz(-pi) q[2];
rz(-1.5124973) q[3];
sx q[3];
rz(-2.9115136) q[3];
sx q[3];
rz(-1.0844025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.094940946) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(0.60232317) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(0.93835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8896821) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(-2.67814) q[0];
rz(-0.015333029) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-2.8212246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.244359) q[0];
sx q[0];
rz(-2.9787052) q[0];
sx q[0];
rz(-0.49400671) q[0];
rz(-pi) q[1];
rz(-1.3213083) q[2];
sx q[2];
rz(-1.87566) q[2];
sx q[2];
rz(-2.9314205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7586472) q[1];
sx q[1];
rz(-1.0898814) q[1];
sx q[1];
rz(-2.2702899) q[1];
rz(-pi) q[2];
rz(1.2804561) q[3];
sx q[3];
rz(-1.2429894) q[3];
sx q[3];
rz(2.7442055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.028712332) q[2];
sx q[2];
rz(-0.72673231) q[2];
sx q[2];
rz(1.0090562) q[2];
rz(-2.3446337) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(0.33826452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0491895) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(3.0567723) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(-1.3176805) q[2];
sx q[2];
rz(-2.8921078) q[2];
sx q[2];
rz(0.38550384) q[2];
rz(-2.6081035) q[3];
sx q[3];
rz(-2.7937952) q[3];
sx q[3];
rz(-1.4813012) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
