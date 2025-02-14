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
rz(-2.6925447) q[0];
sx q[0];
rz(-1.1392051) q[0];
sx q[0];
rz(2.6989302) q[0];
rz(-2.9709385) q[1];
sx q[1];
rz(-0.79167384) q[1];
sx q[1];
rz(2.4125374) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1281888) q[0];
sx q[0];
rz(-0.81522757) q[0];
sx q[0];
rz(-2.5772153) q[0];
rz(-pi) q[1];
rz(0.38972028) q[2];
sx q[2];
rz(-1.5716296) q[2];
sx q[2];
rz(-0.69609797) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9022791) q[1];
sx q[1];
rz(-2.2444176) q[1];
sx q[1];
rz(0.31142733) q[1];
x q[2];
rz(-2.1990875) q[3];
sx q[3];
rz(-1.9071336) q[3];
sx q[3];
rz(-2.9021341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2810104) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(2.3561467) q[2];
rz(-0.9663409) q[3];
sx q[3];
rz(-2.0735425) q[3];
sx q[3];
rz(2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6099089) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(-2.5081778) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(-2.8107218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71459717) q[0];
sx q[0];
rz(-1.9007555) q[0];
sx q[0];
rz(-0.64100463) q[0];
rz(-3.1149338) q[2];
sx q[2];
rz(-2.0377198) q[2];
sx q[2];
rz(-2.9012762) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0111661) q[1];
sx q[1];
rz(-0.067421801) q[1];
sx q[1];
rz(1.0481337) q[1];
rz(-1.5832354) q[3];
sx q[3];
rz(-1.4861938) q[3];
sx q[3];
rz(2.0432664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50515455) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(-0.46112296) q[2];
rz(-0.72502208) q[3];
sx q[3];
rz(-2.6889763) q[3];
sx q[3];
rz(-0.66864526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649699) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(-0.53832501) q[0];
rz(3.1321943) q[1];
sx q[1];
rz(-0.63176578) q[1];
sx q[1];
rz(1.1993154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45505782) q[0];
sx q[0];
rz(-2.0780586) q[0];
sx q[0];
rz(-3.090409) q[0];
rz(-pi) q[1];
rz(-1.4298158) q[2];
sx q[2];
rz(-0.89242291) q[2];
sx q[2];
rz(-1.8110665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4863094) q[1];
sx q[1];
rz(-1.9260672) q[1];
sx q[1];
rz(-0.97092295) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8546507) q[3];
sx q[3];
rz(-0.7470658) q[3];
sx q[3];
rz(-3.0892854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7678396) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(-0.72244942) q[2];
rz(2.5433698) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(-1.9237062) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5825321) q[0];
sx q[0];
rz(-1.2578955) q[0];
sx q[0];
rz(0.019056888) q[0];
rz(1.6825153) q[1];
sx q[1];
rz(-2.9190813) q[1];
sx q[1];
rz(2.6313475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2697578) q[0];
sx q[0];
rz(-1.3063916) q[0];
sx q[0];
rz(0.94839903) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4195477) q[2];
sx q[2];
rz(-1.8068442) q[2];
sx q[2];
rz(1.8473193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45106798) q[1];
sx q[1];
rz(-0.96123403) q[1];
sx q[1];
rz(1.068348) q[1];
rz(1.7137547) q[3];
sx q[3];
rz(-1.1605439) q[3];
sx q[3];
rz(0.15907915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9973008) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(-2.9193381) q[2];
rz(-0.079515919) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(-2.5047746) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6154489) q[0];
sx q[0];
rz(-1.0960824) q[0];
sx q[0];
rz(1.3413606) q[0];
rz(2.6306131) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(2.708639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0315597) q[0];
sx q[0];
rz(-1.827914) q[0];
sx q[0];
rz(1.6297618) q[0];
rz(-1.1673981) q[2];
sx q[2];
rz(-1.3737384) q[2];
sx q[2];
rz(-2.6464406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6420685) q[1];
sx q[1];
rz(-2.7244901) q[1];
sx q[1];
rz(0.92315556) q[1];
x q[2];
rz(-2.8540887) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(0.35121976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13010919) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(2.1437272) q[2];
rz(2.4326371) q[3];
sx q[3];
rz(-2.7537789) q[3];
sx q[3];
rz(-0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82584941) q[0];
sx q[0];
rz(-0.5873) q[0];
sx q[0];
rz(-2.24776) q[0];
rz(-1.029344) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(3.0267874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.288088) q[0];
sx q[0];
rz(-1.3279606) q[0];
sx q[0];
rz(0.22260861) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2132145) q[2];
sx q[2];
rz(-0.96065694) q[2];
sx q[2];
rz(-1.6865735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0521212) q[1];
sx q[1];
rz(-1.8473746) q[1];
sx q[1];
rz(2.8514423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28088681) q[3];
sx q[3];
rz(-1.6320737) q[3];
sx q[3];
rz(2.3063502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43742314) q[2];
sx q[2];
rz(-2.2544474) q[2];
sx q[2];
rz(-2.9278921) q[2];
rz(-2.5050488) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32695025) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-0.11006926) q[0];
rz(-0.07668177) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(-1.1383879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5611794) q[0];
sx q[0];
rz(-1.3512159) q[0];
sx q[0];
rz(-2.7835569) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78779548) q[2];
sx q[2];
rz(-1.5681364) q[2];
sx q[2];
rz(0.12622368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1092905) q[1];
sx q[1];
rz(-1.7545953) q[1];
sx q[1];
rz(1.079783) q[1];
rz(-pi) q[2];
rz(0.24647572) q[3];
sx q[3];
rz(-1.3630056) q[3];
sx q[3];
rz(-2.7026619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8818605) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(0.38121769) q[2];
rz(-3.0242053) q[3];
sx q[3];
rz(-0.50506794) q[3];
sx q[3];
rz(0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7748902) q[0];
sx q[0];
rz(-0.77228868) q[0];
sx q[0];
rz(-2.6430687) q[0];
rz(-2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(-3.0438429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691188) q[0];
sx q[0];
rz(-1.1398672) q[0];
sx q[0];
rz(2.8043967) q[0];
x q[1];
rz(2.7404064) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(-0.40415472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6548582) q[1];
sx q[1];
rz(-0.86198258) q[1];
sx q[1];
rz(0.72388078) q[1];
rz(-0.0093383647) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(-3.0967979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8349614) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(-2.5508896) q[2];
rz(-0.090204209) q[3];
sx q[3];
rz(-2.7019751) q[3];
sx q[3];
rz(0.91501594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10839323) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(-2.3868308) q[0];
rz(-0.33590487) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(2.1051443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901609) q[0];
sx q[0];
rz(-2.2519172) q[0];
sx q[0];
rz(1.7220884) q[0];
x q[1];
rz(-3.1024643) q[2];
sx q[2];
rz(-0.48155287) q[2];
sx q[2];
rz(1.0796384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8120665) q[1];
sx q[1];
rz(-0.5181075) q[1];
sx q[1];
rz(0.078703367) q[1];
rz(-1.5124973) q[3];
sx q[3];
rz(-0.23007904) q[3];
sx q[3];
rz(-2.0571902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.094940946) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(2.6851192) q[2];
rz(-0.60232317) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(2.203234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25191054) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(-2.67814) q[0];
rz(0.015333029) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-0.32036805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25523089) q[0];
sx q[0];
rz(-1.7140653) q[0];
sx q[0];
rz(-1.6485639) q[0];
rz(-pi) q[1];
rz(-1.8202844) q[2];
sx q[2];
rz(-1.2659327) q[2];
sx q[2];
rz(0.21017212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55904222) q[1];
sx q[1];
rz(-0.96331412) q[1];
sx q[1];
rz(-0.59848018) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69971754) q[3];
sx q[3];
rz(-2.7072002) q[3];
sx q[3];
rz(2.7909129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1128803) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(2.1325364) q[2];
rz(2.3446337) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(2.8033281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
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
rz(3.0778733) q[2];
sx q[2];
rz(-1.812171) q[2];
sx q[2];
rz(-2.4952427) q[2];
rz(1.3884801) q[3];
sx q[3];
rz(-1.8686466) q[3];
sx q[3];
rz(2.2212088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
