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
rz(4.2807978) q[0];
sx q[0];
rz(9.8674404) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(-2.3499188) q[1];
sx q[1];
rz(0.72905529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9937482) q[0];
sx q[0];
rz(-1.9707075) q[0];
sx q[0];
rz(2.4104959) q[0];
rz(-1.5698956) q[2];
sx q[2];
rz(-1.1810762) q[2];
sx q[2];
rz(-2.2665521) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6119365) q[1];
sx q[1];
rz(-1.3289598) q[1];
sx q[1];
rz(-2.2685758) q[1];
rz(0.9425052) q[3];
sx q[3];
rz(-1.9071336) q[3];
sx q[3];
rz(-2.9021341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86058229) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(0.78544593) q[2];
rz(2.1752518) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(-2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53168374) q[0];
sx q[0];
rz(-2.5737679) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(-2.5081778) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(-2.8107218) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2658871) q[0];
sx q[0];
rz(-0.71015753) q[0];
sx q[0];
rz(0.52010923) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.10373) q[2];
sx q[2];
rz(-1.5946009) q[2];
sx q[2];
rz(-1.8231152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4875191) q[1];
sx q[1];
rz(-1.6292058) q[1];
sx q[1];
rz(3.1078981) q[1];
x q[2];
rz(-1.5832354) q[3];
sx q[3];
rz(-1.4861938) q[3];
sx q[3];
rz(-1.0983262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6364381) q[2];
sx q[2];
rz(-0.81710368) q[2];
sx q[2];
rz(0.46112296) q[2];
rz(-0.72502208) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(-2.4729474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649699) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(2.6032676) q[0];
rz(-0.0093983924) q[1];
sx q[1];
rz(-0.63176578) q[1];
sx q[1];
rz(-1.9422772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56012404) q[0];
sx q[0];
rz(-2.631979) q[0];
sx q[0];
rz(1.4790003) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.968987) q[2];
sx q[2];
rz(-2.4510018) q[2];
sx q[2];
rz(-2.0334854) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9925132) q[1];
sx q[1];
rz(-1.01304) q[1];
sx q[1];
rz(2.7191619) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84405293) q[3];
sx q[3];
rz(-1.7622602) q[3];
sx q[3];
rz(-1.307631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37375307) q[2];
sx q[2];
rz(-1.6497296) q[2];
sx q[2];
rz(0.72244942) q[2];
rz(-0.59822285) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(1.2178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5825321) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(-3.1225358) q[0];
rz(1.4590774) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(-0.51024514) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8718349) q[0];
sx q[0];
rz(-1.8352011) q[0];
sx q[0];
rz(2.1931936) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8809659) q[2];
sx q[2];
rz(-2.2686743) q[2];
sx q[2];
rz(3.0681899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8277792) q[1];
sx q[1];
rz(-2.372532) q[1];
sx q[1];
rz(0.60375795) q[1];
x q[2];
rz(-1.4278379) q[3];
sx q[3];
rz(-1.9810487) q[3];
sx q[3];
rz(2.9825135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9973008) q[2];
sx q[2];
rz(-1.4520626) q[2];
sx q[2];
rz(0.22225456) q[2];
rz(0.079515919) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52614373) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(1.3413606) q[0];
rz(-2.6306131) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(-2.708639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110033) q[0];
sx q[0];
rz(-1.3136787) q[0];
sx q[0];
rz(1.6297618) q[0];
rz(-1.1673981) q[2];
sx q[2];
rz(-1.7678542) q[2];
sx q[2];
rz(-0.49515206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6420685) q[1];
sx q[1];
rz(-2.7244901) q[1];
sx q[1];
rz(-2.2184371) q[1];
rz(-pi) q[2];
rz(-2.8540887) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(0.35121976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13010919) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(2.1437272) q[2];
rz(2.4326371) q[3];
sx q[3];
rz(-2.7537789) q[3];
sx q[3];
rz(2.1628105) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82584941) q[0];
sx q[0];
rz(-0.5873) q[0];
sx q[0];
rz(2.24776) q[0];
rz(1.029344) q[1];
sx q[1];
rz(-0.7411595) q[1];
sx q[1];
rz(-0.11480521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66291565) q[0];
sx q[0];
rz(-1.3548242) q[0];
sx q[0];
rz(-1.3220644) q[0];
x q[1];
rz(2.6774339) q[2];
sx q[2];
rz(-2.4460829) q[2];
sx q[2];
rz(2.0328731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.25947523) q[1];
sx q[1];
rz(-2.7434556) q[1];
sx q[1];
rz(2.3601301) q[1];
x q[2];
rz(-1.6345665) q[3];
sx q[3];
rz(-1.290451) q[3];
sx q[3];
rz(2.3883723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7041695) q[2];
sx q[2];
rz(-2.2544474) q[2];
sx q[2];
rz(2.9278921) q[2];
rz(-0.63654381) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(-2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32695025) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(-3.0649109) q[1];
sx q[1];
rz(-2.2814467) q[1];
sx q[1];
rz(-1.1383879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46318836) q[0];
sx q[0];
rz(-2.7240755) q[0];
sx q[0];
rz(0.56708401) q[0];
rz(-pi) q[1];
rz(-0.78779548) q[2];
sx q[2];
rz(-1.5734563) q[2];
sx q[2];
rz(3.015369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3507211) q[1];
sx q[1];
rz(-2.6199503) q[1];
sx q[1];
rz(-1.9463368) q[1];
rz(-1.7848694) q[3];
sx q[3];
rz(-1.3297289) q[3];
sx q[3];
rz(-1.1837219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8818605) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(2.760375) q[2];
rz(0.11738736) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(3.05486) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7748902) q[0];
sx q[0];
rz(-0.77228868) q[0];
sx q[0];
rz(-2.6430687) q[0];
rz(2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(3.0438429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17247385) q[0];
sx q[0];
rz(-1.1398672) q[0];
sx q[0];
rz(-2.8043967) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7404064) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(0.40415472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5619288) q[1];
sx q[1];
rz(-1.0439928) q[1];
sx q[1];
rz(0.71814037) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1322543) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(3.0967979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8349614) q[2];
sx q[2];
rz(-1.1502879) q[2];
sx q[2];
rz(-2.5508896) q[2];
rz(0.090204209) q[3];
sx q[3];
rz(-2.7019751) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10839323) q[0];
sx q[0];
rz(-0.84302253) q[0];
sx q[0];
rz(2.3868308) q[0];
rz(0.33590487) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(-2.1051443) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6265246) q[0];
sx q[0];
rz(-1.6881516) q[0];
sx q[0];
rz(2.4548454) q[0];
rz(0.039128379) q[2];
sx q[2];
rz(-0.48155287) q[2];
sx q[2];
rz(1.0796384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9687313) q[1];
sx q[1];
rz(-1.5318499) q[1];
sx q[1];
rz(0.51677468) q[1];
x q[2];
rz(1.5124973) q[3];
sx q[3];
rz(-0.23007904) q[3];
sx q[3];
rz(-1.0844025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(-2.5392695) q[3];
sx q[3];
rz(-0.18776247) q[3];
sx q[3];
rz(-0.93835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8896821) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(-0.46345261) q[0];
rz(3.1262596) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-2.8212246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8863618) q[0];
sx q[0];
rz(-1.4275274) q[0];
sx q[0];
rz(1.6485639) q[0];
rz(0.3139852) q[2];
sx q[2];
rz(-1.8085524) q[2];
sx q[2];
rz(1.2842922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5825504) q[1];
sx q[1];
rz(-0.96331412) q[1];
sx q[1];
rz(2.5431125) q[1];
rz(2.4418751) q[3];
sx q[3];
rz(-0.43439242) q[3];
sx q[3];
rz(2.7909129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1128803) q[2];
sx q[2];
rz(-0.72673231) q[2];
sx q[2];
rz(1.0090562) q[2];
rz(-0.79695898) q[3];
sx q[3];
rz(-1.2538486) q[3];
sx q[3];
rz(0.33826452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0924031) q[0];
sx q[0];
rz(-1.5143464) q[0];
sx q[0];
rz(-1.2595246) q[0];
rz(-0.084820329) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(-1.3289497) q[2];
sx q[2];
rz(-1.6326661) q[2];
sx q[2];
rz(2.2018955) q[2];
rz(2.6081035) q[3];
sx q[3];
rz(-0.34779741) q[3];
sx q[3];
rz(1.6602914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
