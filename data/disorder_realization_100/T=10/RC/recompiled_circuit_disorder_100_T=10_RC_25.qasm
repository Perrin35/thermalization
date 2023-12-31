OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(0.64414135) q[1];
sx q[1];
rz(7.9692366) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(0.88168983) q[0];
x q[1];
rz(-2.2725265) q[2];
sx q[2];
rz(-1.7742187) q[2];
sx q[2];
rz(0.8806526) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75016025) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(2.0425914) q[1];
rz(-pi) q[2];
rz(-1.5815758) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(-2.3673494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(3.0018905) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.2765983) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.93719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4780087) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(-2.708486) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18560974) q[2];
sx q[2];
rz(-2.2679272) q[2];
sx q[2];
rz(2.6999161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(-1.4644535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2180087) q[3];
sx q[3];
rz(-0.70397607) q[3];
sx q[3];
rz(-2.9670027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.191167) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1026099) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(-0.15977504) q[0];
rz(-pi) q[1];
rz(-0.28352719) q[2];
sx q[2];
rz(-1.3361738) q[2];
sx q[2];
rz(1.7262176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.46398677) q[1];
sx q[1];
rz(-1.7772487) q[1];
sx q[1];
rz(-0.21519214) q[1];
rz(-1.1690087) q[3];
sx q[3];
rz(-2.0766633) q[3];
sx q[3];
rz(1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(2.1276786) q[2];
rz(-1.7845456) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2134593) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(-0.34780985) q[0];
x q[1];
rz(-1.5825908) q[2];
sx q[2];
rz(-1.6147385) q[2];
sx q[2];
rz(2.5269788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1624565) q[1];
sx q[1];
rz(-0.41185954) q[1];
sx q[1];
rz(-1.3084175) q[1];
rz(-1.5971848) q[3];
sx q[3];
rz(-2.4370586) q[3];
sx q[3];
rz(1.1772616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61280695) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.7808419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4454942) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(-1.7646837) q[0];
rz(-pi) q[1];
rz(-2.3809793) q[2];
sx q[2];
rz(-1.0906272) q[2];
sx q[2];
rz(1.3590571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9779224) q[1];
sx q[1];
rz(-2.0241892) q[1];
sx q[1];
rz(-1.8015253) q[1];
rz(-pi) q[2];
rz(0.97134437) q[3];
sx q[3];
rz(-1.506862) q[3];
sx q[3];
rz(-1.5435227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3551066) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.8967569) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48588612) q[0];
sx q[0];
rz(-0.61300346) q[0];
sx q[0];
rz(-2.1711151) q[0];
rz(-pi) q[1];
x q[1];
rz(2.198344) q[2];
sx q[2];
rz(-2.3970282) q[2];
sx q[2];
rz(-2.9668691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5092768) q[1];
sx q[1];
rz(-0.61513153) q[1];
sx q[1];
rz(2.5110911) q[1];
rz(3.0589468) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.28073072) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(-1.77805) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-0.71969676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5170383) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(-2.4863003) q[0];
rz(-pi) q[1];
rz(-2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(-0.072349116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4622725) q[1];
sx q[1];
rz(-1.6151036) q[1];
sx q[1];
rz(2.5965152) q[1];
x q[2];
rz(0.68088716) q[3];
sx q[3];
rz(-1.3586992) q[3];
sx q[3];
rz(-1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.482243) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(0.18653175) q[0];
rz(2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.6465181) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(2.5748475) q[0];
x q[1];
rz(-1.375884) q[2];
sx q[2];
rz(-2.2576828) q[2];
sx q[2];
rz(1.9538823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.013156803) q[1];
sx q[1];
rz(-1.6746579) q[1];
sx q[1];
rz(0.26058148) q[1];
rz(-pi) q[2];
rz(2.5008194) q[3];
sx q[3];
rz(-0.91859222) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(3.126826) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(-1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(0.87798464) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.7810129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36688766) q[0];
sx q[0];
rz(-2.5686712) q[0];
sx q[0];
rz(1.5013055) q[0];
rz(-pi) q[1];
rz(1.6352064) q[2];
sx q[2];
rz(-0.62811479) q[2];
sx q[2];
rz(1.8154085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1189551) q[1];
sx q[1];
rz(-1.9070909) q[1];
sx q[1];
rz(-0.35204661) q[1];
x q[2];
rz(0.15501539) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(-0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(2.2488135) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2980625) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(2.4216901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(-1.8887595) q[0];
rz(-pi) q[1];
rz(-0.30883046) q[2];
sx q[2];
rz(-1.575982) q[2];
sx q[2];
rz(-1.3877102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9626999) q[1];
sx q[1];
rz(-1.8176515) q[1];
sx q[1];
rz(0.015951338) q[1];
rz(0.73086892) q[3];
sx q[3];
rz(-1.0167529) q[3];
sx q[3];
rz(-0.73885599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(-2.1574216) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
