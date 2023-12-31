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
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(0.64414135) q[1];
sx q[1];
rz(7.9692366) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830084) q[0];
sx q[0];
rz(-1.9134247) q[0];
sx q[0];
rz(-1.1230099) q[0];
x q[1];
rz(2.2725265) q[2];
sx q[2];
rz(-1.7742187) q[2];
sx q[2];
rz(-0.8806526) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75016025) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(-2.0425914) q[1];
rz(-pi) q[2];
rz(0.047921501) q[3];
sx q[3];
rz(-1.5815634) q[3];
sx q[3];
rz(-0.79706942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(3.0018905) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(-2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4780087) q[0];
sx q[0];
rz(-1.9933797) q[0];
sx q[0];
rz(2.708486) q[0];
rz(-pi) q[1];
rz(-2.9559829) q[2];
sx q[2];
rz(-0.87366548) q[2];
sx q[2];
rz(2.6999161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.9351164) q[1];
sx q[1];
rz(-1.4644535) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9235839) q[3];
sx q[3];
rz(-0.70397607) q[3];
sx q[3];
rz(-2.9670027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(-2.5235126) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.9504257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(-1.114203) q[0];
x q[1];
rz(2.4345258) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(-2.6235839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.7813492) q[1];
sx q[1];
rz(1.7819808) q[1];
x q[2];
rz(1.1690087) q[3];
sx q[3];
rz(-2.0766633) q[3];
sx q[3];
rz(1.5798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8106666) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(1.3775685) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(0.92299443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92813338) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(2.7937828) q[0];
rz(-pi) q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.6147385) q[2];
sx q[2];
rz(2.5269788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7912485) q[1];
sx q[1];
rz(-1.4667759) q[1];
sx q[1];
rz(-1.9700325) q[1];
x q[2];
rz(-0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5287857) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.3607508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795792) q[0];
sx q[0];
rz(-0.29173457) q[0];
sx q[0];
rz(2.4289262) q[0];
x q[1];
rz(0.76061337) q[2];
sx q[2];
rz(-1.0906272) q[2];
sx q[2];
rz(-1.7825356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67140019) q[1];
sx q[1];
rz(-0.50506401) q[1];
sx q[1];
rz(2.7027674) q[1];
x q[2];
rz(-0.077386463) q[3];
sx q[3];
rz(-0.97273982) q[3];
sx q[3];
rz(-0.070904562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(-2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(0.27522603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21101418) q[0];
sx q[0];
rz(-1.0761346) q[0];
sx q[0];
rz(2.7633694) q[0];
rz(0.49595828) q[2];
sx q[2];
rz(-2.1514116) q[2];
sx q[2];
rz(-0.95326391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77995283) q[1];
sx q[1];
rz(-2.0556903) q[1];
sx q[1];
rz(1.1760902) q[1];
rz(-pi) q[2];
rz(1.6861378) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(-1.3184402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.747921) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(-0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-2.4218959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450711) q[0];
sx q[0];
rz(-2.4356027) q[0];
sx q[0];
rz(0.4476053) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(0.072349116) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3230348) q[1];
sx q[1];
rz(-2.5948988) q[1];
sx q[1];
rz(0.085303765) q[1];
x q[2];
rz(2.4607055) q[3];
sx q[3];
rz(-1.3586992) q[3];
sx q[3];
rz(1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(0.68816319) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(-2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.6465181) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197501) q[0];
sx q[0];
rz(-1.9103721) q[0];
sx q[0];
rz(1.3440488) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4453156) q[2];
sx q[2];
rz(-1.4204724) q[2];
sx q[2];
rz(-0.50762774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1284359) q[1];
sx q[1];
rz(-1.4669347) q[1];
sx q[1];
rz(2.8810112) q[1];
rz(-pi) q[2];
rz(2.2349615) q[3];
sx q[3];
rz(-0.88007054) q[3];
sx q[3];
rz(-3.0674792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(-1.7285255) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.7810129) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6920647) q[0];
sx q[0];
rz(-2.1421616) q[0];
sx q[0];
rz(3.0968303) q[0];
rz(-pi) q[1];
rz(-3.0948823) q[2];
sx q[2];
rz(-2.197406) q[2];
sx q[2];
rz(-1.4057297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4728002) q[1];
sx q[1];
rz(-1.2392514) q[1];
sx q[1];
rz(-1.9272976) q[1];
rz(-pi) q[2];
rz(-0.15501539) q[3];
sx q[3];
rz(-1.985637) q[3];
sx q[3];
rz(-0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90912949) q[0];
sx q[0];
rz(-1.4913519) q[0];
sx q[0];
rz(1.8151054) q[0];
rz(-0.017059762) q[2];
sx q[2];
rz(-0.3088726) q[2];
sx q[2];
rz(-0.19933867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9626999) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(3.1256413) q[1];
rz(-pi) q[2];
rz(-0.87749691) q[3];
sx q[3];
rz(-2.174456) q[3];
sx q[3];
rz(-2.7503777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7883041) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(-0.2164671) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-2.4344865) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(0.032565928) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
