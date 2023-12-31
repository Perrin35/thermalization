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
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(-2.2599028) q[0];
x q[1];
rz(2.2725265) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(-2.2609401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.213233) q[1];
sx q[1];
rz(-1.302049) q[1];
sx q[1];
rz(1.0006204) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22110181) q[3];
sx q[3];
rz(-3.0924774) q[3];
sx q[3];
rz(0.55288314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(3.0018905) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(-2.2712767) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(-1.93719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5088168) q[0];
sx q[0];
rz(-2.5460089) q[0];
sx q[0];
rz(-2.3217208) q[0];
x q[1];
rz(1.3538829) q[2];
sx q[2];
rz(-0.7173983) q[2];
sx q[2];
rz(2.4153828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8783427) q[1];
sx q[1];
rz(-1.471457) q[1];
sx q[1];
rz(-2.7753825) q[1];
x q[2];
rz(-2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(-1.1225835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(-0.016050054) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.191167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(-1.114203) q[0];
rz(-pi) q[1];
rz(-0.70706681) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(-2.6235839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9900073) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(-1.3596119) q[1];
rz(2.5268764) q[3];
sx q[3];
rz(-0.63496548) q[3];
sx q[3];
rz(0.8599417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3337293) q[0];
sx q[0];
rz(-2.6834052) q[0];
sx q[0];
rz(-2.3966167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.5268541) q[2];
sx q[2];
rz(0.61461385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.1715602) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.022425671) q[3];
sx q[3];
rz(-0.86655819) q[3];
sx q[3];
rz(-1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(0.64309684) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5287857) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(-1.3607508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0589941) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(2.9181705) q[0];
rz(-pi) q[1];
rz(-0.76061337) q[2];
sx q[2];
rz(-1.0906272) q[2];
sx q[2];
rz(-1.3590571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6319316) q[1];
sx q[1];
rz(-1.7778548) q[1];
sx q[1];
rz(-0.46405554) q[1];
rz(-pi) q[2];
rz(-0.077386463) q[3];
sx q[3];
rz(-2.1688528) q[3];
sx q[3];
rz(0.070904562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(3.0267267) q[2];
rz(0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(2.8663666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5462287) q[0];
sx q[0];
rz(-1.2397791) q[0];
sx q[0];
rz(-2.0966895) q[0];
rz(-pi) q[1];
rz(0.92991021) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(0.90604679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6323159) q[1];
sx q[1];
rz(-0.61513153) q[1];
sx q[1];
rz(-2.5110911) q[1];
x q[2];
rz(3.0589468) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(-1.1766528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5966406) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(2.4218959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73368209) q[0];
sx q[0];
rz(-0.94607293) q[0];
sx q[0];
rz(1.9242994) q[0];
rz(2.3709488) q[2];
sx q[2];
rz(-1.0174123) q[2];
sx q[2];
rz(1.9945952) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81855782) q[1];
sx q[1];
rz(-0.54669387) q[1];
sx q[1];
rz(-3.0562889) q[1];
x q[2];
rz(2.8119874) q[3];
sx q[3];
rz(-0.70809396) q[3];
sx q[3];
rz(-0.13347382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(2.8614614) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.4950745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4743487) q[0];
sx q[0];
rz(-1.357204) q[0];
sx q[0];
rz(-0.347802) q[0];
rz(-pi) q[1];
rz(1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(-1.1877103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5852768) q[1];
sx q[1];
rz(-1.3116515) q[1];
sx q[1];
rz(1.4633333) q[1];
rz(-pi) q[2];
rz(-2.5008194) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(-0.014766679) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(-1.3605798) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44952794) q[0];
sx q[0];
rz(-2.1421616) q[0];
sx q[0];
rz(0.044762386) q[0];
rz(0.94366818) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(2.949122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9578743) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(2.3493489) q[1];
rz(-pi) q[2];
rz(1.9901048) q[3];
sx q[3];
rz(-1.7125704) q[3];
sx q[3];
rz(-1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(-2.2488135) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(-0.7199026) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4601446) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(-0.081865099) q[0];
rz(-1.5762395) q[2];
sx q[2];
rz(-1.2619702) q[2];
sx q[2];
rz(-2.960161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1137005) q[1];
sx q[1];
rz(-2.8942331) q[1];
sx q[1];
rz(-1.6340096) q[1];
rz(-2.4107237) q[3];
sx q[3];
rz(-1.0167529) q[3];
sx q[3];
rz(2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.8855689) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532886) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
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
rz(2.1574216) q[3];
sx q[3];
rz(-1.543672) q[3];
sx q[3];
rz(-2.6437222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
