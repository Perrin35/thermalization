OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.202012) q[0];
sx q[0];
rz(-2.7913845) q[0];
sx q[0];
rz(0.36663088) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830084) q[0];
sx q[0];
rz(-1.2281679) q[0];
sx q[0];
rz(-1.1230099) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8777962) q[2];
sx q[2];
rz(-0.88636878) q[2];
sx q[2];
rz(-0.85927187) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3914324) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(2.0425914) q[1];
rz(1.5815758) q[3];
sx q[3];
rz(-1.5228776) q[3];
sx q[3];
rz(-2.3673494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(-2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.93719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28019529) q[0];
sx q[0];
rz(-1.178) q[0];
sx q[0];
rz(-1.1108206) q[0];
x q[1];
rz(0.86512489) q[2];
sx q[2];
rz(-1.428831) q[2];
sx q[2];
rz(-2.1324468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26324998) q[1];
sx q[1];
rz(-1.6701356) q[1];
sx q[1];
rz(0.36621014) q[1];
x q[2];
rz(-1.9235839) q[3];
sx q[3];
rz(-2.4376166) q[3];
sx q[3];
rz(2.9670027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(2.5235126) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.191167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-2.7911721) q[0];
sx q[0];
rz(1.114203) q[0];
rz(-pi) q[1];
rz(0.28352719) q[2];
sx q[2];
rz(-1.3361738) q[2];
sx q[2];
rz(-1.7262176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.7813492) q[1];
sx q[1];
rz(-1.3596119) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5268764) q[3];
sx q[3];
rz(-0.63496548) q[3];
sx q[3];
rz(2.281651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-1.0323662) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6072236) q[0];
sx q[0];
rz(-1.2396221) q[0];
sx q[0];
rz(-1.2481199) q[0];
rz(0.26206215) q[2];
sx q[2];
rz(-3.0960961) q[2];
sx q[2];
rz(2.2646575) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35034414) q[1];
sx q[1];
rz(-1.4667759) q[1];
sx q[1];
rz(1.1715602) q[1];
rz(1.5444078) q[3];
sx q[3];
rz(-0.70453405) q[3];
sx q[3];
rz(1.964331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-0.91941419) q[0];
rz(-2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.7808419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960984) q[0];
sx q[0];
rz(-1.7901663) q[0];
sx q[0];
rz(-1.7646837) q[0];
x q[1];
rz(0.76061337) q[2];
sx q[2];
rz(-2.0509655) q[2];
sx q[2];
rz(1.7825356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.67140019) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(0.43882521) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1702483) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(-1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(3.0267267) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-2.8663666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48588612) q[0];
sx q[0];
rz(-0.61300346) q[0];
sx q[0];
rz(0.97047752) q[0];
rz(-pi) q[1];
rz(-2.2116824) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(-2.2355459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77995283) q[1];
sx q[1];
rz(-2.0556903) q[1];
sx q[1];
rz(-1.9655025) q[1];
x q[2];
rz(0.082645881) q[3];
sx q[3];
rz(-0.9517037) q[3];
sx q[3];
rz(-1.1766528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-2.6829524) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(-2.5512364) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.77805) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(-0.71969676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62455432) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(0.65529234) q[0];
rz(-2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(3.0692435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67932018) q[1];
sx q[1];
rz(-1.6151036) q[1];
sx q[1];
rz(-0.54507749) q[1];
rz(-2.4607055) q[3];
sx q[3];
rz(-1.3586992) q[3];
sx q[3];
rz(-1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(0.68816319) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625531) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(-0.56674515) q[0];
x q[1];
rz(-2.9096793) q[2];
sx q[2];
rz(-2.4319318) q[2];
sx q[2];
rz(-2.2556925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.013156803) q[1];
sx q[1];
rz(-1.4669347) q[1];
sx q[1];
rz(-2.8810112) q[1];
rz(-pi) q[2];
rz(-2.3319578) q[3];
sx q[3];
rz(-1.0757043) q[3];
sx q[3];
rz(1.9593057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(-1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(-2.263608) q[0];
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
rz(-0.44952794) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(-0.044762386) q[0];
rz(2.1979245) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(2.949122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9578743) q[1];
sx q[1];
rz(-0.4819594) q[1];
sx q[1];
rz(-0.79224371) q[1];
rz(0.15501539) q[3];
sx q[3];
rz(-1.985637) q[3];
sx q[3];
rz(-2.993194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(-0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(-1.8887595) q[0];
x q[1];
rz(3.1245329) q[2];
sx q[2];
rz(-2.8327201) q[2];
sx q[2];
rz(-2.942254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7457909) q[1];
sx q[1];
rz(-1.5862641) q[1];
sx q[1];
rz(-1.8176816) q[1];
x q[2];
rz(0.87749691) q[3];
sx q[3];
rz(-2.174456) q[3];
sx q[3];
rz(-0.391215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(-1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532886) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(-1.6197694) q[3];
sx q[3];
rz(-2.5544142) q[3];
sx q[3];
rz(2.1094473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];