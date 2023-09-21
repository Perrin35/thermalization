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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(-2.2599028) q[0];
x q[1];
rz(1.2615112) q[2];
sx q[2];
rz(-0.72578428) q[2];
sx q[2];
rz(-2.6860565) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.213233) q[1];
sx q[1];
rz(-1.302049) q[1];
sx q[1];
rz(2.1409722) q[1];
rz(1.5600169) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(0.39164266) q[2];
rz(3.0018905) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4166819) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(-2.2712767) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63277584) q[0];
sx q[0];
rz(-2.5460089) q[0];
sx q[0];
rz(0.81987183) q[0];
rz(-pi) q[1];
rz(-1.7877098) q[2];
sx q[2];
rz(-2.4241944) q[2];
sx q[2];
rz(0.72620981) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26324998) q[1];
sx q[1];
rz(-1.471457) q[1];
sx q[1];
rz(-2.7753825) q[1];
rz(-1.2180087) q[3];
sx q[3];
rz(-2.4376166) q[3];
sx q[3];
rz(0.17458992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59042674) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.191167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58142692) q[0];
sx q[0];
rz(-1.7227356) q[0];
sx q[0];
rz(1.8877958) q[0];
rz(-0.28352719) q[2];
sx q[2];
rz(-1.3361738) q[2];
sx q[2];
rz(-1.4153751) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.788113) q[1];
sx q[1];
rz(-0.2971) q[1];
sx q[1];
rz(-0.77570813) q[1];
x q[2];
rz(1.972584) q[3];
sx q[3];
rz(-1.0649293) q[3];
sx q[3];
rz(-1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-2.9372835) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-0.92299443) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2134593) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(2.7937828) q[0];
rz(-pi) q[1];
rz(-1.5590018) q[2];
sx q[2];
rz(-1.5268541) q[2];
sx q[2];
rz(-0.61461385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.877362) q[1];
sx q[1];
rz(-1.1738395) q[1];
sx q[1];
rz(3.0287659) q[1];
x q[2];
rz(-2.2751585) q[3];
sx q[3];
rz(-1.5878864) q[3];
sx q[3];
rz(2.7279502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-0.50951177) q[2];
rz(-2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(-0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(-1.3607508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0589941) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(0.22342213) q[0];
x q[1];
rz(-2.4945716) q[2];
sx q[2];
rz(-0.87304742) q[2];
sx q[2];
rz(-2.9014212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67140019) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(2.7027674) q[1];
rz(0.077386463) q[3];
sx q[3];
rz(-0.97273982) q[3];
sx q[3];
rz(-3.0706881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78648606) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(2.8663666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.595364) q[0];
sx q[0];
rz(-1.2397791) q[0];
sx q[0];
rz(-1.0449032) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.198344) q[2];
sx q[2];
rz(-2.3970282) q[2];
sx q[2];
rz(-0.17472357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77995283) q[1];
sx q[1];
rz(-2.0556903) q[1];
sx q[1];
rz(1.1760902) q[1];
rz(-pi) q[2];
rz(0.082645881) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(1.1766528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(-2.4300872) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
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
rz(1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-0.71969676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62455432) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(2.4863003) q[0];
rz(-pi) q[1];
rz(0.85992203) q[2];
sx q[2];
rz(-2.2051174) q[2];
sx q[2];
rz(0.89564043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2232659) q[1];
sx q[1];
rz(-2.1152788) q[1];
sx q[1];
rz(-1.5189927) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68088716) q[3];
sx q[3];
rz(-1.3586992) q[3];
sx q[3];
rz(1.4500446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-2.4534295) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(-0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(0.18653175) q[0];
rz(-0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.4950745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218425) q[0];
sx q[0];
rz(-1.9103721) q[0];
sx q[0];
rz(-1.7975438) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.375884) q[2];
sx q[2];
rz(-2.2576828) q[2];
sx q[2];
rz(1.9538823) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1868134) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(0.38444744) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64077326) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(-1.4130672) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(-0.87798464) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(-1.7810129) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6920647) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(0.044762386) q[0];
rz(-pi) q[1];
rz(-3.0948823) q[2];
sx q[2];
rz(-2.197406) q[2];
sx q[2];
rz(1.7358629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9578743) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(-0.79224371) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2336041) q[3];
sx q[3];
rz(-2.7003151) q[3];
sx q[3];
rz(-0.518276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(-3.0950586) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(-0.7199026) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35346183) q[0];
sx q[0];
rz(-2.8849368) q[0];
sx q[0];
rz(1.8887595) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5762395) q[2];
sx q[2];
rz(-1.2619702) q[2];
sx q[2];
rz(-0.18143166) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1788927) q[1];
sx q[1];
rz(-1.8176515) q[1];
sx q[1];
rz(0.015951338) q[1];
rz(-pi) q[2];
rz(2.3941444) q[3];
sx q[3];
rz(-0.8851074) q[3];
sx q[3];
rz(-1.3626584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7883041) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(0.36454501) q[2];
sx q[2];
rz(-0.74082965) q[2];
sx q[2];
rz(-0.34168591) q[2];
rz(-1.5218232) q[3];
sx q[3];
rz(-0.58717849) q[3];
sx q[3];
rz(-1.0321454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
