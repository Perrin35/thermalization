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
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.769387) q[0];
sx q[0];
rz(-1.9908449) q[0];
sx q[0];
rz(-0.37680349) q[0];
x q[1];
rz(-1.2615112) q[2];
sx q[2];
rz(-2.4158084) q[2];
sx q[2];
rz(-2.6860565) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75016025) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(-1.0990012) q[1];
rz(-0.047921501) q[3];
sx q[3];
rz(-1.5600292) q[3];
sx q[3];
rz(-0.79706942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1047487) q[2];
sx q[2];
rz(-1.8061183) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(-0.87031594) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6635839) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(0.43310662) q[0];
rz(-1.3538829) q[2];
sx q[2];
rz(-2.4241944) q[2];
sx q[2];
rz(2.4153828) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.9351164) q[1];
sx q[1];
rz(1.6771392) q[1];
x q[2];
rz(0.89800091) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(2.0190092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.191167) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1026099) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(0.15977504) q[0];
rz(-2.8580655) q[2];
sx q[2];
rz(-1.8054188) q[2];
sx q[2];
rz(-1.4153751) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6776059) q[1];
sx q[1];
rz(-1.7772487) q[1];
sx q[1];
rz(0.21519214) q[1];
x q[2];
rz(2.5997945) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(1.013914) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(-1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(2.9372835) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(2.2185982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3337293) q[0];
sx q[0];
rz(-2.6834052) q[0];
sx q[0];
rz(0.74497594) q[0];
rz(-pi) q[1];
rz(-0.26206215) q[2];
sx q[2];
rz(-0.045496551) q[2];
sx q[2];
rz(2.2646575) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1624565) q[1];
sx q[1];
rz(-0.41185954) q[1];
sx q[1];
rz(-1.3084175) q[1];
rz(-0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(-0.91941419) q[0];
rz(-2.1557504) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(1.3607508) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4454942) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(-1.7646837) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76061337) q[2];
sx q[2];
rz(-1.0906272) q[2];
sx q[2];
rz(-1.3590571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6319316) q[1];
sx q[1];
rz(-1.7778548) q[1];
sx q[1];
rz(-0.46405554) q[1];
rz(-0.077386463) q[3];
sx q[3];
rz(-0.97273982) q[3];
sx q[3];
rz(3.0706881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(-0.45587513) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.8967569) q[0];
rz(-2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(0.27522603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48588612) q[0];
sx q[0];
rz(-0.61300346) q[0];
sx q[0];
rz(-2.1711151) q[0];
x q[1];
rz(-0.49595828) q[2];
sx q[2];
rz(-2.1514116) q[2];
sx q[2];
rz(0.95326391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6323159) q[1];
sx q[1];
rz(-2.5264611) q[1];
sx q[1];
rz(-2.5110911) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95008696) q[3];
sx q[3];
rz(-1.5035149) q[3];
sx q[3];
rz(-2.7954807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-0.45864027) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450711) q[0];
sx q[0];
rz(-0.70598999) q[0];
sx q[0];
rz(0.4476053) q[0];
x q[1];
rz(2.3709488) q[2];
sx q[2];
rz(-2.1241803) q[2];
sx q[2];
rz(-1.9945952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81855782) q[1];
sx q[1];
rz(-0.54669387) q[1];
sx q[1];
rz(0.085303765) q[1];
x q[2];
rz(1.8411438) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(-0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625531) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(-2.5748475) q[0];
x q[1];
rz(0.23191339) q[2];
sx q[2];
rz(-2.4319318) q[2];
sx q[2];
rz(-2.2556925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1868134) q[1];
sx q[1];
rz(-0.28007945) q[1];
sx q[1];
rz(0.38444744) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90663119) q[3];
sx q[3];
rz(-0.88007054) q[3];
sx q[3];
rz(-0.074113473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8994393) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(3.126826) q[2];
rz(1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(-1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(0.87798464) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(-1.3605798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44952794) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(0.044762386) q[0];
rz(0.94366818) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(0.19247069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0226375) q[1];
sx q[1];
rz(-1.2345018) q[1];
sx q[1];
rz(2.789546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9079886) q[3];
sx q[3];
rz(-0.4412776) q[3];
sx q[3];
rz(0.518276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-2.2488135) q[2];
rz(-0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601446) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(-3.0597276) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5762395) q[2];
sx q[2];
rz(-1.8796225) q[2];
sx q[2];
rz(2.960161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9626999) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(-3.1256413) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2640957) q[3];
sx q[3];
rz(-2.174456) q[3];
sx q[3];
rz(-2.7503777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.8855689) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-0.2164671) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(-3.1090267) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
