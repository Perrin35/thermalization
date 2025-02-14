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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(5.5362267) q[1];
sx q[1];
rz(9.8510392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881271) q[0];
sx q[0];
rz(-2.3511887) q[0];
sx q[0];
rz(1.3270204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6381016) q[2];
sx q[2];
rz(-1.3174743) q[2];
sx q[2];
rz(-1.5510786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1805686) q[1];
sx q[1];
rz(-2.5371309) q[1];
sx q[1];
rz(3.040041) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0996885) q[3];
sx q[3];
rz(-2.503241) q[3];
sx q[3];
rz(-1.9698576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8969741) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(-0.2429602) q[2];
rz(0.36863676) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.31545562) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(2.7742703) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(0.64250362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8272907) q[0];
sx q[0];
rz(-1.1513984) q[0];
sx q[0];
rz(-2.6227996) q[0];
rz(-pi) q[1];
rz(-1.1682061) q[2];
sx q[2];
rz(-1.0559096) q[2];
sx q[2];
rz(1.3001315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52740188) q[1];
sx q[1];
rz(-2.3878015) q[1];
sx q[1];
rz(0.68282737) q[1];
x q[2];
rz(-1.8090463) q[3];
sx q[3];
rz(-0.49360156) q[3];
sx q[3];
rz(0.77405888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6392886) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(1.3703692) q[2];
rz(-0.227452) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(-1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74533904) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(-1.1981717) q[0];
rz(-1.0802957) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(-0.094873039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8137275) q[0];
sx q[0];
rz(-1.828232) q[0];
sx q[0];
rz(-1.0125005) q[0];
x q[1];
rz(0.91135613) q[2];
sx q[2];
rz(-1.240864) q[2];
sx q[2];
rz(1.0418237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.171189) q[1];
sx q[1];
rz(-0.53817777) q[1];
sx q[1];
rz(-1.490905) q[1];
x q[2];
rz(-1.9547988) q[3];
sx q[3];
rz(-1.8285164) q[3];
sx q[3];
rz(-0.4674165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0723116) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(1.9629924) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.9788205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5793295) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(0.33600268) q[0];
rz(-2.621189) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(3.0373108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1395499) q[0];
sx q[0];
rz(-2.6536396) q[0];
sx q[0];
rz(-2.3298708) q[0];
x q[1];
rz(2.9851172) q[2];
sx q[2];
rz(-1.5847407) q[2];
sx q[2];
rz(0.61323159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6426864) q[1];
sx q[1];
rz(-2.3489174) q[1];
sx q[1];
rz(2.098987) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6882668) q[3];
sx q[3];
rz(-1.4651872) q[3];
sx q[3];
rz(-2.3345514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5229554) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.1445507) q[2];
rz(0.54730225) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8892141) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(0.55437535) q[0];
rz(-0.91122183) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(-2.8401781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2269991) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(0.27497681) q[0];
x q[1];
rz(-2.7533135) q[2];
sx q[2];
rz(-1.7634321) q[2];
sx q[2];
rz(-1.9658058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0423066) q[1];
sx q[1];
rz(-1.9507244) q[1];
sx q[1];
rz(-2.9126588) q[1];
x q[2];
rz(1.9456995) q[3];
sx q[3];
rz(-1.8640567) q[3];
sx q[3];
rz(2.9305262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(0.0069228355) q[2];
rz(-2.210468) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865006) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(1.5337926) q[0];
rz(2.4389229) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-2.839397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5836499) q[0];
sx q[0];
rz(-1.4690225) q[0];
sx q[0];
rz(-0.90910615) q[0];
x q[1];
rz(2.7415761) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(0.053002593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8540878) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(-1.8604398) q[1];
rz(-pi) q[2];
rz(-1.4190361) q[3];
sx q[3];
rz(-1.6223063) q[3];
sx q[3];
rz(-0.53799483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89221421) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(1.2146436) q[3];
sx q[3];
rz(-2.6483783) q[3];
sx q[3];
rz(-0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(2.0330698) q[0];
rz(0.33379894) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(2.0557859) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54667771) q[0];
sx q[0];
rz(-1.6502066) q[0];
sx q[0];
rz(-2.9360442) q[0];
x q[1];
rz(1.9721617) q[2];
sx q[2];
rz(-1.4039206) q[2];
sx q[2];
rz(2.3292975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8657596) q[1];
sx q[1];
rz(-1.0414904) q[1];
sx q[1];
rz(2.5342026) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0039075) q[3];
sx q[3];
rz(-0.9387278) q[3];
sx q[3];
rz(1.1732303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(1.9773352) q[2];
rz(-2.8734251) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(-2.3837762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6044354) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(1.1489768) q[0];
rz(-2.8474498) q[1];
sx q[1];
rz(-1.5584757) q[1];
sx q[1];
rz(1.0424967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71871686) q[0];
sx q[0];
rz(-1.9008753) q[0];
sx q[0];
rz(0.80083682) q[0];
x q[1];
rz(0.22886054) q[2];
sx q[2];
rz(-2.3386152) q[2];
sx q[2];
rz(-0.91509089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62952215) q[1];
sx q[1];
rz(-2.7334556) q[1];
sx q[1];
rz(-0.44993181) q[1];
rz(-pi) q[2];
rz(1.5339666) q[3];
sx q[3];
rz(-0.66434089) q[3];
sx q[3];
rz(-0.0078594154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59757549) q[2];
sx q[2];
rz(-2.47051) q[2];
sx q[2];
rz(-1.3647122) q[2];
rz(0.65258604) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.515601) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(3.1248202) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(-2.9877072) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(-0.31799728) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1915386) q[2];
sx q[2];
rz(-2.1851106) q[2];
sx q[2];
rz(1.4221734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2528884) q[1];
sx q[1];
rz(-1.9333464) q[1];
sx q[1];
rz(-1.8860555) q[1];
rz(-1.8762053) q[3];
sx q[3];
rz(-3.0620259) q[3];
sx q[3];
rz(-0.1480535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1880356) q[2];
sx q[2];
rz(-0.81622684) q[2];
sx q[2];
rz(-0.44450644) q[2];
rz(2.9514173) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1403777) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(2.0709399) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(0.79107034) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2915724) q[0];
sx q[0];
rz(-2.794461) q[0];
sx q[0];
rz(1.8308543) q[0];
rz(2.5612381) q[2];
sx q[2];
rz(-0.46541801) q[2];
sx q[2];
rz(2.4632955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.715383) q[1];
sx q[1];
rz(-1.2150303) q[1];
sx q[1];
rz(-0.4224311) q[1];
rz(1.3338575) q[3];
sx q[3];
rz(-2.624369) q[3];
sx q[3];
rz(0.57843971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4967686) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(1.6688639) q[2];
rz(1.5276927) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(-1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.938217) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(2.6381941) q[2];
sx q[2];
rz(-1.361327) q[2];
sx q[2];
rz(0.8045902) q[2];
rz(-0.41889965) q[3];
sx q[3];
rz(-2.2875026) q[3];
sx q[3];
rz(0.84411375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
