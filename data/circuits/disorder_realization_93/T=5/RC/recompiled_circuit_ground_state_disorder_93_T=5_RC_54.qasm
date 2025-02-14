OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8146347) q[0];
sx q[0];
rz(3.1231413) q[0];
sx q[0];
rz(9.2508247) q[0];
rz(1.8040909) q[1];
sx q[1];
rz(4.017158) q[1];
sx q[1];
rz(12.234966) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1179292) q[0];
sx q[0];
rz(-0.34043306) q[0];
sx q[0];
rz(2.6871919) q[0];
x q[1];
rz(0.5919257) q[2];
sx q[2];
rz(-2.3560696) q[2];
sx q[2];
rz(-0.80410081) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9921205) q[1];
sx q[1];
rz(-1.7594584) q[1];
sx q[1];
rz(0.68575286) q[1];
rz(-2.9339497) q[3];
sx q[3];
rz(-2.9831726) q[3];
sx q[3];
rz(-2.0655565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6203674) q[2];
sx q[2];
rz(-1.8769033) q[2];
sx q[2];
rz(2.550726) q[2];
rz(1.2938195) q[3];
sx q[3];
rz(-1.7655244) q[3];
sx q[3];
rz(-2.6025492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6227601) q[0];
sx q[0];
rz(-2.9157214) q[0];
sx q[0];
rz(0.87509218) q[0];
rz(-0.092763364) q[1];
sx q[1];
rz(-2.0848672) q[1];
sx q[1];
rz(1.0088395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1160285) q[0];
sx q[0];
rz(-0.43034986) q[0];
sx q[0];
rz(-1.5034281) q[0];
rz(-pi) q[1];
rz(-2.5549803) q[2];
sx q[2];
rz(-1.5924708) q[2];
sx q[2];
rz(2.7119111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0663026) q[1];
sx q[1];
rz(-1.7199893) q[1];
sx q[1];
rz(0.91764024) q[1];
rz(-pi) q[2];
rz(-3.1280531) q[3];
sx q[3];
rz(-1.6289181) q[3];
sx q[3];
rz(-2.6839369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2729317) q[2];
sx q[2];
rz(-1.2596143) q[2];
sx q[2];
rz(-1.3966365) q[2];
rz(2.7331288) q[3];
sx q[3];
rz(-0.68578458) q[3];
sx q[3];
rz(2.7062611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8946266) q[0];
sx q[0];
rz(-2.8530687) q[0];
sx q[0];
rz(0.95046473) q[0];
rz(-1.8679999) q[1];
sx q[1];
rz(-1.7968977) q[1];
sx q[1];
rz(1.7705932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0799373) q[0];
sx q[0];
rz(-2.2191471) q[0];
sx q[0];
rz(-1.7440304) q[0];
x q[1];
rz(-1.6865963) q[2];
sx q[2];
rz(-2.0017923) q[2];
sx q[2];
rz(0.64799448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1493686) q[1];
sx q[1];
rz(-2.0724587) q[1];
sx q[1];
rz(-2.9976974) q[1];
rz(1.4617861) q[3];
sx q[3];
rz(-0.9338769) q[3];
sx q[3];
rz(-3.1240841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3889968) q[2];
sx q[2];
rz(-1.4138556) q[2];
sx q[2];
rz(-1.625212) q[2];
rz(-0.73836941) q[3];
sx q[3];
rz(-1.6298031) q[3];
sx q[3];
rz(-0.71877688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3589288) q[0];
sx q[0];
rz(-1.6343225) q[0];
sx q[0];
rz(-1.9184817) q[0];
rz(-2.3861528) q[1];
sx q[1];
rz(-2.0016045) q[1];
sx q[1];
rz(3.0208407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31407315) q[0];
sx q[0];
rz(-1.733092) q[0];
sx q[0];
rz(1.4050964) q[0];
rz(0.37379548) q[2];
sx q[2];
rz(-0.5875365) q[2];
sx q[2];
rz(-0.85174207) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9452136) q[1];
sx q[1];
rz(-2.2506013) q[1];
sx q[1];
rz(0.86974135) q[1];
rz(-pi) q[2];
rz(1.5232682) q[3];
sx q[3];
rz(-1.589425) q[3];
sx q[3];
rz(2.9513074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52189031) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(1.5187368) q[2];
rz(1.3269199) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-2.0324619) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0066852) q[0];
sx q[0];
rz(-1.4445855) q[0];
sx q[0];
rz(2.4476442) q[0];
rz(0.42849439) q[1];
sx q[1];
rz(-0.65281147) q[1];
sx q[1];
rz(-1.8983715) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4994298) q[0];
sx q[0];
rz(-1.3754396) q[0];
sx q[0];
rz(3.0950154) q[0];
x q[1];
rz(-0.17036338) q[2];
sx q[2];
rz(-1.7936094) q[2];
sx q[2];
rz(2.6210268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7299594) q[1];
sx q[1];
rz(-0.45242971) q[1];
sx q[1];
rz(-2.1001454) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6596131) q[3];
sx q[3];
rz(-1.6631806) q[3];
sx q[3];
rz(2.5661812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1474233) q[2];
sx q[2];
rz(-1.3012412) q[2];
sx q[2];
rz(0.48747882) q[2];
rz(-0.068988919) q[3];
sx q[3];
rz(-2.4853849) q[3];
sx q[3];
rz(2.5330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92031205) q[0];
sx q[0];
rz(-2.1423036) q[0];
sx q[0];
rz(2.6699303) q[0];
rz(-0.77700514) q[1];
sx q[1];
rz(-1.124137) q[1];
sx q[1];
rz(1.5583001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9474424) q[0];
sx q[0];
rz(-1.5628038) q[0];
sx q[0];
rz(-3.1119909) q[0];
rz(-pi) q[1];
rz(-0.92783096) q[2];
sx q[2];
rz(-0.6428679) q[2];
sx q[2];
rz(-0.57807482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1123429) q[1];
sx q[1];
rz(-1.3549374) q[1];
sx q[1];
rz(2.6964746) q[1];
x q[2];
rz(-0.72649148) q[3];
sx q[3];
rz(-1.3950751) q[3];
sx q[3];
rz(-0.54783152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9787057) q[2];
sx q[2];
rz(-1.6521613) q[2];
sx q[2];
rz(-1.9784652) q[2];
rz(0.85824054) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(-0.8391909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100886) q[0];
sx q[0];
rz(-1.6781582) q[0];
sx q[0];
rz(0.32291821) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-1.1406621) q[1];
sx q[1];
rz(1.742977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271287) q[0];
sx q[0];
rz(-1.4076978) q[0];
sx q[0];
rz(-2.968279) q[0];
x q[1];
rz(-0.9594907) q[2];
sx q[2];
rz(-1.2388907) q[2];
sx q[2];
rz(1.6417981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7306644) q[1];
sx q[1];
rz(-2.321744) q[1];
sx q[1];
rz(-1.5265205) q[1];
x q[2];
rz(-1.4780462) q[3];
sx q[3];
rz(-1.6997011) q[3];
sx q[3];
rz(-0.034667548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6814331) q[2];
sx q[2];
rz(-0.56542772) q[2];
sx q[2];
rz(-2.4088755) q[2];
rz(3.1198464) q[3];
sx q[3];
rz(-1.6201868) q[3];
sx q[3];
rz(0.17797962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1499504) q[0];
sx q[0];
rz(-3.0198779) q[0];
sx q[0];
rz(-0.44418401) q[0];
rz(1.3581879) q[1];
sx q[1];
rz(-1.5779747) q[1];
sx q[1];
rz(-1.1202728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60591423) q[0];
sx q[0];
rz(-1.6736503) q[0];
sx q[0];
rz(0.82983183) q[0];
x q[1];
rz(0.48768576) q[2];
sx q[2];
rz(-1.2151866) q[2];
sx q[2];
rz(-1.6297225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7071735) q[1];
sx q[1];
rz(-1.4212233) q[1];
sx q[1];
rz(-2.892848) q[1];
rz(-pi) q[2];
rz(-1.9838263) q[3];
sx q[3];
rz(-1.1404697) q[3];
sx q[3];
rz(-2.4347507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.292395) q[2];
sx q[2];
rz(-0.85453832) q[2];
sx q[2];
rz(-1.9067859) q[2];
rz(-0.97709996) q[3];
sx q[3];
rz(-1.2950725) q[3];
sx q[3];
rz(-0.60404122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6848171) q[0];
sx q[0];
rz(-0.026739459) q[0];
sx q[0];
rz(-0.34715125) q[0];
rz(-0.30255643) q[1];
sx q[1];
rz(-1.1056933) q[1];
sx q[1];
rz(0.0060630719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0971507) q[0];
sx q[0];
rz(-1.7527075) q[0];
sx q[0];
rz(3.1245435) q[0];
rz(-pi) q[1];
rz(-1.9337186) q[2];
sx q[2];
rz(-2.4293682) q[2];
sx q[2];
rz(2.1729529) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4655759) q[1];
sx q[1];
rz(-2.1853059) q[1];
sx q[1];
rz(0.37218233) q[1];
x q[2];
rz(-2.361287) q[3];
sx q[3];
rz(-1.5374476) q[3];
sx q[3];
rz(2.1005971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89846316) q[2];
sx q[2];
rz(-0.96776882) q[2];
sx q[2];
rz(1.751162) q[2];
rz(-0.27680382) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(-2.0738585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7368363) q[0];
sx q[0];
rz(-0.75045776) q[0];
sx q[0];
rz(0.87338895) q[0];
rz(0.5492754) q[1];
sx q[1];
rz(-0.6540238) q[1];
sx q[1];
rz(-0.80680791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0334369) q[0];
sx q[0];
rz(-1.8315541) q[0];
sx q[0];
rz(-2.5013862) q[0];
rz(-pi) q[1];
rz(-0.64555706) q[2];
sx q[2];
rz(-2.9462313) q[2];
sx q[2];
rz(2.8091009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6178231) q[1];
sx q[1];
rz(-1.4549991) q[1];
sx q[1];
rz(-1.3208742) q[1];
rz(-0.30183123) q[3];
sx q[3];
rz(-2.1595397) q[3];
sx q[3];
rz(1.6993239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6640545) q[2];
sx q[2];
rz(-1.3128023) q[2];
sx q[2];
rz(0.99267268) q[2];
rz(2.4382675) q[3];
sx q[3];
rz(-2.0086918) q[3];
sx q[3];
rz(0.68737427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138388) q[0];
sx q[0];
rz(-2.3532372) q[0];
sx q[0];
rz(-1.1113356) q[0];
rz(-2.5913024) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(0.69725488) q[2];
sx q[2];
rz(-1.6959126) q[2];
sx q[2];
rz(-0.87903862) q[2];
rz(0.46389085) q[3];
sx q[3];
rz(-1.9474467) q[3];
sx q[3];
rz(1.8578702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
