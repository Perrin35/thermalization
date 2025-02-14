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
rz(2.6142081) q[0];
sx q[0];
rz(-2.5684147) q[0];
sx q[0];
rz(-0.61668026) q[0];
rz(-1.0083899) q[1];
sx q[1];
rz(3.8808793) q[1];
sx q[1];
rz(9.0864658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36114281) q[0];
sx q[0];
rz(-0.2451714) q[0];
sx q[0];
rz(2.992379) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7912999) q[2];
sx q[2];
rz(-1.2630579) q[2];
sx q[2];
rz(-1.6215633) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.68012) q[1];
sx q[1];
rz(-1.0306038) q[1];
sx q[1];
rz(1.4137181) q[1];
x q[2];
rz(-1.7565207) q[3];
sx q[3];
rz(-0.76558569) q[3];
sx q[3];
rz(-3.035745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7484635) q[2];
sx q[2];
rz(-1.7270361) q[2];
sx q[2];
rz(0.65790042) q[2];
rz(-0.45514485) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(-1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667792) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(-0.28208062) q[0];
rz(2.8981949) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(2.3928221) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3385612) q[0];
sx q[0];
rz(-2.411532) q[0];
sx q[0];
rz(0.29166136) q[0];
rz(-pi) q[1];
rz(1.44913) q[2];
sx q[2];
rz(-1.5783383) q[2];
sx q[2];
rz(-2.6040524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48909125) q[1];
sx q[1];
rz(-1.4785119) q[1];
sx q[1];
rz(1.0670061) q[1];
rz(-pi) q[2];
rz(1.5010819) q[3];
sx q[3];
rz(-1.3718714) q[3];
sx q[3];
rz(2.6903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28027174) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(-1.8005499) q[2];
rz(1.1567814) q[3];
sx q[3];
rz(-0.78514922) q[3];
sx q[3];
rz(2.3780499) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84083104) q[0];
sx q[0];
rz(-0.16591993) q[0];
sx q[0];
rz(-0.65735835) q[0];
rz(1.9961458) q[1];
sx q[1];
rz(-2.1871388) q[1];
sx q[1];
rz(-0.81833902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79761222) q[0];
sx q[0];
rz(-1.3320141) q[0];
sx q[0];
rz(0.35399951) q[0];
rz(-pi) q[1];
rz(-1.4253699) q[2];
sx q[2];
rz(-0.31842318) q[2];
sx q[2];
rz(0.91998902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0967674) q[1];
sx q[1];
rz(-0.92508537) q[1];
sx q[1];
rz(-0.94618465) q[1];
x q[2];
rz(1.7210427) q[3];
sx q[3];
rz(-0.45805061) q[3];
sx q[3];
rz(-0.86298215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98693097) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(1.7884802) q[2];
rz(-0.96495572) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.2009784) q[0];
rz(3.1383842) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(2.9046955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85962109) q[0];
sx q[0];
rz(-2.1016723) q[0];
sx q[0];
rz(2.6174699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7014933) q[2];
sx q[2];
rz(-2.2039737) q[2];
sx q[2];
rz(0.53772949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5665) q[1];
sx q[1];
rz(-1.3761576) q[1];
sx q[1];
rz(0.95928488) q[1];
rz(-pi) q[2];
rz(0.8221738) q[3];
sx q[3];
rz(-2.195925) q[3];
sx q[3];
rz(-1.353598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0995471) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(2.715204) q[2];
rz(1.03553) q[3];
sx q[3];
rz(-1.4617045) q[3];
sx q[3];
rz(0.6216014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7972888) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(1.9328851) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(-0.0016317687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372555) q[0];
sx q[0];
rz(-1.4044824) q[0];
sx q[0];
rz(-2.8939899) q[0];
rz(-pi) q[1];
rz(-3.0157005) q[2];
sx q[2];
rz(-2.841904) q[2];
sx q[2];
rz(2.3131049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.092068191) q[1];
sx q[1];
rz(-0.18316575) q[1];
sx q[1];
rz(1.3080255) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21065055) q[3];
sx q[3];
rz(-1.306476) q[3];
sx q[3];
rz(1.5791062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95167595) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(-2.9803989) q[2];
rz(0.4392043) q[3];
sx q[3];
rz(-0.74857155) q[3];
sx q[3];
rz(0.85328931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31140232) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(0.79469529) q[0];
rz(-1.7757802) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(2.6672003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8891539) q[0];
sx q[0];
rz(-1.948101) q[0];
sx q[0];
rz(-1.0214367) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2016868) q[2];
sx q[2];
rz(-2.0966736) q[2];
sx q[2];
rz(2.3591218) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5833334) q[1];
sx q[1];
rz(-2.4129902) q[1];
sx q[1];
rz(0.63459756) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7571449) q[3];
sx q[3];
rz(-1.0214503) q[3];
sx q[3];
rz(2.5540496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5279493) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(-0.5274241) q[2];
rz(2.534965) q[3];
sx q[3];
rz(-2.9546723) q[3];
sx q[3];
rz(-2.0691779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793176) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(0.79757565) q[0];
rz(2.2480615) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(-2.8614047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2036616) q[0];
sx q[0];
rz(-1.9660006) q[0];
sx q[0];
rz(0.18919887) q[0];
x q[1];
rz(0.65164874) q[2];
sx q[2];
rz(-2.0927395) q[2];
sx q[2];
rz(2.3892982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1692062) q[1];
sx q[1];
rz(-1.5582339) q[1];
sx q[1];
rz(2.0533086) q[1];
rz(2.7817621) q[3];
sx q[3];
rz(-1.0092648) q[3];
sx q[3];
rz(2.7240745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4013275) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(1.9291482) q[2];
rz(-2.9017743) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(0.56813204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78385335) q[0];
sx q[0];
rz(-1.724406) q[0];
sx q[0];
rz(1.3215815) q[0];
rz(2.9603738) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(0.063057335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1940799) q[0];
sx q[0];
rz(-0.79160684) q[0];
sx q[0];
rz(0.20157728) q[0];
rz(-pi) q[1];
rz(3.0538043) q[2];
sx q[2];
rz(-1.286418) q[2];
sx q[2];
rz(-1.2149917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2109697) q[1];
sx q[1];
rz(-1.9038868) q[1];
sx q[1];
rz(1.3107783) q[1];
rz(-pi) q[2];
rz(-0.0029011528) q[3];
sx q[3];
rz(-1.408268) q[3];
sx q[3];
rz(1.3292461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41398373) q[2];
sx q[2];
rz(-1.0980282) q[2];
sx q[2];
rz(1.6723527) q[2];
rz(-1.7478878) q[3];
sx q[3];
rz(-1.7017476) q[3];
sx q[3];
rz(-2.9318455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6032228) q[0];
sx q[0];
rz(-0.61036888) q[0];
sx q[0];
rz(0.99697733) q[0];
rz(0.793055) q[1];
sx q[1];
rz(-1.4184364) q[1];
sx q[1];
rz(1.1096032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5376268) q[0];
sx q[0];
rz(-1.0111799) q[0];
sx q[0];
rz(-1.2252612) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47596653) q[2];
sx q[2];
rz(-0.78078166) q[2];
sx q[2];
rz(2.4300574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3396757) q[1];
sx q[1];
rz(-2.0387421) q[1];
sx q[1];
rz(0.94593327) q[1];
rz(2.505198) q[3];
sx q[3];
rz(-1.706746) q[3];
sx q[3];
rz(-3.1042447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7290466) q[2];
sx q[2];
rz(-1.46526) q[2];
sx q[2];
rz(1.868978) q[2];
rz(2.0160969) q[3];
sx q[3];
rz(-0.19386217) q[3];
sx q[3];
rz(0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37128714) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(3.0433997) q[0];
rz(-1.6432317) q[1];
sx q[1];
rz(-1.9941092) q[1];
sx q[1];
rz(2.3032545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9591374) q[0];
sx q[0];
rz(-2.3172624) q[0];
sx q[0];
rz(-1.8126382) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25213253) q[2];
sx q[2];
rz(-1.6246535) q[2];
sx q[2];
rz(-2.5359652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3915836) q[1];
sx q[1];
rz(-2.6152088) q[1];
sx q[1];
rz(-0.98441225) q[1];
rz(-2.7097852) q[3];
sx q[3];
rz(-1.3696967) q[3];
sx q[3];
rz(-0.74936282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3457634) q[2];
sx q[2];
rz(-2.2553307) q[2];
sx q[2];
rz(-0.32540992) q[2];
rz(-2.365153) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(2.5344892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008739) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(2.8554032) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(0.4109702) q[2];
sx q[2];
rz(-2.3373418) q[2];
sx q[2];
rz(1.892754) q[2];
rz(-2.552916) q[3];
sx q[3];
rz(-1.2856144) q[3];
sx q[3];
rz(-1.9348617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
