OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(-1.1450333) q[0];
sx q[0];
rz(1.0166919) q[0];
rz(2.297205) q[1];
sx q[1];
rz(-2.3250186) q[1];
sx q[1];
rz(-0.56632298) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53439683) q[0];
sx q[0];
rz(-1.6677987) q[0];
sx q[0];
rz(2.1792381) q[0];
rz(-1.3728605) q[2];
sx q[2];
rz(-1.4948442) q[2];
sx q[2];
rz(-1.4966485) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6366292) q[1];
sx q[1];
rz(-0.43049059) q[1];
sx q[1];
rz(0.19383793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3641961) q[3];
sx q[3];
rz(-0.99541621) q[3];
sx q[3];
rz(-2.8386338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6301959) q[2];
sx q[2];
rz(-0.87367311) q[2];
sx q[2];
rz(1.1279747) q[2];
rz(1.5026211) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(-2.716841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2750435) q[0];
sx q[0];
rz(-2.6847222) q[0];
sx q[0];
rz(-3.1006815) q[0];
rz(2.5919137) q[1];
sx q[1];
rz(-0.37893852) q[1];
sx q[1];
rz(0.080370195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59118862) q[0];
sx q[0];
rz(-2.4192717) q[0];
sx q[0];
rz(0.0063615464) q[0];
rz(-0.88826142) q[2];
sx q[2];
rz(-2.3739034) q[2];
sx q[2];
rz(1.7861799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.077420635) q[1];
sx q[1];
rz(-1.5158093) q[1];
sx q[1];
rz(1.085678) q[1];
x q[2];
rz(-1.188813) q[3];
sx q[3];
rz(-0.63900164) q[3];
sx q[3];
rz(1.5572214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72767672) q[2];
sx q[2];
rz(-0.84299403) q[2];
sx q[2];
rz(-0.28398871) q[2];
rz(1.5208288) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(0.76333299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3182227) q[0];
sx q[0];
rz(-0.1544054) q[0];
sx q[0];
rz(-2.7247317) q[0];
rz(0.93337026) q[1];
sx q[1];
rz(-2.2552172) q[1];
sx q[1];
rz(-2.0534168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4517072) q[0];
sx q[0];
rz(-2.5390115) q[0];
sx q[0];
rz(-2.0116647) q[0];
rz(-pi) q[1];
rz(-2.5809885) q[2];
sx q[2];
rz(-1.5365155) q[2];
sx q[2];
rz(-2.3115668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58055604) q[1];
sx q[1];
rz(-1.8449748) q[1];
sx q[1];
rz(-0.31117532) q[1];
rz(-pi) q[2];
rz(-2.9253694) q[3];
sx q[3];
rz(-0.29400405) q[3];
sx q[3];
rz(0.8651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0205959) q[2];
sx q[2];
rz(-2.4211297) q[2];
sx q[2];
rz(-0.32568112) q[2];
rz(2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(0.71780786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242103) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(3.000946) q[0];
rz(-2.7463101) q[1];
sx q[1];
rz(-1.5970767) q[1];
sx q[1];
rz(-2.7093754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89807876) q[0];
sx q[0];
rz(-1.167341) q[0];
sx q[0];
rz(-2.0824758) q[0];
x q[1];
rz(3.0555058) q[2];
sx q[2];
rz(-1.0660835) q[2];
sx q[2];
rz(1.1881669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.144369) q[1];
sx q[1];
rz(-2.8565116) q[1];
sx q[1];
rz(-2.8530733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81813502) q[3];
sx q[3];
rz(-0.77518089) q[3];
sx q[3];
rz(1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29493368) q[2];
sx q[2];
rz(-1.1120956) q[2];
sx q[2];
rz(2.6726932) q[2];
rz(2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(1.2676574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0676607) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(-0.78039783) q[0];
rz(-2.228915) q[1];
sx q[1];
rz(-2.3527805) q[1];
sx q[1];
rz(1.3891634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71043832) q[0];
sx q[0];
rz(-1.9214993) q[0];
sx q[0];
rz(1.2855269) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58163403) q[2];
sx q[2];
rz(-1.4624782) q[2];
sx q[2];
rz(1.1896287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4215736) q[1];
sx q[1];
rz(-1.4602082) q[1];
sx q[1];
rz(-2.5441267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.463988) q[3];
sx q[3];
rz(-2.2011122) q[3];
sx q[3];
rz(-2.2614071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7954365) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(-0.1795086) q[2];
rz(1.6620212) q[3];
sx q[3];
rz(-2.447465) q[3];
sx q[3];
rz(-0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96106225) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.4638715) q[0];
rz(3.1278817) q[1];
sx q[1];
rz(-1.4669712) q[1];
sx q[1];
rz(-2.1965068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4714519) q[0];
sx q[0];
rz(-1.9197122) q[0];
sx q[0];
rz(-0.80100153) q[0];
x q[1];
rz(-1.7795947) q[2];
sx q[2];
rz(-1.8755302) q[2];
sx q[2];
rz(-1.7727675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42691222) q[1];
sx q[1];
rz(-1.987793) q[1];
sx q[1];
rz(2.6257573) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2919214) q[3];
sx q[3];
rz(-0.73007089) q[3];
sx q[3];
rz(1.7770313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1549687) q[2];
sx q[2];
rz(-1.0633609) q[2];
sx q[2];
rz(-1.0490136) q[2];
rz(1.5341885) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(-1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5904215) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(-2.3792939) q[0];
rz(-2.6679692) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(-0.90829888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0001091) q[0];
sx q[0];
rz(-2.3628652) q[0];
sx q[0];
rz(-1.5634006) q[0];
rz(-0.17356603) q[2];
sx q[2];
rz(-1.3527292) q[2];
sx q[2];
rz(-1.9517) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.087203793) q[1];
sx q[1];
rz(-1.3339195) q[1];
sx q[1];
rz(-0.014976784) q[1];
rz(-1.5987464) q[3];
sx q[3];
rz(-1.5833491) q[3];
sx q[3];
rz(0.2994699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.121754) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(-0.961595) q[2];
rz(2.614295) q[3];
sx q[3];
rz(-1.3798102) q[3];
sx q[3];
rz(0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2335662) q[0];
sx q[0];
rz(-0.59995025) q[0];
sx q[0];
rz(0.34616923) q[0];
rz(-1.8442122) q[1];
sx q[1];
rz(-2.4751414) q[1];
sx q[1];
rz(0.60874879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93701836) q[0];
sx q[0];
rz(-1.4868344) q[0];
sx q[0];
rz(0.89421009) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4791711) q[2];
sx q[2];
rz(-2.8067744) q[2];
sx q[2];
rz(1.3864558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94683719) q[1];
sx q[1];
rz(-1.5018592) q[1];
sx q[1];
rz(-2.8762455) q[1];
x q[2];
rz(-0.56383384) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(0.13217029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.346647) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(2.4018304) q[3];
sx q[3];
rz(-1.0659822) q[3];
sx q[3];
rz(0.8968001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5935434) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(-1.8400486) q[1];
sx q[1];
rz(-1.5751244) q[1];
sx q[1];
rz(-2.8155933) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4833925) q[0];
sx q[0];
rz(-2.3478386) q[0];
sx q[0];
rz(-1.159159) q[0];
x q[1];
rz(0.65970274) q[2];
sx q[2];
rz(-0.40605011) q[2];
sx q[2];
rz(-2.7999807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1797267) q[1];
sx q[1];
rz(-0.71007198) q[1];
sx q[1];
rz(-1.9741139) q[1];
rz(1.4526618) q[3];
sx q[3];
rz(-1.6836327) q[3];
sx q[3];
rz(-1.090534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43805435) q[2];
sx q[2];
rz(-1.2039098) q[2];
sx q[2];
rz(-0.34029141) q[2];
rz(-1.641364) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-2.5435508) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21094766) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(3.0979544) q[0];
rz(-0.71815193) q[1];
sx q[1];
rz(-1.8225881) q[1];
sx q[1];
rz(1.4290379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008416) q[0];
sx q[0];
rz(-1.740983) q[0];
sx q[0];
rz(-2.9508181) q[0];
rz(-pi) q[1];
rz(1.7068091) q[2];
sx q[2];
rz(-2.2439661) q[2];
sx q[2];
rz(-1.4006795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2344115) q[1];
sx q[1];
rz(-0.49691191) q[1];
sx q[1];
rz(0.18053825) q[1];
x q[2];
rz(-1.9472856) q[3];
sx q[3];
rz(-1.912279) q[3];
sx q[3];
rz(-1.7870513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4622197) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(2.098162) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(-1.2354318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048024561) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(-1.7780766) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(0.72259283) q[2];
sx q[2];
rz(-1.9836139) q[2];
sx q[2];
rz(0.2000533) q[2];
rz(-1.0092953) q[3];
sx q[3];
rz(-2.3875368) q[3];
sx q[3];
rz(-0.33215678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
