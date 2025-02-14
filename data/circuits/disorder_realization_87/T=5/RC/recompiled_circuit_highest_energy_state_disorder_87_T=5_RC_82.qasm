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
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(0.33831212) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7871246) q[0];
sx q[0];
rz(-1.5347052) q[0];
sx q[0];
rz(0.24255328) q[0];
x q[1];
rz(1.9954483) q[2];
sx q[2];
rz(-0.82595794) q[2];
sx q[2];
rz(-2.7934157) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9788614) q[1];
sx q[1];
rz(-2.581209) q[1];
sx q[1];
rz(2.8864229) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1756535) q[3];
sx q[3];
rz(-2.3199816) q[3];
sx q[3];
rz(-0.3607817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3931291) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(-0.65790042) q[2];
rz(-2.6864478) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748134) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(0.28208062) q[0];
rz(-2.8981949) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(0.74877053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42001777) q[0];
sx q[0];
rz(-0.87792464) q[0];
sx q[0];
rz(-1.3189032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6924627) q[2];
sx q[2];
rz(-1.5632544) q[2];
sx q[2];
rz(0.53754025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48909125) q[1];
sx q[1];
rz(-1.6630807) q[1];
sx q[1];
rz(-1.0670061) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5010819) q[3];
sx q[3];
rz(-1.3718714) q[3];
sx q[3];
rz(-0.45123842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28027174) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(1.8005499) q[2];
rz(-1.9848112) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(0.7635428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007616) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(0.65735835) q[0];
rz(1.1454469) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(-0.81833902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3439804) q[0];
sx q[0];
rz(-1.3320141) q[0];
sx q[0];
rz(-0.35399951) q[0];
rz(-pi) q[1];
rz(3.0938593) q[2];
sx q[2];
rz(-1.2558508) q[2];
sx q[2];
rz(-2.3746109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1166621) q[1];
sx q[1];
rz(-2.0567472) q[1];
sx q[1];
rz(0.74851739) q[1];
rz(-pi) q[2];
rz(-1.7210427) q[3];
sx q[3];
rz(-0.45805061) q[3];
sx q[3];
rz(0.86298215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98693097) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(1.7884802) q[2];
rz(-2.1766369) q[3];
sx q[3];
rz(-1.302364) q[3];
sx q[3];
rz(-0.42199782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.2009784) q[0];
rz(-0.0032084223) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(2.9046955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99588441) q[0];
sx q[0];
rz(-2.0170324) q[0];
sx q[0];
rz(-0.97490411) q[0];
rz(0.17572524) q[2];
sx q[2];
rz(-2.496886) q[2];
sx q[2];
rz(2.3852661) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.264852) q[1];
sx q[1];
rz(-0.63794604) q[1];
sx q[1];
rz(-1.9015584) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3637216) q[3];
sx q[3];
rz(-0.98613769) q[3];
sx q[3];
rz(-0.28077048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34430382) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(-0.34991831) q[0];
rz(1.9328851) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(3.1399609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9043372) q[0];
sx q[0];
rz(-1.4044824) q[0];
sx q[0];
rz(0.24760274) q[0];
x q[1];
rz(2.8441378) q[2];
sx q[2];
rz(-1.607873) q[2];
sx q[2];
rz(2.2789479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0495245) q[1];
sx q[1];
rz(-0.18316575) q[1];
sx q[1];
rz(1.3080255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2285813) q[3];
sx q[3];
rz(-2.8051441) q[3];
sx q[3];
rz(-0.89321619) q[3];
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
rz(2.7023884) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(-2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8301903) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(0.79469529) q[0];
rz(-1.7757802) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(2.6672003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2524388) q[0];
sx q[0];
rz(-1.1934917) q[0];
sx q[0];
rz(1.0214367) q[0];
rz(-pi) q[1];
rz(0.55616711) q[2];
sx q[2];
rz(-0.63236299) q[2];
sx q[2];
rz(1.4390505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3620356) q[1];
sx q[1];
rz(-1.0048702) q[1];
sx q[1];
rz(2.0574244) q[1];
rz(-pi) q[2];
rz(-1.0211518) q[3];
sx q[3];
rz(-2.4826038) q[3];
sx q[3];
rz(-0.07168183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61364335) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(-2.6141686) q[2];
rz(0.6066277) q[3];
sx q[3];
rz(-2.9546723) q[3];
sx q[3];
rz(-1.0724148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793176) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(-2.344017) q[0];
rz(0.89353117) q[1];
sx q[1];
rz(-0.83502665) q[1];
sx q[1];
rz(-2.8614047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55927568) q[0];
sx q[0];
rz(-1.3963376) q[0];
sx q[0];
rz(-1.1691537) q[0];
x q[1];
rz(-2.3827219) q[2];
sx q[2];
rz(-2.3312097) q[2];
sx q[2];
rz(0.23960613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57443212) q[1];
sx q[1];
rz(-0.48266294) q[1];
sx q[1];
rz(1.5978651) q[1];
x q[2];
rz(-2.7817621) q[3];
sx q[3];
rz(-1.0092648) q[3];
sx q[3];
rz(0.4175182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4013275) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(-1.2124445) q[2];
rz(2.9017743) q[3];
sx q[3];
rz(-0.66767728) q[3];
sx q[3];
rz(0.56813204) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78385335) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(-1.3215815) q[0];
rz(2.9603738) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(0.063057335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6609333) q[0];
sx q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(2.3602135) q[0];
x q[1];
rz(-3.0538043) q[2];
sx q[2];
rz(-1.8551747) q[2];
sx q[2];
rz(1.926601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6136421) q[1];
sx q[1];
rz(-2.7220352) q[1];
sx q[1];
rz(-0.63906868) q[1];
rz(-pi) q[2];
rz(-3.1386915) q[3];
sx q[3];
rz(-1.7333247) q[3];
sx q[3];
rz(-1.8123466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(-1.6723527) q[2];
rz(-1.7478878) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(-0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6032228) q[0];
sx q[0];
rz(-0.61036888) q[0];
sx q[0];
rz(0.99697733) q[0];
rz(-0.793055) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(-2.0319895) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.1446321) q[2];
sx q[2];
rz(-0.89489102) q[2];
sx q[2];
rz(-0.083783178) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.931568) q[1];
sx q[1];
rz(-0.76142516) q[1];
sx q[1];
rz(-0.85827338) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7392735) q[3];
sx q[3];
rz(-2.2003897) q[3];
sx q[3];
rz(-1.633267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41254607) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(1.868978) q[2];
rz(-1.1254958) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(-0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7703055) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(3.0433997) q[0];
rz(1.498361) q[1];
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
x q[1];
rz(-1.5151843) q[2];
sx q[2];
rz(-1.8225553) q[2];
sx q[2];
rz(2.1625569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0467207) q[1];
sx q[1];
rz(-2.0025684) q[1];
sx q[1];
rz(2.8304965) q[1];
x q[2];
rz(1.7915947) q[3];
sx q[3];
rz(-1.9933369) q[3];
sx q[3];
rz(0.72964668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7958293) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(2.8161827) q[2];
rz(-2.365153) q[3];
sx q[3];
rz(-2.1263945) q[3];
sx q[3];
rz(0.6071035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14071874) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(2.8554032) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(-2.7306225) q[2];
sx q[2];
rz(-2.3373418) q[2];
sx q[2];
rz(1.892754) q[2];
rz(2.552916) q[3];
sx q[3];
rz(-1.8559783) q[3];
sx q[3];
rz(1.2067309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
