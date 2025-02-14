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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(0.80817428) q[0];
rz(1.9975245) q[1];
sx q[1];
rz(3.9132325) q[1];
sx q[1];
rz(11.003718) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842386) q[0];
sx q[0];
rz(-2.1942217) q[0];
sx q[0];
rz(2.4592072) q[0];
rz(1.8380661) q[2];
sx q[2];
rz(-1.7198945) q[2];
sx q[2];
rz(-1.2318512) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6526323) q[1];
sx q[1];
rz(-0.51981407) q[1];
sx q[1];
rz(1.7955154) q[1];
rz(-pi) q[2];
rz(-0.64227958) q[3];
sx q[3];
rz(-0.056996973) q[3];
sx q[3];
rz(1.9798724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19848862) q[2];
sx q[2];
rz(-0.21185943) q[2];
sx q[2];
rz(-1.0342106) q[2];
rz(-2.4487623) q[3];
sx q[3];
rz(-1.0537909) q[3];
sx q[3];
rz(1.0606162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3028054) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(2.5651108) q[0];
rz(-0.020596404) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(-2.1131262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3555455) q[0];
sx q[0];
rz(-1.78616) q[0];
sx q[0];
rz(2.8785365) q[0];
rz(-pi) q[1];
rz(0.307085) q[2];
sx q[2];
rz(-1.2132309) q[2];
sx q[2];
rz(-1.7890695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2814863) q[1];
sx q[1];
rz(-0.53412837) q[1];
sx q[1];
rz(-2.2661282) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9628434) q[3];
sx q[3];
rz(-2.1806751) q[3];
sx q[3];
rz(-1.6926241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2241406) q[2];
sx q[2];
rz(-0.96106207) q[2];
sx q[2];
rz(-1.2646328) q[2];
rz(0.97359109) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(2.8784331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6384386) q[0];
sx q[0];
rz(-1.8373024) q[0];
sx q[0];
rz(0.85357443) q[0];
rz(2.0897934) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(0.015017088) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6081469) q[0];
sx q[0];
rz(-2.2813873) q[0];
sx q[0];
rz(-0.011562499) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9237346) q[2];
sx q[2];
rz(-0.098420489) q[2];
sx q[2];
rz(0.86505666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7370809) q[1];
sx q[1];
rz(-0.85857262) q[1];
sx q[1];
rz(-1.4732857) q[1];
rz(0.28983966) q[3];
sx q[3];
rz(-1.5095169) q[3];
sx q[3];
rz(-1.0672399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1198279) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(2.8601698) q[2];
rz(1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948831) q[0];
sx q[0];
rz(-1.895772) q[0];
sx q[0];
rz(-2.967714) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.4215819) q[1];
sx q[1];
rz(1.7132267) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9857432) q[0];
sx q[0];
rz(-1.5850889) q[0];
sx q[0];
rz(1.3784268) q[0];
rz(-pi) q[1];
rz(-2.7816781) q[2];
sx q[2];
rz(-1.6891306) q[2];
sx q[2];
rz(1.182076) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.169379) q[1];
sx q[1];
rz(-1.5323109) q[1];
sx q[1];
rz(2.4840218) q[1];
rz(-0.36767204) q[3];
sx q[3];
rz(-1.9353364) q[3];
sx q[3];
rz(0.035797771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4398769) q[2];
sx q[2];
rz(-0.82650799) q[2];
sx q[2];
rz(-2.7244549) q[2];
rz(0.16962984) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(-2.4197742) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7222662) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(2.3720692) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(1.6228898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8474591) q[0];
sx q[0];
rz(-1.8222229) q[0];
sx q[0];
rz(-2.2537116) q[0];
rz(2.2686636) q[2];
sx q[2];
rz(-2.180122) q[2];
sx q[2];
rz(0.84963679) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5656083) q[1];
sx q[1];
rz(-0.19153015) q[1];
sx q[1];
rz(1.4148786) q[1];
rz(1.4368016) q[3];
sx q[3];
rz(-1.6686642) q[3];
sx q[3];
rz(1.1214453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7421444) q[2];
sx q[2];
rz(-2.1592906) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(-2.4230912) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(-2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28498483) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(-1.5307776) q[0];
rz(-1.4467622) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(1.4379427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3510144) q[0];
sx q[0];
rz(-1.5503128) q[0];
sx q[0];
rz(1.162942) q[0];
x q[1];
rz(0.5912663) q[2];
sx q[2];
rz(-1.8857393) q[2];
sx q[2];
rz(0.43063088) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9077545) q[1];
sx q[1];
rz(-1.4871039) q[1];
sx q[1];
rz(0.5096883) q[1];
rz(1.5678649) q[3];
sx q[3];
rz(-0.81122196) q[3];
sx q[3];
rz(2.7377759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1630254) q[2];
sx q[2];
rz(-1.7820396) q[2];
sx q[2];
rz(1.2320409) q[2];
rz(1.1466522) q[3];
sx q[3];
rz(-1.3403284) q[3];
sx q[3];
rz(0.31057096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6335886) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(1.1401796) q[0];
rz(-0.89705244) q[1];
sx q[1];
rz(-1.3373809) q[1];
sx q[1];
rz(0.030489771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59647467) q[0];
sx q[0];
rz(-1.4545355) q[0];
sx q[0];
rz(1.3482987) q[0];
rz(-2.4989481) q[2];
sx q[2];
rz(-1.8360999) q[2];
sx q[2];
rz(0.20390262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8463609) q[1];
sx q[1];
rz(-1.8627916) q[1];
sx q[1];
rz(0.83004119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.072149) q[3];
sx q[3];
rz(-2.3460918) q[3];
sx q[3];
rz(-1.3971869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19411479) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(1.8428295) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(1.8208767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3461935) q[0];
sx q[0];
rz(-0.29618707) q[0];
sx q[0];
rz(-1.7561703) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(-0.48929712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4327571) q[0];
sx q[0];
rz(-2.2422195) q[0];
sx q[0];
rz(1.5455724) q[0];
rz(-pi) q[1];
rz(-1.2345612) q[2];
sx q[2];
rz(-1.5362985) q[2];
sx q[2];
rz(-0.86755841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.341574) q[1];
sx q[1];
rz(-1.4765146) q[1];
sx q[1];
rz(-1.0977919) q[1];
rz(1.321526) q[3];
sx q[3];
rz(-1.6695654) q[3];
sx q[3];
rz(2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74799246) q[2];
sx q[2];
rz(-1.1425428) q[2];
sx q[2];
rz(0.46621123) q[2];
rz(-1.5318058) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(2.8209749) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4894678) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(3.0925282) q[0];
rz(-0.24066726) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(3.1191471) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4197893) q[0];
sx q[0];
rz(-1.5947486) q[0];
sx q[0];
rz(2.1502058) q[0];
x q[1];
rz(-2.4149553) q[2];
sx q[2];
rz(-0.49984806) q[2];
sx q[2];
rz(-2.8195087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3901523) q[1];
sx q[1];
rz(-0.94822394) q[1];
sx q[1];
rz(0.97761131) q[1];
rz(2.6050354) q[3];
sx q[3];
rz(-0.76013764) q[3];
sx q[3];
rz(0.11790568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0040794) q[2];
sx q[2];
rz(-1.9065607) q[2];
sx q[2];
rz(0.9683041) q[2];
rz(0.83640313) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4150998) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(-2.7628164) q[0];
rz(-0.0060161034) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(-2.5078497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597719) q[0];
sx q[0];
rz(-2.2359747) q[0];
sx q[0];
rz(-2.0048672) q[0];
x q[1];
rz(-1.815755) q[2];
sx q[2];
rz(-0.58734054) q[2];
sx q[2];
rz(-2.7224685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5744575) q[1];
sx q[1];
rz(-2.0525041) q[1];
sx q[1];
rz(2.5888799) q[1];
x q[2];
rz(-0.59468693) q[3];
sx q[3];
rz(-0.98238889) q[3];
sx q[3];
rz(-2.7168517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4358431) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(0.61326927) q[2];
rz(0.67522007) q[3];
sx q[3];
rz(-0.37874159) q[3];
sx q[3];
rz(-2.3238382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262065) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(0.040493852) q[2];
sx q[2];
rz(-2.6815985) q[2];
sx q[2];
rz(0.56618377) q[2];
rz(-1.6982895) q[3];
sx q[3];
rz(-2.5155407) q[3];
sx q[3];
rz(-1.2609353) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
