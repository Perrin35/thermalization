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
rz(-0.4001652) q[0];
sx q[0];
rz(-0.53181177) q[0];
sx q[0];
rz(2.3398633) q[0];
rz(-2.5126558) q[1];
sx q[1];
rz(-1.3087654) q[1];
sx q[1];
rz(-0.3892678) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4621534) q[0];
sx q[0];
rz(-0.91374699) q[0];
sx q[0];
rz(1.916196) q[0];
x q[1];
rz(-2.8966337) q[2];
sx q[2];
rz(-0.18671045) q[2];
sx q[2];
rz(2.1706131) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3013346) q[1];
sx q[1];
rz(-1.0691027) q[1];
sx q[1];
rz(0.69712104) q[1];
rz(1.5262768) q[3];
sx q[3];
rz(-2.4137437) q[3];
sx q[3];
rz(-2.0522709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6486711) q[2];
sx q[2];
rz(-2.6914983) q[2];
sx q[2];
rz(-2.8802803) q[2];
rz(2.785545) q[3];
sx q[3];
rz(-2.0331523) q[3];
sx q[3];
rz(1.0935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68591958) q[0];
sx q[0];
rz(-2.5643667) q[0];
sx q[0];
rz(-2.8976231) q[0];
rz(-2.7120554) q[1];
sx q[1];
rz(-1.669084) q[1];
sx q[1];
rz(2.4006749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572621) q[0];
sx q[0];
rz(-1.9004915) q[0];
sx q[0];
rz(1.9614205) q[0];
x q[1];
rz(-1.6119611) q[2];
sx q[2];
rz(-1.9518234) q[2];
sx q[2];
rz(1.5059901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0216995) q[1];
sx q[1];
rz(-0.60696536) q[1];
sx q[1];
rz(0.40268437) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71535625) q[3];
sx q[3];
rz(-1.1872059) q[3];
sx q[3];
rz(2.5959248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7168768) q[2];
sx q[2];
rz(-1.6191142) q[2];
sx q[2];
rz(0.87872046) q[2];
rz(2.1667513) q[3];
sx q[3];
rz(-1.2474371) q[3];
sx q[3];
rz(-1.6802906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.2002624) q[0];
sx q[0];
rz(-0.55568475) q[0];
sx q[0];
rz(2.9479807) q[0];
rz(-0.10784736) q[1];
sx q[1];
rz(-0.71433181) q[1];
sx q[1];
rz(2.329619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0064751) q[0];
sx q[0];
rz(-1.1163981) q[0];
sx q[0];
rz(3.019089) q[0];
rz(-pi) q[1];
rz(0.8167972) q[2];
sx q[2];
rz(-1.9765761) q[2];
sx q[2];
rz(-2.7614373) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0206611) q[1];
sx q[1];
rz(-1.963208) q[1];
sx q[1];
rz(0.25041469) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8609686) q[3];
sx q[3];
rz(-0.45553744) q[3];
sx q[3];
rz(-2.303567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0632443) q[2];
sx q[2];
rz(-1.2816659) q[2];
sx q[2];
rz(-0.0035302103) q[2];
rz(-2.6053536) q[3];
sx q[3];
rz(-0.25501525) q[3];
sx q[3];
rz(-2.6030538) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2980767) q[0];
sx q[0];
rz(-1.902782) q[0];
sx q[0];
rz(0.82999825) q[0];
rz(1.7371381) q[1];
sx q[1];
rz(-0.91157118) q[1];
sx q[1];
rz(-2.8703168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7260482) q[0];
sx q[0];
rz(-1.7689287) q[0];
sx q[0];
rz(1.8285486) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4936594) q[2];
sx q[2];
rz(-2.234314) q[2];
sx q[2];
rz(-1.5761988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3437561) q[1];
sx q[1];
rz(-2.7833354) q[1];
sx q[1];
rz(-0.078726493) q[1];
rz(-0.15437281) q[3];
sx q[3];
rz(-0.8027146) q[3];
sx q[3];
rz(1.5824535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5156877) q[2];
sx q[2];
rz(-0.38334623) q[2];
sx q[2];
rz(2.3670727) q[2];
rz(-1.0033311) q[3];
sx q[3];
rz(-2.7480875) q[3];
sx q[3];
rz(0.5459319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14458732) q[0];
sx q[0];
rz(-1.0969176) q[0];
sx q[0];
rz(-1.0928094) q[0];
rz(1.3358759) q[1];
sx q[1];
rz(-2.3922258) q[1];
sx q[1];
rz(0.54378477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519125) q[0];
sx q[0];
rz(-2.8767715) q[0];
sx q[0];
rz(1.7811243) q[0];
x q[1];
rz(0.67853482) q[2];
sx q[2];
rz(-1.4894052) q[2];
sx q[2];
rz(2.396109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.136506) q[1];
sx q[1];
rz(-0.6906116) q[1];
sx q[1];
rz(0.68271272) q[1];
rz(-pi) q[2];
rz(1.5724935) q[3];
sx q[3];
rz(-0.6203531) q[3];
sx q[3];
rz(-1.9242632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.233923) q[2];
sx q[2];
rz(-0.69207865) q[2];
sx q[2];
rz(-0.27883369) q[2];
rz(2.8558266) q[3];
sx q[3];
rz(-2.2279492) q[3];
sx q[3];
rz(0.41551503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.435442) q[0];
sx q[0];
rz(-0.64318648) q[0];
sx q[0];
rz(1.0253133) q[0];
rz(-3.0030491) q[1];
sx q[1];
rz(-1.5516611) q[1];
sx q[1];
rz(-2.9655546) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4885725) q[0];
sx q[0];
rz(-2.5551642) q[0];
sx q[0];
rz(-0.47393786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7654714) q[2];
sx q[2];
rz(-1.2174213) q[2];
sx q[2];
rz(0.88549685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83505364) q[1];
sx q[1];
rz(-1.7944778) q[1];
sx q[1];
rz(2.064631) q[1];
rz(-pi) q[2];
rz(-1.2071836) q[3];
sx q[3];
rz(-2.6904958) q[3];
sx q[3];
rz(2.0157984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90870086) q[2];
sx q[2];
rz(-0.42753926) q[2];
sx q[2];
rz(-2.3517081) q[2];
rz(0.32957736) q[3];
sx q[3];
rz(-1.3720082) q[3];
sx q[3];
rz(0.65143877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.71094197) q[0];
sx q[0];
rz(-0.21345226) q[0];
sx q[0];
rz(-0.99367225) q[0];
rz(-1.0020533) q[1];
sx q[1];
rz(-0.99169815) q[1];
sx q[1];
rz(1.6755982) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3661909) q[0];
sx q[0];
rz(-2.4537773) q[0];
sx q[0];
rz(0.28190502) q[0];
rz(-pi) q[1];
rz(1.3411677) q[2];
sx q[2];
rz(-2.672019) q[2];
sx q[2];
rz(-1.1498677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0198145) q[1];
sx q[1];
rz(-1.7516859) q[1];
sx q[1];
rz(-1.4223076) q[1];
rz(-0.83169072) q[3];
sx q[3];
rz(-1.2281831) q[3];
sx q[3];
rz(2.6765842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34622908) q[2];
sx q[2];
rz(-2.2987404) q[2];
sx q[2];
rz(2.8818434) q[2];
rz(2.3527457) q[3];
sx q[3];
rz(-2.2987821) q[3];
sx q[3];
rz(-1.3938676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9144834) q[0];
sx q[0];
rz(-1.8552584) q[0];
sx q[0];
rz(0.9911384) q[0];
rz(1.8046851) q[1];
sx q[1];
rz(-2.3146345) q[1];
sx q[1];
rz(-1.4088668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1240272) q[0];
sx q[0];
rz(-2.6691737) q[0];
sx q[0];
rz(2.8851465) q[0];
x q[1];
rz(-0.089167909) q[2];
sx q[2];
rz(-2.1074416) q[2];
sx q[2];
rz(0.23039699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0947021) q[1];
sx q[1];
rz(-0.45744236) q[1];
sx q[1];
rz(0.66384683) q[1];
x q[2];
rz(2.9314133) q[3];
sx q[3];
rz(-2.4553799) q[3];
sx q[3];
rz(0.33542882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1651429) q[2];
sx q[2];
rz(-1.4850478) q[2];
sx q[2];
rz(2.5531947) q[2];
rz(2.8277561) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(-0.61217827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2345851) q[0];
sx q[0];
rz(-1.122965) q[0];
sx q[0];
rz(3.0871952) q[0];
rz(0.76039487) q[1];
sx q[1];
rz(-0.94579983) q[1];
sx q[1];
rz(-2.0423582) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1791819) q[0];
sx q[0];
rz(-2.1405309) q[0];
sx q[0];
rz(0.197844) q[0];
rz(-pi) q[1];
rz(0.80087687) q[2];
sx q[2];
rz(-2.184424) q[2];
sx q[2];
rz(-0.18979812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.526872) q[1];
sx q[1];
rz(-0.81382591) q[1];
sx q[1];
rz(0.24095638) q[1];
rz(-pi) q[2];
rz(-1.4599007) q[3];
sx q[3];
rz(-1.520021) q[3];
sx q[3];
rz(-2.5230356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2796563) q[2];
sx q[2];
rz(-1.6013689) q[2];
sx q[2];
rz(0.81735617) q[2];
rz(-3.127408) q[3];
sx q[3];
rz(-2.7458906) q[3];
sx q[3];
rz(-2.0524249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8297183) q[0];
sx q[0];
rz(-1.051396) q[0];
sx q[0];
rz(0.22148393) q[0];
rz(2.7102176) q[1];
sx q[1];
rz(-2.2384877) q[1];
sx q[1];
rz(1.0198786) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1843625) q[0];
sx q[0];
rz(-1.5505322) q[0];
sx q[0];
rz(-0.046837687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20163147) q[2];
sx q[2];
rz(-1.4678737) q[2];
sx q[2];
rz(2.0992744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74032593) q[1];
sx q[1];
rz(-1.6049084) q[1];
sx q[1];
rz(-1.9032065) q[1];
rz(-pi) q[2];
rz(-0.13615578) q[3];
sx q[3];
rz(-1.7014967) q[3];
sx q[3];
rz(-3.1112373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7027616) q[2];
sx q[2];
rz(-1.2063382) q[2];
sx q[2];
rz(1.8863401) q[2];
rz(-3.1158279) q[3];
sx q[3];
rz(-1.8119101) q[3];
sx q[3];
rz(0.51938957) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41853607) q[0];
sx q[0];
rz(-2.3522455) q[0];
sx q[0];
rz(1.1242207) q[0];
rz(-0.33456805) q[1];
sx q[1];
rz(-2.6523013) q[1];
sx q[1];
rz(2.1025067) q[1];
rz(0.46229565) q[2];
sx q[2];
rz(-1.8097327) q[2];
sx q[2];
rz(2.3859522) q[2];
rz(0.48503899) q[3];
sx q[3];
rz(-2.1342604) q[3];
sx q[3];
rz(0.66935183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
