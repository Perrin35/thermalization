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
rz(1.963653) q[0];
sx q[0];
rz(-0.093129245) q[0];
sx q[0];
rz(-2.4094474) q[0];
rz(1.572345) q[1];
sx q[1];
rz(-2.2557926) q[1];
sx q[1];
rz(-0.43876171) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2193091) q[0];
sx q[0];
rz(-1.5153335) q[0];
sx q[0];
rz(1.903141) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83005367) q[2];
sx q[2];
rz(-2.7657653) q[2];
sx q[2];
rz(2.2245882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4784412) q[1];
sx q[1];
rz(-1.8290797) q[1];
sx q[1];
rz(-1.5228935) q[1];
rz(-1.6192299) q[3];
sx q[3];
rz(-2.0287073) q[3];
sx q[3];
rz(1.4428136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0438805) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(-0.91862339) q[2];
rz(-1.0527323) q[3];
sx q[3];
rz(-0.91679263) q[3];
sx q[3];
rz(-0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87616462) q[0];
sx q[0];
rz(-2.3227203) q[0];
sx q[0];
rz(-2.2692666) q[0];
rz(-1.0126975) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-2.8025119) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.942303) q[0];
sx q[0];
rz(-0.45594076) q[0];
sx q[0];
rz(0.29642563) q[0];
rz(-pi) q[1];
rz(-1.048597) q[2];
sx q[2];
rz(-1.9777918) q[2];
sx q[2];
rz(-2.5841449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2660632) q[1];
sx q[1];
rz(-1.5782981) q[1];
sx q[1];
rz(3.0999567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9267155) q[3];
sx q[3];
rz(-0.3444852) q[3];
sx q[3];
rz(3.0657299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8104711) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(2.5362711) q[2];
rz(-0.33453861) q[3];
sx q[3];
rz(-0.3722705) q[3];
sx q[3];
rz(-0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89638585) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(-1.4728004) q[0];
rz(2.8145166) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(2.8715141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6657406) q[0];
sx q[0];
rz(-2.088463) q[0];
sx q[0];
rz(-1.3435908) q[0];
rz(1.142092) q[2];
sx q[2];
rz(-1.639591) q[2];
sx q[2];
rz(2.4479933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0097937219) q[1];
sx q[1];
rz(-2.0094618) q[1];
sx q[1];
rz(-0.75097221) q[1];
x q[2];
rz(-3.0614254) q[3];
sx q[3];
rz(-1.5669549) q[3];
sx q[3];
rz(0.14417917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0373056) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(1.2874862) q[2];
rz(1.9262975) q[3];
sx q[3];
rz(-1.7561965) q[3];
sx q[3];
rz(-0.17770411) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6215068) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(-0.31975123) q[0];
rz(0.41140914) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(0.080332669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7112996) q[0];
sx q[0];
rz(-1.7691433) q[0];
sx q[0];
rz(1.2567149) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9991987) q[2];
sx q[2];
rz(-2.3300397) q[2];
sx q[2];
rz(2.591557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6957533) q[1];
sx q[1];
rz(-1.8421122) q[1];
sx q[1];
rz(-0.81636565) q[1];
x q[2];
rz(1.8378958) q[3];
sx q[3];
rz(-2.1648277) q[3];
sx q[3];
rz(2.475762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2692928) q[2];
sx q[2];
rz(-0.23739693) q[2];
sx q[2];
rz(0.044895127) q[2];
rz(2.6063555) q[3];
sx q[3];
rz(-1.6343296) q[3];
sx q[3];
rz(0.11307344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1503247) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(-0.87565652) q[0];
rz(2.9278897) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(-1.5692086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.050135) q[0];
sx q[0];
rz(-0.3803072) q[0];
sx q[0];
rz(2.1652392) q[0];
rz(-pi) q[1];
rz(-0.43691152) q[2];
sx q[2];
rz(-1.3624117) q[2];
sx q[2];
rz(0.85308272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45156577) q[1];
sx q[1];
rz(-2.4358349) q[1];
sx q[1];
rz(0.16132055) q[1];
rz(0.27690378) q[3];
sx q[3];
rz(-0.33280643) q[3];
sx q[3];
rz(1.8148418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.956942) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(-1.0682028) q[2];
rz(-2.5351561) q[3];
sx q[3];
rz(-1.9304099) q[3];
sx q[3];
rz(-3.1312805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(1.7506208) q[0];
rz(0.5237611) q[1];
sx q[1];
rz(-2.5660089) q[1];
sx q[1];
rz(-1.1837122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42126781) q[0];
sx q[0];
rz(-0.29365942) q[0];
sx q[0];
rz(1.378242) q[0];
x q[1];
rz(-0.95518388) q[2];
sx q[2];
rz(-0.092530017) q[2];
sx q[2];
rz(-1.1902155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27734807) q[1];
sx q[1];
rz(-1.7226761) q[1];
sx q[1];
rz(0.95452301) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11012245) q[3];
sx q[3];
rz(-1.7043132) q[3];
sx q[3];
rz(0.27329521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70048731) q[2];
sx q[2];
rz(-0.59868559) q[2];
sx q[2];
rz(1.0318476) q[2];
rz(1.4619689) q[3];
sx q[3];
rz(-2.4470191) q[3];
sx q[3];
rz(1.7803378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4782891) q[0];
sx q[0];
rz(-2.7155868) q[0];
sx q[0];
rz(-1.7167094) q[0];
rz(-0.58319631) q[1];
sx q[1];
rz(-1.9164663) q[1];
sx q[1];
rz(-3.084175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17744495) q[0];
sx q[0];
rz(-2.3741407) q[0];
sx q[0];
rz(1.5758118) q[0];
x q[1];
rz(-2.4671747) q[2];
sx q[2];
rz(-2.4092374) q[2];
sx q[2];
rz(0.44841097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.098328248) q[1];
sx q[1];
rz(-0.66126138) q[1];
sx q[1];
rz(0.28820451) q[1];
rz(1.972265) q[3];
sx q[3];
rz(-0.9801995) q[3];
sx q[3];
rz(0.99311738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8251557) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(-1.3081029) q[2];
rz(-1.3206652) q[3];
sx q[3];
rz(-0.84669176) q[3];
sx q[3];
rz(0.86110419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945781) q[0];
sx q[0];
rz(-2.4569643) q[0];
sx q[0];
rz(1.8735877) q[0];
rz(-2.0199203) q[1];
sx q[1];
rz(-1.8662165) q[1];
sx q[1];
rz(0.65006382) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9509264) q[0];
sx q[0];
rz(-1.1703849) q[0];
sx q[0];
rz(2.1413442) q[0];
rz(-1.0380657) q[2];
sx q[2];
rz(-2.1830171) q[2];
sx q[2];
rz(0.74712268) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6942581) q[1];
sx q[1];
rz(-1.9617394) q[1];
sx q[1];
rz(2.5503134) q[1];
rz(-2.1459747) q[3];
sx q[3];
rz(-1.1523968) q[3];
sx q[3];
rz(-1.0632689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0012297) q[2];
sx q[2];
rz(-1.659212) q[2];
sx q[2];
rz(-0.39230997) q[2];
rz(2.9366734) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1392764) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(-1.683715) q[0];
rz(2.814759) q[1];
sx q[1];
rz(-1.146233) q[1];
sx q[1];
rz(1.0439509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4503339) q[0];
sx q[0];
rz(-0.46713167) q[0];
sx q[0];
rz(-1.1410261) q[0];
x q[1];
rz(1.4307666) q[2];
sx q[2];
rz(-1.1159889) q[2];
sx q[2];
rz(0.32250139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9731703) q[1];
sx q[1];
rz(-2.7379824) q[1];
sx q[1];
rz(1.7803867) q[1];
x q[2];
rz(0.92071988) q[3];
sx q[3];
rz(-1.8333922) q[3];
sx q[3];
rz(-0.85185862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30847183) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(-0.92588818) q[2];
rz(1.067767) q[3];
sx q[3];
rz(-2.0177149) q[3];
sx q[3];
rz(1.7759751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85497102) q[0];
sx q[0];
rz(-1.4737031) q[0];
sx q[0];
rz(0.58706748) q[0];
rz(-0.58285561) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(-2.6606182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60763393) q[0];
sx q[0];
rz(-1.1635499) q[0];
sx q[0];
rz(1.3429705) q[0];
rz(-pi) q[1];
rz(2.0430668) q[2];
sx q[2];
rz(-2.741217) q[2];
sx q[2];
rz(-3.056293) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96718699) q[1];
sx q[1];
rz(-0.71651006) q[1];
sx q[1];
rz(2.4661676) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39980657) q[3];
sx q[3];
rz(-1.7553864) q[3];
sx q[3];
rz(-1.0291444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9440072) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(2.3594989) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(-1.9935002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999577) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(1.6089454) q[1];
sx q[1];
rz(-1.1679222) q[1];
sx q[1];
rz(2.4293778) q[1];
rz(0.92171348) q[2];
sx q[2];
rz(-2.1940008) q[2];
sx q[2];
rz(-3.1354207) q[2];
rz(-0.28859517) q[3];
sx q[3];
rz(-0.51641083) q[3];
sx q[3];
rz(-0.33752059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
