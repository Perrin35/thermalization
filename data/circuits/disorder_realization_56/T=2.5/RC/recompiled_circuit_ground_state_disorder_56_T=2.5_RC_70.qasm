OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(-0.30012497) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(4.6680968) q[1];
sx q[1];
rz(11.065281) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34632698) q[0];
sx q[0];
rz(-3.1259968) q[0];
sx q[0];
rz(2.8898284) q[0];
x q[1];
rz(1.1741224) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(-1.3371181) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.056880205) q[1];
sx q[1];
rz(-0.69082971) q[1];
sx q[1];
rz(2.7954196) q[1];
rz(-0.0048801076) q[3];
sx q[3];
rz(-0.80880755) q[3];
sx q[3];
rz(-2.463256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(-0.70297757) q[2];
rz(2.5718555) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(-1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1341781) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(1.9084357) q[0];
rz(0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(2.0436683) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63609517) q[0];
sx q[0];
rz(-1.5042802) q[0];
sx q[0];
rz(1.7719222) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64867257) q[2];
sx q[2];
rz(-2.5498516) q[2];
sx q[2];
rz(-1.6206738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.532053) q[1];
sx q[1];
rz(-1.5960448) q[1];
sx q[1];
rz(-1.4933426) q[1];
x q[2];
rz(-1.6637832) q[3];
sx q[3];
rz(-2.3602897) q[3];
sx q[3];
rz(2.6546991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42830959) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(1.606288) q[2];
rz(0.71896583) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83800256) q[0];
sx q[0];
rz(-2.328023) q[0];
sx q[0];
rz(-0.59233061) q[0];
rz(-1.8968286) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2687362) q[0];
sx q[0];
rz(-1.4352192) q[0];
sx q[0];
rz(2.6536634) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2592741) q[2];
sx q[2];
rz(-1.6711418) q[2];
sx q[2];
rz(1.060263) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9551202) q[1];
sx q[1];
rz(-1.8545517) q[1];
sx q[1];
rz(-2.1452745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0592732) q[3];
sx q[3];
rz(-1.8104189) q[3];
sx q[3];
rz(-2.3828196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(-1.2325475) q[2];
rz(2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(0.48503748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1211014) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(-0.030601587) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(-1.0928924) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26390938) q[0];
sx q[0];
rz(-1.0315686) q[0];
sx q[0];
rz(2.3506021) q[0];
rz(-0.061528607) q[2];
sx q[2];
rz(-2.3752651) q[2];
sx q[2];
rz(1.4166703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92387256) q[1];
sx q[1];
rz(-1.3663379) q[1];
sx q[1];
rz(-2.2108498) q[1];
x q[2];
rz(2.4678342) q[3];
sx q[3];
rz(-2.1762037) q[3];
sx q[3];
rz(0.76402367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0577724) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(0.27285451) q[2];
rz(1.7504292) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(0.951989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4680173) q[0];
sx q[0];
rz(-1.3382358) q[0];
sx q[0];
rz(1.8433628) q[0];
rz(-2.4507554) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(1.615049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012155449) q[0];
sx q[0];
rz(-0.94506663) q[0];
sx q[0];
rz(2.246494) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66730209) q[2];
sx q[2];
rz(-2.7030253) q[2];
sx q[2];
rz(1.2325317) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5954819) q[1];
sx q[1];
rz(-0.87015663) q[1];
sx q[1];
rz(0.033820669) q[1];
rz(2.4590309) q[3];
sx q[3];
rz(-2.6877211) q[3];
sx q[3];
rz(1.0938494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4096628) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(-0.16732495) q[2];
rz(1.7512199) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(0.65756857) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996027) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(-1.2784736) q[0];
rz(1.4356042) q[1];
sx q[1];
rz(-1.4878788) q[1];
sx q[1];
rz(-0.37567589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4091051) q[0];
sx q[0];
rz(-1.8408436) q[0];
sx q[0];
rz(-2.9069515) q[0];
x q[1];
rz(-0.85083811) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(-1.6825324) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76498079) q[1];
sx q[1];
rz(-1.3479479) q[1];
sx q[1];
rz(1.2443903) q[1];
x q[2];
rz(1.3270946) q[3];
sx q[3];
rz(-0.82984314) q[3];
sx q[3];
rz(-2.0654701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(-2.2806878) q[2];
rz(1.143645) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(-0.27826571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025573108) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(0.75468165) q[0];
rz(-0.93983752) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(1.2831203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1719753) q[0];
sx q[0];
rz(-1.2640868) q[0];
sx q[0];
rz(0.0071010751) q[0];
rz(-pi) q[1];
rz(-2.0286328) q[2];
sx q[2];
rz(-1.9001868) q[2];
sx q[2];
rz(3.0731346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.922328) q[1];
sx q[1];
rz(-1.1891051) q[1];
sx q[1];
rz(-1.2569409) q[1];
x q[2];
rz(0.8831034) q[3];
sx q[3];
rz(-1.1576443) q[3];
sx q[3];
rz(-2.2108111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54520404) q[2];
sx q[2];
rz(-1.117492) q[2];
sx q[2];
rz(1.0835353) q[2];
rz(2.3986473) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899406) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(-0.75491828) q[0];
rz(-0.49916357) q[1];
sx q[1];
rz(-0.9639591) q[1];
sx q[1];
rz(-2.4430433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4071199) q[0];
sx q[0];
rz(-1.880058) q[0];
sx q[0];
rz(0.63068181) q[0];
rz(1.4261943) q[2];
sx q[2];
rz(-1.1100262) q[2];
sx q[2];
rz(-0.20587155) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10741988) q[1];
sx q[1];
rz(-1.0386969) q[1];
sx q[1];
rz(-3.1183395) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1404834) q[3];
sx q[3];
rz(-1.9898763) q[3];
sx q[3];
rz(0.97930479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90427202) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(0.80238706) q[2];
rz(0.55195105) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(-0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9629843) q[0];
sx q[0];
rz(-1.4405788) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(1.7604609) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(-1.6345056) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7716146) q[0];
sx q[0];
rz(-2.7214337) q[0];
sx q[0];
rz(-2.7922947) q[0];
x q[1];
rz(-1.1610031) q[2];
sx q[2];
rz(-1.2772182) q[2];
sx q[2];
rz(-0.30337203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.448882) q[1];
sx q[1];
rz(-1.7160985) q[1];
sx q[1];
rz(2.7753745) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3146995) q[3];
sx q[3];
rz(-1.5103467) q[3];
sx q[3];
rz(1.353144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(-1.8692807) q[2];
rz(1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24499527) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-0.63648883) q[1];
sx q[1];
rz(0.12558118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.784362) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(0.016655075) q[0];
rz(0.56350817) q[2];
sx q[2];
rz(-2.3803408) q[2];
sx q[2];
rz(-0.13680563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6148147) q[1];
sx q[1];
rz(-1.1452951) q[1];
sx q[1];
rz(-2.6761901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93070458) q[3];
sx q[3];
rz(-0.32880201) q[3];
sx q[3];
rz(-2.9931675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-2.3385091) q[2];
rz(-0.17656365) q[3];
sx q[3];
rz(-1.3253788) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9618027) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(-1.8190307) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(-2.8925469) q[2];
sx q[2];
rz(-2.8916044) q[2];
sx q[2];
rz(1.2492346) q[2];
rz(-1.6311036) q[3];
sx q[3];
rz(-1.2338439) q[3];
sx q[3];
rz(-3.0046786) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
