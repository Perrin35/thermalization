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
rz(-0.99875206) q[0];
sx q[0];
rz(-0.771703) q[0];
sx q[0];
rz(-1.2326711) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(-0.86644679) q[1];
sx q[1];
rz(0.37860695) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0888935) q[0];
sx q[0];
rz(-1.3406154) q[0];
sx q[0];
rz(-0.038904638) q[0];
rz(-pi) q[1];
rz(-0.98388715) q[2];
sx q[2];
rz(-2.4574033) q[2];
sx q[2];
rz(3.0923051) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5224865) q[1];
sx q[1];
rz(-2.2402856) q[1];
sx q[1];
rz(-1.2944024) q[1];
rz(-pi) q[2];
rz(2.5687417) q[3];
sx q[3];
rz(-1.2490954) q[3];
sx q[3];
rz(-0.44627809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0014701) q[2];
sx q[2];
rz(-1.1727611) q[2];
sx q[2];
rz(0.44974652) q[2];
rz(-2.5887515) q[3];
sx q[3];
rz(-0.55838412) q[3];
sx q[3];
rz(0.21137485) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1028035) q[0];
sx q[0];
rz(-1.9961822) q[0];
sx q[0];
rz(-0.75622028) q[0];
rz(2.4839632) q[1];
sx q[1];
rz(-2.2537474) q[1];
sx q[1];
rz(-0.58849803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257826) q[0];
sx q[0];
rz(-1.9414158) q[0];
sx q[0];
rz(-0.70612208) q[0];
x q[1];
rz(1.334085) q[2];
sx q[2];
rz(-0.76067096) q[2];
sx q[2];
rz(1.0631281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5246702) q[1];
sx q[1];
rz(-1.2782885) q[1];
sx q[1];
rz(-3.133417) q[1];
rz(-pi) q[2];
rz(-1.9733834) q[3];
sx q[3];
rz(-2.5036771) q[3];
sx q[3];
rz(2.0762553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6077891) q[2];
sx q[2];
rz(-1.6705931) q[2];
sx q[2];
rz(3.0652453) q[2];
rz(2.7401183) q[3];
sx q[3];
rz(-2.6186826) q[3];
sx q[3];
rz(2.0126066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56979316) q[0];
sx q[0];
rz(-2.5558668) q[0];
sx q[0];
rz(-2.8049923) q[0];
rz(-1.0353237) q[1];
sx q[1];
rz(-0.88408771) q[1];
sx q[1];
rz(-0.050315637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6487728) q[0];
sx q[0];
rz(-1.1002212) q[0];
sx q[0];
rz(-0.64071606) q[0];
x q[1];
rz(-0.59479721) q[2];
sx q[2];
rz(-1.7414879) q[2];
sx q[2];
rz(-1.6015944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5118647) q[1];
sx q[1];
rz(-1.9105163) q[1];
sx q[1];
rz(1.7930255) q[1];
rz(1.2462685) q[3];
sx q[3];
rz(-0.59610329) q[3];
sx q[3];
rz(-2.3625797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3745554) q[2];
sx q[2];
rz(-2.2646246) q[2];
sx q[2];
rz(2.1602737) q[2];
rz(1.0513002) q[3];
sx q[3];
rz(-1.6853354) q[3];
sx q[3];
rz(0.013896996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0392847) q[0];
sx q[0];
rz(-1.4816875) q[0];
sx q[0];
rz(-0.977595) q[0];
rz(0.54999894) q[1];
sx q[1];
rz(-0.98098749) q[1];
sx q[1];
rz(-0.94211334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18888488) q[0];
sx q[0];
rz(-1.8994944) q[0];
sx q[0];
rz(1.230775) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6245234) q[2];
sx q[2];
rz(-1.9280528) q[2];
sx q[2];
rz(-2.3600649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4244183) q[1];
sx q[1];
rz(-2.5150329) q[1];
sx q[1];
rz(-1.8612629) q[1];
rz(-pi) q[2];
rz(2.3650424) q[3];
sx q[3];
rz(-0.7012434) q[3];
sx q[3];
rz(1.4724191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.00086870988) q[2];
sx q[2];
rz(-2.7694747) q[2];
sx q[2];
rz(-1.5127399) q[2];
rz(-0.77423972) q[3];
sx q[3];
rz(-2.1858678) q[3];
sx q[3];
rz(2.9775508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0416439) q[0];
sx q[0];
rz(-2.994717) q[0];
sx q[0];
rz(0.80648333) q[0];
rz(2.3355314) q[1];
sx q[1];
rz(-1.1763923) q[1];
sx q[1];
rz(2.3890266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4680639) q[0];
sx q[0];
rz(-1.6861334) q[0];
sx q[0];
rz(0.88729422) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7474598) q[2];
sx q[2];
rz(-1.6051636) q[2];
sx q[2];
rz(-3.0913196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8991685) q[1];
sx q[1];
rz(-0.43546477) q[1];
sx q[1];
rz(-0.240761) q[1];
rz(-pi) q[2];
rz(-2.7769776) q[3];
sx q[3];
rz(-0.43496736) q[3];
sx q[3];
rz(1.9221969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8807184) q[2];
sx q[2];
rz(-1.4148834) q[2];
sx q[2];
rz(-0.073337642) q[2];
rz(0.65555769) q[3];
sx q[3];
rz(-0.95584241) q[3];
sx q[3];
rz(2.5478794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23621479) q[0];
sx q[0];
rz(-1.0991993) q[0];
sx q[0];
rz(-0.69860506) q[0];
rz(1.7504494) q[1];
sx q[1];
rz(-2.6222484) q[1];
sx q[1];
rz(3.0379675) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7983652) q[0];
sx q[0];
rz(-1.5916887) q[0];
sx q[0];
rz(0.038488009) q[0];
x q[1];
rz(1.6334055) q[2];
sx q[2];
rz(-0.52982989) q[2];
sx q[2];
rz(1.7018715) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5840533) q[1];
sx q[1];
rz(-1.7615093) q[1];
sx q[1];
rz(0.10277842) q[1];
rz(-0.68714924) q[3];
sx q[3];
rz(-2.8212929) q[3];
sx q[3];
rz(-0.57673467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43332064) q[2];
sx q[2];
rz(-0.89858133) q[2];
sx q[2];
rz(-1.1517322) q[2];
rz(3.057462) q[3];
sx q[3];
rz(-0.79456544) q[3];
sx q[3];
rz(-0.35552037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5750835) q[0];
sx q[0];
rz(-0.36324781) q[0];
sx q[0];
rz(-1.890924) q[0];
rz(-1.7364712) q[1];
sx q[1];
rz(-0.57650081) q[1];
sx q[1];
rz(-1.0692495) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.273928) q[0];
sx q[0];
rz(-1.8262415) q[0];
sx q[0];
rz(-0.61966611) q[0];
x q[1];
rz(-1.092717) q[2];
sx q[2];
rz(-1.5826477) q[2];
sx q[2];
rz(0.48323378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9907367) q[1];
sx q[1];
rz(-1.4889608) q[1];
sx q[1];
rz(0.86151716) q[1];
rz(-1.2573383) q[3];
sx q[3];
rz(-2.8953585) q[3];
sx q[3];
rz(1.1321958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30026597) q[2];
sx q[2];
rz(-1.0043251) q[2];
sx q[2];
rz(2.0036073) q[2];
rz(2.581572) q[3];
sx q[3];
rz(-2.0798101) q[3];
sx q[3];
rz(-1.2573857) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1244125) q[0];
sx q[0];
rz(-0.37473285) q[0];
sx q[0];
rz(0.66456932) q[0];
rz(-2.7542704) q[1];
sx q[1];
rz(-2.1699984) q[1];
sx q[1];
rz(-1.9688781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.634501) q[0];
sx q[0];
rz(-0.74247737) q[0];
sx q[0];
rz(0.93284705) q[0];
rz(-pi) q[1];
rz(0.46229273) q[2];
sx q[2];
rz(-1.3707118) q[2];
sx q[2];
rz(-2.878396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.091306578) q[1];
sx q[1];
rz(-1.0027409) q[1];
sx q[1];
rz(2.594669) q[1];
x q[2];
rz(0.4807289) q[3];
sx q[3];
rz(-0.58739118) q[3];
sx q[3];
rz(1.4241854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49501255) q[2];
sx q[2];
rz(-0.4979555) q[2];
sx q[2];
rz(0.23769561) q[2];
rz(2.2869535) q[3];
sx q[3];
rz(-1.5187289) q[3];
sx q[3];
rz(-2.2304227) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89331996) q[0];
sx q[0];
rz(-1.8017636) q[0];
sx q[0];
rz(2.7082537) q[0];
rz(-1.3029107) q[1];
sx q[1];
rz(-1.358526) q[1];
sx q[1];
rz(-0.059159577) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24640326) q[0];
sx q[0];
rz(-1.8120017) q[0];
sx q[0];
rz(-3.1144322) q[0];
rz(-pi) q[1];
rz(-2.2325977) q[2];
sx q[2];
rz(-1.7285936) q[2];
sx q[2];
rz(2.1775624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98461565) q[1];
sx q[1];
rz(-1.8530775) q[1];
sx q[1];
rz(-2.8033189) q[1];
x q[2];
rz(0.99265679) q[3];
sx q[3];
rz(-1.3380782) q[3];
sx q[3];
rz(2.6925399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2032418) q[2];
sx q[2];
rz(-1.4138736) q[2];
sx q[2];
rz(-2.4436277) q[2];
rz(0.62567726) q[3];
sx q[3];
rz(-2.8901849) q[3];
sx q[3];
rz(-1.3231369) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980001) q[0];
sx q[0];
rz(-0.55857825) q[0];
sx q[0];
rz(-2.887605) q[0];
rz(2.1506073) q[1];
sx q[1];
rz(-1.6041944) q[1];
sx q[1];
rz(-0.00016798642) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22826057) q[0];
sx q[0];
rz(-2.0096632) q[0];
sx q[0];
rz(2.069838) q[0];
rz(-pi) q[1];
rz(-0.047880574) q[2];
sx q[2];
rz(-1.6872395) q[2];
sx q[2];
rz(2.6021985) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7533299) q[1];
sx q[1];
rz(-2.7811353) q[1];
sx q[1];
rz(-0.75812104) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17642085) q[3];
sx q[3];
rz(-1.9690445) q[3];
sx q[3];
rz(0.12052025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1508472) q[2];
sx q[2];
rz(-1.3517697) q[2];
sx q[2];
rz(-3.0899437) q[2];
rz(1.026574) q[3];
sx q[3];
rz(-2.5995422) q[3];
sx q[3];
rz(0.76977229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6386665) q[0];
sx q[0];
rz(-0.8994871) q[0];
sx q[0];
rz(0.58216397) q[0];
rz(0.40099405) q[1];
sx q[1];
rz(-2.4759226) q[1];
sx q[1];
rz(-0.84025875) q[1];
rz(-2.1037965) q[2];
sx q[2];
rz(-2.0280975) q[2];
sx q[2];
rz(1.0191647) q[2];
rz(-1.6760679) q[3];
sx q[3];
rz(-1.7240763) q[3];
sx q[3];
rz(-0.45818466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
