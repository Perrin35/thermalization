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
rz(1.9089215) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(2.2751459) q[1];
sx q[1];
rz(12.187764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0526992) q[0];
sx q[0];
rz(-1.8009773) q[0];
sx q[0];
rz(3.102688) q[0];
rz(-pi) q[1];
rz(0.9742175) q[2];
sx q[2];
rz(-1.2132036) q[2];
sx q[2];
rz(1.9973988) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12593658) q[1];
sx q[1];
rz(-1.7864461) q[1];
sx q[1];
rz(0.68839781) q[1];
rz(1.1932329) q[3];
sx q[3];
rz(-1.0306944) q[3];
sx q[3];
rz(-1.3256955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0014701) q[2];
sx q[2];
rz(-1.1727611) q[2];
sx q[2];
rz(-0.44974652) q[2];
rz(-0.55284119) q[3];
sx q[3];
rz(-0.55838412) q[3];
sx q[3];
rz(2.9302178) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387892) q[0];
sx q[0];
rz(-1.1454104) q[0];
sx q[0];
rz(2.3853724) q[0];
rz(2.4839632) q[1];
sx q[1];
rz(-2.2537474) q[1];
sx q[1];
rz(-0.58849803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158101) q[0];
sx q[0];
rz(-1.2001769) q[0];
sx q[0];
rz(-2.4354706) q[0];
rz(-pi) q[1];
rz(-2.3173554) q[2];
sx q[2];
rz(-1.7331799) q[2];
sx q[2];
rz(0.33467143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1900764) q[1];
sx q[1];
rz(-1.562968) q[1];
sx q[1];
rz(-1.8633134) q[1];
x q[2];
rz(-1.1682092) q[3];
sx q[3];
rz(-2.5036771) q[3];
sx q[3];
rz(1.0653374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6077891) q[2];
sx q[2];
rz(-1.6705931) q[2];
sx q[2];
rz(-3.0652453) q[2];
rz(-2.7401183) q[3];
sx q[3];
rz(-2.6186826) q[3];
sx q[3];
rz(-2.0126066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5717995) q[0];
sx q[0];
rz(-0.58572584) q[0];
sx q[0];
rz(-2.8049923) q[0];
rz(1.0353237) q[1];
sx q[1];
rz(-2.2575049) q[1];
sx q[1];
rz(3.091277) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6487728) q[0];
sx q[0];
rz(-2.0413715) q[0];
sx q[0];
rz(-2.5008766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59479721) q[2];
sx q[2];
rz(-1.4001048) q[2];
sx q[2];
rz(-1.5399982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2256098) q[1];
sx q[1];
rz(-0.40357737) q[1];
sx q[1];
rz(-0.55761375) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21304275) q[3];
sx q[3];
rz(-1.0096978) q[3];
sx q[3];
rz(-1.9764501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3745554) q[2];
sx q[2];
rz(-2.2646246) q[2];
sx q[2];
rz(-0.98131895) q[2];
rz(1.0513002) q[3];
sx q[3];
rz(-1.4562573) q[3];
sx q[3];
rz(3.1276957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0392847) q[0];
sx q[0];
rz(-1.6599052) q[0];
sx q[0];
rz(2.1639977) q[0];
rz(2.5915937) q[1];
sx q[1];
rz(-0.98098749) q[1];
sx q[1];
rz(0.94211334) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6459759) q[0];
sx q[0];
rz(-1.8919404) q[0];
sx q[0];
rz(-0.34713388) q[0];
rz(0.14288516) q[2];
sx q[2];
rz(-2.7804903) q[2];
sx q[2];
rz(2.2074769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3636605) q[1];
sx q[1];
rz(-0.97423755) q[1];
sx q[1];
rz(-2.9371757) q[1];
rz(-pi) q[2];
rz(1.0364384) q[3];
sx q[3];
rz(-2.0490408) q[3];
sx q[3];
rz(-0.75936067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.00086870988) q[2];
sx q[2];
rz(-0.372118) q[2];
sx q[2];
rz(-1.5127399) q[2];
rz(-0.77423972) q[3];
sx q[3];
rz(-0.95572487) q[3];
sx q[3];
rz(0.16404185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0999488) q[0];
sx q[0];
rz(-0.14687563) q[0];
sx q[0];
rz(-2.3351093) q[0];
rz(-2.3355314) q[1];
sx q[1];
rz(-1.1763923) q[1];
sx q[1];
rz(0.75256601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8984932) q[0];
sx q[0];
rz(-2.4499735) q[0];
sx q[0];
rz(-1.3893632) q[0];
x q[1];
rz(-3.0522992) q[2];
sx q[2];
rz(-2.7460423) q[2];
sx q[2];
rz(-1.7035005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8991685) q[1];
sx q[1];
rz(-2.7061279) q[1];
sx q[1];
rz(0.240761) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7349929) q[3];
sx q[3];
rz(-1.9754275) q[3];
sx q[3];
rz(-1.5238289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8807184) q[2];
sx q[2];
rz(-1.4148834) q[2];
sx q[2];
rz(-3.068255) q[2];
rz(-0.65555769) q[3];
sx q[3];
rz(-2.1857502) q[3];
sx q[3];
rz(2.5478794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9053779) q[0];
sx q[0];
rz(-1.0991993) q[0];
sx q[0];
rz(0.69860506) q[0];
rz(1.7504494) q[1];
sx q[1];
rz(-0.51934424) q[1];
sx q[1];
rz(-3.0379675) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2267644) q[0];
sx q[0];
rz(-1.5323167) q[0];
sx q[0];
rz(-1.5498885) q[0];
rz(-pi) q[1];
rz(1.6334055) q[2];
sx q[2];
rz(-0.52982989) q[2];
sx q[2];
rz(-1.4397211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0857398) q[1];
sx q[1];
rz(-0.21634783) q[1];
sx q[1];
rz(2.0592709) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68714924) q[3];
sx q[3];
rz(-0.32029974) q[3];
sx q[3];
rz(-0.57673467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.708272) q[2];
sx q[2];
rz(-0.89858133) q[2];
sx q[2];
rz(1.1517322) q[2];
rz(3.057462) q[3];
sx q[3];
rz(-2.3470272) q[3];
sx q[3];
rz(-2.7860723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5750835) q[0];
sx q[0];
rz(-0.36324781) q[0];
sx q[0];
rz(1.890924) q[0];
rz(-1.4051215) q[1];
sx q[1];
rz(-0.57650081) q[1];
sx q[1];
rz(1.0692495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0435283) q[0];
sx q[0];
rz(-2.4777924) q[0];
sx q[0];
rz(2.7190156) q[0];
rz(-pi) q[1];
rz(-1.092717) q[2];
sx q[2];
rz(-1.5826477) q[2];
sx q[2];
rz(0.48323378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6266963) q[1];
sx q[1];
rz(-2.428423) q[1];
sx q[1];
rz(-1.4455224) q[1];
rz(1.2573383) q[3];
sx q[3];
rz(-2.8953585) q[3];
sx q[3];
rz(-1.1321958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8413267) q[2];
sx q[2];
rz(-2.1372676) q[2];
sx q[2];
rz(2.0036073) q[2];
rz(0.56002069) q[3];
sx q[3];
rz(-2.0798101) q[3];
sx q[3];
rz(-1.8842069) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244125) q[0];
sx q[0];
rz(-0.37473285) q[0];
sx q[0];
rz(-0.66456932) q[0];
rz(-2.7542704) q[1];
sx q[1];
rz(-0.97159425) q[1];
sx q[1];
rz(1.9688781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50709162) q[0];
sx q[0];
rz(-0.74247737) q[0];
sx q[0];
rz(-2.2087456) q[0];
x q[1];
rz(0.46229273) q[2];
sx q[2];
rz(-1.3707118) q[2];
sx q[2];
rz(-2.878396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9786611) q[1];
sx q[1];
rz(-1.1169945) q[1];
sx q[1];
rz(2.2125208) q[1];
x q[2];
rz(-2.6608638) q[3];
sx q[3];
rz(-0.58739118) q[3];
sx q[3];
rz(-1.7174073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49501255) q[2];
sx q[2];
rz(-0.4979555) q[2];
sx q[2];
rz(-2.903897) q[2];
rz(-0.85463917) q[3];
sx q[3];
rz(-1.5187289) q[3];
sx q[3];
rz(0.91116992) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89331996) q[0];
sx q[0];
rz(-1.3398291) q[0];
sx q[0];
rz(-2.7082537) q[0];
rz(-1.8386819) q[1];
sx q[1];
rz(-1.7830667) q[1];
sx q[1];
rz(-0.059159577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8951894) q[0];
sx q[0];
rz(-1.8120017) q[0];
sx q[0];
rz(3.1144322) q[0];
rz(-pi) q[1];
rz(-2.9425609) q[2];
sx q[2];
rz(-0.91863721) q[2];
sx q[2];
rz(2.6566344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48850075) q[1];
sx q[1];
rz(-1.2464123) q[1];
sx q[1];
rz(1.2725194) q[1];
x q[2];
rz(2.1489359) q[3];
sx q[3];
rz(-1.8035144) q[3];
sx q[3];
rz(2.6925399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2032418) q[2];
sx q[2];
rz(-1.7277191) q[2];
sx q[2];
rz(0.69796491) q[2];
rz(-2.5159154) q[3];
sx q[3];
rz(-2.8901849) q[3];
sx q[3];
rz(-1.3231369) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043592602) q[0];
sx q[0];
rz(-0.55857825) q[0];
sx q[0];
rz(0.25398764) q[0];
rz(2.1506073) q[1];
sx q[1];
rz(-1.5373983) q[1];
sx q[1];
rz(-3.1414247) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68037409) q[0];
sx q[0];
rz(-0.6520642) q[0];
sx q[0];
rz(-0.79508938) q[0];
rz(-1.1824047) q[2];
sx q[2];
rz(-3.0157308) q[2];
sx q[2];
rz(2.2110155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7385369) q[1];
sx q[1];
rz(-1.3118032) q[1];
sx q[1];
rz(1.8243755) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1669257) q[3];
sx q[3];
rz(-1.40831) q[3];
sx q[3];
rz(-1.5193017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1508472) q[2];
sx q[2];
rz(-1.3517697) q[2];
sx q[2];
rz(0.051648971) q[2];
rz(1.026574) q[3];
sx q[3];
rz(-0.54205042) q[3];
sx q[3];
rz(2.3718204) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3401716) q[2];
sx q[2];
rz(-2.4540171) q[2];
sx q[2];
rz(-3.050741) q[2];
rz(2.9874728) q[3];
sx q[3];
rz(-1.4667635) q[3];
sx q[3];
rz(1.0964805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
