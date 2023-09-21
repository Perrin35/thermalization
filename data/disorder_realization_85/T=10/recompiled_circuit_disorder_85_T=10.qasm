OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.207333624362946) q[0];
sx q[0];
rz(3.73195579846437) q[0];
sx q[0];
rz(9.05376020669147) q[0];
rz(-0.381292164325714) q[1];
sx q[1];
rz(2.5420904477411) q[1];
sx q[1];
rz(11.1903435945432) q[1];
cx q[1],q[0];
rz(0.57173079252243) q[0];
sx q[0];
rz(4.10802713234956) q[0];
sx q[0];
rz(10.4532739877622) q[0];
rz(0.491012930870056) q[2];
sx q[2];
rz(4.05836639006669) q[2];
sx q[2];
rz(10.1354425907056) q[2];
cx q[2],q[1];
rz(0.103780969977379) q[1];
sx q[1];
rz(4.56072130997712) q[1];
sx q[1];
rz(9.75592822431728) q[1];
rz(-0.0561532825231552) q[3];
sx q[3];
rz(3.89382109244401) q[3];
sx q[3];
rz(9.97741750477954) q[3];
cx q[3],q[2];
rz(1.0842889547348) q[2];
sx q[2];
rz(3.54192710121209) q[2];
sx q[2];
rz(10.4136960863988) q[2];
rz(0.752542853355408) q[3];
sx q[3];
rz(4.28745201428468) q[3];
sx q[3];
rz(10.2555489897649) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.207240208983421) q[0];
sx q[0];
rz(3.25385686208541) q[0];
sx q[0];
rz(10.6047091245572) q[0];
rz(-0.997693657875061) q[1];
sx q[1];
rz(4.42479899724061) q[1];
sx q[1];
rz(10.1490906834523) q[1];
cx q[1],q[0];
rz(0.225310102105141) q[0];
sx q[0];
rz(3.90123841364915) q[0];
sx q[0];
rz(9.88877902030154) q[0];
rz(-1.17643141746521) q[2];
sx q[2];
rz(2.9260202815109) q[2];
sx q[2];
rz(12.8992316484372) q[2];
cx q[2],q[1];
rz(0.380729407072067) q[1];
sx q[1];
rz(3.55337572296197) q[1];
sx q[1];
rz(9.95983991622134) q[1];
rz(0.879942417144775) q[3];
sx q[3];
rz(2.83755311568315) q[3];
sx q[3];
rz(9.3344731092374) q[3];
cx q[3],q[2];
rz(-0.905368447303772) q[2];
sx q[2];
rz(4.95277443726594) q[2];
sx q[2];
rz(9.78666359781429) q[2];
rz(0.136067166924477) q[3];
sx q[3];
rz(3.69729486306245) q[3];
sx q[3];
rz(9.47040213494703) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.14194619655609) q[0];
sx q[0];
rz(4.23522678216035) q[0];
sx q[0];
rz(10.8202116250913) q[0];
rz(0.462292581796646) q[1];
sx q[1];
rz(2.71700829465921) q[1];
sx q[1];
rz(10.6438634157102) q[1];
cx q[1],q[0];
rz(0.314194709062576) q[0];
sx q[0];
rz(3.58726069529588) q[0];
sx q[0];
rz(10.0721755981366) q[0];
rz(1.27499186992645) q[2];
sx q[2];
rz(4.03265527089173) q[2];
sx q[2];
rz(9.97543088196918) q[2];
cx q[2],q[1];
rz(-0.298750072717667) q[1];
sx q[1];
rz(4.18395796616609) q[1];
sx q[1];
rz(10.715190267555) q[1];
rz(-0.0562682300806046) q[3];
sx q[3];
rz(3.41000309784944) q[3];
sx q[3];
rz(10.4997886180799) q[3];
cx q[3],q[2];
rz(0.421577125787735) q[2];
sx q[2];
rz(3.88183936675126) q[2];
sx q[2];
rz(7.97921619414493) q[2];
rz(0.56882631778717) q[3];
sx q[3];
rz(3.99131909211213) q[3];
sx q[3];
rz(9.53529494851037) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.223594263195992) q[0];
sx q[0];
rz(3.58979129989678) q[0];
sx q[0];
rz(9.94303790330097) q[0];
rz(0.715424299240112) q[1];
sx q[1];
rz(4.25787845452363) q[1];
sx q[1];
rz(10.25153449773) q[1];
cx q[1],q[0];
rz(0.0872388333082199) q[0];
sx q[0];
rz(4.2103363593393) q[0];
sx q[0];
rz(9.89616603254482) q[0];
rz(0.134379431605339) q[2];
sx q[2];
rz(3.94925669034059) q[2];
sx q[2];
rz(9.09623629450008) q[2];
cx q[2],q[1];
rz(-0.0668478533625603) q[1];
sx q[1];
rz(4.42758611043031) q[1];
sx q[1];
rz(10.5975917339246) q[1];
rz(0.84639984369278) q[3];
sx q[3];
rz(1.80762341816957) q[3];
sx q[3];
rz(9.95737866162463) q[3];
cx q[3],q[2];
rz(0.152390792965889) q[2];
sx q[2];
rz(3.34057168860967) q[2];
sx q[2];
rz(10.8036901712339) q[2];
rz(0.0723230242729187) q[3];
sx q[3];
rz(3.95402643282945) q[3];
sx q[3];
rz(11.0701622724454) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.826926589012146) q[0];
sx q[0];
rz(3.76171991427476) q[0];
sx q[0];
rz(10.5506653547208) q[0];
rz(0.902448832988739) q[1];
sx q[1];
rz(4.11548909743363) q[1];
sx q[1];
rz(9.70994765161678) q[1];
cx q[1],q[0];
rz(1.06411588191986) q[0];
sx q[0];
rz(2.66490939457948) q[0];
sx q[0];
rz(9.517301259927) q[0];
rz(1.08918142318726) q[2];
sx q[2];
rz(2.41846052010591) q[2];
sx q[2];
rz(9.11899573206111) q[2];
cx q[2],q[1];
rz(0.83851820230484) q[1];
sx q[1];
rz(3.45464867551858) q[1];
sx q[1];
rz(9.42341895344808) q[1];
rz(1.16915917396545) q[3];
sx q[3];
rz(4.62321534951264) q[3];
sx q[3];
rz(9.55562198757335) q[3];
cx q[3],q[2];
rz(0.29331648349762) q[2];
sx q[2];
rz(2.63603761990602) q[2];
sx q[2];
rz(9.96423128842517) q[2];
rz(0.306824922561646) q[3];
sx q[3];
rz(4.02611205180223) q[3];
sx q[3];
rz(8.97056514619991) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.258131116628647) q[0];
sx q[0];
rz(3.8920243700319) q[0];
sx q[0];
rz(9.58179088532134) q[0];
rz(-0.693333864212036) q[1];
sx q[1];
rz(2.26088288624818) q[1];
sx q[1];
rz(10.791789507858) q[1];
cx q[1],q[0];
rz(-0.0889286622405052) q[0];
sx q[0];
rz(3.21731445391709) q[0];
sx q[0];
rz(8.97353554367229) q[0];
rz(0.603468298912048) q[2];
sx q[2];
rz(3.99788585503633) q[2];
sx q[2];
rz(9.13467562793895) q[2];
cx q[2],q[1];
rz(-0.659466028213501) q[1];
sx q[1];
rz(3.44218602974946) q[1];
sx q[1];
rz(9.14868206381007) q[1];
rz(-0.103619866073132) q[3];
sx q[3];
rz(4.36115971406037) q[3];
sx q[3];
rz(10.2988617181699) q[3];
cx q[3],q[2];
rz(-0.296290159225464) q[2];
sx q[2];
rz(3.99160370429093) q[2];
sx q[2];
rz(9.82824549674197) q[2];
rz(0.481633126735687) q[3];
sx q[3];
rz(4.21375218232209) q[3];
sx q[3];
rz(9.94401512145206) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.242131739854813) q[0];
sx q[0];
rz(4.02488246758515) q[0];
sx q[0];
rz(10.2986407637517) q[0];
rz(-0.447723984718323) q[1];
sx q[1];
rz(3.88060030539567) q[1];
sx q[1];
rz(11.3956774234693) q[1];
cx q[1],q[0];
rz(1.31603515148163) q[0];
sx q[0];
rz(3.16543277923996) q[0];
sx q[0];
rz(10.7339695453565) q[0];
rz(-1.15877294540405) q[2];
sx q[2];
rz(3.26831931074197) q[2];
sx q[2];
rz(11.8088760137479) q[2];
cx q[2],q[1];
rz(1.2504186630249) q[1];
sx q[1];
rz(4.10112831194932) q[1];
sx q[1];
rz(10.1659834742467) q[1];
rz(1.0366370677948) q[3];
sx q[3];
rz(4.44934466679628) q[3];
sx q[3];
rz(10.5402817487638) q[3];
cx q[3],q[2];
rz(2.04470229148865) q[2];
sx q[2];
rz(3.7108001430803) q[2];
sx q[2];
rz(8.81402382849857) q[2];
rz(0.475105345249176) q[3];
sx q[3];
rz(4.23215440114076) q[3];
sx q[3];
rz(10.3525191903035) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.241918131709099) q[0];
sx q[0];
rz(3.02452589024837) q[0];
sx q[0];
rz(9.72190389632388) q[0];
rz(1.74694526195526) q[1];
sx q[1];
rz(4.29257586796815) q[1];
sx q[1];
rz(10.0709100723188) q[1];
cx q[1],q[0];
rz(0.329833656549454) q[0];
sx q[0];
rz(2.91148971219594) q[0];
sx q[0];
rz(8.96685416101619) q[0];
rz(-1.24936759471893) q[2];
sx q[2];
rz(2.64769634802873) q[2];
sx q[2];
rz(11.4981660604398) q[2];
cx q[2],q[1];
rz(0.384780019521713) q[1];
sx q[1];
rz(5.054467471438) q[1];
sx q[1];
rz(10.4635635375898) q[1];
rz(0.633873879909515) q[3];
sx q[3];
rz(4.00968209107453) q[3];
sx q[3];
rz(9.40932618583694) q[3];
cx q[3],q[2];
rz(0.0646765530109406) q[2];
sx q[2];
rz(4.09464588959748) q[2];
sx q[2];
rz(9.87961540221378) q[2];
rz(-0.701397657394409) q[3];
sx q[3];
rz(4.1720062812143) q[3];
sx q[3];
rz(10.558802342407) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.699445366859436) q[0];
sx q[0];
rz(3.73064139683778) q[0];
sx q[0];
rz(10.2222848892133) q[0];
rz(0.517567157745361) q[1];
sx q[1];
rz(3.95420995553071) q[1];
sx q[1];
rz(9.53319562076732) q[1];
cx q[1],q[0];
rz(0.781288206577301) q[0];
sx q[0];
rz(3.80786046584184) q[0];
sx q[0];
rz(9.51831234841748) q[0];
rz(-0.252656847238541) q[2];
sx q[2];
rz(3.32600337465341) q[2];
sx q[2];
rz(10.2694743037145) q[2];
cx q[2],q[1];
rz(0.1154959872365) q[1];
sx q[1];
rz(2.86710876424844) q[1];
sx q[1];
rz(10.0848985075872) q[1];
rz(1.75767457485199) q[3];
sx q[3];
rz(4.01675322850282) q[3];
sx q[3];
rz(9.8076213657777) q[3];
cx q[3],q[2];
rz(2.12915277481079) q[2];
sx q[2];
rz(4.48845043976838) q[2];
sx q[2];
rz(9.76473111509486) q[2];
rz(0.418395042419434) q[3];
sx q[3];
rz(2.54516259034211) q[3];
sx q[3];
rz(10.1503692030828) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.53239345550537) q[0];
sx q[0];
rz(3.53549981315667) q[0];
sx q[0];
rz(8.74595836400195) q[0];
rz(-0.364183068275452) q[1];
sx q[1];
rz(4.58576026757295) q[1];
sx q[1];
rz(9.47993674724504) q[1];
cx q[1],q[0];
rz(0.776277422904968) q[0];
sx q[0];
rz(2.75898295839364) q[0];
sx q[0];
rz(11.1695016384046) q[0];
rz(1.01110219955444) q[2];
sx q[2];
rz(4.53373065789277) q[2];
sx q[2];
rz(9.87389022707149) q[2];
cx q[2],q[1];
rz(0.897964715957642) q[1];
sx q[1];
rz(4.12635979254777) q[1];
sx q[1];
rz(8.7577690243642) q[1];
rz(0.0746598690748215) q[3];
sx q[3];
rz(3.59870308836038) q[3];
sx q[3];
rz(9.74648374914333) q[3];
cx q[3],q[2];
rz(0.983833432197571) q[2];
sx q[2];
rz(4.16301360924775) q[2];
sx q[2];
rz(8.93460161089107) q[2];
rz(0.137520775198936) q[3];
sx q[3];
rz(4.24118080933625) q[3];
sx q[3];
rz(10.3628655433576) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.0482234954834) q[0];
sx q[0];
rz(4.07405588229234) q[0];
sx q[0];
rz(9.82654569148227) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.276412755250931) q[1];
sx q[1];
rz(3.59051889379556) q[1];
sx q[1];
rz(10.9345360755841) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.570307612419128) q[2];
sx q[2];
rz(3.48932284315164) q[2];
sx q[2];
rz(9.69851506351634) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.398756206035614) q[3];
sx q[3];
rz(3.41036808689172) q[3];
sx q[3];
rz(9.26025398670837) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];