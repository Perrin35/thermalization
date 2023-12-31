OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(1.4899878) q[0];
sx q[0];
rz(8.4943354) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(1.1073444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410755) q[0];
sx q[0];
rz(-2.374875) q[0];
sx q[0];
rz(-0.47857743) q[0];
rz(-pi) q[1];
rz(0.47714699) q[2];
sx q[2];
rz(-2.2331182) q[2];
sx q[2];
rz(-1.1610111) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2335637) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(-1.8506552) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9355448) q[3];
sx q[3];
rz(-1.8555292) q[3];
sx q[3];
rz(-2.2176544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.8135653) q[2];
rz(-0.32087457) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-0.11238614) q[0];
sx q[0];
rz(0.88062084) q[0];
rz(1.2940787) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-0.81726384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016804455) q[0];
sx q[0];
rz(-1.4298555) q[0];
sx q[0];
rz(1.9190448) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0613282) q[2];
sx q[2];
rz(-1.763952) q[2];
sx q[2];
rz(-0.9888538) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17727795) q[1];
sx q[1];
rz(-0.69523584) q[1];
sx q[1];
rz(-2.7362105) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2137787) q[3];
sx q[3];
rz(-2.6414053) q[3];
sx q[3];
rz(-0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(0.98637995) q[3];
sx q[3];
rz(-1.7247518) q[3];
sx q[3];
rz(1.4646437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(2.1752775) q[0];
rz(-0.19293383) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(1.4470709) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149949) q[0];
sx q[0];
rz(-1.5883755) q[0];
sx q[0];
rz(0.21501644) q[0];
x q[1];
rz(-1.3005199) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(-1.9384055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24070534) q[1];
sx q[1];
rz(-1.1966238) q[1];
sx q[1];
rz(-0.64970533) q[1];
rz(-pi) q[2];
rz(-2.8964642) q[3];
sx q[3];
rz(-1.9760625) q[3];
sx q[3];
rz(0.95642904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43859279) q[2];
sx q[2];
rz(-2.6575228) q[2];
sx q[2];
rz(-1.1052216) q[2];
rz(0.74622074) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(-2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352017) q[0];
sx q[0];
rz(-0.52440301) q[0];
sx q[0];
rz(-1.6756469) q[0];
rz(0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-0.38898653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60880946) q[0];
sx q[0];
rz(-2.7739077) q[0];
sx q[0];
rz(0.68433783) q[0];
x q[1];
rz(-1.2126036) q[2];
sx q[2];
rz(-0.51414031) q[2];
sx q[2];
rz(-2.6280623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0491838) q[1];
sx q[1];
rz(-2.1158943) q[1];
sx q[1];
rz(2.0712907) q[1];
rz(-pi) q[2];
rz(2.3936478) q[3];
sx q[3];
rz(-1.7224632) q[3];
sx q[3];
rz(2.9880854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7230364) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(-0.45670613) q[2];
rz(-1.5151954) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(-2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(-0.49450758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3798824) q[0];
sx q[0];
rz(-2.8344791) q[0];
sx q[0];
rz(-2.7193927) q[0];
rz(0.039426609) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(2.6089422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33489409) q[1];
sx q[1];
rz(-0.94167275) q[1];
sx q[1];
rz(1.6449528) q[1];
rz(-3.0881808) q[3];
sx q[3];
rz(-2.6969389) q[3];
sx q[3];
rz(-0.85944552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(-2.8430856) q[2];
rz(-0.31202894) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(-0.19590713) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(-1.235199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2081731) q[0];
sx q[0];
rz(-2.6801077) q[0];
sx q[0];
rz(2.3955406) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2321129) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(-2.0896926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8630353) q[1];
sx q[1];
rz(-0.99579358) q[1];
sx q[1];
rz(1.0689736) q[1];
x q[2];
rz(-0.33552334) q[3];
sx q[3];
rz(-2.8918859) q[3];
sx q[3];
rz(-1.5229131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-0.43169272) q[2];
rz(-1.7290944) q[3];
sx q[3];
rz(-2.4192211) q[3];
sx q[3];
rz(-3.0055962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7504808) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(0.11211638) q[0];
rz(-2.926459) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047785) q[0];
sx q[0];
rz(-1.830173) q[0];
sx q[0];
rz(-0.39774261) q[0];
rz(-pi) q[1];
rz(-2.8261975) q[2];
sx q[2];
rz(-1.3239667) q[2];
sx q[2];
rz(2.18404) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73247319) q[1];
sx q[1];
rz(-1.2670244) q[1];
sx q[1];
rz(0.71883808) q[1];
rz(-pi) q[2];
rz(-0.21344276) q[3];
sx q[3];
rz(-1.1323954) q[3];
sx q[3];
rz(1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.001174288) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.4814608) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(1.9326899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-0.27420843) q[0];
sx q[0];
rz(1.2742395) q[0];
rz(2.3102343) q[2];
sx q[2];
rz(-2.0908329) q[2];
sx q[2];
rz(-1.3607197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0749082) q[1];
sx q[1];
rz(-1.9004596) q[1];
sx q[1];
rz(-0.7633022) q[1];
x q[2];
rz(-2.2015757) q[3];
sx q[3];
rz(-1.7490083) q[3];
sx q[3];
rz(2.7959787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37844354) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086583) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(1.2110442) q[1];
sx q[1];
rz(-0.90463224) q[1];
sx q[1];
rz(-0.39168721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2101333) q[0];
sx q[0];
rz(-1.1362695) q[0];
sx q[0];
rz(-1.5772485) q[0];
rz(-0.45231818) q[2];
sx q[2];
rz(-1.7512133) q[2];
sx q[2];
rz(2.7197321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.081025) q[1];
sx q[1];
rz(-2.276366) q[1];
sx q[1];
rz(2.702436) q[1];
rz(-pi) q[2];
rz(3.016032) q[3];
sx q[3];
rz(-0.52913044) q[3];
sx q[3];
rz(0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(2.3796066) q[3];
sx q[3];
rz(-0.31969324) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-15/(14*pi)) q[0];
sx q[0];
rz(-2.8388192) q[0];
sx q[0];
rz(-2.5706932) q[0];
rz(-1.7123429) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(2.9796519) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93962651) q[0];
sx q[0];
rz(-1.8585748) q[0];
sx q[0];
rz(2.1956325) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0091153092) q[2];
sx q[2];
rz(-1.3822123) q[2];
sx q[2];
rz(-1.7878704) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3871799) q[1];
sx q[1];
rz(-1.3993235) q[1];
sx q[1];
rz(-3.0110554) q[1];
rz(-pi) q[2];
rz(1.5374244) q[3];
sx q[3];
rz(-2.0936692) q[3];
sx q[3];
rz(-0.066730412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(-1.9421633) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(-1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7006871) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.3700925) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(0.8926819) q[2];
sx q[2];
rz(-2.953139) q[2];
sx q[2];
rz(2.1072731) q[2];
rz(1.3508425) q[3];
sx q[3];
rz(-1.150288) q[3];
sx q[3];
rz(1.8899959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
