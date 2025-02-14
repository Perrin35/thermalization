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
rz(0.14035913) q[0];
sx q[0];
rz(-1.6816513) q[0];
sx q[0];
rz(0.85293823) q[0];
rz(1.9536904) q[1];
sx q[1];
rz(-2.8953084) q[1];
sx q[1];
rz(-0.33382094) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66957314) q[0];
sx q[0];
rz(-1.8209576) q[0];
sx q[0];
rz(2.0894368) q[0];
rz(2.5025057) q[2];
sx q[2];
rz(-0.74153343) q[2];
sx q[2];
rz(1.4394176) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1934564) q[1];
sx q[1];
rz(-1.4206593) q[1];
sx q[1];
rz(0.18530042) q[1];
rz(-pi) q[2];
rz(1.6556486) q[3];
sx q[3];
rz(-1.4699664) q[3];
sx q[3];
rz(0.75750297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21981257) q[2];
sx q[2];
rz(-0.87679902) q[2];
sx q[2];
rz(2.5019257) q[2];
rz(-2.2031247) q[3];
sx q[3];
rz(-2.7191021) q[3];
sx q[3];
rz(1.2708906) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6476145) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(-1.6139503) q[0];
rz(0.24213067) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(-0.72227824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76748192) q[0];
sx q[0];
rz(-2.6038873) q[0];
sx q[0];
rz(-2.2925823) q[0];
rz(-pi) q[1];
rz(2.3467335) q[2];
sx q[2];
rz(-1.7599918) q[2];
sx q[2];
rz(-0.64293381) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5879155) q[1];
sx q[1];
rz(-1.0727296) q[1];
sx q[1];
rz(0.95446511) q[1];
x q[2];
rz(1.8037968) q[3];
sx q[3];
rz(-0.63574857) q[3];
sx q[3];
rz(2.014117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0904514) q[2];
sx q[2];
rz(-0.34153667) q[2];
sx q[2];
rz(0.086645834) q[2];
rz(-2.5610949) q[3];
sx q[3];
rz(-1.2248421) q[3];
sx q[3];
rz(0.95180029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314303) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(-1.6828368) q[0];
rz(0.1156062) q[1];
sx q[1];
rz(-1.1980201) q[1];
sx q[1];
rz(0.16673949) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78555303) q[0];
sx q[0];
rz(-2.5305439) q[0];
sx q[0];
rz(-1.9370228) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7202994) q[2];
sx q[2];
rz(-1.2308106) q[2];
sx q[2];
rz(-0.30145633) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3003242) q[1];
sx q[1];
rz(-1.027123) q[1];
sx q[1];
rz(-2.6946696) q[1];
rz(-pi) q[2];
rz(0.7451959) q[3];
sx q[3];
rz(-1.4663854) q[3];
sx q[3];
rz(-0.94693434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4358257) q[2];
sx q[2];
rz(-2.9220118) q[2];
sx q[2];
rz(-1.1337918) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.8688801) q[3];
sx q[3];
rz(2.4165966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.13919203) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(-0.27714002) q[0];
rz(2.3587522) q[1];
sx q[1];
rz(-1.6354086) q[1];
sx q[1];
rz(-0.35071075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6576516) q[0];
sx q[0];
rz(-0.707905) q[0];
sx q[0];
rz(-1.4483676) q[0];
rz(-pi) q[1];
x q[1];
rz(0.088557505) q[2];
sx q[2];
rz(-0.32677256) q[2];
sx q[2];
rz(-1.8984924) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33138946) q[1];
sx q[1];
rz(-0.94600979) q[1];
sx q[1];
rz(0.58831711) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1923741) q[3];
sx q[3];
rz(-0.11648341) q[3];
sx q[3];
rz(2.0476215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1002525) q[2];
sx q[2];
rz(-1.2790054) q[2];
sx q[2];
rz(-1.1502385) q[2];
rz(0.1746812) q[3];
sx q[3];
rz(-2.8476069) q[3];
sx q[3];
rz(3.0512419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1133076) q[0];
sx q[0];
rz(-0.57357967) q[0];
sx q[0];
rz(-1.9453402) q[0];
rz(2.664227) q[1];
sx q[1];
rz(-0.82672516) q[1];
sx q[1];
rz(0.54944077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73348896) q[0];
sx q[0];
rz(-0.86015742) q[0];
sx q[0];
rz(-0.77581328) q[0];
x q[1];
rz(1.3515768) q[2];
sx q[2];
rz(-2.0984841) q[2];
sx q[2];
rz(1.4702733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6345646) q[1];
sx q[1];
rz(-1.1643049) q[1];
sx q[1];
rz(-0.1312934) q[1];
x q[2];
rz(0.82990928) q[3];
sx q[3];
rz(-0.58778712) q[3];
sx q[3];
rz(0.15443328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4598733) q[2];
sx q[2];
rz(-0.37178603) q[2];
sx q[2];
rz(-2.8968774) q[2];
rz(-1.6230445) q[3];
sx q[3];
rz(-2.0270429) q[3];
sx q[3];
rz(-0.092930704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2731584) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(0.42700818) q[0];
rz(1.931841) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(-3.1337813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88290993) q[0];
sx q[0];
rz(-0.83246233) q[0];
sx q[0];
rz(-0.29405221) q[0];
x q[1];
rz(-1.023667) q[2];
sx q[2];
rz(-1.5627699) q[2];
sx q[2];
rz(3.0723177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1949948) q[1];
sx q[1];
rz(-1.5203069) q[1];
sx q[1];
rz(-2.4908134) q[1];
x q[2];
rz(0.2711556) q[3];
sx q[3];
rz(-2.249345) q[3];
sx q[3];
rz(-0.0038879768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4012287) q[2];
sx q[2];
rz(-2.2682891) q[2];
sx q[2];
rz(1.1391696) q[2];
rz(0.27979699) q[3];
sx q[3];
rz(-1.8179025) q[3];
sx q[3];
rz(1.6789702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4918268) q[0];
sx q[0];
rz(-0.38182807) q[0];
sx q[0];
rz(-2.5873798) q[0];
rz(0.062049374) q[1];
sx q[1];
rz(-0.65826145) q[1];
sx q[1];
rz(2.2668692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85200441) q[0];
sx q[0];
rz(-1.6920493) q[0];
sx q[0];
rz(-0.79033388) q[0];
rz(-3.1360097) q[2];
sx q[2];
rz(-1.5783596) q[2];
sx q[2];
rz(0.39952229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8043912) q[1];
sx q[1];
rz(-0.60742765) q[1];
sx q[1];
rz(-1.1056221) q[1];
x q[2];
rz(-0.90662065) q[3];
sx q[3];
rz(-1.8795089) q[3];
sx q[3];
rz(1.3332092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41500652) q[2];
sx q[2];
rz(-1.0174624) q[2];
sx q[2];
rz(-0.93878186) q[2];
rz(-0.69432652) q[3];
sx q[3];
rz(-0.66803473) q[3];
sx q[3];
rz(-0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071455) q[0];
sx q[0];
rz(-2.5818765) q[0];
sx q[0];
rz(-2.8330084) q[0];
rz(1.6584819) q[1];
sx q[1];
rz(-1.7959271) q[1];
sx q[1];
rz(-2.9187091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4321992) q[0];
sx q[0];
rz(-2.319992) q[0];
sx q[0];
rz(0.51986583) q[0];
x q[1];
rz(-0.22391386) q[2];
sx q[2];
rz(-1.6708879) q[2];
sx q[2];
rz(1.8318286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7208689) q[1];
sx q[1];
rz(-1.6096186) q[1];
sx q[1];
rz(1.8160519) q[1];
rz(-pi) q[2];
rz(-2.127366) q[3];
sx q[3];
rz(-0.96650079) q[3];
sx q[3];
rz(0.44318553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8181204) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.9984261) q[3];
sx q[3];
rz(-1.6616471) q[3];
sx q[3];
rz(1.0955048) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42300647) q[0];
sx q[0];
rz(-1.1439332) q[0];
sx q[0];
rz(-3.123172) q[0];
rz(-0.1704692) q[1];
sx q[1];
rz(-1.2455995) q[1];
sx q[1];
rz(0.19759321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391846) q[0];
sx q[0];
rz(-1.5614334) q[0];
sx q[0];
rz(1.739517) q[0];
x q[1];
rz(1.7066585) q[2];
sx q[2];
rz(-2.7721375) q[2];
sx q[2];
rz(3.0970705) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3090299) q[1];
sx q[1];
rz(-1.7062999) q[1];
sx q[1];
rz(1.2700646) q[1];
rz(-pi) q[2];
rz(0.39579795) q[3];
sx q[3];
rz(-1.5185818) q[3];
sx q[3];
rz(-0.26480103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1175055) q[2];
sx q[2];
rz(-1.0988289) q[2];
sx q[2];
rz(0.39829028) q[2];
rz(2.6309218) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(2.2834856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9061822) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(-0.19943516) q[0];
rz(0.27345744) q[1];
sx q[1];
rz(-0.43379915) q[1];
sx q[1];
rz(2.1308897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11226296) q[0];
sx q[0];
rz(-2.4284017) q[0];
sx q[0];
rz(-2.7351456) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083551466) q[2];
sx q[2];
rz(-2.9054135) q[2];
sx q[2];
rz(-2.186508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6124552) q[1];
sx q[1];
rz(-2.3061064) q[1];
sx q[1];
rz(1.6409954) q[1];
rz(-pi) q[2];
rz(1.8685249) q[3];
sx q[3];
rz(-1.862251) q[3];
sx q[3];
rz(3.1133661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79591862) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(2.3460713) q[2];
rz(-2.7703721) q[3];
sx q[3];
rz(-1.3859387) q[3];
sx q[3];
rz(0.25446874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(0.026074792) q[0];
sx q[0];
rz(-1.7255029) q[0];
sx q[0];
rz(1.2988476) q[0];
rz(2.5904291) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(-2.416257) q[2];
sx q[2];
rz(-1.4015522) q[2];
sx q[2];
rz(1.3649469) q[2];
rz(0.24088151) q[3];
sx q[3];
rz(-1.1987562) q[3];
sx q[3];
rz(-0.22507122) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
