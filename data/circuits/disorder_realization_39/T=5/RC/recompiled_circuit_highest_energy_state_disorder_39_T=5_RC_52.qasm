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
rz(-3.0012335) q[0];
sx q[0];
rz(-1.4599414) q[0];
sx q[0];
rz(2.2886544) q[0];
rz(1.9536904) q[1];
sx q[1];
rz(-2.8953084) q[1];
sx q[1];
rz(-0.33382094) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1000017) q[0];
sx q[0];
rz(-1.0698338) q[0];
sx q[0];
rz(2.8554625) q[0];
rz(-pi) q[1];
rz(-1.0707955) q[2];
sx q[2];
rz(-0.99783932) q[2];
sx q[2];
rz(2.2286602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1934564) q[1];
sx q[1];
rz(-1.4206593) q[1];
sx q[1];
rz(-2.9562922) q[1];
x q[2];
rz(0.6972972) q[3];
sx q[3];
rz(-0.13168959) q[3];
sx q[3];
rz(3.0856709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9217801) q[2];
sx q[2];
rz(-2.2647936) q[2];
sx q[2];
rz(0.63966695) q[2];
rz(2.2031247) q[3];
sx q[3];
rz(-2.7191021) q[3];
sx q[3];
rz(1.870702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6476145) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(1.5276424) q[0];
rz(2.899462) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(-2.4193144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3741107) q[0];
sx q[0];
rz(-0.53770533) q[0];
sx q[0];
rz(0.84901036) q[0];
rz(-pi) q[1];
rz(-0.79485915) q[2];
sx q[2];
rz(-1.3816009) q[2];
sx q[2];
rz(0.64293381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4237883) q[1];
sx q[1];
rz(-2.3700938) q[1];
sx q[1];
rz(2.3256734) q[1];
x q[2];
rz(1.3377959) q[3];
sx q[3];
rz(-0.63574857) q[3];
sx q[3];
rz(-2.014117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0904514) q[2];
sx q[2];
rz(-2.800056) q[2];
sx q[2];
rz(3.0549468) q[2];
rz(2.5610949) q[3];
sx q[3];
rz(-1.9167506) q[3];
sx q[3];
rz(-2.1897924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314303) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(-1.4587559) q[0];
rz(3.0259865) q[1];
sx q[1];
rz(-1.1980201) q[1];
sx q[1];
rz(-0.16673949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3560396) q[0];
sx q[0];
rz(-0.61104873) q[0];
sx q[0];
rz(1.9370228) q[0];
x q[1];
rz(0.42129321) q[2];
sx q[2];
rz(-1.910782) q[2];
sx q[2];
rz(-0.30145633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6550811) q[1];
sx q[1];
rz(-1.1919199) q[1];
sx q[1];
rz(-0.98030555) q[1];
rz(-pi) q[2];
rz(-1.7124246) q[3];
sx q[3];
rz(-0.83061355) q[3];
sx q[3];
rz(0.52799559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4358257) q[2];
sx q[2];
rz(-0.21958084) q[2];
sx q[2];
rz(1.1337918) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.8688801) q[3];
sx q[3];
rz(2.4165966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.13919203) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(-2.8644526) q[0];
rz(2.3587522) q[1];
sx q[1];
rz(-1.506184) q[1];
sx q[1];
rz(0.35071075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6445054) q[0];
sx q[0];
rz(-0.86928029) q[0];
sx q[0];
rz(3.0374466) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8160118) q[2];
sx q[2];
rz(-1.5424038) q[2];
sx q[2];
rz(2.7300026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8102032) q[1];
sx q[1];
rz(-0.94600979) q[1];
sx q[1];
rz(-0.58831711) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0983885) q[3];
sx q[3];
rz(-1.4625878) q[3];
sx q[3];
rz(-1.6668591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1002525) q[2];
sx q[2];
rz(-1.8625872) q[2];
sx q[2];
rz(1.9913541) q[2];
rz(-0.1746812) q[3];
sx q[3];
rz(-0.29398578) q[3];
sx q[3];
rz(-0.090350769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133076) q[0];
sx q[0];
rz(-2.568013) q[0];
sx q[0];
rz(-1.9453402) q[0];
rz(2.664227) q[1];
sx q[1];
rz(-0.82672516) q[1];
sx q[1];
rz(0.54944077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8910946) q[0];
sx q[0];
rz(-2.142494) q[0];
sx q[0];
rz(0.88776161) q[0];
rz(-pi) q[1];
rz(2.7844593) q[2];
sx q[2];
rz(-2.5741842) q[2];
sx q[2];
rz(-1.0536989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0116033) q[1];
sx q[1];
rz(-1.450256) q[1];
sx q[1];
rz(1.9804363) q[1];
rz(-1.1138238) q[3];
sx q[3];
rz(-1.9544056) q[3];
sx q[3];
rz(2.3758604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68171936) q[2];
sx q[2];
rz(-2.7698066) q[2];
sx q[2];
rz(-0.2447153) q[2];
rz(-1.6230445) q[3];
sx q[3];
rz(-1.1145498) q[3];
sx q[3];
rz(0.092930704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86843425) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(2.7145845) q[0];
rz(1.931841) q[1];
sx q[1];
rz(-1.5245583) q[1];
sx q[1];
rz(3.1337813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2526379) q[0];
sx q[0];
rz(-1.7868306) q[0];
sx q[0];
rz(-2.3310069) q[0];
x q[1];
rz(-1.023667) q[2];
sx q[2];
rz(-1.5788227) q[2];
sx q[2];
rz(-3.0723177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1949948) q[1];
sx q[1];
rz(-1.5203069) q[1];
sx q[1];
rz(-2.4908134) q[1];
rz(-pi) q[2];
rz(0.2711556) q[3];
sx q[3];
rz(-0.89224766) q[3];
sx q[3];
rz(-3.1377047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4012287) q[2];
sx q[2];
rz(-2.2682891) q[2];
sx q[2];
rz(-1.1391696) q[2];
rz(2.8617957) q[3];
sx q[3];
rz(-1.8179025) q[3];
sx q[3];
rz(1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6497659) q[0];
sx q[0];
rz(-0.38182807) q[0];
sx q[0];
rz(2.5873798) q[0];
rz(-0.062049374) q[1];
sx q[1];
rz(-2.4833312) q[1];
sx q[1];
rz(-0.87472349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84034568) q[0];
sx q[0];
rz(-2.3537427) q[0];
sx q[0];
rz(1.3993174) q[0];
rz(-1.563233) q[2];
sx q[2];
rz(-1.5652135) q[2];
sx q[2];
rz(1.9703608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8043912) q[1];
sx q[1];
rz(-2.534165) q[1];
sx q[1];
rz(-2.0359705) q[1];
rz(-0.90662065) q[3];
sx q[3];
rz(-1.2620838) q[3];
sx q[3];
rz(-1.3332092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7265861) q[2];
sx q[2];
rz(-2.1241302) q[2];
sx q[2];
rz(-2.2028108) q[2];
rz(-2.4472661) q[3];
sx q[3];
rz(-0.66803473) q[3];
sx q[3];
rz(-2.9850849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071455) q[0];
sx q[0];
rz(-0.55971611) q[0];
sx q[0];
rz(2.8330084) q[0];
rz(1.4831108) q[1];
sx q[1];
rz(-1.7959271) q[1];
sx q[1];
rz(2.9187091) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5102986) q[0];
sx q[0];
rz(-1.9430854) q[0];
sx q[0];
rz(-2.3907803) q[0];
rz(-pi) q[1];
rz(-2.7168458) q[2];
sx q[2];
rz(-0.24492376) q[2];
sx q[2];
rz(0.67451678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8376779) q[1];
sx q[1];
rz(-2.8933446) q[1];
sx q[1];
rz(1.7294238) q[1];
x q[2];
rz(0.68273516) q[3];
sx q[3];
rz(-2.020524) q[3];
sx q[3];
rz(-2.3537824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3234723) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.1431665) q[3];
sx q[3];
rz(-1.4799456) q[3];
sx q[3];
rz(-2.0460879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42300647) q[0];
sx q[0];
rz(-1.1439332) q[0];
sx q[0];
rz(-0.018420694) q[0];
rz(0.1704692) q[1];
sx q[1];
rz(-1.2455995) q[1];
sx q[1];
rz(-0.19759321) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.402408) q[0];
sx q[0];
rz(-1.5801593) q[0];
sx q[0];
rz(1.4020756) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.204448) q[2];
sx q[2];
rz(-1.6197259) q[2];
sx q[2];
rz(1.6530703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8325628) q[1];
sx q[1];
rz(-1.7062999) q[1];
sx q[1];
rz(-1.2700646) q[1];
x q[2];
rz(-1.5142158) q[3];
sx q[3];
rz(-1.9660249) q[3];
sx q[3];
rz(-1.8574018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0240872) q[2];
sx q[2];
rz(-1.0988289) q[2];
sx q[2];
rz(-0.39829028) q[2];
rz(-0.5106709) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(-0.85810703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9061822) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(-2.9421575) q[0];
rz(-0.27345744) q[1];
sx q[1];
rz(-0.43379915) q[1];
sx q[1];
rz(-2.1308897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0293297) q[0];
sx q[0];
rz(-0.71319095) q[0];
sx q[0];
rz(-0.40644706) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083551466) q[2];
sx q[2];
rz(-0.23617911) q[2];
sx q[2];
rz(-0.95508466) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5291374) q[1];
sx q[1];
rz(-2.3061064) q[1];
sx q[1];
rz(-1.5005972) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2730678) q[3];
sx q[3];
rz(-1.2793417) q[3];
sx q[3];
rz(-3.1133661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.345674) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(0.79552135) q[2];
rz(2.7703721) q[3];
sx q[3];
rz(-1.755654) q[3];
sx q[3];
rz(0.25446874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1155179) q[0];
sx q[0];
rz(-1.7255029) q[0];
sx q[0];
rz(1.2988476) q[0];
rz(-0.55116354) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(-1.3462832) q[2];
sx q[2];
rz(-2.2835352) q[2];
sx q[2];
rz(2.78751) q[2];
rz(-0.24088151) q[3];
sx q[3];
rz(-1.9428365) q[3];
sx q[3];
rz(2.9165214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
