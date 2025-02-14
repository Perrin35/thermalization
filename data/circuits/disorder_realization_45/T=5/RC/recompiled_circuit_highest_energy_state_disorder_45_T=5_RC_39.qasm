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
rz(-1.4650605) q[0];
sx q[0];
rz(-0.2413916) q[0];
sx q[0];
rz(0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(-0.018996039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762704) q[0];
sx q[0];
rz(-0.92060584) q[0];
sx q[0];
rz(-2.8702535) q[0];
rz(-pi) q[1];
x q[1];
rz(0.094474205) q[2];
sx q[2];
rz(-1.3850537) q[2];
sx q[2];
rz(0.31982012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87352301) q[1];
sx q[1];
rz(-1.8933081) q[1];
sx q[1];
rz(-1.3871865) q[1];
rz(-pi) q[2];
rz(2.9844445) q[3];
sx q[3];
rz(-1.8007319) q[3];
sx q[3];
rz(0.63369753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0211109) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(0.79180229) q[2];
rz(2.0558489) q[3];
sx q[3];
rz(-1.9224242) q[3];
sx q[3];
rz(2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.87078142) q[0];
sx q[0];
rz(-0.15853156) q[0];
sx q[0];
rz(0.32919163) q[0];
rz(-1.5022494) q[1];
sx q[1];
rz(-0.90198016) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5188816) q[0];
sx q[0];
rz(-0.71463138) q[0];
sx q[0];
rz(2.463752) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8453526) q[2];
sx q[2];
rz(-2.5736817) q[2];
sx q[2];
rz(2.8366249) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93041673) q[1];
sx q[1];
rz(-2.8143344) q[1];
sx q[1];
rz(0.14463592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6738191) q[3];
sx q[3];
rz(-0.80108445) q[3];
sx q[3];
rz(-2.5689606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4210356) q[2];
sx q[2];
rz(-1.9562419) q[2];
sx q[2];
rz(2.0015008) q[2];
rz(-0.0044048443) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(-0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36188257) q[0];
sx q[0];
rz(-0.10040586) q[0];
sx q[0];
rz(-0.34565872) q[0];
rz(1.0935621) q[1];
sx q[1];
rz(-2.3193017) q[1];
sx q[1];
rz(1.0543157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39798007) q[0];
sx q[0];
rz(-0.52310399) q[0];
sx q[0];
rz(-2.4260957) q[0];
x q[1];
rz(3.0385706) q[2];
sx q[2];
rz(-0.81562519) q[2];
sx q[2];
rz(1.4519004) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2002751) q[1];
sx q[1];
rz(-1.6345342) q[1];
sx q[1];
rz(2.6885808) q[1];
rz(-pi) q[2];
rz(2.2369713) q[3];
sx q[3];
rz(-1.5890749) q[3];
sx q[3];
rz(1.3410904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67636079) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(-0.54222995) q[2];
rz(-0.17255653) q[3];
sx q[3];
rz(-1.6514643) q[3];
sx q[3];
rz(-2.0358613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365823) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(1.2611058) q[0];
rz(2.005596) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(-0.13066185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4379556) q[0];
sx q[0];
rz(-1.2017631) q[0];
sx q[0];
rz(-2.8519408) q[0];
x q[1];
rz(-0.50640743) q[2];
sx q[2];
rz(-1.3739488) q[2];
sx q[2];
rz(2.0692689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9798292) q[1];
sx q[1];
rz(-1.6024622) q[1];
sx q[1];
rz(1.0624485) q[1];
x q[2];
rz(2.4774136) q[3];
sx q[3];
rz(-1.2746433) q[3];
sx q[3];
rz(0.68063762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7589492) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(0.22155133) q[2];
rz(2.4394636) q[3];
sx q[3];
rz(-1.5789072) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57109433) q[0];
sx q[0];
rz(-1.8966738) q[0];
sx q[0];
rz(-3.0076497) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.1993273) q[1];
sx q[1];
rz(-1.3901002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54467359) q[0];
sx q[0];
rz(-1.8159165) q[0];
sx q[0];
rz(-2.0952203) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2657792) q[2];
sx q[2];
rz(-1.3172564) q[2];
sx q[2];
rz(0.60088742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83727969) q[1];
sx q[1];
rz(-1.8472478) q[1];
sx q[1];
rz(0.56378341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36365328) q[3];
sx q[3];
rz(-1.8580274) q[3];
sx q[3];
rz(-1.938323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(1.585539) q[2];
rz(2.8499917) q[3];
sx q[3];
rz(-2.6652938) q[3];
sx q[3];
rz(2.8628023) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6911102) q[0];
sx q[0];
rz(-2.1463558) q[0];
sx q[0];
rz(-2.6158748) q[0];
rz(0.37824962) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(2.8674616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12488406) q[0];
sx q[0];
rz(-1.5965466) q[0];
sx q[0];
rz(-2.5514249) q[0];
rz(2.3304105) q[2];
sx q[2];
rz(-2.5497782) q[2];
sx q[2];
rz(-0.20840157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8822599) q[1];
sx q[1];
rz(-1.3162398) q[1];
sx q[1];
rz(1.1544636) q[1];
rz(2.0977705) q[3];
sx q[3];
rz(-0.29114215) q[3];
sx q[3];
rz(3.1028735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5452925) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(-0.25823414) q[2];
rz(-1.5958512) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(-0.405092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3426568) q[0];
sx q[0];
rz(-2.6234143) q[0];
sx q[0];
rz(0.10547353) q[0];
rz(0.27490973) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(2.8555433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28562322) q[0];
sx q[0];
rz(-0.21846314) q[0];
sx q[0];
rz(-0.56665786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83467612) q[2];
sx q[2];
rz(-1.7913941) q[2];
sx q[2];
rz(0.59360628) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5040999) q[1];
sx q[1];
rz(-1.5063725) q[1];
sx q[1];
rz(1.9570051) q[1];
rz(0.55283847) q[3];
sx q[3];
rz(-1.1854608) q[3];
sx q[3];
rz(0.60147731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1470571) q[2];
sx q[2];
rz(-0.68238443) q[2];
sx q[2];
rz(2.6915754) q[2];
rz(-2.5543645) q[3];
sx q[3];
rz(-1.3415895) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29379544) q[0];
sx q[0];
rz(-1.6814517) q[0];
sx q[0];
rz(-1.1867123) q[0];
rz(2.8222491) q[1];
sx q[1];
rz(-2.1198544) q[1];
sx q[1];
rz(-2.879338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569993) q[0];
sx q[0];
rz(-2.5525064) q[0];
sx q[0];
rz(2.7688945) q[0];
rz(-pi) q[1];
rz(3.0056001) q[2];
sx q[2];
rz(-1.1609224) q[2];
sx q[2];
rz(-0.52025822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7062807) q[1];
sx q[1];
rz(-1.4888933) q[1];
sx q[1];
rz(2.1141073) q[1];
x q[2];
rz(-2.3877034) q[3];
sx q[3];
rz(-2.5544832) q[3];
sx q[3];
rz(2.052149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.78476) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(2.812815) q[2];
rz(-1.5152991) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78557712) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(0.27895862) q[0];
rz(-0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(-1.4252211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12869975) q[0];
sx q[0];
rz(-0.94966799) q[0];
sx q[0];
rz(1.5170044) q[0];
x q[1];
rz(0.80323146) q[2];
sx q[2];
rz(-2.0798426) q[2];
sx q[2];
rz(-2.5521297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9659075) q[1];
sx q[1];
rz(-1.07465) q[1];
sx q[1];
rz(-1.7860852) q[1];
x q[2];
rz(1.2012977) q[3];
sx q[3];
rz(-0.69108057) q[3];
sx q[3];
rz(-1.406351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7384537) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(-2.9816755) q[2];
rz(-0.81965172) q[3];
sx q[3];
rz(-1.9194226) q[3];
sx q[3];
rz(-2.405449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4495471) q[0];
sx q[0];
rz(-1.291438) q[0];
sx q[0];
rz(-2.8458169) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(1.6326509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7180824) q[0];
sx q[0];
rz(-2.5959466) q[0];
sx q[0];
rz(2.3033124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9231173) q[2];
sx q[2];
rz(-1.6197471) q[2];
sx q[2];
rz(-1.6214329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9450284) q[1];
sx q[1];
rz(-1.6383645) q[1];
sx q[1];
rz(-0.089971926) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4911777) q[3];
sx q[3];
rz(-2.5456508) q[3];
sx q[3];
rz(0.015344674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8136924) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(-2.5076765) q[2];
rz(-3.0564803) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(-0.77272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91372981) q[0];
sx q[0];
rz(-1.5416523) q[0];
sx q[0];
rz(-0.1097485) q[0];
rz(1.9204503) q[1];
sx q[1];
rz(-0.84538645) q[1];
sx q[1];
rz(0.56199817) q[1];
rz(-0.22460266) q[2];
sx q[2];
rz(-1.9997659) q[2];
sx q[2];
rz(-1.275296) q[2];
rz(1.3391277) q[3];
sx q[3];
rz(-1.5944407) q[3];
sx q[3];
rz(-3.1138314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
