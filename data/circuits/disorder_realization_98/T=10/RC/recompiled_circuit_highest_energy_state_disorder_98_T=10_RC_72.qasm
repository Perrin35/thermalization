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
rz(0.75818169) q[0];
sx q[0];
rz(-2.7739006) q[0];
sx q[0];
rz(-2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(2.9089622) q[1];
sx q[1];
rz(11.46041) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6152412) q[0];
sx q[0];
rz(-0.59068524) q[0];
sx q[0];
rz(2.3870941) q[0];
rz(-0.33943098) q[2];
sx q[2];
rz(-0.28305211) q[2];
sx q[2];
rz(1.7189738) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76274058) q[1];
sx q[1];
rz(-1.5334852) q[1];
sx q[1];
rz(-0.0040405063) q[1];
x q[2];
rz(-1.5638142) q[3];
sx q[3];
rz(-1.8551702) q[3];
sx q[3];
rz(2.6740536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37630633) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(0.37962309) q[2];
rz(1.6342573) q[3];
sx q[3];
rz(-2.0416656) q[3];
sx q[3];
rz(-2.1716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64604243) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(0.90091339) q[0];
rz(-3.0100929) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(2.6995755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501502) q[0];
sx q[0];
rz(-1.8953865) q[0];
sx q[0];
rz(-1.6571655) q[0];
rz(-0.0046185812) q[2];
sx q[2];
rz(-1.5717662) q[2];
sx q[2];
rz(1.0790973) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1902679) q[1];
sx q[1];
rz(-0.76808483) q[1];
sx q[1];
rz(2.4332341) q[1];
x q[2];
rz(1.6824746) q[3];
sx q[3];
rz(-1.1726716) q[3];
sx q[3];
rz(-2.3071837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76571959) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(-2.5326552) q[2];
rz(0.40306148) q[3];
sx q[3];
rz(-1.1761222) q[3];
sx q[3];
rz(2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32560638) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(-1.1380648) q[0];
rz(1.552938) q[1];
sx q[1];
rz(-2.373003) q[1];
sx q[1];
rz(-2.3345711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0332542) q[0];
sx q[0];
rz(-1.6699808) q[0];
sx q[0];
rz(-1.5981111) q[0];
rz(1.7337695) q[2];
sx q[2];
rz(-2.0103243) q[2];
sx q[2];
rz(-1.3541612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5182636) q[1];
sx q[1];
rz(-2.9179472) q[1];
sx q[1];
rz(-0.56648751) q[1];
rz(-pi) q[2];
rz(0.93385796) q[3];
sx q[3];
rz(-0.45972201) q[3];
sx q[3];
rz(0.62593725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1300065) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(1.9278795) q[2];
rz(3.0774434) q[3];
sx q[3];
rz(-1.4870653) q[3];
sx q[3];
rz(0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5478058) q[0];
sx q[0];
rz(-1.2603899) q[0];
sx q[0];
rz(-2.2547145) q[0];
rz(-2.9828494) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(1.278272) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072479362) q[0];
sx q[0];
rz(-2.2769391) q[0];
sx q[0];
rz(2.1581277) q[0];
rz(-1.099894) q[2];
sx q[2];
rz(-2.0713901) q[2];
sx q[2];
rz(2.4125227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5090829) q[1];
sx q[1];
rz(-0.69493587) q[1];
sx q[1];
rz(-1.4768697) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.015993) q[3];
sx q[3];
rz(-1.4420866) q[3];
sx q[3];
rz(0.57818613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13901237) q[2];
sx q[2];
rz(-2.2860892) q[2];
sx q[2];
rz(2.1878237) q[2];
rz(-1.9468797) q[3];
sx q[3];
rz(-0.84269968) q[3];
sx q[3];
rz(0.4755303) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139528) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(0.59999505) q[0];
rz(0.68296877) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(2.8111828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4709188) q[0];
sx q[0];
rz(-0.74711159) q[0];
sx q[0];
rz(1.8404191) q[0];
x q[1];
rz(-0.94769623) q[2];
sx q[2];
rz(-2.8665339) q[2];
sx q[2];
rz(0.81240053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8176387) q[1];
sx q[1];
rz(-1.3778701) q[1];
sx q[1];
rz(-1.0798389) q[1];
rz(2.7269269) q[3];
sx q[3];
rz(-2.2288481) q[3];
sx q[3];
rz(1.0533223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7052475) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(-0.62502965) q[2];
rz(-2.446512) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(-0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92796749) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(1.3881418) q[1];
sx q[1];
rz(-0.87166798) q[1];
sx q[1];
rz(-1.6067827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4563303) q[0];
sx q[0];
rz(-1.4233411) q[0];
sx q[0];
rz(-0.035650226) q[0];
rz(-pi) q[1];
rz(2.0892145) q[2];
sx q[2];
rz(-2.7847342) q[2];
sx q[2];
rz(-1.8722347) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88253792) q[1];
sx q[1];
rz(-1.1657622) q[1];
sx q[1];
rz(2.4292913) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2535048) q[3];
sx q[3];
rz(-2.7219238) q[3];
sx q[3];
rz(0.19007006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2024978) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(-1.6004174) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.3991791) q[3];
sx q[3];
rz(0.30103621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0114667) q[0];
sx q[0];
rz(-3.0045894) q[0];
sx q[0];
rz(-0.89114183) q[0];
rz(2.4961684) q[1];
sx q[1];
rz(-1.7452469) q[1];
sx q[1];
rz(-0.83271629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1306097) q[0];
sx q[0];
rz(-1.1101221) q[0];
sx q[0];
rz(-1.7230665) q[0];
rz(1.4786167) q[2];
sx q[2];
rz(-2.5155009) q[2];
sx q[2];
rz(2.7220059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5834076) q[1];
sx q[1];
rz(-1.3761569) q[1];
sx q[1];
rz(0.50048142) q[1];
rz(0.54033684) q[3];
sx q[3];
rz(-1.9129921) q[3];
sx q[3];
rz(-2.1006753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5523395) q[2];
sx q[2];
rz(-0.66698843) q[2];
sx q[2];
rz(1.5935295) q[2];
rz(2.7175236) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(-2.8467395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1770723) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(0.82343423) q[0];
rz(-1.1973165) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(1.4607325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13393672) q[0];
sx q[0];
rz(-1.2526576) q[0];
sx q[0];
rz(0.11548059) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1637971) q[2];
sx q[2];
rz(-1.7109181) q[2];
sx q[2];
rz(1.3205547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6951235) q[1];
sx q[1];
rz(-0.90988509) q[1];
sx q[1];
rz(1.5170694) q[1];
x q[2];
rz(0.36007215) q[3];
sx q[3];
rz(-1.2692034) q[3];
sx q[3];
rz(2.2693279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0921649) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(3.0180422) q[2];
rz(-2.3030247) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(-0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11414828) q[0];
sx q[0];
rz(-1.0076948) q[0];
sx q[0];
rz(1.4134407) q[0];
rz(2.24276) q[1];
sx q[1];
rz(-2.2347968) q[1];
sx q[1];
rz(-2.5487505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4384189) q[0];
sx q[0];
rz(-0.98383622) q[0];
sx q[0];
rz(-2.274375) q[0];
x q[1];
rz(0.18211629) q[2];
sx q[2];
rz(-1.4859293) q[2];
sx q[2];
rz(3.0013468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1072126) q[1];
sx q[1];
rz(-1.1606367) q[1];
sx q[1];
rz(2.7849547) q[1];
rz(0.96372202) q[3];
sx q[3];
rz(-2.1603909) q[3];
sx q[3];
rz(0.88966767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8574519) q[2];
sx q[2];
rz(-0.74791869) q[2];
sx q[2];
rz(0.22135529) q[2];
rz(-1.0434693) q[3];
sx q[3];
rz(-2.2011493) q[3];
sx q[3];
rz(-0.38844696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8429883) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(2.2931732) q[0];
rz(3.0449955) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(0.75327795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3101075) q[0];
sx q[0];
rz(-1.3436755) q[0];
sx q[0];
rz(0.102553) q[0];
rz(3.1084849) q[2];
sx q[2];
rz(-2.3121142) q[2];
sx q[2];
rz(0.25698369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92434525) q[1];
sx q[1];
rz(-2.0567523) q[1];
sx q[1];
rz(0.43083453) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6274458) q[3];
sx q[3];
rz(-2.3580868) q[3];
sx q[3];
rz(2.099382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0640556) q[2];
sx q[2];
rz(-0.83751837) q[2];
sx q[2];
rz(0.45292863) q[2];
rz(-0.60106599) q[3];
sx q[3];
rz(-1.4671003) q[3];
sx q[3];
rz(0.99630228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8375028) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(2.9651463) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(1.5701736) q[2];
sx q[2];
rz(-1.9194308) q[2];
sx q[2];
rz(0.44375833) q[2];
rz(-0.64519683) q[3];
sx q[3];
rz(-0.83958902) q[3];
sx q[3];
rz(-2.9506172) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
