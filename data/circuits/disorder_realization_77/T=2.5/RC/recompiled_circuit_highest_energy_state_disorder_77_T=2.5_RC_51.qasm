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
rz(-0.4975118) q[0];
sx q[0];
rz(4.4805718) q[0];
sx q[0];
rz(6.5988402) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(0.60400909) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16306388) q[0];
sx q[0];
rz(-1.8730358) q[0];
sx q[0];
rz(0.053044293) q[0];
rz(-pi) q[1];
rz(2.2411795) q[2];
sx q[2];
rz(-1.1657823) q[2];
sx q[2];
rz(-0.23460282) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0993201) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(2.7541942) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8798001) q[3];
sx q[3];
rz(-1.4674392) q[3];
sx q[3];
rz(0.37783937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(-1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58594054) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(1.9441388) q[0];
rz(-0.50152957) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(-1.7053568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2510808) q[0];
sx q[0];
rz(-0.62411004) q[0];
sx q[0];
rz(1.9735543) q[0];
rz(-1.4714437) q[2];
sx q[2];
rz(-0.87026419) q[2];
sx q[2];
rz(1.5347028) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2982218) q[1];
sx q[1];
rz(-2.9195115) q[1];
sx q[1];
rz(-2.6574617) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3615492) q[3];
sx q[3];
rz(-1.1236785) q[3];
sx q[3];
rz(-2.1802878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88840914) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-0.02296981) q[2];
rz(2.2394771) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70812923) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(-2.1500812) q[0];
rz(2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69312364) q[0];
sx q[0];
rz(-2.4938994) q[0];
sx q[0];
rz(2.7231587) q[0];
rz(-0.63629432) q[2];
sx q[2];
rz(-0.71517309) q[2];
sx q[2];
rz(1.6805122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4146862) q[1];
sx q[1];
rz(-0.79830805) q[1];
sx q[1];
rz(1.9195791) q[1];
rz(-pi) q[2];
rz(-2.1622873) q[3];
sx q[3];
rz(-2.0404173) q[3];
sx q[3];
rz(-2.8036433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.027355) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(1.2973971) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39795136) q[0];
sx q[0];
rz(-2.209111) q[0];
sx q[0];
rz(-2.3241296) q[0];
rz(-0.3262597) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(0.78525966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9049587) q[0];
sx q[0];
rz(-1.3044622) q[0];
sx q[0];
rz(1.8840709) q[0];
x q[1];
rz(2.8815124) q[2];
sx q[2];
rz(-0.65979119) q[2];
sx q[2];
rz(1.5957956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1568415) q[1];
sx q[1];
rz(-1.9923216) q[1];
sx q[1];
rz(-2.0578029) q[1];
x q[2];
rz(-1.6503851) q[3];
sx q[3];
rz(-1.9246939) q[3];
sx q[3];
rz(0.45512629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(2.2480887) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(-1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23984443) q[0];
sx q[0];
rz(-2.8185066) q[0];
sx q[0];
rz(-1.595994) q[0];
rz(1.0218989) q[1];
sx q[1];
rz(-2.0183759) q[1];
sx q[1];
rz(-2.9603069) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820532) q[0];
sx q[0];
rz(-1.7632369) q[0];
sx q[0];
rz(1.0715225) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6659545) q[2];
sx q[2];
rz(-1.5283094) q[2];
sx q[2];
rz(-1.3459537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0376772) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(1.7930536) q[1];
x q[2];
rz(2.6986609) q[3];
sx q[3];
rz(-1.7128403) q[3];
sx q[3];
rz(-1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1229317) q[2];
sx q[2];
rz(-2.4666726) q[2];
sx q[2];
rz(-1.8956511) q[2];
rz(2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-0.30884185) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(-3.0016532) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788743) q[0];
sx q[0];
rz(-2.1669183) q[0];
sx q[0];
rz(-2.8403502) q[0];
x q[1];
rz(1.762319) q[2];
sx q[2];
rz(-0.46326783) q[2];
sx q[2];
rz(-0.5199711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4637101) q[1];
sx q[1];
rz(-1.2695128) q[1];
sx q[1];
rz(0.82656411) q[1];
rz(-pi) q[2];
rz(2.841137) q[3];
sx q[3];
rz(-2.2078035) q[3];
sx q[3];
rz(-2.9954122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-0.24197401) q[2];
rz(2.3954929) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.022843) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(1.7005327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25172397) q[0];
sx q[0];
rz(-1.8419203) q[0];
sx q[0];
rz(0.95377484) q[0];
rz(2.7602642) q[2];
sx q[2];
rz(-1.6994119) q[2];
sx q[2];
rz(0.53024697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59380925) q[1];
sx q[1];
rz(-1.9873706) q[1];
sx q[1];
rz(-1.3273456) q[1];
rz(1.3541) q[3];
sx q[3];
rz(-0.85734493) q[3];
sx q[3];
rz(2.5503623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.8249576) q[2];
rz(2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.104326) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-0.88687801) q[0];
rz(-2.9414224) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(-1.0221457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048269317) q[0];
sx q[0];
rz(-1.5962068) q[0];
sx q[0];
rz(1.7019468) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3273507) q[2];
sx q[2];
rz(-0.98746429) q[2];
sx q[2];
rz(-2.5098206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31412582) q[1];
sx q[1];
rz(-0.74605251) q[1];
sx q[1];
rz(1.57342) q[1];
rz(-pi) q[2];
rz(2.6273055) q[3];
sx q[3];
rz(-2.5011241) q[3];
sx q[3];
rz(-0.91292229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6518121) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(-2.1583648) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(-2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(1.5552833) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(2.3568025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12778388) q[0];
sx q[0];
rz(-1.133636) q[0];
sx q[0];
rz(-2.5378835) q[0];
rz(0.92918877) q[2];
sx q[2];
rz(-0.25114533) q[2];
sx q[2];
rz(-0.26640688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5604374) q[1];
sx q[1];
rz(-2.0056917) q[1];
sx q[1];
rz(-1.5337318) q[1];
rz(0.611245) q[3];
sx q[3];
rz(-1.9211624) q[3];
sx q[3];
rz(-0.65920748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(1.2189216) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(-1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(-0.75540677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9616325) q[0];
sx q[0];
rz(-0.20997071) q[0];
sx q[0];
rz(1.5487973) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.5546397) q[2];
sx q[2];
rz(1.1262058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81138071) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(2.5274171) q[1];
rz(-pi) q[2];
rz(2.741119) q[3];
sx q[3];
rz(-2.4254834) q[3];
sx q[3];
rz(-3.0587089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5648254) q[2];
sx q[2];
rz(-2.8456523) q[2];
sx q[2];
rz(0.15288615) q[2];
rz(-1.8704002) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(0.66711867) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-1.1350994) q[2];
sx q[2];
rz(-1.0183327) q[2];
sx q[2];
rz(-0.20295126) q[2];
rz(0.85930227) q[3];
sx q[3];
rz(-0.51325428) q[3];
sx q[3];
rz(-0.2556066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
