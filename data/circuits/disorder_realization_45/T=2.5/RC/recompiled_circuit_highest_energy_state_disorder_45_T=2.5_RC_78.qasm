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
rz(-0.51195872) q[0];
sx q[0];
rz(6.6867642) q[0];
sx q[0];
rz(9.4288958) q[0];
rz(0.68139684) q[1];
sx q[1];
rz(-1.2702785) q[1];
sx q[1];
rz(-1.7008002) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3128634) q[0];
sx q[0];
rz(-1.848319) q[0];
sx q[0];
rz(0.1779495) q[0];
rz(2.7782562) q[2];
sx q[2];
rz(-0.78509313) q[2];
sx q[2];
rz(1.8232283) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3057249) q[1];
sx q[1];
rz(-2.222371) q[1];
sx q[1];
rz(2.0209794) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1295092) q[3];
sx q[3];
rz(-1.7444495) q[3];
sx q[3];
rz(-2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9422841) q[2];
sx q[2];
rz(-2.0950623) q[2];
sx q[2];
rz(-0.47041565) q[2];
rz(0.89620245) q[3];
sx q[3];
rz(-1.6737409) q[3];
sx q[3];
rz(1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50297058) q[0];
sx q[0];
rz(-0.48415411) q[0];
sx q[0];
rz(-2.7726987) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.8577441) q[1];
sx q[1];
rz(0.2643815) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6579012) q[0];
sx q[0];
rz(-2.2122447) q[0];
sx q[0];
rz(0.19199706) q[0];
rz(-1.7272378) q[2];
sx q[2];
rz(-1.3231734) q[2];
sx q[2];
rz(-0.81849097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7355613) q[1];
sx q[1];
rz(-1.2792865) q[1];
sx q[1];
rz(-0.31707615) q[1];
rz(-pi) q[2];
rz(3.0204078) q[3];
sx q[3];
rz(-0.78599343) q[3];
sx q[3];
rz(-1.3898897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5276864) q[2];
sx q[2];
rz(-2.9759585) q[2];
sx q[2];
rz(1.3045093) q[2];
rz(-0.3197318) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(0.50362292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9673135) q[0];
sx q[0];
rz(-2.9496084) q[0];
sx q[0];
rz(2.018003) q[0];
rz(-2.7536821) q[1];
sx q[1];
rz(-1.8759517) q[1];
sx q[1];
rz(2.8111615) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39836568) q[0];
sx q[0];
rz(-1.8332229) q[0];
sx q[0];
rz(0.18375476) q[0];
rz(-pi) q[1];
rz(-1.6638043) q[2];
sx q[2];
rz(-0.96538097) q[2];
sx q[2];
rz(-0.85697848) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8113241) q[1];
sx q[1];
rz(-1.2172592) q[1];
sx q[1];
rz(-1.2297022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5009242) q[3];
sx q[3];
rz(-1.0721777) q[3];
sx q[3];
rz(-2.9151268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72184163) q[2];
sx q[2];
rz(-1.1417049) q[2];
sx q[2];
rz(1.9640131) q[2];
rz(-1.2490595) q[3];
sx q[3];
rz(-1.306059) q[3];
sx q[3];
rz(-0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053442001) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(-3.0528659) q[0];
rz(-0.42452043) q[1];
sx q[1];
rz(-2.4219234) q[1];
sx q[1];
rz(-1.4094062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30520327) q[0];
sx q[0];
rz(-0.38581784) q[0];
sx q[0];
rz(-1.057339) q[0];
rz(-pi) q[1];
rz(-2.2327044) q[2];
sx q[2];
rz(-0.7134921) q[2];
sx q[2];
rz(-0.59407083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.030050412) q[1];
sx q[1];
rz(-0.71978986) q[1];
sx q[1];
rz(2.5538302) q[1];
rz(2.0982101) q[3];
sx q[3];
rz(-2.5549915) q[3];
sx q[3];
rz(2.9709904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7523664) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-2.7092095) q[2];
rz(-0.45807517) q[3];
sx q[3];
rz(-1.4832486) q[3];
sx q[3];
rz(-2.1336011) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3573989) q[0];
sx q[0];
rz(-0.14506871) q[0];
sx q[0];
rz(-2.8302637) q[0];
rz(-1.1773102) q[1];
sx q[1];
rz(-1.177634) q[1];
sx q[1];
rz(-2.6572773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949849) q[0];
sx q[0];
rz(-1.4123045) q[0];
sx q[0];
rz(-0.096317795) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3179146) q[2];
sx q[2];
rz(-0.30159471) q[2];
sx q[2];
rz(-0.67267928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2023485) q[1];
sx q[1];
rz(-0.40196291) q[1];
sx q[1];
rz(-1.9415929) q[1];
x q[2];
rz(-0.14967646) q[3];
sx q[3];
rz(-1.6972741) q[3];
sx q[3];
rz(1.5095131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8208661) q[2];
sx q[2];
rz(-1.7118688) q[2];
sx q[2];
rz(0.42496625) q[2];
rz(-0.10415569) q[3];
sx q[3];
rz(-0.36915532) q[3];
sx q[3];
rz(-0.66515508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7729618) q[0];
sx q[0];
rz(-0.63987982) q[0];
sx q[0];
rz(-2.9848918) q[0];
rz(2.8796097) q[1];
sx q[1];
rz(-1.1527088) q[1];
sx q[1];
rz(-1.6798457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77685415) q[0];
sx q[0];
rz(-2.1663453) q[0];
sx q[0];
rz(0.2100581) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.842672) q[2];
sx q[2];
rz(-2.4266234) q[2];
sx q[2];
rz(0.78503099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9380381) q[1];
sx q[1];
rz(-0.21760908) q[1];
sx q[1];
rz(-0.59534351) q[1];
rz(-0.42243345) q[3];
sx q[3];
rz(-2.3452873) q[3];
sx q[3];
rz(-0.56047201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4713952) q[2];
sx q[2];
rz(-1.4963701) q[2];
sx q[2];
rz(0.71693286) q[2];
rz(-0.21643058) q[3];
sx q[3];
rz(-1.8559034) q[3];
sx q[3];
rz(1.9560248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4099429) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(1.2662079) q[0];
rz(2.3833497) q[1];
sx q[1];
rz(-1.9634602) q[1];
sx q[1];
rz(-2.0936802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9564932) q[0];
sx q[0];
rz(-2.4742715) q[0];
sx q[0];
rz(1.0769847) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0133466) q[2];
sx q[2];
rz(-1.1801022) q[2];
sx q[2];
rz(1.3732131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5541682) q[1];
sx q[1];
rz(-1.5012073) q[1];
sx q[1];
rz(1.2072366) q[1];
x q[2];
rz(-0.45067421) q[3];
sx q[3];
rz(-1.3503805) q[3];
sx q[3];
rz(0.31639755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8083814) q[2];
sx q[2];
rz(-1.0071249) q[2];
sx q[2];
rz(0.0086616596) q[2];
rz(-2.1374785) q[3];
sx q[3];
rz(-0.4929556) q[3];
sx q[3];
rz(1.0021817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.778331) q[0];
sx q[0];
rz(-2.9982428) q[0];
sx q[0];
rz(2.1891201) q[0];
rz(-1.5337503) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(2.105377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0909611) q[0];
sx q[0];
rz(-2.1078175) q[0];
sx q[0];
rz(-1.0616674) q[0];
x q[1];
rz(0.63663738) q[2];
sx q[2];
rz(-1.271406) q[2];
sx q[2];
rz(0.68975291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1055323) q[1];
sx q[1];
rz(-2.1478473) q[1];
sx q[1];
rz(0.60571508) q[1];
rz(-0.50004543) q[3];
sx q[3];
rz(-0.79179057) q[3];
sx q[3];
rz(-2.4337976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2967534) q[2];
sx q[2];
rz(-2.4256458) q[2];
sx q[2];
rz(-2.1412264) q[2];
rz(-1.865546) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(-2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412096) q[0];
sx q[0];
rz(-2.2333133) q[0];
sx q[0];
rz(-0.25417438) q[0];
rz(1.3379478) q[1];
sx q[1];
rz(-0.86054069) q[1];
sx q[1];
rz(-2.3635704) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8786273) q[0];
sx q[0];
rz(-0.74889442) q[0];
sx q[0];
rz(-0.48567943) q[0];
x q[1];
rz(-1.4073952) q[2];
sx q[2];
rz(-1.3721352) q[2];
sx q[2];
rz(0.29115788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7791479) q[1];
sx q[1];
rz(-1.0082642) q[1];
sx q[1];
rz(1.0845119) q[1];
rz(1.3555525) q[3];
sx q[3];
rz(-1.4741033) q[3];
sx q[3];
rz(1.7413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.079387) q[2];
sx q[2];
rz(-1.6239245) q[2];
sx q[2];
rz(0.13266955) q[2];
rz(-1.1072055) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.4440822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4001813) q[0];
sx q[0];
rz(-1.1962698) q[0];
sx q[0];
rz(-2.3924526) q[0];
rz(-0.49508849) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(-1.7503768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0780236) q[0];
sx q[0];
rz(-1.1207523) q[0];
sx q[0];
rz(2.4343632) q[0];
rz(-pi) q[1];
x q[1];
rz(1.130213) q[2];
sx q[2];
rz(-1.2017219) q[2];
sx q[2];
rz(-2.1598494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0799082) q[1];
sx q[1];
rz(-2.7322953) q[1];
sx q[1];
rz(3.0319935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0416609) q[3];
sx q[3];
rz(-1.1253998) q[3];
sx q[3];
rz(-0.22280773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5007925) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(-2.2643209) q[2];
rz(-0.52214617) q[3];
sx q[3];
rz(-2.3418661) q[3];
sx q[3];
rz(-1.4570025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397639) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(-1.8632035) q[1];
sx q[1];
rz(-2.8291193) q[1];
sx q[1];
rz(3.1183174) q[1];
rz(-1.6017492) q[2];
sx q[2];
rz(-0.94018117) q[2];
sx q[2];
rz(-1.2311819) q[2];
rz(-0.69576453) q[3];
sx q[3];
rz(-0.56752612) q[3];
sx q[3];
rz(-0.023434536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
