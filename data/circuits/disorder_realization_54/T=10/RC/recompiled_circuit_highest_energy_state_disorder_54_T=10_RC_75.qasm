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
rz(-0.71159166) q[0];
sx q[0];
rz(-2.0366259) q[0];
sx q[0];
rz(1.1397064) q[0];
rz(-3.9330339) q[1];
sx q[1];
rz(7.3242261) q[1];
sx q[1];
rz(8.2952226) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0985467) q[0];
sx q[0];
rz(-0.73972964) q[0];
sx q[0];
rz(-1.84124) q[0];
rz(2.4784577) q[2];
sx q[2];
rz(-1.5485933) q[2];
sx q[2];
rz(2.210169) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0186386) q[1];
sx q[1];
rz(-1.2771416) q[1];
sx q[1];
rz(-0.83033009) q[1];
rz(-1.6057062) q[3];
sx q[3];
rz(-1.2381609) q[3];
sx q[3];
rz(1.2768318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1019885) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(2.2339036) q[2];
rz(3.0229819) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-0.25092009) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(3.0372341) q[0];
rz(1.1907499) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(2.6844535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9741544) q[0];
sx q[0];
rz(-1.3773241) q[0];
sx q[0];
rz(-0.03137145) q[0];
rz(-pi) q[1];
rz(2.1188583) q[2];
sx q[2];
rz(-0.75621683) q[2];
sx q[2];
rz(0.05847419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7434843) q[1];
sx q[1];
rz(-2.8116075) q[1];
sx q[1];
rz(1.8673926) q[1];
rz(-pi) q[2];
rz(1.3966772) q[3];
sx q[3];
rz(-2.4655364) q[3];
sx q[3];
rz(0.0043178252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(0.15698329) q[2];
rz(-0.5213151) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(-0.014001525) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-2.795862) q[0];
rz(-1.6861457) q[1];
sx q[1];
rz(-1.2724178) q[1];
sx q[1];
rz(1.5706496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81996213) q[0];
sx q[0];
rz(-1.4819711) q[0];
sx q[0];
rz(-2.1169244) q[0];
rz(-pi) q[1];
rz(0.0035946058) q[2];
sx q[2];
rz(-2.1934641) q[2];
sx q[2];
rz(1.1284356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6863043) q[1];
sx q[1];
rz(-1.521811) q[1];
sx q[1];
rz(-2.9279885) q[1];
rz(-1.8351052) q[3];
sx q[3];
rz(-0.66625957) q[3];
sx q[3];
rz(0.245417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7741144) q[2];
sx q[2];
rz(-2.9134637) q[2];
sx q[2];
rz(-2.1841614) q[2];
rz(-2.0188324) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289309) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(1.3900025) q[1];
sx q[1];
rz(-1.8252239) q[1];
sx q[1];
rz(1.4768627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6679803) q[0];
sx q[0];
rz(-1.851871) q[0];
sx q[0];
rz(2.5949508) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9042911) q[2];
sx q[2];
rz(-2.8899404) q[2];
sx q[2];
rz(-0.66019765) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37472607) q[1];
sx q[1];
rz(-2.0927168) q[1];
sx q[1];
rz(-2.7609227) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.847193) q[3];
sx q[3];
rz(-1.2528462) q[3];
sx q[3];
rz(-2.2459787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59820286) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(-0.99679917) q[2];
rz(-1.8761084) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-2.8640981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0811512) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(1.0720217) q[0];
rz(2.187166) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(1.6486453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38781751) q[0];
sx q[0];
rz(-2.3507581) q[0];
sx q[0];
rz(-2.8641939) q[0];
rz(-pi) q[1];
rz(0.49333879) q[2];
sx q[2];
rz(-2.3559921) q[2];
sx q[2];
rz(2.0480905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4665597) q[1];
sx q[1];
rz(-1.4283115) q[1];
sx q[1];
rz(0.33054466) q[1];
x q[2];
rz(-2.8478283) q[3];
sx q[3];
rz(-0.33393327) q[3];
sx q[3];
rz(-1.6373843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91168779) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(0.041042717) q[3];
sx q[3];
rz(-1.5464455) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(1.1596229) q[0];
rz(-1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(-1.9083091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430968) q[0];
sx q[0];
rz(-2.2944106) q[0];
sx q[0];
rz(1.9886524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3911925) q[2];
sx q[2];
rz(-0.078243464) q[2];
sx q[2];
rz(0.83759826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6813876) q[1];
sx q[1];
rz(-2.3355995) q[1];
sx q[1];
rz(0.56350033) q[1];
rz(-pi) q[2];
rz(-2.7090577) q[3];
sx q[3];
rz(-2.5721241) q[3];
sx q[3];
rz(1.3506952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50501862) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(-1.9912857) q[2];
rz(0.73927528) q[3];
sx q[3];
rz(-1.247568) q[3];
sx q[3];
rz(-2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4742541) q[0];
sx q[0];
rz(-0.69304729) q[0];
sx q[0];
rz(0.44878238) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-2.0346784) q[1];
sx q[1];
rz(-1.1610228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338116) q[0];
sx q[0];
rz(-0.49612653) q[0];
sx q[0];
rz(2.5563942) q[0];
rz(-pi) q[1];
rz(1.0801187) q[2];
sx q[2];
rz(-2.0188234) q[2];
sx q[2];
rz(2.6971044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22299448) q[1];
sx q[1];
rz(-1.9004993) q[1];
sx q[1];
rz(-1.8818284) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1584985) q[3];
sx q[3];
rz(-0.59368836) q[3];
sx q[3];
rz(-2.3883461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9728969) q[2];
sx q[2];
rz(-1.8099512) q[2];
sx q[2];
rz(-1.9608344) q[2];
rz(-2.0598038) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(-2.7150174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5931554) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(2.4539808) q[0];
rz(-1.9718735) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(-2.7489472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9967742) q[0];
sx q[0];
rz(-0.4985362) q[0];
sx q[0];
rz(1.2720435) q[0];
rz(-pi) q[1];
rz(-2.8296109) q[2];
sx q[2];
rz(-1.7034966) q[2];
sx q[2];
rz(-3.0614982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4325786) q[1];
sx q[1];
rz(-2.6582632) q[1];
sx q[1];
rz(0.18928115) q[1];
rz(-0.18889922) q[3];
sx q[3];
rz(-1.4903063) q[3];
sx q[3];
rz(0.98505011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9795867) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(1.4957734) q[2];
rz(1.5227854) q[3];
sx q[3];
rz(-0.77097547) q[3];
sx q[3];
rz(2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9252121) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(1.9981617) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(-2.3480031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0045354768) q[0];
sx q[0];
rz(-1.5001778) q[0];
sx q[0];
rz(0.15839346) q[0];
x q[1];
rz(-3.0751444) q[2];
sx q[2];
rz(-2.1796103) q[2];
sx q[2];
rz(-0.43283909) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55793437) q[1];
sx q[1];
rz(-1.0685295) q[1];
sx q[1];
rz(-1.7947547) q[1];
rz(-pi) q[2];
rz(1.3511485) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.3212412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6952343) q[2];
sx q[2];
rz(-1.0506722) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(-0.71600437) q[3];
sx q[3];
rz(-2.1456238) q[3];
sx q[3];
rz(2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-0.85549131) q[0];
sx q[0];
rz(0.18095428) q[0];
rz(-1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-2.1612371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0553833) q[0];
sx q[0];
rz(-0.50859857) q[0];
sx q[0];
rz(0.17091708) q[0];
x q[1];
rz(1.1372526) q[2];
sx q[2];
rz(-1.989365) q[2];
sx q[2];
rz(-2.7035382) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36114039) q[1];
sx q[1];
rz(-1.7516416) q[1];
sx q[1];
rz(0.7493345) q[1];
rz(-0.47360955) q[3];
sx q[3];
rz(-0.43821143) q[3];
sx q[3];
rz(-1.1058863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.02562) q[2];
sx q[2];
rz(-2.9774234) q[2];
sx q[2];
rz(-1.3244965) q[2];
rz(-0.65274158) q[3];
sx q[3];
rz(-0.96941152) q[3];
sx q[3];
rz(0.038221922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7272335) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(0.090567055) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(-2.3694585) q[3];
sx q[3];
rz(-1.6498389) q[3];
sx q[3];
rz(-2.6966358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
