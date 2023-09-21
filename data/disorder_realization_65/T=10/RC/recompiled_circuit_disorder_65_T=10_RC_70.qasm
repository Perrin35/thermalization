OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(2.0895045) q[0];
sx q[0];
rz(10.917561) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5045835) q[0];
sx q[0];
rz(-1.5890117) q[0];
sx q[0];
rz(-1.5655184) q[0];
rz(1.7311814) q[2];
sx q[2];
rz(-2.3540902) q[2];
sx q[2];
rz(3.0410142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97779146) q[1];
sx q[1];
rz(-1.4492387) q[1];
sx q[1];
rz(1.4896643) q[1];
rz(0.53332897) q[3];
sx q[3];
rz(-1.1733857) q[3];
sx q[3];
rz(2.1189789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.5011903) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879724) q[0];
sx q[0];
rz(-3.0164218) q[0];
sx q[0];
rz(0.20010389) q[0];
rz(1.1019727) q[2];
sx q[2];
rz(-1.2906244) q[2];
sx q[2];
rz(2.3902992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.048763927) q[1];
sx q[1];
rz(-1.3991014) q[1];
sx q[1];
rz(2.462639) q[1];
x q[2];
rz(1.5562015) q[3];
sx q[3];
rz(-0.39108927) q[3];
sx q[3];
rz(-2.1434243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(0.93018053) q[0];
rz(2.1554699) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(0.090678064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40239247) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(-1.5261569) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4012666) q[3];
sx q[3];
rz(-2.4781514) q[3];
sx q[3];
rz(1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79608941) q[0];
sx q[0];
rz(-2.8755099) q[0];
sx q[0];
rz(-0.43568133) q[0];
x q[1];
rz(0.97781424) q[2];
sx q[2];
rz(-1.9359509) q[2];
sx q[2];
rz(-2.1821373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9149575) q[1];
sx q[1];
rz(-1.8640222) q[1];
sx q[1];
rz(-2.3073767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3991688) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(-2.1267668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4887061) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(2.2573684) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852172) q[0];
sx q[0];
rz(-1.3803122) q[0];
sx q[0];
rz(2.9604244) q[0];
rz(-pi) q[1];
rz(-0.40712936) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(-1.4332353) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3430749) q[1];
sx q[1];
rz(-1.6371562) q[1];
sx q[1];
rz(-2.9970478) q[1];
rz(-pi) q[2];
rz(2.1518425) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(2.1637672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(-2.3905335) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9861525) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-0.040779671) q[0];
x q[1];
rz(-1.6386119) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(-0.80578795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25989306) q[1];
sx q[1];
rz(-1.1105781) q[1];
sx q[1];
rz(-2.0395181) q[1];
x q[2];
rz(-1.045741) q[3];
sx q[3];
rz(-1.9805555) q[3];
sx q[3];
rz(-2.4214782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(1.4356027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5933701) q[0];
sx q[0];
rz(-0.75836327) q[0];
sx q[0];
rz(-0.81292696) q[0];
x q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(-0.4617304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3064733) q[1];
sx q[1];
rz(-1.2700348) q[1];
sx q[1];
rz(1.1919828) q[1];
rz(-pi) q[2];
rz(-0.88019754) q[3];
sx q[3];
rz(-0.62303632) q[3];
sx q[3];
rz(-2.8043889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(-1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85132861) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.3652323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(-1.4343065) q[0];
x q[1];
rz(1.6985967) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(-2.4335361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9284968) q[1];
sx q[1];
rz(-1.0398541) q[1];
sx q[1];
rz(-1.5099768) q[1];
x q[2];
rz(0.39799277) q[3];
sx q[3];
rz(-2.6905305) q[3];
sx q[3];
rz(0.034702452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(-0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7565219) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(-2.4247652) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5588435) q[2];
sx q[2];
rz(-2.0593615) q[2];
sx q[2];
rz(1.919463) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82211557) q[1];
sx q[1];
rz(-1.9648106) q[1];
sx q[1];
rz(2.134321) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68848916) q[3];
sx q[3];
rz(-0.8430891) q[3];
sx q[3];
rz(-1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.753153) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5441226) q[0];
sx q[0];
rz(2.1795576) q[0];
rz(-pi) q[1];
rz(1.1881234) q[2];
sx q[2];
rz(-1.6560935) q[2];
sx q[2];
rz(-1.6558937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2529777) q[1];
sx q[1];
rz(-1.1974317) q[1];
sx q[1];
rz(0.044815973) q[1];
rz(-pi) q[2];
rz(0.90922728) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(2.6565068) q[2];
sx q[2];
rz(-0.1242287) q[2];
sx q[2];
rz(2.6835174) q[2];
rz(0.86919541) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];