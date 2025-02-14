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
rz(-1.1779397) q[0];
sx q[0];
rz(-3.0484634) q[0];
sx q[0];
rz(2.4094474) q[0];
rz(1.572345) q[1];
sx q[1];
rz(0.88580004) q[1];
sx q[1];
rz(9.8635397) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33235655) q[0];
sx q[0];
rz(-1.9026103) q[0];
sx q[0];
rz(0.058666243) q[0];
rz(0.83005367) q[2];
sx q[2];
rz(-2.7657653) q[2];
sx q[2];
rz(-2.2245882) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4784412) q[1];
sx q[1];
rz(-1.8290797) q[1];
sx q[1];
rz(1.5228935) q[1];
x q[2];
rz(2.6832163) q[3];
sx q[3];
rz(-1.5273558) q[3];
sx q[3];
rz(2.9921852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0977122) q[2];
sx q[2];
rz(-3.1011797) q[2];
sx q[2];
rz(-2.2229693) q[2];
rz(2.0888603) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.87616462) q[0];
sx q[0];
rz(-0.81887236) q[0];
sx q[0];
rz(-0.87232605) q[0];
rz(2.1288952) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-2.8025119) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144086) q[0];
sx q[0];
rz(-2.0054584) q[0];
sx q[0];
rz(-1.7130769) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85810272) q[2];
sx q[2];
rz(-0.65026186) q[2];
sx q[2];
rz(1.6158783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69557938) q[1];
sx q[1];
rz(-1.6124311) q[1];
sx q[1];
rz(1.5783046) q[1];
x q[2];
rz(-1.2463739) q[3];
sx q[3];
rz(-1.4528465) q[3];
sx q[3];
rz(1.1583386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8104711) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(0.60532153) q[2];
rz(2.807054) q[3];
sx q[3];
rz(-0.3722705) q[3];
sx q[3];
rz(-0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89638585) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(1.4728004) q[0];
rz(-0.32707602) q[1];
sx q[1];
rz(-1.8388803) q[1];
sx q[1];
rz(-2.8715141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9327346) q[0];
sx q[0];
rz(-1.3737824) q[0];
sx q[0];
rz(0.52893649) q[0];
rz(-pi) q[1];
rz(1.4065341) q[2];
sx q[2];
rz(-2.7077419) q[2];
sx q[2];
rz(1.0264772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9880319) q[1];
sx q[1];
rz(-0.84756339) q[1];
sx q[1];
rz(2.5392697) q[1];
rz(-0.047932054) q[3];
sx q[3];
rz(-0.080259003) q[3];
sx q[3];
rz(1.3788392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.104287) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(-1.8541065) q[2];
rz(1.2152952) q[3];
sx q[3];
rz(-1.3853962) q[3];
sx q[3];
rz(2.9638885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6215068) q[0];
sx q[0];
rz(-2.6984213) q[0];
sx q[0];
rz(0.31975123) q[0];
rz(-2.7301835) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(-3.06126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7112996) q[0];
sx q[0];
rz(-1.3724494) q[0];
sx q[0];
rz(1.2567149) q[0];
rz(-pi) q[1];
rz(-0.80647237) q[2];
sx q[2];
rz(-1.6739168) q[2];
sx q[2];
rz(-1.1191302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44583933) q[1];
sx q[1];
rz(-1.8421122) q[1];
sx q[1];
rz(-0.81636565) q[1];
rz(-pi) q[2];
rz(-1.3036968) q[3];
sx q[3];
rz(-2.1648277) q[3];
sx q[3];
rz(-0.66583064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8722998) q[2];
sx q[2];
rz(-2.9041957) q[2];
sx q[2];
rz(3.0966975) q[2];
rz(2.6063555) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.991268) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(2.2659361) q[0];
rz(0.21370299) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(1.5692086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0914577) q[0];
sx q[0];
rz(-2.7612855) q[0];
sx q[0];
rz(2.1652392) q[0];
rz(-pi) q[1];
rz(-2.7046811) q[2];
sx q[2];
rz(-1.3624117) q[2];
sx q[2];
rz(2.2885099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8991291) q[1];
sx q[1];
rz(-1.4664259) q[1];
sx q[1];
rz(2.4422798) q[1];
rz(-pi) q[2];
rz(-0.27690378) q[3];
sx q[3];
rz(-2.8087862) q[3];
sx q[3];
rz(-1.3267508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18465061) q[2];
sx q[2];
rz(-0.1447548) q[2];
sx q[2];
rz(2.0733898) q[2];
rz(0.60643658) q[3];
sx q[3];
rz(-1.9304099) q[3];
sx q[3];
rz(-3.1312805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40474263) q[0];
sx q[0];
rz(-2.8917942) q[0];
sx q[0];
rz(-1.7506208) q[0];
rz(2.6178316) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(-1.1837122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9212657) q[0];
sx q[0];
rz(-1.8588716) q[0];
sx q[0];
rz(0.057805268) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1864088) q[2];
sx q[2];
rz(-0.092530017) q[2];
sx q[2];
rz(1.1902155) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6377012) q[1];
sx q[1];
rz(-2.5092485) q[1];
sx q[1];
rz(1.8296627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88504412) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(2.1751753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.70048731) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(2.109745) q[2];
rz(-1.4619689) q[3];
sx q[3];
rz(-0.69457355) q[3];
sx q[3];
rz(-1.3612548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66330355) q[0];
sx q[0];
rz(-2.7155868) q[0];
sx q[0];
rz(1.4248832) q[0];
rz(2.5583963) q[1];
sx q[1];
rz(-1.9164663) q[1];
sx q[1];
rz(0.057417631) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7518508) q[0];
sx q[0];
rz(-1.5673141) q[0];
sx q[0];
rz(-2.338242) q[0];
rz(-pi) q[1];
rz(0.67441794) q[2];
sx q[2];
rz(-2.4092374) q[2];
sx q[2];
rz(-2.6931817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8805931) q[1];
sx q[1];
rz(-0.94124244) q[1];
sx q[1];
rz(-1.3531319) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52773169) q[3];
sx q[3];
rz(-0.70042983) q[3];
sx q[3];
rz(2.79984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3164369) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(-1.3081029) q[2];
rz(1.3206652) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(0.86110419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84701455) q[0];
sx q[0];
rz(-2.4569643) q[0];
sx q[0];
rz(-1.8735877) q[0];
rz(1.1216724) q[1];
sx q[1];
rz(-1.8662165) q[1];
sx q[1];
rz(-2.4915288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92597678) q[0];
sx q[0];
rz(-2.4575666) q[0];
sx q[0];
rz(-0.90606545) q[0];
rz(-pi) q[1];
rz(0.68393884) q[2];
sx q[2];
rz(-1.9994447) q[2];
sx q[2];
rz(0.49698439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2686058) q[1];
sx q[1];
rz(-1.0293616) q[1];
sx q[1];
rz(-1.1100015) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88553377) q[3];
sx q[3];
rz(-0.69708744) q[3];
sx q[3];
rz(3.0892885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0012297) q[2];
sx q[2];
rz(-1.659212) q[2];
sx q[2];
rz(2.7492827) q[2];
rz(2.9366734) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1392764) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(-0.32683364) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(-1.0439509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4503339) q[0];
sx q[0];
rz(-0.46713167) q[0];
sx q[0];
rz(-1.1410261) q[0];
rz(-pi) q[1];
rz(-2.6828917) q[2];
sx q[2];
rz(-1.6965116) q[2];
sx q[2];
rz(1.1864551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3957246) q[1];
sx q[1];
rz(-1.9650794) q[1];
sx q[1];
rz(0.08862011) q[1];
rz(-pi) q[2];
rz(1.9887505) q[3];
sx q[3];
rz(-2.4476789) q[3];
sx q[3];
rz(-1.0479463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.30847183) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(0.92588818) q[2];
rz(2.0738257) q[3];
sx q[3];
rz(-2.0177149) q[3];
sx q[3];
rz(-1.7759751) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85497102) q[0];
sx q[0];
rz(-1.6678896) q[0];
sx q[0];
rz(-2.5545252) q[0];
rz(2.558737) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(-2.6606182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0635446) q[0];
sx q[0];
rz(-0.46351156) q[0];
sx q[0];
rz(0.48234756) q[0];
rz(-pi) q[1];
rz(1.9312385) q[2];
sx q[2];
rz(-1.7490462) q[2];
sx q[2];
rz(2.0958063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1744057) q[1];
sx q[1];
rz(-2.4250826) q[1];
sx q[1];
rz(-2.4661676) q[1];
rz(-pi) q[2];
x q[2];
rz(1.770788) q[3];
sx q[3];
rz(-1.963435) q[3];
sx q[3];
rz(-0.61905608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19758548) q[2];
sx q[2];
rz(-0.53637594) q[2];
sx q[2];
rz(-1.8120922) q[2];
rz(-2.3594989) q[3];
sx q[3];
rz(-2.8160281) q[3];
sx q[3];
rz(-1.9935002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.44163497) q[0];
sx q[0];
rz(-1.7229719) q[0];
sx q[0];
rz(-1.4650719) q[0];
rz(1.5326473) q[1];
sx q[1];
rz(-1.9736704) q[1];
sx q[1];
rz(-0.71221487) q[1];
rz(-2.4423626) q[2];
sx q[2];
rz(-2.2742827) q[2];
sx q[2];
rz(-0.90894528) q[2];
rz(1.731012) q[3];
sx q[3];
rz(-1.0777149) q[3];
sx q[3];
rz(3.1330681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
