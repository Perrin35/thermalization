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
rz(-0.72379273) q[0];
sx q[0];
rz(-2.0279217) q[0];
sx q[0];
rz(0.96473515) q[0];
rz(3.0678897) q[1];
sx q[1];
rz(-2.5001723) q[1];
sx q[1];
rz(1.4944271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49046966) q[0];
sx q[0];
rz(-2.2283366) q[0];
sx q[0];
rz(-1.0762908) q[0];
rz(-pi) q[1];
rz(-1.2531075) q[2];
sx q[2];
rz(-1.6699546) q[2];
sx q[2];
rz(2.3929208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19077483) q[1];
sx q[1];
rz(-0.72682805) q[1];
sx q[1];
rz(0.72587691) q[1];
x q[2];
rz(1.4919287) q[3];
sx q[3];
rz(-0.53073629) q[3];
sx q[3];
rz(2.4473878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73704314) q[2];
sx q[2];
rz(-1.6208384) q[2];
sx q[2];
rz(-0.21657319) q[2];
rz(0.51566044) q[3];
sx q[3];
rz(-1.0786846) q[3];
sx q[3];
rz(-3.0576341) q[3];
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
rz(-0.11230042) q[0];
sx q[0];
rz(-0.071542112) q[0];
sx q[0];
rz(-1.6610425) q[0];
rz(-2.5484565) q[1];
sx q[1];
rz(-2.1430404) q[1];
sx q[1];
rz(2.8592529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4078699) q[0];
sx q[0];
rz(-1.5998915) q[0];
sx q[0];
rz(-0.12761527) q[0];
x q[1];
rz(3.00782) q[2];
sx q[2];
rz(-1.3743128) q[2];
sx q[2];
rz(1.5143732) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39277252) q[1];
sx q[1];
rz(-1.7085679) q[1];
sx q[1];
rz(-2.5974136) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2766791) q[3];
sx q[3];
rz(-2.5967715) q[3];
sx q[3];
rz(-0.13710216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7181519) q[2];
sx q[2];
rz(-2.580692) q[2];
sx q[2];
rz(-2.1924428) q[2];
rz(-0.078350457) q[3];
sx q[3];
rz(-1.6716985) q[3];
sx q[3];
rz(1.9115492) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80701989) q[0];
sx q[0];
rz(-0.0581352) q[0];
sx q[0];
rz(2.9491501) q[0];
rz(-0.9737393) q[1];
sx q[1];
rz(-2.1111919) q[1];
sx q[1];
rz(-1.2724426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9285043) q[0];
sx q[0];
rz(-0.72505378) q[0];
sx q[0];
rz(0.37037767) q[0];
x q[1];
rz(0.22120938) q[2];
sx q[2];
rz(-0.84748805) q[2];
sx q[2];
rz(2.4257367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99062186) q[1];
sx q[1];
rz(-0.34268296) q[1];
sx q[1];
rz(-2.3938378) q[1];
rz(1.4951823) q[3];
sx q[3];
rz(-2.6845884) q[3];
sx q[3];
rz(1.5357032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56866139) q[2];
sx q[2];
rz(-0.91721407) q[2];
sx q[2];
rz(2.5407963) q[2];
rz(-2.3066547) q[3];
sx q[3];
rz(-0.91885126) q[3];
sx q[3];
rz(-2.6981573) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3804974) q[0];
sx q[0];
rz(-1.1325855) q[0];
sx q[0];
rz(-1.7219211) q[0];
rz(-0.28600606) q[1];
sx q[1];
rz(-2.0068469) q[1];
sx q[1];
rz(-2.6223415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7560067) q[0];
sx q[0];
rz(-2.950188) q[0];
sx q[0];
rz(0.49218224) q[0];
rz(-2.6250766) q[2];
sx q[2];
rz(-0.20843796) q[2];
sx q[2];
rz(2.1303653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1545461) q[1];
sx q[1];
rz(-1.0602762) q[1];
sx q[1];
rz(-1.6493919) q[1];
x q[2];
rz(-0.54597609) q[3];
sx q[3];
rz(-1.5432976) q[3];
sx q[3];
rz(-3.0101903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1661735) q[2];
sx q[2];
rz(-0.68048802) q[2];
sx q[2];
rz(0.55230459) q[2];
rz(0.78260261) q[3];
sx q[3];
rz(-1.3026594) q[3];
sx q[3];
rz(0.18463126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85823378) q[0];
sx q[0];
rz(-0.62223804) q[0];
sx q[0];
rz(1.2828113) q[0];
rz(1.917631) q[1];
sx q[1];
rz(-1.5354278) q[1];
sx q[1];
rz(-0.70972365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388572) q[0];
sx q[0];
rz(-1.5488911) q[0];
sx q[0];
rz(-0.028907348) q[0];
x q[1];
rz(-0.46825306) q[2];
sx q[2];
rz(-1.5522458) q[2];
sx q[2];
rz(-1.5520688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1695474) q[1];
sx q[1];
rz(-1.2575713) q[1];
sx q[1];
rz(-3.0671228) q[1];
rz(-pi) q[2];
rz(2.9622795) q[3];
sx q[3];
rz(-1.9783522) q[3];
sx q[3];
rz(2.5236764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84511406) q[2];
sx q[2];
rz(-2.5402386) q[2];
sx q[2];
rz(-0.55310407) q[2];
rz(2.3002355) q[3];
sx q[3];
rz(-2.0935121) q[3];
sx q[3];
rz(-2.0060189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3784921) q[0];
sx q[0];
rz(-2.0846413) q[0];
sx q[0];
rz(-0.77504778) q[0];
rz(-1.6686324) q[1];
sx q[1];
rz(-2.5170363) q[1];
sx q[1];
rz(-0.83736173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178029) q[0];
sx q[0];
rz(-2.508381) q[0];
sx q[0];
rz(1.7342315) q[0];
rz(-pi) q[1];
rz(0.35995324) q[2];
sx q[2];
rz(-0.48243603) q[2];
sx q[2];
rz(-3.0878518) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36951276) q[1];
sx q[1];
rz(-1.1489778) q[1];
sx q[1];
rz(-2.0095411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9894137) q[3];
sx q[3];
rz(-1.8004724) q[3];
sx q[3];
rz(0.21462378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2396635) q[2];
sx q[2];
rz(-1.9749125) q[2];
sx q[2];
rz(0.76266328) q[2];
rz(0.20721063) q[3];
sx q[3];
rz(-1.2615633) q[3];
sx q[3];
rz(-2.5337849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36829456) q[0];
sx q[0];
rz(-2.2070856) q[0];
sx q[0];
rz(-1.6873129) q[0];
rz(-1.2794718) q[1];
sx q[1];
rz(-2.2984633) q[1];
sx q[1];
rz(1.4131193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9521546) q[0];
sx q[0];
rz(-1.3690152) q[0];
sx q[0];
rz(-0.39315572) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3033875) q[2];
sx q[2];
rz(-2.0186485) q[2];
sx q[2];
rz(0.034016646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.6225902) q[1];
sx q[1];
rz(-1.1347924) q[1];
sx q[1];
rz(1.9080129) q[1];
rz(2.787049) q[3];
sx q[3];
rz(-1.156329) q[3];
sx q[3];
rz(-2.8855151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9965245) q[2];
sx q[2];
rz(-0.86864305) q[2];
sx q[2];
rz(-0.16150148) q[2];
rz(0.7343556) q[3];
sx q[3];
rz(-2.2752094) q[3];
sx q[3];
rz(1.8374779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2187918) q[0];
sx q[0];
rz(-2.6113593) q[0];
sx q[0];
rz(-2.9042322) q[0];
rz(2.3573719) q[1];
sx q[1];
rz(-1.3619245) q[1];
sx q[1];
rz(1.721419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99731555) q[0];
sx q[0];
rz(-2.0610556) q[0];
sx q[0];
rz(1.9107633) q[0];
x q[1];
rz(-0.88241045) q[2];
sx q[2];
rz(-1.7441445) q[2];
sx q[2];
rz(1.8515996) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7569468) q[1];
sx q[1];
rz(-1.8638041) q[1];
sx q[1];
rz(-2.5811367) q[1];
rz(-pi) q[2];
rz(-0.79717199) q[3];
sx q[3];
rz(-0.67250055) q[3];
sx q[3];
rz(-0.029880015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64804849) q[2];
sx q[2];
rz(-1.250114) q[2];
sx q[2];
rz(-1.2250712) q[2];
rz(1.8152292) q[3];
sx q[3];
rz(-2.1882961) q[3];
sx q[3];
rz(-1.1939322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5305283) q[0];
sx q[0];
rz(-0.079212991) q[0];
sx q[0];
rz(-2.5031669) q[0];
rz(0.05052677) q[1];
sx q[1];
rz(-1.2082929) q[1];
sx q[1];
rz(2.5343177) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0078871) q[0];
sx q[0];
rz(-1.2528208) q[0];
sx q[0];
rz(1.1668276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88596099) q[2];
sx q[2];
rz(-1.6185624) q[2];
sx q[2];
rz(-1.1807962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.133195) q[1];
sx q[1];
rz(-1.6411157) q[1];
sx q[1];
rz(1.2446872) q[1];
rz(-pi) q[2];
rz(-1.3783941) q[3];
sx q[3];
rz(-1.7524613) q[3];
sx q[3];
rz(-2.3013129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4813469) q[2];
sx q[2];
rz(-1.9504184) q[2];
sx q[2];
rz(2.353239) q[2];
rz(-2.3039019) q[3];
sx q[3];
rz(-2.3205784) q[3];
sx q[3];
rz(-1.2825509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7343219) q[0];
sx q[0];
rz(-0.73198524) q[0];
sx q[0];
rz(-2.6328073) q[0];
rz(1.5599627) q[1];
sx q[1];
rz(-2.1211233) q[1];
sx q[1];
rz(-0.4090974) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1412774) q[0];
sx q[0];
rz(-2.8594236) q[0];
sx q[0];
rz(-0.98833584) q[0];
rz(2.6652135) q[2];
sx q[2];
rz(-1.1217818) q[2];
sx q[2];
rz(0.19461122) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.075114014) q[1];
sx q[1];
rz(-1.7914322) q[1];
sx q[1];
rz(1.5815939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1942098) q[3];
sx q[3];
rz(-1.9344182) q[3];
sx q[3];
rz(-1.0417494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6009377) q[2];
sx q[2];
rz(-1.6161852) q[2];
sx q[2];
rz(-1.8589004) q[2];
rz(3.1278074) q[3];
sx q[3];
rz(-1.1252334) q[3];
sx q[3];
rz(-1.3330601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8448821) q[0];
sx q[0];
rz(-1.7322576) q[0];
sx q[0];
rz(-2.5175293) q[0];
rz(-2.2804672) q[1];
sx q[1];
rz(-0.54665165) q[1];
sx q[1];
rz(0.13269592) q[1];
rz(1.672198) q[2];
sx q[2];
rz(-1.8894926) q[2];
sx q[2];
rz(2.8429902) q[2];
rz(-0.80879296) q[3];
sx q[3];
rz(-0.8142796) q[3];
sx q[3];
rz(0.7249226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
