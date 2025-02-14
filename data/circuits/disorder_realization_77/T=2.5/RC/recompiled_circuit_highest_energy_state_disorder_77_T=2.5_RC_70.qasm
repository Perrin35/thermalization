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
rz(-1.9714126) q[1];
sx q[1];
rz(-2.7538731) q[1];
sx q[1];
rz(2.5375836) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9785288) q[0];
sx q[0];
rz(-1.8730358) q[0];
sx q[0];
rz(3.0885484) q[0];
rz(-pi) q[1];
rz(-0.9004132) q[2];
sx q[2];
rz(-1.1657823) q[2];
sx q[2];
rz(2.9069898) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4166222) q[1];
sx q[1];
rz(-1.23471) q[1];
sx q[1];
rz(-2.1137456) q[1];
x q[2];
rz(1.2617926) q[3];
sx q[3];
rz(-1.6741535) q[3];
sx q[3];
rz(-0.37783937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.151256) q[2];
sx q[2];
rz(-1.1948816) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.6319298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5556521) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(-1.1974539) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(-1.7053568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89051188) q[0];
sx q[0];
rz(-2.5174826) q[0];
sx q[0];
rz(1.1680383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6701489) q[2];
sx q[2];
rz(-2.2713285) q[2];
sx q[2];
rz(1.6068899) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9402071) q[1];
sx q[1];
rz(-1.6734945) q[1];
sx q[1];
rz(-2.9443355) q[1];
x q[2];
rz(-1.7800434) q[3];
sx q[3];
rz(-1.1236785) q[3];
sx q[3];
rz(0.9613049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70812923) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(-2.1500812) q[0];
rz(-2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(-1.2352157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69312364) q[0];
sx q[0];
rz(-2.4938994) q[0];
sx q[0];
rz(0.41843398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0472544) q[2];
sx q[2];
rz(-2.1263577) q[2];
sx q[2];
rz(0.68651344) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0491525) q[1];
sx q[1];
rz(-1.3235281) q[1];
sx q[1];
rz(0.80353261) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54872437) q[3];
sx q[3];
rz(-1.0503891) q[3];
sx q[3];
rz(-1.6136374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11423763) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.39795136) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(-0.81746307) q[0];
rz(-0.3262597) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(-0.78525966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8924849) q[0];
sx q[0];
rz(-1.2689225) q[0];
sx q[0];
rz(2.8623146) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3738858) q[2];
sx q[2];
rz(-0.93683883) q[2];
sx q[2];
rz(-1.9206573) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98475115) q[1];
sx q[1];
rz(-1.1492711) q[1];
sx q[1];
rz(2.0578029) q[1];
rz(-pi) q[2];
rz(-2.7866632) q[3];
sx q[3];
rz(-1.4961492) q[3];
sx q[3];
rz(-1.9982893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(2.5784967) q[2];
rz(0.8935039) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(-1.9068498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017482) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(1.5455986) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(2.9603069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90732987) q[0];
sx q[0];
rz(-1.0815623) q[0];
sx q[0];
rz(0.2184043) q[0];
rz(0.092575707) q[2];
sx q[2];
rz(-2.6642054) q[2];
sx q[2];
rz(2.8344748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1039155) q[1];
sx q[1];
rz(-1.5988776) q[1];
sx q[1];
rz(-1.3485391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32205228) q[3];
sx q[3];
rz(-0.46370927) q[3];
sx q[3];
rz(-0.85953322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.01866092) q[2];
sx q[2];
rz(-2.4666726) q[2];
sx q[2];
rz(-1.8956511) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6084006) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(0.30884185) q[0];
rz(0.32288512) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(3.0016532) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2353781) q[0];
sx q[0];
rz(-1.3227292) q[0];
sx q[0];
rz(2.1885022) q[0];
x q[1];
rz(1.3792737) q[2];
sx q[2];
rz(-2.6783248) q[2];
sx q[2];
rz(-0.5199711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9819239) q[1];
sx q[1];
rz(-0.86729151) q[1];
sx q[1];
rz(2.7419006) q[1];
x q[2];
rz(-2.2298584) q[3];
sx q[3];
rz(-1.3305802) q[3];
sx q[3];
rz(-1.6068589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(0.24197401) q[2];
rz(-0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.411875) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(1.2710849) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(-1.44106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1313167) q[0];
sx q[0];
rz(-2.1621341) q[0];
sx q[0];
rz(2.813126) q[0];
x q[1];
rz(0.3344603) q[2];
sx q[2];
rz(-2.7401667) q[2];
sx q[2];
rz(-2.4106467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.064441) q[1];
sx q[1];
rz(-1.7930601) q[1];
sx q[1];
rz(2.7138316) q[1];
rz(-pi) q[2];
rz(-0.7251803) q[3];
sx q[3];
rz(-1.4075052) q[3];
sx q[3];
rz(2.0189328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21961221) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(1.316635) q[2];
rz(-2.5804139) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(-1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.104326) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(0.88687801) q[0];
rz(-2.9414224) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(2.1194469) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048269317) q[0];
sx q[0];
rz(-1.5962068) q[0];
sx q[0];
rz(1.4396458) q[0];
rz(-pi) q[1];
rz(0.80506246) q[2];
sx q[2];
rz(-2.2230122) q[2];
sx q[2];
rz(2.7307939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8238942) q[1];
sx q[1];
rz(-2.3168457) q[1];
sx q[1];
rz(-0.0024248799) q[1];
rz(2.5660146) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(0.23250599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6518121) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(2.32302) q[2];
rz(-2.1583648) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(-2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.1061358) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(1.5863093) q[0];
rz(-0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(0.78479016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1589543) q[0];
sx q[0];
rz(-2.1110015) q[0];
sx q[0];
rz(2.0870952) q[0];
rz(2.2124039) q[2];
sx q[2];
rz(-0.25114533) q[2];
sx q[2];
rz(0.26640688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1363298) q[1];
sx q[1];
rz(-1.6044093) q[1];
sx q[1];
rz(-0.43515794) q[1];
rz(-0.611245) q[3];
sx q[3];
rz(-1.2204303) q[3];
sx q[3];
rz(2.4823852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(-2.3308241) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588381) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(2.3639823) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(-0.75540677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574673) q[0];
sx q[0];
rz(-1.3608772) q[0];
sx q[0];
rz(-0.0046878417) q[0];
rz(-pi) q[1];
rz(1.5532812) q[2];
sx q[2];
rz(-1.1743769) q[2];
sx q[2];
rz(2.6902386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0499784) q[1];
sx q[1];
rz(-2.0777933) q[1];
sx q[1];
rz(1.9732287) q[1];
rz(-pi) q[2];
rz(2.741119) q[3];
sx q[3];
rz(-2.4254834) q[3];
sx q[3];
rz(-3.0587089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(-0.15288615) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818903) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(-1.9181171) q[1];
sx q[1];
rz(-1.9203095) q[1];
sx q[1];
rz(-0.65345678) q[1];
rz(0.6003004) q[2];
sx q[2];
rz(-2.4523198) q[2];
sx q[2];
rz(-0.92858989) q[2];
rz(-2.7889403) q[3];
sx q[3];
rz(-1.1897539) q[3];
sx q[3];
rz(-2.6170058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
