OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6141619) q[0];
sx q[0];
rz(-0.80978137) q[0];
sx q[0];
rz(2.6101987) q[0];
rz(-2.9040789) q[1];
sx q[1];
rz(-1.3637435) q[1];
sx q[1];
rz(-1.9385424) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703585) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(-1.7869851) q[0];
rz(-1.7877903) q[2];
sx q[2];
rz(-1.8197682) q[2];
sx q[2];
rz(-0.72460246) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0002999) q[1];
sx q[1];
rz(-1.3360026) q[1];
sx q[1];
rz(-1.6293807) q[1];
rz(-pi) q[2];
rz(-2.8486512) q[3];
sx q[3];
rz(-1.6062692) q[3];
sx q[3];
rz(2.5778511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.6137326) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(-2.3221827) q[0];
rz(0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0430543) q[0];
sx q[0];
rz(-0.96999723) q[0];
sx q[0];
rz(-2.0354802) q[0];
rz(-pi) q[1];
rz(-1.9752713) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(-1.5420367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5310865) q[1];
sx q[1];
rz(-2.5889581) q[1];
sx q[1];
rz(0.83341877) q[1];
x q[2];
rz(1.1191788) q[3];
sx q[3];
rz(-1.041881) q[3];
sx q[3];
rz(0.12714735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1566029) q[2];
sx q[2];
rz(-1.9282324) q[2];
sx q[2];
rz(1.9796237) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6067628) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(1.3820232) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(0.64750013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8436463) q[0];
sx q[0];
rz(-2.2046411) q[0];
sx q[0];
rz(1.514939) q[0];
rz(-2.4261549) q[2];
sx q[2];
rz(-1.1416661) q[2];
sx q[2];
rz(-1.5769584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69933575) q[1];
sx q[1];
rz(-1.934821) q[1];
sx q[1];
rz(2.060021) q[1];
x q[2];
rz(-0.89214274) q[3];
sx q[3];
rz(-2.3193079) q[3];
sx q[3];
rz(-3.126614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2993762) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(2.1598699) q[3];
sx q[3];
rz(-1.2795307) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340387) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(2.2696944) q[0];
rz(-2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(-0.29104582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440118) q[0];
sx q[0];
rz(-0.62950069) q[0];
sx q[0];
rz(-1.9660216) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.307345) q[2];
sx q[2];
rz(-0.4220037) q[2];
sx q[2];
rz(-0.01854245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0709666) q[1];
sx q[1];
rz(-0.91849209) q[1];
sx q[1];
rz(0.91826203) q[1];
rz(-pi) q[2];
rz(1.8157186) q[3];
sx q[3];
rz(-2.5513253) q[3];
sx q[3];
rz(-0.11412379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3691833) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-2.0777168) q[3];
sx q[3];
rz(-0.86597401) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0020224) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(-2.5158665) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-2.5114139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5298313) q[0];
sx q[0];
rz(-2.1453619) q[0];
sx q[0];
rz(1.4252421) q[0];
rz(-pi) q[1];
rz(-2.3216256) q[2];
sx q[2];
rz(-1.8659288) q[2];
sx q[2];
rz(0.68857869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4634134) q[1];
sx q[1];
rz(-2.5555875) q[1];
sx q[1];
rz(1.9834118) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.583902) q[3];
sx q[3];
rz(-1.4413177) q[3];
sx q[3];
rz(1.6491216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(1.2498614) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(-2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(0.71682799) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(-2.3278918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1610634) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(1.784523) q[0];
rz(0.72197638) q[2];
sx q[2];
rz(-1.1140031) q[2];
sx q[2];
rz(-1.7320088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5972283) q[1];
sx q[1];
rz(-1.9042943) q[1];
sx q[1];
rz(0.79798214) q[1];
rz(-0.30721174) q[3];
sx q[3];
rz(-1.9091515) q[3];
sx q[3];
rz(-0.070778155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-0.6860835) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(0.66326052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30090573) q[0];
sx q[0];
rz(-1.5947043) q[0];
sx q[0];
rz(-2.2211214) q[0];
rz(-pi) q[1];
rz(2.5040607) q[2];
sx q[2];
rz(-1.3946643) q[2];
sx q[2];
rz(0.1917563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85147334) q[1];
sx q[1];
rz(-2.0514601) q[1];
sx q[1];
rz(-0.22766797) q[1];
rz(-pi) q[2];
rz(2.0459941) q[3];
sx q[3];
rz(-1.0605269) q[3];
sx q[3];
rz(-2.1197532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(3.1265124) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052977234) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(-2.836851) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(-2.6838141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25154197) q[0];
sx q[0];
rz(-2.4180782) q[0];
sx q[0];
rz(-1.0735805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0395983) q[2];
sx q[2];
rz(-2.2441412) q[2];
sx q[2];
rz(1.4199867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0486054) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(-0.98577568) q[1];
x q[2];
rz(2.5896766) q[3];
sx q[3];
rz(-1.8822) q[3];
sx q[3];
rz(-0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(1.7101074) q[2];
rz(1.7163904) q[3];
sx q[3];
rz(-2.5148354) q[3];
sx q[3];
rz(1.8553998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15437056) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(0.95348683) q[0];
rz(-1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(-0.58473933) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23823243) q[0];
sx q[0];
rz(-1.0833217) q[0];
sx q[0];
rz(1.4438629) q[0];
rz(1.936475) q[2];
sx q[2];
rz(-2.0846016) q[2];
sx q[2];
rz(2.6363381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19496275) q[1];
sx q[1];
rz(-1.1319036) q[1];
sx q[1];
rz(1.6473325) q[1];
rz(-1.0869914) q[3];
sx q[3];
rz(-2.8414747) q[3];
sx q[3];
rz(-0.3008315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(0.27627036) q[2];
rz(-2.590498) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(-0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(2.4861091) q[0];
rz(-1.9845225) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(-1.6190593) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(1.9615016) q[0];
rz(-pi) q[1];
rz(0.68799512) q[2];
sx q[2];
rz(-1.3032459) q[2];
sx q[2];
rz(-2.5431395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9177502) q[1];
sx q[1];
rz(-0.94864861) q[1];
sx q[1];
rz(-1.7734581) q[1];
x q[2];
rz(1.1679681) q[3];
sx q[3];
rz(-1.1790923) q[3];
sx q[3];
rz(-0.21564461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(2.299451) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-3.1150637) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(-1.7799829) q[2];
sx q[2];
rz(-1.5447415) q[2];
sx q[2];
rz(0.011199657) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
