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
rz(-0.37201878) q[0];
sx q[0];
rz(3.4932669) q[0];
sx q[0];
rz(9.3695661) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(3.0076495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75823821) q[0];
sx q[0];
rz(-1.3516434) q[0];
sx q[0];
rz(1.4683649) q[0];
rz(-pi) q[1];
rz(-1.3710466) q[2];
sx q[2];
rz(-2.1416975) q[2];
sx q[2];
rz(2.8303245) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4856) q[1];
sx q[1];
rz(-2.6645086) q[1];
sx q[1];
rz(2.361195) q[1];
rz(-pi) q[2];
rz(-1.6381959) q[3];
sx q[3];
rz(-0.52636468) q[3];
sx q[3];
rz(2.0215066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.967531) q[2];
sx q[2];
rz(-1.2207737) q[2];
sx q[2];
rz(-2.3270712) q[2];
rz(2.9461765) q[3];
sx q[3];
rz(-0.82868367) q[3];
sx q[3];
rz(1.0403847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8622417) q[0];
sx q[0];
rz(-2.8793654) q[0];
sx q[0];
rz(2.1694515) q[0];
rz(0.60130087) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(1.7960637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20939553) q[0];
sx q[0];
rz(-0.69898134) q[0];
sx q[0];
rz(2.2369034) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53671543) q[2];
sx q[2];
rz(-1.3379352) q[2];
sx q[2];
rz(0.95271207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0474009) q[1];
sx q[1];
rz(-1.7128452) q[1];
sx q[1];
rz(-1.1588276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18109326) q[3];
sx q[3];
rz(-1.5646184) q[3];
sx q[3];
rz(0.97038499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95416516) q[2];
sx q[2];
rz(-1.8623872) q[2];
sx q[2];
rz(2.1698451) q[2];
rz(2.1238972) q[3];
sx q[3];
rz(-1.1757937) q[3];
sx q[3];
rz(3.0903604) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63075066) q[0];
sx q[0];
rz(-0.96931163) q[0];
sx q[0];
rz(-0.72072679) q[0];
rz(-0.12598704) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(-0.40649498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2590721) q[0];
sx q[0];
rz(-1.2693624) q[0];
sx q[0];
rz(0.92275019) q[0];
rz(-0.7438306) q[2];
sx q[2];
rz(-1.1421428) q[2];
sx q[2];
rz(0.33996087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7123901) q[1];
sx q[1];
rz(-0.99502975) q[1];
sx q[1];
rz(-2.5151252) q[1];
rz(1.93531) q[3];
sx q[3];
rz(-1.3008833) q[3];
sx q[3];
rz(-0.70333896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6005818) q[2];
sx q[2];
rz(-1.7966248) q[2];
sx q[2];
rz(2.7715032) q[2];
rz(-3.0626512) q[3];
sx q[3];
rz(-0.97349662) q[3];
sx q[3];
rz(2.6551042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.9321891) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(-0.47602794) q[0];
rz(1.1141106) q[1];
sx q[1];
rz(-2.1962491) q[1];
sx q[1];
rz(-1.6260737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0335122) q[0];
sx q[0];
rz(-0.42244222) q[0];
sx q[0];
rz(-2.6080934) q[0];
x q[1];
rz(2.1944207) q[2];
sx q[2];
rz(-2.2786254) q[2];
sx q[2];
rz(0.90367452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9263401) q[1];
sx q[1];
rz(-0.36141342) q[1];
sx q[1];
rz(0.24140668) q[1];
rz(-pi) q[2];
rz(-1.2021121) q[3];
sx q[3];
rz(-0.75206176) q[3];
sx q[3];
rz(0.78898174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6473306) q[2];
sx q[2];
rz(-1.5212955) q[2];
sx q[2];
rz(0.94804135) q[2];
rz(0.4387795) q[3];
sx q[3];
rz(-0.56883562) q[3];
sx q[3];
rz(-2.2426864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(2.7137252) q[0];
rz(0.95871344) q[1];
sx q[1];
rz(-1.2370647) q[1];
sx q[1];
rz(2.0895035) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7822582) q[0];
sx q[0];
rz(-2.3365031) q[0];
sx q[0];
rz(-0.87053086) q[0];
x q[1];
rz(-2.1688227) q[2];
sx q[2];
rz(-2.5199119) q[2];
sx q[2];
rz(1.6046815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0316443) q[1];
sx q[1];
rz(-0.5809902) q[1];
sx q[1];
rz(-1.7419001) q[1];
rz(-pi) q[2];
rz(-2.0922727) q[3];
sx q[3];
rz(-2.0745894) q[3];
sx q[3];
rz(-1.4544482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.170257) q[2];
sx q[2];
rz(-1.8361788) q[2];
sx q[2];
rz(-2.7830284) q[2];
rz(-1.1809008) q[3];
sx q[3];
rz(-0.87555331) q[3];
sx q[3];
rz(-0.53926474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81826687) q[0];
sx q[0];
rz(-1.2610672) q[0];
sx q[0];
rz(2.4928424) q[0];
rz(-2.5541041) q[1];
sx q[1];
rz(-1.3222539) q[1];
sx q[1];
rz(-0.23744753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5295805) q[0];
sx q[0];
rz(-1.6027614) q[0];
sx q[0];
rz(0.096613823) q[0];
x q[1];
rz(2.7858829) q[2];
sx q[2];
rz(-0.28537073) q[2];
sx q[2];
rz(-1.732638) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0157369) q[1];
sx q[1];
rz(-2.3185711) q[1];
sx q[1];
rz(-1.5344844) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1319207) q[3];
sx q[3];
rz(-0.791852) q[3];
sx q[3];
rz(0.71969024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4805523) q[2];
sx q[2];
rz(-0.66142267) q[2];
sx q[2];
rz(-0.95575571) q[2];
rz(1.7620979) q[3];
sx q[3];
rz(-2.5199514) q[3];
sx q[3];
rz(2.9852338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5718004) q[0];
sx q[0];
rz(-0.0026230165) q[0];
sx q[0];
rz(3.0963335) q[0];
rz(2.3472002) q[1];
sx q[1];
rz(-0.88961273) q[1];
sx q[1];
rz(0.20283595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2204593) q[0];
sx q[0];
rz(-2.5510802) q[0];
sx q[0];
rz(-2.8648225) q[0];
rz(-pi) q[1];
rz(0.21059307) q[2];
sx q[2];
rz(-2.3989429) q[2];
sx q[2];
rz(-2.3499678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.07405532) q[1];
sx q[1];
rz(-2.592784) q[1];
sx q[1];
rz(-2.9850053) q[1];
rz(2.2928004) q[3];
sx q[3];
rz(-1.5986852) q[3];
sx q[3];
rz(1.8384118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63991919) q[2];
sx q[2];
rz(-1.7379652) q[2];
sx q[2];
rz(2.5301834) q[2];
rz(1.1008788) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(-0.27975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47073498) q[0];
sx q[0];
rz(-1.1962471) q[0];
sx q[0];
rz(-1.083495) q[0];
rz(0.96393839) q[1];
sx q[1];
rz(-2.0699392) q[1];
sx q[1];
rz(0.49577698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9993047) q[0];
sx q[0];
rz(-1.5919627) q[0];
sx q[0];
rz(-1.643996) q[0];
x q[1];
rz(1.8476358) q[2];
sx q[2];
rz(-1.3949035) q[2];
sx q[2];
rz(-1.7664282) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26357061) q[1];
sx q[1];
rz(-1.1568034) q[1];
sx q[1];
rz(2.1946226) q[1];
rz(-pi) q[2];
rz(-0.029997901) q[3];
sx q[3];
rz(-0.98326937) q[3];
sx q[3];
rz(-0.8647747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.701391) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(2.4526147) q[2];
rz(-1.5135328) q[3];
sx q[3];
rz(-0.83008927) q[3];
sx q[3];
rz(1.0015063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643519) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(2.054457) q[0];
rz(0.76308909) q[1];
sx q[1];
rz(-1.3975846) q[1];
sx q[1];
rz(2.9147002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.034982) q[0];
sx q[0];
rz(-0.7312432) q[0];
sx q[0];
rz(1.2319706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8925324) q[2];
sx q[2];
rz(-1.2870803) q[2];
sx q[2];
rz(2.4987881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1400661) q[1];
sx q[1];
rz(-2.9140009) q[1];
sx q[1];
rz(-2.9796991) q[1];
x q[2];
rz(-0.50548203) q[3];
sx q[3];
rz(-1.0311605) q[3];
sx q[3];
rz(-1.1330549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.521296) q[2];
sx q[2];
rz(-2.2038348) q[2];
sx q[2];
rz(-2.1916981) q[2];
rz(1.0319483) q[3];
sx q[3];
rz(-2.2622006) q[3];
sx q[3];
rz(0.43752813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49143219) q[0];
sx q[0];
rz(-3.037945) q[0];
sx q[0];
rz(1.1520804) q[0];
rz(3.0488455) q[1];
sx q[1];
rz(-0.68789613) q[1];
sx q[1];
rz(0.50382096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4277305) q[0];
sx q[0];
rz(-0.52885884) q[0];
sx q[0];
rz(1.669017) q[0];
x q[1];
rz(1.8015734) q[2];
sx q[2];
rz(-1.3084931) q[2];
sx q[2];
rz(-2.3364146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31695932) q[1];
sx q[1];
rz(-1.0976296) q[1];
sx q[1];
rz(1.6184357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9484048) q[3];
sx q[3];
rz(-1.8657627) q[3];
sx q[3];
rz(1.1136536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6759701) q[2];
sx q[2];
rz(-1.114782) q[2];
sx q[2];
rz(0.62391227) q[2];
rz(2.5850249) q[3];
sx q[3];
rz(-0.68816319) q[3];
sx q[3];
rz(-2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9482166) q[0];
sx q[0];
rz(-0.93127903) q[0];
sx q[0];
rz(-2.2122526) q[0];
rz(-2.9931862) q[1];
sx q[1];
rz(-1.2444617) q[1];
sx q[1];
rz(-0.1958227) q[1];
rz(2.6779867) q[2];
sx q[2];
rz(-1.4624034) q[2];
sx q[2];
rz(0.69093888) q[2];
rz(-1.2986167) q[3];
sx q[3];
rz(-1.4930354) q[3];
sx q[3];
rz(-2.2477575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
