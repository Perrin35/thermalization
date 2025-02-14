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
rz(2.4438357) q[0];
sx q[0];
rz(4.3351686) q[0];
sx q[0];
rz(9.6293443) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(-1.3890356) q[1];
sx q[1];
rz(1.3995481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7273063) q[0];
sx q[0];
rz(-2.901863) q[0];
sx q[0];
rz(1.9637462) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84487652) q[2];
sx q[2];
rz(-2.5683218) q[2];
sx q[2];
rz(2.5115761) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5198361) q[1];
sx q[1];
rz(-2.7108828) q[1];
sx q[1];
rz(-0.21060305) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70943009) q[3];
sx q[3];
rz(-1.7755847) q[3];
sx q[3];
rz(3.0258609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43510398) q[2];
sx q[2];
rz(-3.1092643) q[2];
sx q[2];
rz(-0.48933634) q[2];
rz(-2.1183744) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(1.1255012) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045687549) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(-0.80279654) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(-0.057770483) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36299713) q[0];
sx q[0];
rz(-2.5311806) q[0];
sx q[0];
rz(-1.1108062) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28654307) q[2];
sx q[2];
rz(-1.6725958) q[2];
sx q[2];
rz(-1.941178) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80177414) q[1];
sx q[1];
rz(-1.8027824) q[1];
sx q[1];
rz(1.07527) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3996734) q[3];
sx q[3];
rz(-2.4226396) q[3];
sx q[3];
rz(-2.3477682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(-1.2954953) q[2];
rz(-2.1748491) q[3];
sx q[3];
rz(-2.3714378) q[3];
sx q[3];
rz(-0.75743341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029723786) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(1.7080074) q[0];
rz(-0.068610527) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(2.5624018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1807208) q[0];
sx q[0];
rz(-1.5128883) q[0];
sx q[0];
rz(1.6625893) q[0];
x q[1];
rz(-3.1379329) q[2];
sx q[2];
rz(-0.47572593) q[2];
sx q[2];
rz(-1.9329247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3855558) q[1];
sx q[1];
rz(-1.7969404) q[1];
sx q[1];
rz(1.534034) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3287401) q[3];
sx q[3];
rz(-0.62004706) q[3];
sx q[3];
rz(2.4923785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27089831) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(-0.18343055) q[2];
rz(-0.80426788) q[3];
sx q[3];
rz(-2.1427514) q[3];
sx q[3];
rz(2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(-0.59952366) q[0];
rz(0.41246688) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-2.1309158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0507999) q[0];
sx q[0];
rz(-2.7677892) q[0];
sx q[0];
rz(2.7888377) q[0];
rz(-1.5013736) q[2];
sx q[2];
rz(-0.63328082) q[2];
sx q[2];
rz(0.4236003) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49941555) q[1];
sx q[1];
rz(-0.66153971) q[1];
sx q[1];
rz(1.1206349) q[1];
x q[2];
rz(1.7092199) q[3];
sx q[3];
rz(-2.3218621) q[3];
sx q[3];
rz(-0.26811308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3707054) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(2.5192449) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(2.5573825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7790826) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(-1.0269748) q[0];
rz(-0.55038553) q[1];
sx q[1];
rz(-0.06409476) q[1];
sx q[1];
rz(1.2170353) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37792045) q[0];
sx q[0];
rz(-0.91912133) q[0];
sx q[0];
rz(2.4125189) q[0];
rz(-1.6486859) q[2];
sx q[2];
rz(-1.0337892) q[2];
sx q[2];
rz(-2.3324049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0266307) q[1];
sx q[1];
rz(-1.7424045) q[1];
sx q[1];
rz(-2.2725355) q[1];
x q[2];
rz(-1.1547791) q[3];
sx q[3];
rz(-2.3868594) q[3];
sx q[3];
rz(0.55547914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43651849) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-2.3684033) q[2];
rz(0.1117205) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(2.2022061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6620827) q[0];
sx q[0];
rz(-2.918512) q[0];
sx q[0];
rz(2.7146085) q[0];
rz(0.93049479) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(2.6771136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1982959) q[0];
sx q[0];
rz(-1.7664599) q[0];
sx q[0];
rz(0.29384675) q[0];
rz(-pi) q[1];
rz(0.97067483) q[2];
sx q[2];
rz(-1.2874741) q[2];
sx q[2];
rz(-1.7194058) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.528094) q[1];
sx q[1];
rz(-1.6069176) q[1];
sx q[1];
rz(-2.7725459) q[1];
rz(-1.9950946) q[3];
sx q[3];
rz(-1.5992377) q[3];
sx q[3];
rz(0.44876307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6102607) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-2.8533234) q[2];
rz(-1.0432976) q[3];
sx q[3];
rz(-0.61560029) q[3];
sx q[3];
rz(0.73879761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814602) q[0];
sx q[0];
rz(-1.8539424) q[0];
sx q[0];
rz(0.75755358) q[0];
rz(3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(3.0933948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91201997) q[0];
sx q[0];
rz(-2.1549112) q[0];
sx q[0];
rz(-0.035365625) q[0];
rz(-pi) q[1];
rz(2.9075895) q[2];
sx q[2];
rz(-2.1562139) q[2];
sx q[2];
rz(-0.3096748) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99302247) q[1];
sx q[1];
rz(-2.3766368) q[1];
sx q[1];
rz(-0.92924849) q[1];
rz(-pi) q[2];
rz(-2.0406575) q[3];
sx q[3];
rz(-1.5239232) q[3];
sx q[3];
rz(2.8702503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4797719) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(3.0721967) q[2];
rz(-1.5875459) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(-2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7540392) q[0];
sx q[0];
rz(-1.0559005) q[0];
sx q[0];
rz(1.4132502) q[0];
rz(-0.82855254) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(0.54263306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0135745) q[0];
sx q[0];
rz(-1.4262222) q[0];
sx q[0];
rz(-1.4506863) q[0];
rz(-1.2130402) q[2];
sx q[2];
rz(-1.925549) q[2];
sx q[2];
rz(-2.8194129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60195758) q[1];
sx q[1];
rz(-1.4550721) q[1];
sx q[1];
rz(-1.6249815) q[1];
rz(2.7067634) q[3];
sx q[3];
rz(-1.2246338) q[3];
sx q[3];
rz(-1.0258254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1594306) q[2];
sx q[2];
rz(-0.40951481) q[2];
sx q[2];
rz(2.8816667) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-0.25769886) q[3];
sx q[3];
rz(2.3475588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(-1.6625241) q[0];
rz(1.6089449) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(0.74302465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4585658) q[0];
sx q[0];
rz(-1.4560149) q[0];
sx q[0];
rz(-2.0857138) q[0];
x q[1];
rz(-2.7588927) q[2];
sx q[2];
rz(-2.0819132) q[2];
sx q[2];
rz(1.9267043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0155269) q[1];
sx q[1];
rz(-3.0558476) q[1];
sx q[1];
rz(-0.87460204) q[1];
rz(-pi) q[2];
rz(0.10492341) q[3];
sx q[3];
rz(-1.5399974) q[3];
sx q[3];
rz(-1.2237751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4674025) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(-0.80545938) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(2.3384371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886803) q[0];
sx q[0];
rz(-2.5997933) q[0];
sx q[0];
rz(0.75206494) q[0];
rz(1.9883142) q[1];
sx q[1];
rz(-0.88429943) q[1];
sx q[1];
rz(2.8582252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3842938) q[0];
sx q[0];
rz(-1.6597676) q[0];
sx q[0];
rz(-2.6020537) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4128017) q[2];
sx q[2];
rz(-0.758095) q[2];
sx q[2];
rz(-1.5604863) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9787406) q[1];
sx q[1];
rz(-1.2417485) q[1];
sx q[1];
rz(0.63905893) q[1];
x q[2];
rz(-1.180319) q[3];
sx q[3];
rz(-0.9538528) q[3];
sx q[3];
rz(-2.4882567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2470384) q[2];
sx q[2];
rz(-3.0587695) q[2];
sx q[2];
rz(1.438633) q[2];
rz(-0.29397193) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1080078) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(0.53957466) q[1];
sx q[1];
rz(-0.78782606) q[1];
sx q[1];
rz(0.16001564) q[1];
rz(-1.1926959) q[2];
sx q[2];
rz(-0.38496502) q[2];
sx q[2];
rz(-2.5733583) q[2];
rz(-2.7593437) q[3];
sx q[3];
rz(-1.3180238) q[3];
sx q[3];
rz(-2.870306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
