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
rz(-2.383411) q[0];
sx q[0];
rz(-0.36769205) q[0];
sx q[0];
rz(2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(1.1059603) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6152412) q[0];
sx q[0];
rz(-0.59068524) q[0];
sx q[0];
rz(2.3870941) q[0];
rz(0.26768406) q[2];
sx q[2];
rz(-1.6639198) q[2];
sx q[2];
rz(0.47506079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3788521) q[1];
sx q[1];
rz(-1.6081075) q[1];
sx q[1];
rz(-0.0040405063) q[1];
rz(-pi) q[2];
rz(-2.8572122) q[3];
sx q[3];
rz(-1.577498) q[3];
sx q[3];
rz(2.0363765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7652863) q[2];
sx q[2];
rz(-2.6708965) q[2];
sx q[2];
rz(2.7619696) q[2];
rz(1.5073353) q[3];
sx q[3];
rz(-1.0999271) q[3];
sx q[3];
rz(-2.1716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.64604243) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(-0.90091339) q[0];
rz(0.13149978) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(-0.44201717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64009066) q[0];
sx q[0];
rz(-1.2462062) q[0];
sx q[0];
rz(1.4844271) q[0];
rz(-pi) q[1];
rz(-3.1369741) q[2];
sx q[2];
rz(-1.5698264) q[2];
sx q[2];
rz(-2.0624954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0626559) q[1];
sx q[1];
rz(-1.0149983) q[1];
sx q[1];
rz(-1.009726) q[1];
rz(-0.40036277) q[3];
sx q[3];
rz(-1.6737079) q[3];
sx q[3];
rz(-2.3617552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3758731) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(-2.5326552) q[2];
rz(2.7385312) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8159863) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(2.0035279) q[0];
rz(-1.552938) q[1];
sx q[1];
rz(-2.373003) q[1];
sx q[1];
rz(2.3345711) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1083384) q[0];
sx q[0];
rz(-1.4716118) q[0];
sx q[0];
rz(1.5434815) q[0];
rz(-pi) q[1];
rz(2.6969063) q[2];
sx q[2];
rz(-1.7181601) q[2];
sx q[2];
rz(0.28648892) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5182636) q[1];
sx q[1];
rz(-2.9179472) q[1];
sx q[1];
rz(0.56648751) q[1];
x q[2];
rz(-1.1919954) q[3];
sx q[3];
rz(-1.8378432) q[3];
sx q[3];
rz(2.7821845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1300065) q[2];
sx q[2];
rz(-2.5277972) q[2];
sx q[2];
rz(-1.2137132) q[2];
rz(-0.064149292) q[3];
sx q[3];
rz(-1.4870653) q[3];
sx q[3];
rz(0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5478058) q[0];
sx q[0];
rz(-1.2603899) q[0];
sx q[0];
rz(-2.2547145) q[0];
rz(-2.9828494) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(-1.8633206) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87065164) q[0];
sx q[0];
rz(-0.88490153) q[0];
sx q[0];
rz(-2.565388) q[0];
rz(-pi) q[1];
rz(2.5910225) q[2];
sx q[2];
rz(-1.9801557) q[2];
sx q[2];
rz(1.0813724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63250971) q[1];
sx q[1];
rz(-2.4466568) q[1];
sx q[1];
rz(1.4768697) q[1];
x q[2];
rz(-2.9904891) q[3];
sx q[3];
rz(-1.0211111) q[3];
sx q[3];
rz(1.0719887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0025803) q[2];
sx q[2];
rz(-2.2860892) q[2];
sx q[2];
rz(-2.1878237) q[2];
rz(1.9468797) q[3];
sx q[3];
rz(-0.84269968) q[3];
sx q[3];
rz(-0.4755303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139528) q[0];
sx q[0];
rz(-2.4251745) q[0];
sx q[0];
rz(-2.5415976) q[0];
rz(2.4586239) q[1];
sx q[1];
rz(-0.78790793) q[1];
sx q[1];
rz(-0.33040985) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4415057) q[0];
sx q[0];
rz(-1.7528025) q[0];
sx q[0];
rz(0.84199961) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7960828) q[2];
sx q[2];
rz(-1.4116294) q[2];
sx q[2];
rz(2.988254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7327795) q[1];
sx q[1];
rz(-2.6169852) q[1];
sx q[1];
rz(-1.9636159) q[1];
x q[2];
rz(-2.0512625) q[3];
sx q[3];
rz(-0.76101979) q[3];
sx q[3];
rz(-2.7121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7052475) q[2];
sx q[2];
rz(-2.0630431) q[2];
sx q[2];
rz(2.516563) q[2];
rz(0.69508067) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(-0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92796749) q[0];
sx q[0];
rz(-1.9897505) q[0];
sx q[0];
rz(-1.3731765) q[0];
rz(-1.3881418) q[1];
sx q[1];
rz(-2.2699247) q[1];
sx q[1];
rz(-1.6067827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2508188) q[0];
sx q[0];
rz(-1.5355331) q[0];
sx q[0];
rz(-1.4232487) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18267554) q[2];
sx q[2];
rz(-1.8790882) q[2];
sx q[2];
rz(-1.8162372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0252367) q[1];
sx q[1];
rz(-0.80151171) q[1];
sx q[1];
rz(2.5610183) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0032773) q[3];
sx q[3];
rz(-1.9683016) q[3];
sx q[3];
rz(0.53526263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2024978) q[2];
sx q[2];
rz(-1.1700583) q[2];
sx q[2];
rz(-1.6004174) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(2.8405564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.130126) q[0];
sx q[0];
rz(-0.13700329) q[0];
sx q[0];
rz(-2.2504508) q[0];
rz(-2.4961684) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(2.3088764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37207212) q[0];
sx q[0];
rz(-1.7070878) q[0];
sx q[0];
rz(-2.6762677) q[0];
rz(-pi) q[1];
rz(2.1948692) q[2];
sx q[2];
rz(-1.6247617) q[2];
sx q[2];
rz(2.0651511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.32729808) q[1];
sx q[1];
rz(-2.6076064) q[1];
sx q[1];
rz(0.3898062) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1770958) q[3];
sx q[3];
rz(-1.0648921) q[3];
sx q[3];
rz(-0.33123744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58925313) q[2];
sx q[2];
rz(-0.66698843) q[2];
sx q[2];
rz(-1.5480631) q[2];
rz(-0.42406905) q[3];
sx q[3];
rz(-1.721902) q[3];
sx q[3];
rz(-0.29485318) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9645204) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(-0.82343423) q[0];
rz(1.1973165) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(1.6808602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4890503) q[0];
sx q[0];
rz(-0.33777896) q[0];
sx q[0];
rz(1.2342288) q[0];
rz(-pi) q[1];
rz(-1.8180188) q[2];
sx q[2];
rz(-2.5342016) q[2];
sx q[2];
rz(3.0957019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6951235) q[1];
sx q[1];
rz(-2.2317076) q[1];
sx q[1];
rz(1.6245232) q[1];
x q[2];
rz(-2.7815205) q[3];
sx q[3];
rz(-1.8723893) q[3];
sx q[3];
rz(-2.2693279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0921649) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(0.1235505) q[2];
rz(-0.83856797) q[3];
sx q[3];
rz(-1.1127915) q[3];
sx q[3];
rz(-0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(3.0274444) q[0];
sx q[0];
rz(-2.1338978) q[0];
sx q[0];
rz(-1.4134407) q[0];
rz(0.89883262) q[1];
sx q[1];
rz(-0.90679589) q[1];
sx q[1];
rz(0.59284219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956447) q[0];
sx q[0];
rz(-2.2586825) q[0];
sx q[0];
rz(-2.3700729) q[0];
rz(-1.4845092) q[2];
sx q[2];
rz(-1.3893428) q[2];
sx q[2];
rz(1.726651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.034380091) q[1];
sx q[1];
rz(-1.980956) q[1];
sx q[1];
rz(0.35663794) q[1];
x q[2];
rz(2.4580703) q[3];
sx q[3];
rz(-2.0647979) q[3];
sx q[3];
rz(-2.0919111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8574519) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(-2.9202374) q[2];
rz(1.0434693) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(2.7531457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8429883) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(-0.84841949) q[0];
rz(3.0449955) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(-2.3883147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83148513) q[0];
sx q[0];
rz(-1.3436755) q[0];
sx q[0];
rz(0.102553) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3123872) q[2];
sx q[2];
rz(-1.5952139) q[2];
sx q[2];
rz(1.8054192) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9946361) q[1];
sx q[1];
rz(-2.5038233) q[1];
sx q[1];
rz(-0.90177782) q[1];
x q[2];
rz(-1.1152335) q[3];
sx q[3];
rz(-2.2325745) q[3];
sx q[3];
rz(-1.7152804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0640556) q[2];
sx q[2];
rz(-2.3040743) q[2];
sx q[2];
rz(2.688664) q[2];
rz(2.5405267) q[3];
sx q[3];
rz(-1.4671003) q[3];
sx q[3];
rz(0.99630228) q[3];
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
rz(-0.30408981) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(-2.9651463) q[1];
sx q[1];
rz(-0.23527589) q[1];
sx q[1];
rz(-1.489524) q[1];
rz(0.34863451) q[2];
sx q[2];
rz(-1.5713816) q[2];
sx q[2];
rz(-1.1272507) q[2];
rz(2.1613359) q[3];
sx q[3];
rz(-0.93385812) q[3];
sx q[3];
rz(1.0360624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
