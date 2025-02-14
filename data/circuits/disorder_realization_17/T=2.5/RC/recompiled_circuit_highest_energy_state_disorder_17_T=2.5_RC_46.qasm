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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(0.59803522) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(-0.66507566) q[1];
sx q[1];
rz(0.15288487) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4476112) q[0];
sx q[0];
rz(-1.6504399) q[0];
sx q[0];
rz(-1.682974) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48798765) q[2];
sx q[2];
rz(-1.5190912) q[2];
sx q[2];
rz(2.0590797) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26273706) q[1];
sx q[1];
rz(-1.3595743) q[1];
sx q[1];
rz(0.12843522) q[1];
rz(-pi) q[2];
rz(-1.5678207) q[3];
sx q[3];
rz(-1.6167534) q[3];
sx q[3];
rz(0.19356914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87024706) q[2];
sx q[2];
rz(-2.953244) q[2];
sx q[2];
rz(1.1790454) q[2];
rz(2.7315268) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(2.2470391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1882741) q[0];
sx q[0];
rz(-2.4725547) q[0];
sx q[0];
rz(-2.8642995) q[0];
rz(0.20031985) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(0.083981363) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7070933) q[0];
sx q[0];
rz(-1.7696119) q[0];
sx q[0];
rz(2.0484224) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72767392) q[2];
sx q[2];
rz(-1.6817589) q[2];
sx q[2];
rz(-0.79906711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0711956) q[1];
sx q[1];
rz(-2.4807811) q[1];
sx q[1];
rz(0.7699716) q[1];
rz(-1.3091607) q[3];
sx q[3];
rz(-2.0896834) q[3];
sx q[3];
rz(-0.46552697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1818485) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(2.3954771) q[2];
rz(-0.6212081) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(1.752689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8212432) q[0];
sx q[0];
rz(-2.4140883) q[0];
sx q[0];
rz(-1.2028836) q[0];
rz(-2.7299643) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(-2.0278377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8343415) q[0];
sx q[0];
rz(-1.1528307) q[0];
sx q[0];
rz(-0.78180914) q[0];
rz(2.162777) q[2];
sx q[2];
rz(-2.4360058) q[2];
sx q[2];
rz(-1.6340812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74677709) q[1];
sx q[1];
rz(-1.703843) q[1];
sx q[1];
rz(0.89292553) q[1];
rz(-2.5630234) q[3];
sx q[3];
rz(-2.1820543) q[3];
sx q[3];
rz(-1.8139773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0959629) q[2];
sx q[2];
rz(-2.3489504) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(0.80471936) q[3];
sx q[3];
rz(-1.2730803) q[3];
sx q[3];
rz(1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3941037) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.6266741) q[0];
rz(-0.50846848) q[1];
sx q[1];
rz(-2.4912806) q[1];
sx q[1];
rz(1.633684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8430458) q[0];
sx q[0];
rz(-0.15104476) q[0];
sx q[0];
rz(1.5848978) q[0];
x q[1];
rz(2.1437683) q[2];
sx q[2];
rz(-2.8854239) q[2];
sx q[2];
rz(-1.0333956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3292309) q[1];
sx q[1];
rz(-0.32608247) q[1];
sx q[1];
rz(-1.2124168) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0359457) q[3];
sx q[3];
rz(-0.7300762) q[3];
sx q[3];
rz(-1.4393782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7278829) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(2.3940562) q[2];
rz(1.9109292) q[3];
sx q[3];
rz(-0.80941284) q[3];
sx q[3];
rz(0.49720732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43664765) q[0];
sx q[0];
rz(-2.0132988) q[0];
sx q[0];
rz(1.6823912) q[0];
rz(-2.7875426) q[1];
sx q[1];
rz(-2.470128) q[1];
sx q[1];
rz(-0.57463247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.282772) q[0];
sx q[0];
rz(-3.0631815) q[0];
sx q[0];
rz(2.9693102) q[0];
x q[1];
rz(2.3726497) q[2];
sx q[2];
rz(-1.3731602) q[2];
sx q[2];
rz(-0.96425024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71855356) q[1];
sx q[1];
rz(-1.3438112) q[1];
sx q[1];
rz(-1.5178229) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7754885) q[3];
sx q[3];
rz(-1.4581956) q[3];
sx q[3];
rz(-1.7623869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6358801) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(-0.20225784) q[2];
rz(1.0328736) q[3];
sx q[3];
rz(-0.60690108) q[3];
sx q[3];
rz(-0.066298299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031484) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(0.37145823) q[0];
rz(0.29608852) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(-0.16337005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0609157) q[0];
sx q[0];
rz(-1.0759584) q[0];
sx q[0];
rz(-1.3617152) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67739112) q[2];
sx q[2];
rz(-2.4755726) q[2];
sx q[2];
rz(2.1474554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5886053) q[1];
sx q[1];
rz(-1.7514015) q[1];
sx q[1];
rz(0.44579472) q[1];
rz(1.3903687) q[3];
sx q[3];
rz(-2.437915) q[3];
sx q[3];
rz(1.1168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0801733) q[2];
sx q[2];
rz(-2.8714608) q[2];
sx q[2];
rz(-0.6529676) q[2];
rz(2.669892) q[3];
sx q[3];
rz(-0.74840778) q[3];
sx q[3];
rz(-0.72699839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2111874) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(-2.1042673) q[0];
rz(0.51517454) q[1];
sx q[1];
rz(-1.2778792) q[1];
sx q[1];
rz(1.3264664) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68149306) q[0];
sx q[0];
rz(-1.1837479) q[0];
sx q[0];
rz(2.1429064) q[0];
rz(0.35197978) q[2];
sx q[2];
rz(-1.5328141) q[2];
sx q[2];
rz(0.69881907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1465254) q[1];
sx q[1];
rz(-0.77153782) q[1];
sx q[1];
rz(3.0027215) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8765728) q[3];
sx q[3];
rz(-1.631284) q[3];
sx q[3];
rz(2.4779589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4482164) q[2];
sx q[2];
rz(-2.6770112) q[2];
sx q[2];
rz(-2.8153937) q[2];
rz(-2.163573) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(-3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4433032) q[0];
sx q[0];
rz(-0.25743085) q[0];
sx q[0];
rz(2.8522016) q[0];
rz(-3.1293213) q[1];
sx q[1];
rz(-0.27996501) q[1];
sx q[1];
rz(1.4853005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3607583) q[0];
sx q[0];
rz(-0.88874431) q[0];
sx q[0];
rz(-2.2498796) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2337306) q[2];
sx q[2];
rz(-0.22960358) q[2];
sx q[2];
rz(-0.52403852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0721133) q[1];
sx q[1];
rz(-0.24976191) q[1];
sx q[1];
rz(0.62432557) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9910219) q[3];
sx q[3];
rz(-0.54383792) q[3];
sx q[3];
rz(0.71292927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54899186) q[2];
sx q[2];
rz(-1.5697378) q[2];
sx q[2];
rz(-0.17295095) q[2];
rz(-1.3744099) q[3];
sx q[3];
rz(-1.3783312) q[3];
sx q[3];
rz(-1.4779199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7462815) q[0];
sx q[0];
rz(-0.44773856) q[0];
sx q[0];
rz(-0.81004274) q[0];
rz(0.41656247) q[1];
sx q[1];
rz(-0.77605334) q[1];
sx q[1];
rz(0.3449482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5589011) q[0];
sx q[0];
rz(-1.4479637) q[0];
sx q[0];
rz(-0.71101153) q[0];
rz(-pi) q[1];
rz(2.9760804) q[2];
sx q[2];
rz(-1.8461746) q[2];
sx q[2];
rz(-2.2294758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7899985) q[1];
sx q[1];
rz(-1.8934544) q[1];
sx q[1];
rz(-2.0568454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79699253) q[3];
sx q[3];
rz(-2.9755962) q[3];
sx q[3];
rz(2.5696563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42744669) q[2];
sx q[2];
rz(-0.88627187) q[2];
sx q[2];
rz(0.084058849) q[2];
rz(-2.2345624) q[3];
sx q[3];
rz(-1.4182988) q[3];
sx q[3];
rz(0.9575873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8091938) q[0];
sx q[0];
rz(-3.0539303) q[0];
sx q[0];
rz(0.3122538) q[0];
rz(2.9088083) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(1.8544474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9159661) q[0];
sx q[0];
rz(-0.73333987) q[0];
sx q[0];
rz(-2.649113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5239117) q[2];
sx q[2];
rz(-2.2094036) q[2];
sx q[2];
rz(-2.435911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97032769) q[1];
sx q[1];
rz(-1.1210151) q[1];
sx q[1];
rz(2.8697017) q[1];
rz(-0.49317499) q[3];
sx q[3];
rz(-0.68636591) q[3];
sx q[3];
rz(0.92991352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9707668) q[2];
sx q[2];
rz(-1.8584741) q[2];
sx q[2];
rz(-1.4975366) q[2];
rz(-1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(2.6523377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.004414) q[0];
sx q[0];
rz(-1.351384) q[0];
sx q[0];
rz(2.2030892) q[0];
rz(1.3068403) q[1];
sx q[1];
rz(-2.2521781) q[1];
sx q[1];
rz(2.3440012) q[1];
rz(-2.2597792) q[2];
sx q[2];
rz(-1.6826993) q[2];
sx q[2];
rz(-0.8347389) q[2];
rz(2.9994503) q[3];
sx q[3];
rz(-2.2323043) q[3];
sx q[3];
rz(-1.4122152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
