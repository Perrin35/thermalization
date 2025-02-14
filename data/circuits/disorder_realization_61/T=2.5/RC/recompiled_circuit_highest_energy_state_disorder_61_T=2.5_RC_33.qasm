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
rz(0.50245589) q[0];
sx q[0];
rz(-2.034019) q[0];
sx q[0];
rz(-1.7662319) q[0];
rz(-3.4604685) q[1];
sx q[1];
rz(4.7825216) q[1];
sx q[1];
rz(8.6934269) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3325746) q[0];
sx q[0];
rz(-1.2515731) q[0];
sx q[0];
rz(2.2126424) q[0];
rz(1.574975) q[2];
sx q[2];
rz(-0.90882817) q[2];
sx q[2];
rz(-0.83534345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13818422) q[1];
sx q[1];
rz(-0.78552442) q[1];
sx q[1];
rz(1.5064606) q[1];
x q[2];
rz(1.4764686) q[3];
sx q[3];
rz(-1.7324311) q[3];
sx q[3];
rz(-1.9548655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2607164) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-2.1991849) q[2];
rz(0.85855329) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(-0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-3.0829561) q[0];
rz(-3.0851641) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(1.1172392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65172136) q[0];
sx q[0];
rz(-2.3014223) q[0];
sx q[0];
rz(2.1108759) q[0];
rz(1.1740721) q[2];
sx q[2];
rz(-2.8182321) q[2];
sx q[2];
rz(-2.4692995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5176425) q[1];
sx q[1];
rz(-0.25282598) q[1];
sx q[1];
rz(1.235591) q[1];
rz(-pi) q[2];
rz(1.4499217) q[3];
sx q[3];
rz(-1.0981005) q[3];
sx q[3];
rz(-2.8877088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(2.7351725) q[2];
rz(2.8255919) q[3];
sx q[3];
rz(-1.6812811) q[3];
sx q[3];
rz(2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(-1.8055441) q[0];
rz(0.27851963) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(0.96885243) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1135546) q[0];
sx q[0];
rz(-2.4062706) q[0];
sx q[0];
rz(0.57082973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79040852) q[2];
sx q[2];
rz(-1.9876754) q[2];
sx q[2];
rz(1.1000772) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90889764) q[1];
sx q[1];
rz(-1.8669584) q[1];
sx q[1];
rz(-3.0975268) q[1];
rz(-pi) q[2];
rz(-2.8519451) q[3];
sx q[3];
rz(-2.395219) q[3];
sx q[3];
rz(-2.1044855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99732533) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(-0.17511314) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(-1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-1.0293707) q[0];
rz(-0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(0.79684657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40596698) q[0];
sx q[0];
rz(-1.7628364) q[0];
sx q[0];
rz(-2.6940932) q[0];
rz(-pi) q[1];
rz(1.5810851) q[2];
sx q[2];
rz(-0.52797374) q[2];
sx q[2];
rz(2.1473403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0681035) q[1];
sx q[1];
rz(-1.3655546) q[1];
sx q[1];
rz(1.2935335) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6096973) q[3];
sx q[3];
rz(-2.4060156) q[3];
sx q[3];
rz(-0.07594219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13021079) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(-1.2288564) q[2];
rz(0.43904385) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(-1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-0.80832344) q[0];
rz(-1.7841548) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(-2.435991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3177174) q[0];
sx q[0];
rz(-0.30254811) q[0];
sx q[0];
rz(3.0251193) q[0];
x q[1];
rz(-1.8417435) q[2];
sx q[2];
rz(-0.72110211) q[2];
sx q[2];
rz(-1.9443898) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8562634) q[1];
sx q[1];
rz(-1.8239905) q[1];
sx q[1];
rz(-0.002753792) q[1];
rz(-pi) q[2];
rz(0.6198778) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(-1.0244964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48050532) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(0.60574469) q[2];
rz(-3.0211966) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-0.81407434) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(0.0023728097) q[0];
rz(-2.782605) q[1];
sx q[1];
rz(-1.2009883) q[1];
sx q[1];
rz(0.40828362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9914382) q[0];
sx q[0];
rz(-1.6065734) q[0];
sx q[0];
rz(1.4442025) q[0];
rz(-0.26770182) q[2];
sx q[2];
rz(-2.1609339) q[2];
sx q[2];
rz(-0.4877643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.017458113) q[1];
sx q[1];
rz(-1.6377021) q[1];
sx q[1];
rz(2.9111039) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5950795) q[3];
sx q[3];
rz(-2.236955) q[3];
sx q[3];
rz(0.72197589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4957054) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(3.0118946) q[2];
rz(0.4194704) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(-2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3892589) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(0.93841249) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(-1.2260812) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48213595) q[0];
sx q[0];
rz(-0.99471751) q[0];
sx q[0];
rz(-2.0123737) q[0];
x q[1];
rz(-1.3834882) q[2];
sx q[2];
rz(-0.97895998) q[2];
sx q[2];
rz(-2.1746152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6437373) q[1];
sx q[1];
rz(-1.8561991) q[1];
sx q[1];
rz(1.1021815) q[1];
rz(2.5874373) q[3];
sx q[3];
rz(-0.99923493) q[3];
sx q[3];
rz(0.038427834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40743264) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(-0.70719353) q[2];
rz(-1.4736942) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(-1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4555175) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(1.0651945) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-2.0069897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77662879) q[0];
sx q[0];
rz(-2.7641649) q[0];
sx q[0];
rz(-1.6799742) q[0];
x q[1];
rz(0.61750268) q[2];
sx q[2];
rz(-0.98101014) q[2];
sx q[2];
rz(1.2794354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5091619) q[1];
sx q[1];
rz(-2.1010509) q[1];
sx q[1];
rz(3.0796618) q[1];
rz(-1.9530961) q[3];
sx q[3];
rz(-2.144882) q[3];
sx q[3];
rz(2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(0.17063394) q[3];
sx q[3];
rz(-0.65244397) q[3];
sx q[3];
rz(-1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(2.3353031) q[0];
rz(-3.131648) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(2.8795805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6638787) q[0];
sx q[0];
rz(-2.0933525) q[0];
sx q[0];
rz(-2.9502003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4100894) q[2];
sx q[2];
rz(-0.93148684) q[2];
sx q[2];
rz(-0.29603816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1478473) q[1];
sx q[1];
rz(-2.6585852) q[1];
sx q[1];
rz(1.1149393) q[1];
rz(2.0724107) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(0.25267664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.749873) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(-0.14287512) q[2];
rz(-0.90397942) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-0.8814632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7664465) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(-0.65650702) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(2.5118929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5507767) q[0];
sx q[0];
rz(-1.970298) q[0];
sx q[0];
rz(2.4062064) q[0];
rz(2.0067425) q[2];
sx q[2];
rz(-2.5647047) q[2];
sx q[2];
rz(-2.1132121) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7460664) q[1];
sx q[1];
rz(-1.6159355) q[1];
sx q[1];
rz(3.1137115) q[1];
rz(2.0721758) q[3];
sx q[3];
rz(-2.5540941) q[3];
sx q[3];
rz(-1.0621493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730561) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(0.14280351) q[1];
sx q[1];
rz(-1.0719943) q[1];
sx q[1];
rz(0.73678585) q[1];
rz(-1.8129195) q[2];
sx q[2];
rz(-1.6264781) q[2];
sx q[2];
rz(-1.3774058) q[2];
rz(-0.96137267) q[3];
sx q[3];
rz(-1.6037887) q[3];
sx q[3];
rz(-2.9870839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
