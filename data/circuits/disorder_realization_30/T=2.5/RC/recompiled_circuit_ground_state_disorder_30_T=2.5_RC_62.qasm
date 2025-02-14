OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(-0.77646065) q[0];
sx q[0];
rz(2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(3.4263098) q[1];
sx q[1];
rz(10.813536) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2251988) q[0];
sx q[0];
rz(-1.7728191) q[0];
sx q[0];
rz(1.6388157) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35438673) q[2];
sx q[2];
rz(-0.72438188) q[2];
sx q[2];
rz(-2.5222833) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60290135) q[1];
sx q[1];
rz(-0.60990342) q[1];
sx q[1];
rz(2.2944488) q[1];
rz(-2.0991058) q[3];
sx q[3];
rz(-1.0431223) q[3];
sx q[3];
rz(-1.8535525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(-0.93519768) q[2];
rz(2.8397371) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(0.39366084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1644156) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(-0.49829495) q[0];
rz(1.7138819) q[1];
sx q[1];
rz(-1.2063113) q[1];
sx q[1];
rz(-2.0548342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3918864) q[0];
sx q[0];
rz(-1.3343617) q[0];
sx q[0];
rz(2.4357585) q[0];
rz(-pi) q[1];
rz(-0.94765122) q[2];
sx q[2];
rz(-2.8171879) q[2];
sx q[2];
rz(1.0783006) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77176911) q[1];
sx q[1];
rz(-0.46716094) q[1];
sx q[1];
rz(-0.36231706) q[1];
rz(-pi) q[2];
rz(0.03647904) q[3];
sx q[3];
rz(-0.31841296) q[3];
sx q[3];
rz(2.1977001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9429417) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(1.8093713) q[2];
rz(-1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.3326299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20951095) q[0];
sx q[0];
rz(-1.7387583) q[0];
sx q[0];
rz(2.984356) q[0];
rz(1.9137742) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(-0.42207178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8764719) q[0];
sx q[0];
rz(-2.9419964) q[0];
sx q[0];
rz(1.100698) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7447168) q[2];
sx q[2];
rz(-1.8444159) q[2];
sx q[2];
rz(1.5652986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0737368) q[1];
sx q[1];
rz(-1.4826117) q[1];
sx q[1];
rz(1.9531308) q[1];
rz(-pi) q[2];
rz(-1.8881599) q[3];
sx q[3];
rz(-0.79447132) q[3];
sx q[3];
rz(2.1411856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.813039) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(2.4075244) q[2];
rz(2.0640533) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(-0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828736) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(-0.016481312) q[0];
rz(0.074542848) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(-0.25513908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48038424) q[0];
sx q[0];
rz(-2.0503649) q[0];
sx q[0];
rz(-1.40856) q[0];
rz(-pi) q[1];
rz(-2.071923) q[2];
sx q[2];
rz(-2.5146896) q[2];
sx q[2];
rz(1.7834341) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8447579) q[1];
sx q[1];
rz(-0.95176178) q[1];
sx q[1];
rz(2.1059787) q[1];
rz(0.11394088) q[3];
sx q[3];
rz(-1.9181532) q[3];
sx q[3];
rz(1.243227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61465803) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(-3.0497293) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(-0.45472586) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62234539) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(-1.1118332) q[0];
rz(0.19019292) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(1.1598738) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69731312) q[0];
sx q[0];
rz(-2.796058) q[0];
sx q[0];
rz(-1.4688501) q[0];
rz(-pi) q[1];
rz(-2.2091523) q[2];
sx q[2];
rz(-1.6886504) q[2];
sx q[2];
rz(-1.3502075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9638744) q[1];
sx q[1];
rz(-1.323112) q[1];
sx q[1];
rz(-2.5801093) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0489063) q[3];
sx q[3];
rz(-1.6097704) q[3];
sx q[3];
rz(-0.029595395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5708892) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-1.0620091) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(1.5084069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45288169) q[0];
sx q[0];
rz(-1.0588366) q[0];
sx q[0];
rz(1.3558615) q[0];
rz(0.52976766) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(1.1518325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36042155) q[0];
sx q[0];
rz(-2.7716405) q[0];
sx q[0];
rz(-0.2759224) q[0];
rz(1.7284313) q[2];
sx q[2];
rz(-0.54716483) q[2];
sx q[2];
rz(-2.5131651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63910145) q[1];
sx q[1];
rz(-1.3768973) q[1];
sx q[1];
rz(-0.37096937) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21422503) q[3];
sx q[3];
rz(-1.6856632) q[3];
sx q[3];
rz(1.8095995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6666224) q[2];
sx q[2];
rz(-2.5002561) q[2];
sx q[2];
rz(-0.46710157) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-0.34365383) q[3];
sx q[3];
rz(1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2924627) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(2.9169061) q[0];
rz(-0.16381964) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(0.64635578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689679) q[0];
sx q[0];
rz(-2.5494721) q[0];
sx q[0];
rz(2.0298241) q[0];
x q[1];
rz(1.8936526) q[2];
sx q[2];
rz(-1.3636929) q[2];
sx q[2];
rz(-0.32760179) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19749459) q[1];
sx q[1];
rz(-2.1131599) q[1];
sx q[1];
rz(-1.6587064) q[1];
x q[2];
rz(-0.38742806) q[3];
sx q[3];
rz(-1.6968075) q[3];
sx q[3];
rz(0.84539094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(2.3112467) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(-1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53772563) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(2.895288) q[0];
rz(-2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(2.8776317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4564698) q[0];
sx q[0];
rz(-1.5780266) q[0];
sx q[0];
rz(-1.6662706) q[0];
x q[1];
rz(-1.0954148) q[2];
sx q[2];
rz(-0.84566294) q[2];
sx q[2];
rz(2.1755168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.57264489) q[1];
sx q[1];
rz(-1.8518475) q[1];
sx q[1];
rz(-0.87139178) q[1];
rz(2.5313782) q[3];
sx q[3];
rz(-0.68192484) q[3];
sx q[3];
rz(2.1206926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0366514) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(-0.71511739) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(0.74688545) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320553) q[0];
sx q[0];
rz(-1.6522464) q[0];
sx q[0];
rz(2.706053) q[0];
rz(2.9938193) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(-1.4685644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4184121) q[0];
sx q[0];
rz(-1.7505129) q[0];
sx q[0];
rz(-0.43850684) q[0];
rz(-pi) q[1];
rz(0.69738241) q[2];
sx q[2];
rz(-1.5101523) q[2];
sx q[2];
rz(1.9549119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10994153) q[1];
sx q[1];
rz(-1.9579971) q[1];
sx q[1];
rz(-2.9613858) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6665523) q[3];
sx q[3];
rz(-0.44774017) q[3];
sx q[3];
rz(1.0645176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.396951) q[2];
sx q[2];
rz(-1.7640742) q[2];
sx q[2];
rz(-0.89293876) q[2];
rz(1.8630155) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(2.5355329) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(0.079924718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250076) q[0];
sx q[0];
rz(-1.8226624) q[0];
sx q[0];
rz(0.99117898) q[0];
rz(0.28744185) q[2];
sx q[2];
rz(-0.86311695) q[2];
sx q[2];
rz(1.3112932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.895782) q[1];
sx q[1];
rz(-1.5930098) q[1];
sx q[1];
rz(-2.7345045) q[1];
rz(-pi) q[2];
rz(-0.20584835) q[3];
sx q[3];
rz(-1.3449114) q[3];
sx q[3];
rz(-3.1354648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14110485) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(-2.1263988) q[2];
rz(0.60123932) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64423185) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(-2.5524706) q[1];
sx q[1];
rz(-0.67710572) q[1];
sx q[1];
rz(2.6893375) q[1];
rz(-2.6964006) q[2];
sx q[2];
rz(-1.4839076) q[2];
sx q[2];
rz(-2.1640726) q[2];
rz(2.3169869) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
